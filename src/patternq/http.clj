(ns patternq.http
  "The Pattern Data Commons HTTP API: the default query transport (see
  patternq.db/q), dataset listing, dataset-name -> database-name resolution,
  and measurement matrix file downloads. Uses PATTERNQ_ENDPOINT /
  PATTERNQ_API_KEY like the R, Python and Julia libraries.

  The query service accepts queries only in the JSON form parsed by its
  datalog-json-parser, so Datalog written as Clojure data is converted
  client-side (`query->wire`): symbols become \"?x\" / \"_\" / \"...\" / \"$\"
  strings, keywords \":ns/name\" strings, lists and vectors arrays, pull
  pattern maps objects. Results are converted back so they are identical to
  the Datomic peer's (\":ns/name\" strings -> keywords, find specs
  [?x ...] / ?x . / [?a ?b] reshaped client-side).

  Only POST /query, POST /matrix and GET /api-v1/list are used: read-only."
  (:require [charred.api :as json]
            [clojure.string :as str]
            [clojure.walk :as walk])
  (:import (java.io ByteArrayInputStream InputStream)
           (java.net URI)
           (java.net.http HttpClient HttpRequest HttpRequest$BodyPublishers HttpResponse HttpResponse$BodyHandlers)
           (java.time Duration)
           (java.util.zip GZIPInputStream)))

(set! *warn-on-reflection* true)

(def default-endpoint "https://data-commons.rcrf-dev.org")

(declare endpoint-override)

(defn endpoint []
  (let [env (System/getenv "PATTERNQ_ENDPOINT")]
    (str/replace (or @endpoint-override (if (and env (.startsWith ^String env "http")) env default-endpoint))
                 #"/+$" "")))

(defonce ^:private token-override (atom nil))
(defonce ^:private endpoint-override (atom nil))

(defn set-token!
  "Set the API token for this process (overrides PATTERNQ_API_KEY); nil restores the env var."
  [token]
  (reset! token-override token))

(defn set-endpoint!
  "Set the query service endpoint for this process; nil restores PATTERNQ_ENDPOINT / the default."
  [url]
  (reset! endpoint-override url))

(defn- api-token []
  (or @token-override
      (System/getenv "PATTERNQ_API_KEY")
      (throw (ex-info "No API token: set PATTERNQ_API_KEY in the environment (or patternq.http/set-token!)." {}))))

(defonce ^:private client
  ;; one client for the process: connection (and TLS session) reuse
  (delay (-> (HttpClient/newBuilder)
             (.connectTimeout (Duration/ofSeconds 30))
             (.build))))

(defn- get-json [path]
  (let [^HttpClient client @client
        req (-> (HttpRequest/newBuilder (URI. (str (endpoint) path)))
                (.header "Authorization" (str "Bearer " (api-token)))
                (.header "Accept" "application/json")
                (.timeout (Duration/ofSeconds 60))
                (.GET)
                (.build))
        ^HttpResponse resp (.send client req (HttpResponse$BodyHandlers/ofString))]
    (when-not (= 200 (.statusCode resp))
      (throw (ex-info (str "HTTP " (.statusCode resp) " from " path)
                      {:status (.statusCode resp) :body (.body resp)})))
    (json/read-json (.body resp))))

(defn list-datasets
  "Datasets available to your API key: maps of :dataset, :db (current database
  name), :patient-count, :sample-count, :assays, :tags."
  []
  (->> (get (get-json "/api-v1/list/datasets") "datasets")
       (mapv (fn [d]
               {:dataset (get d "dataset/name")
                :db (get-in d ["dataset/database" "database/name"])
                :patient-count (get d "dataset/patient-count")
                :sample-count (get d "dataset/sample-count")
                :assays (vec (get d "dataset/assays"))
                :tags (vec (get d "dataset/tags"))}))))

(defn resolve-db
  "Dataset name (stable, e.g. \"tcga-uvm\") -> current database name (changes
  on re-import). A database name is passed through."
  [dataset]
  (let [ds (list-datasets)]
    (or (some #(when (= dataset (:dataset %)) (:db %)) ds)
        (some #(when (= dataset (:db %)) (:db %)) ds)
        (throw (ex-info (str "Unknown dataset " dataset) {:dataset dataset})))))

(defn- ^HttpResponse send! [^HttpRequest req handler]
  (.send ^HttpClient @client req handler))

(defn- gunzip-if-needed ^bytes [^bytes bs]
  (if (and (> (alength bs) 1) (= (aget bs 0) (unchecked-byte 0x1f)) (= (aget bs 1) (unchecked-byte 0x8b)))
    (with-open [in (GZIPInputStream. (ByteArrayInputStream. bs))]
      (.readAllBytes in))
    bs))

(defn measurement-matrix
  "Download a measurement matrix (e.g. single-cell counts) by backing-file key
  (see patternq.dataset/measurement-matrices) through the query service; the
  service returns a presigned S3 URL. Returns {:columns [...] :rows [[...]]}
  of the TSV (strings)."
  [db-name matrix-key]
  (let [req (-> (HttpRequest/newBuilder (URI. (str (endpoint) "/matrix/" db-name "/" matrix-key)))
                (.header "Authorization" (str "Bearer " (api-token)))
                (.header "Accept" "text/plain")
                (.header "Content-Type" "application/json")
                (.timeout (Duration/ofSeconds 120))
                (.POST (HttpRequest$BodyPublishers/ofString "{}"))
                (.build))
        ^HttpResponse resp (send! req (HttpResponse$BodyHandlers/ofString))]
    (when-not (= 200 (.statusCode resp))
      (throw (ex-info (str "HTTP " (.statusCode resp) " requesting matrix") {:body (.body resp)})))
    (let [url (str/trim (str (.body resp)))
          ^HttpResponse data (send! (-> (HttpRequest/newBuilder (URI. url)) (.GET) (.build))
                                    (HttpResponse$BodyHandlers/ofByteArray))
          text (String. (gunzip-if-needed ^bytes (.body data)) "UTF-8")
          [header & rows] (json/read-csv text :separator \tab)]
      {:columns (vec header) :rows (vec rows)})))

;; -- queries --

(def ^:dynamic *timeout-ms*
  "Query timeout (ms) sent to the service."
  120000)

(def ^:dynamic *cache*
  "Use the service's S3 result cache (Accept text/plain: presigned URL of the
  gzipped cached result, computed on a miss). false: result inline as JSON,
  cache skipped."
  true)

(def ^:dynamic *refresh-cache*
  "Recompute and overwrite the cached result."
  false)

(defn- wire
  "Clojure Datalog data -> JSON wire form."
  [x]
  (cond (symbol? x) (str x)
        (keyword? x) (str x)
        (map? x) (into {} (map (fn [[k v]] [(if (keyword? k) (str k) (str k)) (wire v)])) x)
        (or (sequential? x) (set? x)) (mapv wire x)
        :else x))

(defn- find-shape
  "[relation-find-elems reshape-fn] for a :find spec."
  [find]
  (cond
    (and (= 2 (count find)) (= '. (second find)))
    [[(first find)] (fn [rows] (ffirst rows))]
    (and (= 1 (count find)) (vector? (first find)) (= '... (last (first find))))
    [[(ffirst find)] (fn [rows] (mapv first rows))]
    (and (= 1 (count find)) (vector? (first find)))
    [(first find) (fn [rows] (first rows))]
    :else [find vec]))

(defn- inline-keyword-args
  "Keyword scalar args (idents) can't be sent as JSON args (they would arrive
  as strings); inline them into the query as literals."
  [{:keys [in] :as query} args]
  (if-not in
    [query args]
    (let [ins (vec in)
          has-db? (= '$ (first ins))
          arg-ins (if has-db? (subvec ins 1) ins)
          pairs (map vector arg-ins args)
          inline (into {} (keep (fn [[v a]] (when (and (symbol? v) (keyword? a)) [v a]))) pairs)
          kept (remove (fn [[v _]] (contains? inline v)) pairs)]
      (if (empty? inline)
        [query args]
        [(-> query
             (assoc :in (vec (concat (when has-db? ['$]) (map first kept))))
             (update :where #(walk/postwalk-replace inline %))
             (update :find #(walk/postwalk-replace inline %)))
         (mapv second kept)]))))

(def ^:private ident-re #"^:[A-Za-z*+!_?-][A-Za-z0-9*+!_?.-]*/[A-Za-z0-9*+!_?.<>=-]+$")

(defn- from-wire
  "JSON result -> peer-shaped Clojure data: \":ns/name\" strings (map keys,
  idents) -> keywords."
  [x]
  (walk/postwalk
    (fn [v]
      (cond (and (string? v) (re-matches ident-re v)) (keyword (subs v 1))
            (map-entry? v) v
            :else v))
    x))

(defn query->wire
  "The JSON body sent to POST /query/<db> for a query map (as from
  patternq.db/find-spec) and args."
  [query args]
  (let [[query args] (inline-keyword-args query args)
        [find-elems _] (find-shape (:find query))
        q (cond-> (assoc query :find find-elems)
            (and (:in query) (not= '$ (first (:in query)))) (update :in #(vec (cons '$ %))))]
    (cond-> {"query" (into {} (map (fn [[k v]] [(str k) (wire v)])) q)
             "timeout" *timeout-ms*}
      (seq args) (assoc "args" (wire (vec args)))
      *refresh-cache* (assoc "refresh-cache" true))))

(defn q
  "Run a query (map form) against an HTTP database handle
  (patternq.db.HttpDb) with args; see patternq.db/q."
  [query db args]
  (let [body (json/write-json-str (query->wire query args))
        [_ reshape] (find-shape (:find query))
        req (-> (HttpRequest/newBuilder (URI. (str (endpoint) "/query/" (:db-name db))))
                (.header "Authorization" (str "Bearer " (api-token)))
                (.header "Accept" (if *cache* "text/plain" "application/json"))
                (.header "Content-Type" "application/json")
                (.timeout (Duration/ofMillis (+ (long *timeout-ms*) 60000)))
                (.POST (HttpRequest$BodyPublishers/ofString body))
                (.build))
        ^HttpResponse resp (send! req (HttpResponse$BodyHandlers/ofString))
        payload (str/trim (str (.body resp)))]
    (when-not (= 200 (.statusCode resp))
      (throw (ex-info (str "Query failed (HTTP " (.statusCode resp) "): " payload)
                      {:status (.statusCode resp) :body payload :db (:db-name db)})))
    (let [res (if (str/starts-with? payload "{")
                (json/read-json payload)
                (let [^HttpResponse data (send! (-> (HttpRequest/newBuilder (URI. payload)) (.GET) (.build))
                                                (HttpResponse$BodyHandlers/ofByteArray))]
                  (json/read-json (String. (gunzip-if-needed ^bytes (.body data)) "UTF-8"))))]
      (when-let [err (get res "error")]
        (throw (ex-info (str "Query error: " err) {:db (:db-name db) :error err})))
      (reset! (:basis-t db) (get res "basis_t"))
      (reshape (from-wire (get res "query_result"))))))
