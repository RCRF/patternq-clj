(ns patternq.db
  "Read-only access to Pattern Data Commons dataset databases.

  Every dataset is its own database. `db` returns a database handle for a
  database name, and `q` runs Datomic Datalog against a handle; every canned
  query in patternq goes through `q`, so all functions work over both
  transports:

  - :http (default): the Pattern Data Commons query service (the same one
    the R, Python and Julia libraries use; see patternq.http). Needs only
    PATTERNQ_API_KEY (and optionally PATTERNQ_ENDPOINT). The handle is a
    patternq.db.HttpDb record.
  - :peer: the Datomic peer, when it is on the classpath (the :peer deps
    alias) and PATTERNQ_DATOMIC_URI names the storage base URI. The handle
    is a Datomic database value.

  Select with the PATTERNQ_TRANSPORT env var (\"http\" or \"peer\") or
  `set-transport!`. Requesting :peer without the peer on the classpath, or
  without a storage URI, is an error (no silent fallback).

  patternq is strictly read-only: it never transacts, creates, renames or
  deletes databases, and never connects to the admin database.")

(set! *warn-on-reflection* true)

(def ^:private admin-db-name "admin-db-1")

;; -- transport selection --

(defonce ^:private transport-override (atom nil))

(defn transport
  "The transport in use: :http (default) or :peer."
  []
  (or @transport-override
      (case (some-> (System/getenv "PATTERNQ_TRANSPORT") clojure.string/lower-case)
        "peer" :peer
        :http)))

(defn set-transport!
  "Set the transport for this process (:http or :peer); nil restores the
  env var / default."
  [t]
  (when-not (contains? #{:http :peer nil} t)
    (throw (ex-info "transport must be :http or :peer" {:transport t})))
  (reset! transport-override t))

;; -- peer (optional) --

(defonce ^:private base-uri-override (atom nil))

(defn base-uri
  "Peer storage base URI in use (PATTERNQ_DATOMIC_URI or set-base-uri!), or
  nil when unset."
  []
  (when-let [uri (or @base-uri-override (System/getenv "PATTERNQ_DATOMIC_URI"))]
    (if (clojure.string/ends-with? uri "/") uri (str uri "/"))))

(defn set-base-uri!
  "Set the peer storage base URI for this process; nil restores the env var."
  [uri]
  (reset! base-uri-override uri))

(defn- peer-fn
  "Resolve a read-only datomic.api fn, with a clear error when the peer is not
  on the classpath."
  [sym]
  (or (try (requiring-resolve (symbol "datomic.api" (name sym)))
           (catch Exception _ nil))
      (throw (ex-info (str "PATTERNQ_TRANSPORT=peer but the Datomic peer is not on the classpath; "
                           "add the :peer alias (clojure -M:peer ...) or use the default HTTP transport.")
                      {:patternq/error :peer-unavailable}))))

(defn peer-available?
  "Is the Datomic peer library on the classpath?"
  []
  (boolean (try (requiring-resolve 'datomic.api/q) (catch Exception _ nil))))

(defn db-uri
  [db-name]
  (when (= admin-db-name db-name)
    (throw (ex-info "patternq does not access the admin database." {:db-name db-name})))
  (if-let [base (base-uri)]
    (str base db-name)
    (throw (ex-info "Peer transport needs the storage base URI: set PATTERNQ_DATOMIC_URI (or patternq.db/set-base-uri!)."
                    {:patternq/error :peer-uri-unset}))))

;; -- HTTP handle --

(defrecord HttpDb [db-name basis-t])

(defn http-db? [x] (instance? HttpDb x))

(defn db
  "Database handle for the dataset database `db-name` (e.g.
  \"H37001-2026-09-21a\") on the current transport."
  [db-name]
  (when (= admin-db-name db-name)
    (throw (ex-info "patternq does not access the admin database." {:db-name db-name})))
  (case (transport)
    :http (->HttpDb db-name (atom nil))
    :peer (let [connect (peer-fn 'connect) dbf (peer-fn 'db)]
            (dbf (connect (db-uri db-name))))))

(defn as-db
  "Accepts a database handle or a database name."
  [db-or-name]
  (if (string? db-or-name) (db db-or-name) db-or-name))

(defn db-name
  "Database name of a handle, when known (HTTP handles and names)."
  [db-or-name]
  (cond (string? db-or-name) db-or-name
        (http-db? db-or-name) (:db-name db-or-name)
        :else nil))

;; -- queries --

(defn- find-spec
  "Normalize a query (vector or map form) into a map with keyword keys."
  [query]
  (if (map? query)
    query
    (loop [[x & more] query k nil out {}]
      (cond (nil? x) (if (seq more) (recur more k out) out)
            (keyword? x) (recur more x (assoc out x []))
            :else (recur more k (update out k conj x))))))

(defn q
  "Run a Datomic query (vector or map form, as for datomic.api/q) against a
  database handle (or name) with args, on either transport. Results have the
  same shape as datomic.api/q: relations are collections of vectors; the
  find specs [?x ...], ?x . and [?a ?b] return collections, scalars and
  tuples."
  [query db-or-name & args]
  (let [db (as-db db-or-name)]
    (if (http-db? db)
      ((requiring-resolve 'patternq.http/q) (find-spec query) db args)
      (apply (peer-fn 'q) query db args))))

(defn pull-many
  "Pull `pattern` for entity ids `eids` (a vector of maps, in no particular
  order for HTTP; same order as eids for the peer)."
  [db-or-name pattern eids]
  (let [db (as-db db-or-name)]
    (if (http-db? db)
      (let [by-id (into {} (map (fn [[m]] [(:db/id m) m]))
                        (q {:find [(list 'pull '?e (vec (cons :db/id pattern)))] :in '[$ [?e ...]]
                            :where '[[?e]]}
                           db (vec eids)))]
        (vec (keep by-id eids)))
      ((peer-fn 'pull-many) db pattern eids))))

(defn pull
  "Pull `pattern` for one entity id or lookup ref ([:attr value])."
  [db-or-name pattern eid]
  (let [db (as-db db-or-name)]
    (if (http-db? db)
      (if (vector? eid)
        (ffirst (q {:find [(list 'pull '?e pattern)] :in '[$ ?v] :where [['?e (first eid) '?v]]}
                   db (second eid)))
        (first (pull-many db pattern [eid])))
      ((peer-fn 'pull) db pattern eid))))

(defn pull-by
  "{value pulled-map} for the entities whose `attr` is one of `values`
  (one query on either transport), e.g. (pull-by db pattern :variant/id ids)."
  [db-or-name pattern attr values]
  (into {} (q {:find ['?v (list 'pull '?e pattern)] :in '[$ [?v ...]] :where [['?e attr '?v]]}
              db-or-name (vec values))))

(defn pull-map
  "{eid pulled-map} for entity ids (one query on either transport)."
  [db-or-name pattern eids]
  (into {} (map (juxt :db/id identity)) (pull-many db-or-name (vec (distinct (cons :db/id pattern))) (vec eids))))

(defn basis-t
  "Basis t of a handle (for HTTP: of the latest query run through it)."
  [db]
  (if (http-db? db) @(:basis-t db) ((peer-fn 'basis-t) db)))

(defn- db-name* [db] (when (http-db? db) (:db-name db)))

(defn provenance
  "Provenance of a query against `db`: database name (if known), basis t,
  timestamp. Attached as metadata to results of patternq canned queries."
  ([db] (provenance db nil))
  ([db db-name]
   {:db (or db-name (db-name* db))
    :basis-t (basis-t db)
    :timestamp (str (java.time.LocalDateTime/now))}))

(defn with-provenance
  [result db db-name]
  (if (instance? clojure.lang.IObj result)
    (with-meta result {:patternq/provenance (provenance db db-name)})
    result))

(defn across-dbs
  "Run `f` (a fn of a database name, returning a seq of maps) against each
  database name and concatenate the results, adding :db to each row. Each
  dataset is its own database, so cohort comparisons run per database."
  [f db-names]
  (vec (mapcat (fn [n] (map #(assoc % :db n) (f n))) db-names)))
