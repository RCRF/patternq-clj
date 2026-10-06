(ns patternq.backpressure-test
  "API back-pressure contract (unify-central dev-docs/API_BACKPRESSURE.md),
  against a local stub server."
  (:require [charred.api :as json]
            [clojure.test :refer [deftest is testing use-fixtures]]
            [patternq.backpressure :as bp]
            [patternq.db :as pdb]
            [patternq.http :as http])
  (:import (com.sun.net.httpserver HttpExchange HttpHandler HttpServer)
           (java.net InetSocketAddress)
           (java.util.concurrent Executors)))

(def problem bp/problem-prefix)

(def ok [200 {"Content-Type" "application/json" "RateLimit" "\"api\";r=99;t=0"}
         {"query_result" [[1]] "basis_t" 7}])

(defn throttled [& {:keys [kind status retry-after retryable]
                    :or {kind "rate-limited" status 429 retryable true}}]
  [status (cond-> {"Content-Type" "application/problem+json"}
            retry-after (assoc "Retry-After" (str retry-after)))
   {"type" (str problem kind) "title" "throttled" "status" status
    "detail" (str kind " detail") "retryable" retryable}])

(def ^:dynamic *sleeps* nil)

(defn stub
  "Serve responses in order (the last one repeats); records paths and the peak
  number of requests in flight."
  [responses & {:keys [delay-ms] :or {delay-ms 0}}]
  (let [queue (atom (vec responses))
        rec (atom {:paths [] :inflight 0 :peak 0})
        server (HttpServer/create (InetSocketAddress. "127.0.0.1" 0) 0)
        pool (Executors/newFixedThreadPool 16)]
    (.createContext server "/"
                    (reify HttpHandler
                      (handle [_ ex]
                        (try
                        (let [^HttpExchange ex ex
                              [status headers body] (ffirst (swap-vals! queue #(if (next %) (subvec % 1) %)))
                              ^bytes data (.getBytes (if (string? body) ^String body (json/write-json-str body)) "UTF-8")]
                          (.readAllBytes (.getRequestBody ex))
                          (swap! rec (fn [r] (let [n (inc (:inflight r))]
                                               (-> r (update :paths conj (str (.getRequestURI ex)))
                                                   (assoc :inflight n :peak (max n (:peak r)))))))
                          (Thread/sleep (long delay-ms))
                          (doseq [[k v] headers] (.add (.getResponseHeaders ex) k v))
                          (.sendResponseHeaders ex status (alength data))
                          (with-open [out (.getResponseBody ex)] (.write out data))
                          (swap! rec update :inflight dec))
                        (catch Throwable e (.printStackTrace e))))))
    (.setExecutor server pool)
    (.start server)
    (http/set-endpoint! (str "http://127.0.0.1:" (.getPort (.getAddress server))))
    (http/set-token! "test-token")
    {:server server :pool pool :rec rec}))

(defn- reset-state! []
  (bp/set-retry-policy! {:max-retries 5 :max-backoff 60 :max-concurrency 4})
  (reset! @#'bp/next-allowed-ms 0)
  (reset! @#'bp/max-timeout bp/default-max-timeout-ms))

(use-fixtures :each
  (fn [t]
    (reset-state!)
    (binding [*sleeps* (atom [])]
      (with-redefs-fn {#'bp/sleep! (fn [ms] (swap! *sleeps* conj ms))}
        (fn []
          (try (t)
               (finally (http/set-endpoint! nil) (http/set-token! nil) (reset-state!))))))))

(def query '{:find [?x] :where [[?x :db/ident]]})

(defn run-query []
  (binding [http/*cache* false]
    (http/q query (pdb/->HttpDb "db1" (atom nil)) [])))

(defmacro with-stub [[sym responses & opts] & body]
  `(let [~sym (stub ~responses ~@opts)]
     (try ~@body
          (finally (.stop ^HttpServer (:server ~sym) 0)
                   (.shutdownNow ^java.util.concurrent.ExecutorService (:pool ~sym))))))

(deftest retries-with-retry-after
  (with-stub [s [(throttled :retry-after 3) (throttled :kind "overloaded" :status 503 :retry-after 1) ok]]
    (is (= [[1]] (run-query)))
    (is (= [3000 1000] @*sleeps*))
    (is (= 3 (count (:paths @(:rec s)))))))

(deftest jittered-backoff-without-retry-after
  (with-stub [_ (concat (repeat 3 [429 {"Content-Type" "text/plain"} "slow down"]) [ok])]
    (bp/set-retry-policy! {:max-backoff 3})
    (run-query)
    (is (= 3 (count @*sleeps*)))
    (is (every? true? (map-indexed (fn [i ms] (<= 0 ms (* 1000 (min 3 (Math/pow 2 (inc i)))))) @*sleeps*)))))

(deftest gives-up-after-max-retries
  (with-stub [s [(throttled :retry-after 1)]]
    (bp/set-retry-policy! {:max-retries 2})
    (let [e (try (run-query) nil (catch clojure.lang.ExceptionInfo e e))]
      (is (:patternq/throttled (ex-data e)))
      (is (= (str problem "rate-limited") (:problem-type (ex-data e))))
      (is (re-find #"rate-limited detail" (ex-message e)))
      (is (re-find #"2 retries" (ex-message e)))
      (is (= 3 (count (:paths @(:rec s))))))))

(deftest not-retried
  (doseq [resp [(throttled :retryable false) [502 {} "bad gateway"] [504 {} "timeout"]
                [400 {"Content-Type" "application/json"} {"error" "bad query"}]]]
    (reset! *sleeps* [])
    (with-stub [s [resp ok]]
      (let [e (try (run-query) nil (catch clojure.lang.ExceptionInfo e e))]
        (is (re-find #"Query failed \(HTTP" (ex-message e)))
        (is (= 1 (count (:paths @(:rec s)))))
        (is (empty? @*sleeps*))))))

(deftest query-timeout-advises-narrowing
  (with-stub [_ [[400 {"Content-Type" "application/json"}
                  {"error" "Query canceled: timeout elapsed" "timeout" true "type" (str problem "query-timeout")}]]]
    (is (thrown-with-msg? clojure.lang.ExceptionInfo #"Narrow the query or page it" (run-query)))))

(deftest query-too-broad-is-surfaced-not-retried
  (with-stub [s [[400 {"Content-Type" "application/json"}
                  {"error" "Query too broad: clause [?e ?a ?v] binds neither attribute nor entity; bind the attribute" "type" (str problem "query-too-broad")
                   "reason" "full-scan" "clause" "[?e ?a ?v]"}] ok]]
    (let [e (try (run-query) nil (catch clojure.lang.ExceptionInfo e e))]
      (is (re-find #"bind the attribute" (ex-message e)))
      (is (= (str problem "query-too-broad") (:problem-type (ex-data e))))
      (is (= 1 (count (:paths @(:rec s)))))
      (is (empty? @*sleeps*)))))

(deftest self-throttles-on-r-zero
  (with-stub [_ [[200 {"Content-Type" "application/json" "RateLimit" "\"api\";r=0;t=1"}
                  {"query_result" [] "basis_t" 1}] ok]]
    (run-query)
    (is (empty? @*sleeps*))
    (run-query)
    (is (= 1 (count @*sleeps*)))
    (is (< 0 (first @*sleeps*) 1001))))

(deftest list-datasets-goes-through-back-pressure
  (with-stub [_ [(throttled :retry-after 1) [200 {"Content-Type" "application/json"} {"datasets" []}]]]
    (is (= [] (http/list-datasets)))
    (is (= [1000] @*sleeps*))))

(deftest concurrency-capped-per-process
  (with-stub [s [ok] :delay-ms 150]
    (bp/set-retry-policy! {:max-concurrency 2})
    (run! deref (doall (repeatedly 6 #(future (run-query)))))
    (is (= 6 (count (:paths @(:rec s)))))
    (is (= 2 (:peak @(:rec s))))))

(deftest learns-the-advertised-timeout-cap
  (with-stub [_ [[200 {"Content-Type" "application/json" "PDC-Query-Max-Timeout-Ms" "60000"}
                  {"query_result" [] "basis_t" 1}]]]
    (run-query)
    (is (= 60000 (bp/max-timeout-ms)))))

(deftest policy-defaults
  (is (= {:max-retries 5 :max-backoff 60 :max-concurrency 4} (bp/retry-policy))))
