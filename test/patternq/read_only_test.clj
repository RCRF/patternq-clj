(ns patternq.read-only-test
  "patternq is strictly read-only on both transports:
  - peer: no code may call (or resolve) Datomic functions that change
    database state, and peer functions are only resolved from a whitelist;
  - HTTP: only the read endpoints of the query service are used."
  (:require [clojure.java.io :as io]
            [clojure.test :refer [deftest is]]))

(defn- sources []
  (for [f (file-seq (io/file "src"))
        :when (.isFile ^java.io.File f)]
    [(.getPath ^java.io.File f) (slurp f)]))

(def forbidden-peer
  #"\b(d|datomic\.api)/(transact|transact-async|create-database|delete-database|rename-database|request-index|gc-storage|sync-excise)\b|\"(transact|transact-async|create-database|delete-database|rename-database|request-index|gc-storage|sync-excise)\"")

(def peer-whitelist #{"q" "connect" "db" "pull" "pull-many" "basis-t"})

(deftest no-state-changing-datomic-calls
  (doseq [[path src] (sources)]
    (is (not (re-find forbidden-peer src))
        (str path " refers to a state-changing Datomic function"))))

(deftest peer-fns-are-whitelisted
  ;; the peer is loaded with requiring-resolve; every resolved datomic.api
  ;; name must be a read-only one
  (doseq [[path src] (sources)
          [_ n] (concat (re-seq #"peer-fn '([a-z-]+)" src)
                        (re-seq #"'datomic\.api/([a-z-]+)" src))]
    (is (contains? peer-whitelist n) (str path " resolves datomic.api/" n)))
  (doseq [[path src] (sources)]
    (is (not (re-find #"\(requiring-resolve \(symbol \"datomic\.api\" (?!\(name sym\))" src))
        (str path " resolves datomic.api dynamically outside peer-fn"))))

(deftest http-read-endpoints-only
  (let [src (slurp (io/file "src/patternq/http.clj"))
        paths (map second (re-seq #"\"(/[a-z0-9/-]+)" src))]
    (is (seq paths))
    (doseq [p paths]
      (is (some #(.startsWith ^String p %) ["/query/" "/matrix/" "/api-v1/list/"])
          (str "http.clj uses endpoint " p)))
    (is (not (re-find #"(?i)transact|\"DELETE\"|\.DELETE|\.PUT" src)))))
