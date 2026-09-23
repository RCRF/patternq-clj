(ns patternq.transport-test
  "Transport selection, the EDN -> JSON query conversion, and (when the peer
  is on the classpath and PATTERNQ_DATOMIC_URI is set) live parity between
  the HTTP and peer transports."
  (:require [clojure.test :refer [deftest is testing]]
            [patternq.clinical :as clinical]
            [patternq.dataset :as pqd]
            [patternq.db :as pdb]
            [patternq.http :as http]
            [patternq.reference :as ref]
            [patternq.variants :as pv]))

(deftest query->wire-conversion
  (let [body (http/query->wire '{:find [?id (pull ?s [* {:sample/subject [:subject/id]}]) (count ?x)]
                                 :in [$ ?imp [?g ...]]
                                 :where [[?s :sample/id ?id]
                                         [(> ?vaf 0.3)]
                                         (not [?v :variant/impact ?imp])
                                         [_ :x/y ?x]]}
                               [:variant.impact/low ["TP53" "KRAS"]])]
    (testing "keyword args are inlined as literals; other args sent"
      (is (= ["$" ["?g" "..."]] (get-in body ["query" ":in"])))
      (is (= [["TP53" "KRAS"]] (get body "args")))
      (is (= ["not" ["?v" ":variant/impact" ":variant.impact/low"]] (get-in body ["query" ":where" 2]))))
    (testing "symbols, keywords, pull maps, expressions"
      (is (= ["?id" ["pull" "?s" ["*" {":sample/subject" [":subject/id"]}]] ["count" "?x"]]
             (get-in body ["query" ":find"])))
      (is (= [[">" "?vaf" 0.3]] (get-in body ["query" ":where" 1])))
      (is (= ["_" ":x/y" "?x"] (get-in body ["query" ":where" 3])))))
  (testing "find specs are sent as relations"
    (is (= ["?x"] (get-in (http/query->wire '{:find [[?x ...]] :where [[_ :a/b ?x]]} []) ["query" ":find"])))
    (is (= ["?x"] (get-in (http/query->wire '{:find [?x .] :where [[_ :a/b ?x]]} []) ["query" ":find"])))
    (is (= "$" (first (get-in (http/query->wire '{:find [?x] :in [?y] :where [[?y :a/b ?x]]} ["v"])
                              ["query" ":in"]))))))

(deftest transport-selection
  (let [orig (pdb/transport)]
    (try
      (pdb/set-transport! :http)
      (is (= :http (pdb/transport)))
      (is (pdb/http-db? (pdb/db "any-db")))
      (is (thrown? clojure.lang.ExceptionInfo (pdb/db "admin-db-1")))
      (is (thrown? clojure.lang.ExceptionInfo (pdb/set-transport! :bogus)))
      (pdb/set-transport! :peer)
      (when-not (pdb/peer-available?)
        (is (thrown-with-msg? clojure.lang.ExceptionInfo #"not on the classpath" (pdb/db "any-db"))))
      (when (pdb/peer-available?)
        (let [uri (pdb/base-uri)]
          (try
            (pdb/set-base-uri! nil)
            (when-not (System/getenv "PATTERNQ_DATOMIC_URI")
              (is (thrown-with-msg? clojure.lang.ExceptionInfo #"PATTERNQ_DATOMIC_URI" (pdb/db "any-db"))))
            (finally (pdb/set-base-uri! uri)))))
      (finally (pdb/set-transport! nil)
               (is (= orig (pdb/transport)))))))

(defn- both
  "Run f under both transports; returns [http-result peer-result]."
  [f]
  (try
    [(do (pdb/set-transport! :http) (f)) (do (pdb/set-transport! :peer) (f))]
    (finally (pdb/set-transport! nil))))

(defn- norm
  "Order-insensitive comparison. :db.type/float values are Floats from the
  peer and Doubles (the float's shortest decimal) over HTTP, so floating
  values are compared at float precision."
  [rows]
  (set (map (fn [r] (update-vals r #(if (float? %) (float %) %))) rows)))

(deftest ^:live http-peer-parity
  (if-not (and (pdb/peer-available?) (pdb/base-uri))
    (println "skipping http-peer-parity: peer not on classpath or PATTERNQ_DATOMIC_URI unset")
    (let [h37001 (http/resolve-db "H37001")
          prince (http/resolve-db "prince-2022")
          same (fn [label f] (let [[a b] (both f)] (is (= (norm a) (norm b)) label)))]
      (same "samples" #(pqd/samples h37001))
      (same "subjects" #(pqd/subjects prince))
      (same "variants" #(pqd/variants h37001))
      (same "gene expression" #(pqd/gene-expression h37001 {:genes ["BAP1" "GNAQ"] :measurement :tpm}))
      (same "measurement types" #(pqd/measurement-types prince "PICI CyTOF Immune Profiling"))
      (same "measurements" #(pqd/measurements prince :percent-of-parent {:measurement-set "PICI CyTOF Immune Profiling"}))
      (same "outcomes" #(clinical/subject-outcomes prince))
      (same "timepoints" #(pqd/timepoints prince))
      (same "variant annotations" #(ref/variant-annotations prince {:genes ["KRAS"]}))
      (same "cnv segments" #(pqd/cnv-segments h37001 {:genes ["BAP1"]}))
      (same "participant variants" #(pv/participant-variants h37001 "H37001" {:exclude-genes #{"OR8U1"}}))
      (let [[a b] (both #(:basis-t (:patternq/provenance (meta (pqd/samples h37001)))))]
        (is (= a b) "basis t")))))
