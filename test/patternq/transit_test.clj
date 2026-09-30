(ns patternq.transit-test
  "Optional transit response formats (patternq.http/*format*): the same
  results as JSON (live, needs PATTERNQ_API_KEY)."
  (:require [clojure.test :refer [deftest is testing]]
            [patternq.dataset :as pqd]
            [patternq.http :as http]))

(defn- close? [a b]
  (if (and (number? a) (number? b))
    (<= (Math/abs (- (double a) (double b))) (* 1e-6 (max 1.0 (Math/abs (double a)))))
    (= a b)))

(defn- same-rows? [xs ys]
  (let [k #(dissoc % :vaf)
        xs (sort-by pr-str (map k xs)) ys (sort-by pr-str (map k ys))]
    (= xs ys)))

(deftest ^:live transit-formats-match-json
  (when (System/getenv "PATTERNQ_API_KEY")
    (let [db (http/resolve-db "tcga-uvm")]
      (doseq [f [pqd/samples pqd/variants]]
        (let [base (binding [http/*cache* false] (vec (f db)))]
          (doseq [fmt [:transit+json :transit+msgpack]]
            (testing (str f " " fmt)
              (let [r (binding [http/*format* fmt] (vec (f db)))]
                (is (= (count base) (count r)))
                (is (same-rows? base r))
                (let [vafs (fn [rows] (sort (keep :vaf rows)))]
                  (is (every? true? (map close? (vafs base) (vafs r)))))))))))))

(deftest unknown-format
  (is (thrown-with-msg? clojure.lang.ExceptionInfo #"\*format\* must be"
                        (binding [http/*format* :edn]
                          (http/q '{:find [?x] :where [[_ :a/b ?x]]} {:db-name "x" :basis-t (atom nil)} [])))))
