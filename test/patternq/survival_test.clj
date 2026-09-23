(ns patternq.survival-test
  (:require [clojure.test :refer [deftest is testing]]
            [patternq.plot :as plot]
            [patternq.survival :as s]))

(defn- approx= [a b tol] (< (Math/abs (- (double a) (double b))) tol))

(defn- rows [times events groups]
  (mapv (fn [t e g] {:time t :event e :group g}) times events groups))

(deftest logrank-matches-r
  (testing "shared cross-language case (R/tests/testthat/test-survival.R)"
    (let [lr (s/logrank-test (rows (range 1 11) [true true false true true true false true false true]
                                   ["a" "b" "a" "a" "b" "a" "b" "b" "a" "b"]))]
      (is (approx= 0.2201257417 (:chisq lr) 1e-9))
      (is (= 1 (:df lr)))))
  (testing "3 groups, ties; reference from R patternq::logrank_test = survival::survdiff"
    (let [lr (s/logrank-test (rows [1 2 3 4 5 6 7 8 9 10 11 12 3 5 7]
                                   [1 0 1 1 0 1 1 1 0 1 1 0 1 1 0]
                                   ["a" "a" "a" "a" "a" "b" "b" "b" "b" "b" "c" "c" "c" "c" "c"]))]
      (is (approx= 5.04930334738172 (:chisq lr) 1e-9))
      (is (approx= 0.0800862040606797 (:p lr) 1e-9))
      (is (= 2 (:df lr)))
      (is (= {"a" 3 "b" 4 "c" 3} (:observed lr)))
      (is (approx= 4.186871 (get-in lr [:expected "b"]) 1e-6))))
  (testing "2 groups"
    (let [lr (s/logrank-test (rows (range 1 11) [1 0 1 1 0 1 1 1 0 1] (concat (repeat 5 "a") (repeat 5 "b"))))]
      (is (approx= 4.9138490041685978 (:chisq lr) 1e-9))
      (is (approx= 0.0266422089217833 (:p lr) 1e-9))))
  (is (nil? (:p (s/logrank-test (rows [1 2] [1 1] ["a" "a"]))))))

(deftest chisq-tail
  (is (approx= 0.05 (s/chisq-upper 3.841458820694124 1) 1e-10))
  (is (approx= 0.05 (s/chisq-upper 5.991464547107979 2) 1e-10))
  (is (approx= 1.0 (s/chisq-upper 0 1) 1e-12)))

(deftest helpers
  (is (= ["low" "high" "high" nil "low"] (s/median-split [1 3 5 nil 2])))
  (is (= ["alive at 1 year" "died within 1 year" nil nil]
         (mapv s/survival-status [{:os 20 :os-event false} {:os 5 :os-event true}
                                  {:os 5 :os-event false} {:os nil :os-event true}])))
  (let [km (s/survival-by-median {"a" 1 "b" 2 "c" 3 "d" 4 "e" nil}
                                 [{:subject-id "a" :os 1 :os-event true} {:subject-id "b" :os 2 :os-event true}
                                  {:subject-id "c" :os 5 :os-event false} {:subject-id "d" :os 6 :os-event true}
                                  {:subject-id "e" :os 3 :os-event true}])]
    (is (= ["low" "low" "high" "high"] (map :group km)))
    (is (number? (:p (s/logrank km)))))
  (let [ch (s/change-from-baseline [{:subject-id "p" :timepoint-id "C1D1" :cell-population "x" :value 2.0}
                                    {:subject-id "p" :timepoint-id "C2D1" :cell-population "x" :value 8.0}
                                    {:subject-id "q" :timepoint-id "C2D1" :cell-population "x" :value 1.0}])]
    (is (= [0.0 2.0] (mapv :change ch))))
  (is (= [0.0 6.0] (mapv :change (s/change-from-baseline [{:subject-id "p" :timepoint-id "C1D1" :value 2.0}
                                                          {:subject-id "p" :timepoint-id "C2D1" :value 8.0}]
                                                         {:method :difference}))))
  (is (= "0.012" (plot/format-p 0.0123))))

(deftest plot-updates
  (let [oc [{:os 1 :os-event true :group "low"} {:os 2 :os-event false :group "high"}
            {:os 3 :os-event true :group "low"} {:os 4 :os-event true :group "high"}]
        spec (plot/survival oc {:group :group :levels ["low" "high"]})]
    (is (= ["low (n=2)" "high (n=2)"] (keep :name (filter :showlegend (map #(assoc % :showlegend (not (false? (:showlegend %)))) (:data spec))))))
    (is (re-find #"log-rank p" (get-in spec [:layout :annotations 0 :text]))))
  (let [spec (plot/heatmap {:row-names ["m1" "m2" "m3"] :col-names ["a" "b" "c"]
                            :values [[1 2 3] [3 2 1] [1 1 2]]}
                           {:scale :row :col-groups {"a" "alive" "b" "died" "c" "alive"}})]
    (is (= 2 (count (:data spec))))
    (is (= "y2" (get-in spec [:data 1 :yaxis])))
    (is (= 2 (count (get-in spec [:layout :annotations])))))
  (let [spec (plot/by-timepoint [{:subject-id "p" :timepoint-id "C1D1" :change 0.0 :g "x"}
                                 {:subject-id "p" :timepoint-id "C2D1" :change 1.0 :g "x"}
                                 {:subject-id "p" :timepoint-id "EOS" :change 1.0 :g "x"}]
                                {:value :change :group :g :lines? true :timepoints ["C1D1" "C2D1"]})]
    (is (= ["C1D1" "C2D1"] (get-in spec [:layout :xaxis :categoryarray])))
    (is (= ["scatter" "box"] (map :type (:data spec))))))
