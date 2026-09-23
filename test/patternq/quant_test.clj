(ns patternq.quant-test
  (:require [clojure.test :refer [deftest is testing]]
            [patternq.quant :as q]))

(defn- approx= [a b] (< (Math/abs (- (double a) (double b))) 1e-9))

(deftest descriptive
  (is (approx= 2.5 (q/mean [1 2 3 4])))
  (is (approx= (/ 5.0 3.0) (q/variance [1 2 3 4])))
  (is (approx= 0.0 (q/variance [7])))
  (is (approx= 3.0 (q/log2 8)))
  (is (approx= 2.5 (q/median [4 1 3 2])))
  (is (approx= 2.5 (q/quantile [1 2 3 4] 0.5))))

(deftest percentile-and-rank
  (is (approx= 50.0 (q/percentile-rank [1 2 3 4] 2.5)))
  (is (approx= 37.5 (q/percentile-rank [1 2 3 4] 2)))   ; (1 + 0.5) / 4
  (is (approx= 0.0 (q/percentile-rank [] 1)))
  (is (= 2 (q/insertion-index [1 3 5 7] 4)))
  (is (= 0 (q/insertion-index [1 3 5 7] 0)))
  (is (= 4 (q/insertion-index [1 3 5 7] 9))))

(deftest z-scores-and-fold-change
  (let [zs (q/z-scores [["A" 1] ["B" 2] ["C" 3]])]
    (is (approx= 0.0 (nth (second zs) 2)))
    (is (approx= (/ 1.0 (Math/sqrt (/ 2.0 3.0))) (nth (last zs) 2))))
  (is (approx= 1.0 (q/z-score [1 2 3] 3)))
  (is (= {"A" 1.0} (update-vals (into {} (q/log-fold-change [["A" 1]] [["A" 3]])) double)))
  (let [[[g m a]] (q/ma-values {"A" 1} {"A" 3})]
    (is (= "A" g)) (is (approx= 1.0 m)) (is (approx= 1.5 a))))

(deftest distances
  (is (approx= 5.0 (q/l2-distance {"a" 3 "b" 4} {})))
  (is (approx= 0.0 (q/cosine-distance {"a" 1 "b" 1} {"a" 2 "b" 2})))
  (is (approx= 1.0 (q/cosine-distance {"a" 1} {"b" 1})))
  (is (approx= 1.0 (q/cosine-distance {"a" 0} {"a" 1}))))

(deftest top-varying-and-ssgsea
  (is (= ["B" "A"] (map :symbol (q/top-varying {"A" [1 1 1] "B" [0 3 15]} 2))))
  (testing "ssGSEA: set at the top of the ranking scores higher than at the bottom"
    (let [sample [["G1" 100] ["G2" 90] ["G3" 80] ["G4" 10] ["G5" 5] ["G6" 1]]]
      (is (> (q/ssgsea-score sample #{"G1" "G2"}) (q/ssgsea-score sample #{"G5" "G6"})))
      (is (nil? (q/ssgsea-score sample #{"NOPE"}))))))

(deftest ssgsea-matches-r-and-python
  ;; reference values computed with the R library's ssgsea_score
  (let [x (map (fn [i] [(str "G" i) (double (- 101 i))]) (range 1 101))]
    (is (< (Math/abs (- 50.0216062883 (q/ssgsea-score x (map #(str "G" %) (range 1 11))))) 1e-8))
    (is (< (Math/abs (- 17.3034501639 (q/ssgsea-score x ["G3" "G17" "G50" "G88"]))) 1e-8))))

(deftest histogram-km-clustering
  (is (= [[[0 1] 3] [[1 2] 1]] (q/->histogram [[0.5 2] [0.9 1] [1.5 1] [5 9]] [0 1 2])))
  (let [km (q/kaplan-meier [{:time 1 :event true} {:time 2 :event false} {:time 3 :event true} {:time 4 :event true}])]
    (is (= [1 2 3 4] (map :time km)))
    (is (approx= 0.75 (:surv (first km))))
    (is (approx= 0.375 (:surv (nth km 2))))
    (is (approx= 0.0 (:surv (last km)))))
  (let [o (q/hclust-order [[0 0] [0 1] [10 10] [11 10]])]
    (is (= #{0 1 2 3} (set o)))
    (is (contains? #{#{0 1} #{2 3}} (set (take 2 o)))))
  (is (= [[-1.0 0.0 1.0]] (q/scale-rows [[1 2 3]]))))
