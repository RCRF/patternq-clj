(ns patternq.results-plot-test
  (:require [clojure.test :refer [deftest is]]
            [patternq.plot :as plot]
            [patternq.results :as res]))

(deftest flatten-pull
  (is (= {:sample-id "S1" :subject-id "P1" :sample-specimen "ffpe"}
         (res/flatten-pull {:sample/id "S1" :db/id 1
                            :sample/subject {:subject/id "P1"}
                            :sample/specimen {:db/ident :sample.specimen/ffpe}})))
  (is (= {:cnv-genes ["A" "B"]} (res/flatten-pull {:cnv/genes [{:gene/hgnc-symbol "A"} {:gene/hgnc-symbol "B"}]})))
  (is (= {:subject-race ["white"]} (res/flatten-pull {:subject/race [{:db/ident :subject.race/white}]}))))

(deftest matrices
  (let [m (res/->matrix [{:sample-id "s1" :g "A" :value 1} {:sample-id "s1" :g "A" :value 3}
                         {:sample-id "s2" :g "B" :value 5}] :g)]
    (is (= ["s1" "s2"] (:row-names m)))
    (is (= [[2.0 nil] [nil 5.0]] (:values m)))
    (is (= 2 (count (res/->long m))))
    (is (= ["B"] (:col-names (res/select-targets m {:exclude ["A"]}))))))

(deftest plot-specs
  (let [vs [{:sample-id "a" :vaf 0.5 :hgnc-symbol "TP53" :impact "high"}
            {:sample-id "b" :vaf 0.2 :hgnc-symbol "TP53"}]]
    (is (= 2 (count (:data (plot/vaf-histogram vs)))))
    (is (= "#2a78d6" (get-in (plot/vaf-histogram vs) [:data 0 :marker :color])))
    (is (= "heatmap" (get-in (plot/mutation-landscape vs) [:data 0 :type]))))
  (let [spec (plot/survival [{:os 1 :os-event true :bor "PR"} {:os 2 :os-event false :bor "PD"}] {:group :bor})]
    (is (= "hv" (get-in spec [:data 0 :line :shape]))))
  (is (= 8 (count (distinct (map (plot/fold-other (map str (range 20))) (map str (range 20))))))))
