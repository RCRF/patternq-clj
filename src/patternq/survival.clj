(ns patternq.survival
  "Survival helpers for outcome-association analyses (as in the PRINCE
  biomarker figures: OS stratified at the median of a baseline biomarker,
  log-rank p-values, landmark survival status), and change from baseline.
  Mirrors R/patternq/R/survival.R. Hand-rolled, no dependencies; the log-rank
  test matches R's patternq::logrank_test (= survival::survdiff).

  Kaplan-Meier estimation itself is patternq.quant/kaplan-meier."
  (:require [patternq.quant :as quant]))

(set! *warn-on-reflection* true)

;; -- chi-square upper tail via the regularized incomplete gamma function --

(defn- log-gamma ^double [^double x]
  ;; Lanczos approximation (g = 7, n = 9)
  (let [c [0.99999999999980993 676.5203681218851 -1259.1392167224028 771.32342877765313
           -176.61502916214059 12.507343278686905 -0.13857109526572012
           9.9843695780195716e-6 1.5056327351493116e-7]]
    (if (< x 0.5)
      (- (Math/log (/ Math/PI (Math/sin (* Math/PI x)))) (log-gamma (- 1.0 x)))
      (let [x (- x 1.0)
            t (+ x 7.5)
            a (reduce (fn [^double acc i] (+ acc (/ (double (c i)) (+ x (double i))))) (double (c 0)) (range 1 9))]
        (+ (* 0.5 (Math/log (* 2 Math/PI))) (* (+ x 0.5) (Math/log t)) (- t) (Math/log a))))))

(defn- gamma-p-series ^double [^double a ^double x]
  (loop [n 1 ap a sum (/ 1.0 a) del (/ 1.0 a)]
    (let [ap (+ ap 1.0) del (* del (/ x ap)) sum (+ sum del)]
      (if (or (< (Math/abs del) (* (Math/abs sum) 1e-15)) (> n 1000))
        (* sum (Math/exp (- (* a (Math/log x)) x (log-gamma a))))
        (recur (inc n) ap sum del)))))

(defn- gamma-q-cf ^double [^double a ^double x]
  ;; Lentz continued fraction for Q(a, x)
  (let [tiny 1e-300]
    (loop [i 1 b (+ x 1.0 (- a)) c (/ 1.0 tiny) d (/ 1.0 (+ x 1.0 (- a))) h (/ 1.0 (+ x 1.0 (- a)))]
      (let [an (* (- i) (- i a))
            b (+ b 2.0)
            d (let [d (+ (* an d) b)] (/ 1.0 (if (< (Math/abs d) tiny) tiny d)))
            c (let [c (+ b (/ an c))] (if (< (Math/abs c) tiny) tiny c))
            del (* d c)
            h (* h del)]
        (if (or (< (Math/abs (- del 1.0)) 1e-15) (> i 1000))
          (* h (Math/exp (- (* a (Math/log x)) x (log-gamma a))))
          (recur (inc i) b c d h))))))

(defn chisq-upper
  "P(X > x) for X ~ chi-square with `df` degrees of freedom."
  ^double [x df]
  (let [x (double x) a (/ (double df) 2.0) y (/ x 2.0)]
    (cond (<= x 0.0) 1.0
          (< y (+ a 1.0)) (- 1.0 (gamma-p-series a y))
          :else (gamma-q-cf a y))))

(defn- solve
  "Solve A x = b (small dense system, Gaussian elimination with pivoting)."
  [a b]
  (let [n (count b)
        m (mapv (fn [row bi] (conj (mapv double row) (double bi))) a b)
        m (reduce (fn [m col]
                    (let [piv (apply max-key #(Math/abs (double (get-in m [% col]))) (range col n))
                          m (assoc m col (m piv) piv (m col))
                          p (double (get-in m [col col]))]
                      (reduce (fn [m r]
                                (if (= r col) m
                                    (let [f (/ (double (get-in m [r col])) p)]
                                      (assoc m r (mapv (fn [x y] (- (double x) (* f (double y)))) (m r) (m col))))))
                              m (range n))))
                  m (range n))]
    (mapv (fn [i] (/ (double (peek (m i))) (double (get-in m [i i])))) (range n))))

(defn logrank-test
  "Mantel-Haenszel log-rank test for a difference in survival between two or
  more groups. rows: seq of maps; opts :time (default :time), :event
  (default :event; boolean or 0/1), :group (default :group). Rows with a nil
  time, event or group are dropped.
  Returns {:chisq :df :p :observed {g n} :expected {g n} :n {g n}}; :p nil
  with fewer than two groups."
  ([rows] (logrank-test rows {}))
  ([rows {:keys [time event group] :or {time :time event :event group :group}}]
   (let [ev? (fn [e] (if (number? e) (pos? (double e)) (boolean e)))
         rows (->> rows
                   (filter #(and (some? (get % time)) (some? (get % event)) (some? (get % group))))
                   (mapv (fn [r] {:t (double (get r time)) :e (ev? (get r event)) :g (str (get r group))})))
         lv (vec (sort (distinct (map :g rows))))
         k (count lv)]
     (if (< k 2)
       {:chisq nil :df 0 :p nil}
       (let [times (sort (distinct (map :t (filter :e rows))))
             zero-v (vec (repeat k (vec (repeat k 0.0))))
             [O E V] (reduce
                       (fn [[O E V] t]
                         (let [at-risk (mapv (fn [g] (count (filter #(and (>= (double (:t %)) (double t)) (= g (:g %))) rows))) lv)
                               d-g (mapv (fn [g] (count (filter #(and (== (double (:t %)) (double t)) (:e %) (= g (:g %))) rows))) lv)
                               n (double (reduce + at-risk))
                               d (double (reduce + d-g))]
                           (if (< n 1.0)
                             [O E V]
                             [(mapv + O d-g)
                              (mapv (fn [e r] (+ (double e) (/ (* d (double r)) n))) E at-risk)
                              (if (> n 1.0)
                                (let [f (/ (* d (- n d)) (* n n (- n 1.0)))]
                                  (vec (for [i (range k)]
                                         (vec (for [j (range k)]
                                                (+ (double (get-in V [i j]))
                                                   (* f (- (if (= i j) (* (double (at-risk i)) n) 0.0)
                                                           (* (double (at-risk i)) (double (at-risk j)))))))))))
                                V)])))
                       [(vec (repeat k 0)) (vec (repeat k 0.0)) zero-v]
                       times)
             idx (range (dec k))
             diff (mapv #(- (double (O %)) (double (E %))) idx)
             vs (mapv (fn [i] (mapv (fn [j] (get-in V [i j])) idx)) idx)
             x (solve vs diff)
             chisq (reduce + (map * diff x))]
         {:chisq chisq :df (dec k) :p (chisq-upper chisq (dec k))
          :observed (zipmap lv O) :expected (zipmap lv E)
          :n (zipmap lv (map (fn [g] (count (filter #(= g (:g %)) rows))) lv))})))))

(defn median-split
  "Label each value \"low\" (below the median) or \"high\" (at or above);
  nil stays nil. opts: :labels [low high]."
  ([xs] (median-split xs {}))
  ([xs {:keys [labels] :or {labels ["low" "high"]}}]
   (let [m (quant/median (remove nil? xs))]
     (mapv #(when (some? %) (if (>= (double %) (double m)) (second labels) (first labels))) xs))))

(defn survival-status
  "Landmark survival status of an outcome row: \"alive at 1 year\" when time
  >= :at, \"died within 1 year\" with an event before :at, nil when censored
  before :at (unknown). opts: :time (:os), :event (:os-event), :at (12,
  months for PRINCE OS), :labels {:alive :died}."
  ([outcome] (survival-status outcome {}))
  ([outcome {:keys [time event at labels]
             :or {time :os event :os-event at 12
                  labels {:alive "alive at 1 year" :died "died within 1 year"}}}]
   (let [t (get outcome time) e (get outcome event)]
     (cond (nil? t) nil
           (>= (double t) (double at)) (:alive labels)
           (and (some? e) (if (number? e) (pos? (double e)) (boolean e))) (:died labels)
           :else nil))))

(defn survival-by-median
  "Split subjects at the median of a subject-level biomarker and test the
  difference in survival (log-rank).

  values: {subject-id value}; outcomes: patternq.clinical/subject-outcomes.
  Returns the outcome rows of subjects with a value, plus :value and :group
  (\"low\"/\"high\"); the log-rank result is in the metadata under
  :patternq/logrank. opts: :time (:os), :event (:os-event)."
  ([values outcomes] (survival-by-median values outcomes {}))
  ([values outcomes {:keys [time event] :or {time :os event :os-event}}]
   (let [values (into {} (remove (comp nil? val)) values)
         rows (filterv #(contains? values (:subject-id %)) outcomes)
         groups (median-split (map #(values (:subject-id %)) rows))
         rows (mapv (fn [r g] (assoc r :value (values (:subject-id r)) :group g)) rows groups)]
     (with-meta rows {:patternq/logrank (logrank-test rows {:time time :event event :group :group})}))))

(defn logrank
  "The log-rank result attached by survival-by-median."
  [rows]
  (:patternq/logrank (meta rows)))

(def ^:private default-series-keys
  [:cell-population :epitope-id :hgnc-symbol :signature :measurement-set :target])

(defn change-from-baseline
  "For each subject (and target, e.g. cell population or protein), the change
  of the value relative to the subject's value at the baseline timepoint.

  rows: maps with :subject-id, :timepoint-id and :value (e.g. from
  patternq.context/add-sample-context). opts:
    :baseline    baseline timepoint id (default \"C1D1\")
    :by          keys identifying a series besides :subject-id (default: the
                 target keys present, e.g. :cell-population)
    :method      :log2-ratio (log2(value / baseline), for frequencies;
                 default), :difference (value - baseline, for values already on
                 a log scale such as Olink NPX) or :ratio
    :pseudocount added to both before a ratio (default 0)
  Returns rows with a baseline value, plus :baseline-value and :change (nil
  when not finite)."
  ([rows] (change-from-baseline rows {}))
  ([rows {:keys [baseline by method pseudocount] :or {baseline "C1D1" method :log2-ratio pseudocount 0}}]
   (let [by (or by (filterv (fn [k] (some #(contains? % k) rows)) default-series-keys))
         key-fn (fn [r] (into [(:subject-id r)] (map #(get r %) by)))
         base (->> rows
                   (filter #(and (= baseline (:timepoint-id %)) (some? (:value %))))
                   (group-by key-fn)
                   (into {} (map (fn [[k rs]] [k (quant/mean (map :value rs))]))))
         pc (double pseudocount)
         change (fn [v b]
                  (let [v (double v) b (double b)
                        c (case method
                            :log2-ratio (quant/log2 (/ (+ v pc) (+ b pc)))
                            :difference (- v b)
                            :ratio (/ (+ v pc) (+ b pc)))]
                    (when (and (not (Double/isNaN c)) (not (Double/isInfinite c))) c)))]
     (vec (for [r rows
                :let [b (base (key-fn r))]
                :when (some? b)]
            (assoc r :baseline-value b :change (when (some? (:value r)) (change (:value r) b))))))))
