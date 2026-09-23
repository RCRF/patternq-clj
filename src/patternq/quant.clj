(ns patternq.quant
  "Descriptive statistics and distances for expression analysis, extracted
  from the unify-central variant forensics analysis (quantitative.clj,
  util.clj). Plain Clojure math with primitive type hints; no dependencies.

  Expression inputs are either seqs of numbers, [[gene value] ...] tuples,
  or {gene value} maps as noted per fn."
  (:require [clojure.set :as set]))

(set! *warn-on-reflection* true)

(defn mean
  ^double [xs]
  (let [n (count xs)]
    (if (zero? n)
      Double/NaN
      (/ (double (reduce (fn [^double acc x] (+ acc (double x))) 0.0 xs)) n))))

(defn variance
  "Sample variance (n - 1 denominator); 0.0 for fewer than 2 values."
  ^double [xs]
  (let [n (count xs)]
    (if (< n 2)
      0.0
      (let [m (mean xs)]
        (/ (double (reduce (fn [^double acc x] (let [d (- (double x) m)] (+ acc (* d d)))) 0.0 xs))
           (dec n))))))

(defn sd ^double [xs] (Math/sqrt (variance xs)))

(defn log2 ^double [x] (/ (Math/log (double x)) (Math/log 2.0)))

(defn median
  [xs]
  (when (seq xs)
    (let [v (vec (sort xs)) n (count v) h (quot n 2)]
      (if (odd? n) (double (nth v h)) (/ (+ (double (nth v (dec h))) (double (nth v h))) 2.0)))))

(defn quantile
  "Linear-interpolated quantile p in [0, 1]."
  [xs p]
  (when (seq xs)
    (let [v (vec (sort xs))
          pos (* (double p) (dec (count v)))
          lo (long (Math/floor pos))
          hi (long (Math/ceil pos))
          frac (- pos lo)]
      (+ (double (nth v lo)) (* frac (- (double (nth v hi)) (double (nth v lo))))))))

(defn percentile-rank
  "Percentile (0-100) of `value` within `data-points`, counting ties as half:
  100 * (count below + 0.5 * count equal) / n."
  ^double [data-points value]
  (let [n (count data-points)
        v (double value)
        below (count (filter #(< (double %) v) data-points))
        equal (count (filter #(== (double %) v) data-points))]
    (if (zero? n)
      0.0
      (max 0.0 (* 100.0 (/ (+ below (* 0.5 equal)) n))))))

(defn insertion-index
  "Index at which `target` would be inserted in ascending `sorted-vec`
  (binary search; for an existing value returns one of its indices)."
  ^long [sorted-vec target]
  (let [sv (vec sorted-vec) t (double target)]
    (loop [low 0 high (dec (count sv))]
      (if (> low high)
        low
        (let [mid (quot (+ low high) 2)
              mv (double (nth sv mid))]
          (cond
            (== mv t) mid
            (< mv t) (recur (inc mid) high)
            :else (recur low (dec mid))))))))

(defn z-scores
  "[[gene value] ...] -> [[gene value z] ...] using the population standard
  deviation of the values."
  [data]
  (let [vals (map (comp double second) data)
        n (count vals)
        m (mean vals)
        s (Math/sqrt (/ (double (reduce + (map #(let [d (- (double %) m)] (* d d)) vals))) n))]
    (mapv (fn [[g v]] [g v (if (zero? s) 0.0 (/ (- (double v) m) s))]) data)))

(defn z-score
  "z of `value` relative to reference `xs` (sample sd)."
  ^double [xs value]
  (let [s (sd xs)]
    (if (zero? s) 0.0 (/ (- (double value) (mean xs)) s))))

(defn log-fold-change
  "Log2 fold change treatment / control with pseudocount 1, over the union of
  genes (missing values count as 0). Inputs: [[gene value] ...] or maps.
  Returns [[gene lfc] ...]."
  [control treatment]
  (let [c (into {} control) t (into {} treatment)]
    (vec (for [g (distinct (concat (keys c) (keys t)))]
           [g (log2 (/ (+ 1.0 (double (get t g 0.0))) (+ 1.0 (double (get c g 0.0)))))]))))

(defn ma-values
  "Two-sample MA values over genes present in either sample: [[gene M A] ...]
  with M = log2((t+1)/(c+1)) and A = mean of log2(x+1)."
  [control treatment]
  (let [c (into {} control) t (into {} treatment)]
    (vec (for [g (distinct (concat (keys c) (keys t)))
               :let [lc (log2 (inc (double (get c g 0.0))))
                     lt (log2 (inc (double (get t g 0.0))))]]
           [g (- lt lc) (/ (+ lt lc) 2.0)]))))

(defn l2-distance
  "Euclidean distance between two {gene value} maps over the genes of the
  first (missing genes in the second count as 0)."
  ^double [gx1 gx2]
  (Math/sqrt (double (reduce-kv (fn [^double acc g v]
                                  (let [d (- (double v) (double (get gx2 g 0.0)))] (+ acc (* d d))))
                                0.0 (into {} gx1)))))

(defn cosine-distance
  "1 - cosine similarity between two {gene value} maps (1.0 if either is all
  zero)."
  ^double [gx1 gx2]
  (let [m1 (into {} gx1) m2 (into {} gx2)
        dot (double (reduce-kv (fn [^double acc g v] (+ acc (* (double v) (double (get m2 g 0.0))))) 0.0 m1))
        n1 (double (reduce + (map #(let [x (double %)] (* x x)) (vals m1))))
        n2 (double (reduce + (map #(let [x (double %)] (* x x)) (vals m2))))]
    (if (or (zero? n1) (zero? n2))
      1.0
      (- 1.0 (/ dot (* (Math/sqrt n1) (Math/sqrt n2)))))))

(defn top-varying
  "Top `n` genes by variance of log2(x + 1).
  gene-data: [{:symbol s :counts [numbers]} ...] (or {gene [values]}).
  Returns [{:symbol :mean :variance} ...] sorted by variance desc."
  [gene-data n]
  (->> (if (map? gene-data) (map (fn [[k v]] {:symbol k :counts v}) gene-data) gene-data)
       (map (fn [{:keys [symbol counts]}]
              (let [lc (map #(log2 (inc (double %))) counts)]
                {:symbol symbol :mean (mean lc) :variance (variance lc)})))
       (sort-by :variance >)
       (take n)
       vec))

(defn ssgsea-score
  "Single-sample GSEA enrichment score (Barbie et al. 2009): genes ordered by
  value (descending); the running sum adds hits weighted by rank^alpha, where
  the highest-expressed gene has the highest rank (n), minus misses, and the
  score integrates the path (sum). Same definition as the R and Python
  libraries. (The unify-central report version weighted by the inverse rank.)
  sample-data: [[gene value] ...] or map; gene-set: coll of symbols;
  alpha default 0.25. Returns nil when no gene (or every gene) is in the set."
  ([sample-data gene-set] (ssgsea-score sample-data gene-set 0.25))
  ([sample-data gene-set alpha]
   (let [gs (set gene-set)
         sorted (sort-by second > (seq sample-data))
         n (count sorted)
         ranked (map-indexed (fn [i [g _]] [g (- n (long i))]) sorted)
         hits (filter (comp gs first) ranked)
         n-h (count hits)
         a (double alpha)
         sum-w (double (reduce + (map #(Math/pow (double (second %)) a) hits)))]
     (if (or (zero? n-h) (= n n-h))
       nil
       (let [miss-step (/ 1.0 (- n n-h))]
         (loop [rs (seq ranked) p-hit 0.0 p-miss 0.0 total 0.0]
           (if-not rs
             total
             (let [[g r] (first rs)
                   hit? (contains? gs g)
                   ph (if hit? (+ p-hit (/ (Math/pow (double r) a) sum-w)) p-hit)
                   pm (if hit? p-miss (+ p-miss miss-step))]
               (recur (next rs) ph pm (+ total (- ph pm)))))))))))

(defn ->histogram
  "Bin [[value count] ...] into half-open intervals from ascending `bins`
  edges: [[[start stop] count] ...]."
  [vals+counts bins]
  (let [intervals (map vector (butlast bins) (rest bins))]
    (mapv (fn [[lo hi]]
            [[lo hi] (reduce + (keep (fn [[v n]] (when (and (>= (double v) (double lo)) (< (double v) (double hi))) n))
                                     vals+counts))])
          intervals)))

(defn kaplan-meier
  "Kaplan-Meier estimate. rows: [{:time t :event bool} ...]. Returns
  [{:time :n-risk :n-event :n-censor :surv} ...] at each distinct time."
  [rows]
  (let [rows (sort-by :time (filter #(some? (:time %)) rows))
        times (distinct (map :time rows))]
    (loop [ts times s 1.0 out []]
      (if-let [t (first ts)]
        (let [at-risk (count (filter #(>= (double (:time %)) (double t)) rows))
              at-t (filter #(== (double (:time %)) (double t)) rows)
              ev (count (filter #(true? (boolean (:event %))) at-t))
              ce (- (count at-t) ev)
              s' (* s (- 1.0 (/ (double ev) at-risk)))]
          (recur (rest ts) s' (conj out {:time t :n-risk at-risk :n-event ev :n-censor ce :surv s'})))
        out))))

;; -- clustering (for heatmaps) --

(defn- euclid ^double [a b]
  (Math/sqrt (double (reduce + (map (fn [x y] (let [d (- (double (or x 0.0)) (double (or y 0.0)))] (* d d))) a b)))))

(defn hclust-order
  "Leaf order of average-linkage agglomerative clustering of `rows`
  (vectors of numbers; nil counts as 0) by Euclidean distance. O(n^3); fine
  for heatmap-sized inputs (hundreds of rows)."
  [rows]
  (let [n (count rows)]
    (if (< n 3)
      (vec (range n))
      (let [rows (vec rows)
            dist (into {} (for [i (range n) j (range (inc i) n)] [[i j] (euclid (rows i) (rows j))]))
            d (fn [i j] (dist (if (< i j) [i j] [j i])))]
        ;; clusters: id -> {:members [...] :order [...]}
        (loop [clusters (into {} (map (fn [i] [i {:members [i] :order [i]}])) (range n))]
          (if (= 1 (count clusters))
            (:order (val (first clusters)))
            (let [ids (vec (keys clusters))
                  [a b] (apply min-key
                               (fn [[a b]]
                                 (let [ma (:members (clusters a)) mb (:members (clusters b))]
                                   (/ (double (reduce + (for [i ma j mb] (d i j)))) (* (count ma) (count mb)))))
                               (for [x (range (count ids)) y (range (inc x) (count ids))] [(ids x) (ids y)]))
                  merged {:members (into (:members (clusters a)) (:members (clusters b)))
                          :order (into (:order (clusters a)) (:order (clusters b)))}]
              (recur (-> clusters (dissoc a b) (assoc (min a b) merged))))))))))

(defn scale-rows
  "z-score each row (vector of numbers, nil kept) by its mean and sample sd."
  [rows]
  (mapv (fn [row]
          (let [xs (remove nil? row)
                m (mean xs)
                s (let [s (sd xs)] (if (or (zero? s) (Double/isNaN s)) 1.0 s))]
            (mapv #(when (some? %) (/ (- (double %) m) s)) row)))
        rows))

(defn overlap
  "Genes present in both collections."
  [a b]
  (set/intersection (set a) (set b)))
