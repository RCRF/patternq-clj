(ns patternq.expression
  "Gene expression helpers for participant vs. reference cohort analysis,
  extracted from the unify-central variant forensics analysis (util.clj,
  quantitative.clj, unify-central.db.plot-queries). Cohorts are other
  dataset databases (e.g. TCGA) passed explicitly by name or value; nothing
  here reads the admin database."
  (:require [patternq.dataset :as pqd]
            [patternq.db :as pdb]
            [patternq.plot :as plot]
            [patternq.quant :as quant]))

(set! *warn-on-reflection* true)

(def default-measurement :measurement/rsem-normalized-count)

(defn- mattr [m]
  (cond (nil? m) default-measurement
        (and (keyword? m) (namespace m)) m
        :else (keyword "measurement" (name m))))

(defn geneset-expression
  "[[hgnc value] ...] for `genes` in one sample (default measurement
  rsem-normalized-count; opts :measurement e.g. :tpm)."
  ([db-or-name genes sample-id] (geneset-expression db-or-name genes sample-id {}))
  ([db-or-name genes sample-id {:keys [measurement]}]
   (vec (pdb/q {:find '[?hgnc ?v]
              :in '[$ ?sample-id [?hgnc ...]]
              :where ['[?g :gene/hgnc-symbol ?hgnc]
                      '[?gp :gene-product/gene ?g]
                      '[?m :measurement/gene-product ?gp]
                      ['?m (mattr measurement) '?v]
                      '[?m :measurement/sample ?s]
                      '[?s :sample/id ?sample-id]]}
             (pdb/as-db db-or-name) sample-id (vec genes)))))

(defn sample->gene-expression
  "All [[hgnc value] ...] for one sample (default rsem-normalized-count)."
  ([db-or-name sample-id] (sample->gene-expression db-or-name sample-id {}))
  ([db-or-name sample-id {:keys [measurement]}]
   (vec (pdb/q {:find '[?hgnc ?v]
              :in '[$ ?sid]
              :where ['[?s :sample/id ?sid]
                      '[?m :measurement/sample ?s]
                      ['?m (mattr measurement) '?v]
                      '[?m :measurement/gene-product ?gp]
                      '[?gp :gene-product/gene ?g]
                      '[?g :gene/hgnc-symbol ?hgnc]]}
             (pdb/as-db db-or-name) sample-id))))

(defn cohort-gene-expression
  "All values of one gene's expression in a (cohort) database."
  ([db-or-name gene] (cohort-gene-expression db-or-name gene {}))
  ([db-or-name gene {:keys [measurement]}]
   (pdb/q {:find '[[?v ...]]
         :with '[?m]
         :in '[$ ?gene]
         :where ['[?g :gene/hgnc-symbol ?gene]
                 '[?gp :gene-product/gene ?g]
                 '[?m :measurement/gene-product ?gp]
                 ['?m (mattr measurement) '?v]]}
        (pdb/as-db db-or-name) gene)))

(defn cohort-batch-expression
  "{gene [values...]} for a batch of genes in a (cohort) database."
  ([db-or-name genes] (cohort-batch-expression db-or-name genes {}))
  ([db-or-name genes {:keys [measurement]}]
   (->> (pdb/q {:find '[?gene ?v ?m]
              :in '[$ [?gene ...]]
              :where ['[?g :gene/hgnc-symbol ?gene]
                      '[?gp :gene-product/gene ?g]
                      '[?m :measurement/gene-product ?gp]
                      ['?m (mattr measurement) '?v]]}
             (pdb/as-db db-or-name) (vec genes))
        (group-by first)
        (into {} (map (fn [[g rows]] [g (mapv second rows)]))))))

(def genex-attr-candidates
  "Gene-expression value attributes in preference order."
  [:measurement/tpm :measurement/fpkm :measurement/rsem-normalized-count
   :measurement/rsem-raw-count :measurement/read-count])

(defn available-genex-attrs
  "Which genex-attr-candidates a database's RNA-seq measurement sets carry,
  in preference order (samples a few measurements per set)."
  [db-or-name]
  (let [db (pdb/as-db db-or-name)
        sets (pdb/q '[:find [?ms ...]
                    :where [?a :assay/technology :assay.technology/RNA-seq] [?a :assay/measurement-sets ?ms]]
                  db)
        present (into #{} (mapcat (fn [ms]
                                    (->> (pdb/q '[:find (sample 25 ?m) . :in $ ?ms
                                                  :where [?ms :measurement-set/measurements ?m]]
                                                db ms)
                                         (pdb/pull-many db '[*]) (mapcat keys))))
                      sets)]
    (filterv present genex-attr-candidates)))

(defn resolve-genex-attrs
  "Pick the expression attribute to compare a participant db with a cohort
  db: a shared attribute when both have one, else each side's best
  (:shared? false: values are not directly comparable)."
  [participant-db cohort-db]
  (let [p (available-genex-attrs participant-db)
        c (available-genex-attrs cohort-db)
        shared (first (filter (set c) p))]
    {:participant-available p :cohort-available c
     :participant-attr (or shared (first p)) :cohort-attr (or shared (first c))
     :shared? (boolean shared)}))

(defn pathway-gene-percentiles
  "Per gene of `genes`: {:hgnc-symbol :percentile :value :cohort-mean} of a
  participant sample within a cohort (genes missing on either side are
  skipped)."
  ([db cohort-db sample-id genes] (pathway-gene-percentiles db cohort-db sample-id genes {}))
  ([db cohort-db sample-id genes opts]
   (let [part (into {} (geneset-expression db genes sample-id opts))
         cohort (cohort-batch-expression cohort-db genes opts)]
     (vec (for [g genes
                :let [c (cohort g) v (part g)]
                :when (and (seq c) (some? v))]
            {:hgnc-symbol g :percentile (quant/percentile-rank c v) :value v :cohort-mean (quant/mean c)})))))

(defn rank-expression
  "Rank of each of a sample's genes within the cohort's sorted values for the
  gene (insertion index), in batches of 100 genes: [{:hgnc-symbol :rank
  :cohort-n :value} ...]."
  ([db cohort-db sample-id] (rank-expression db cohort-db sample-id {}))
  ([db cohort-db sample-id opts]
   (let [part (into {} (sample->gene-expression db sample-id opts))]
     (vec (mapcat (fn [batch]
                    (for [[g vals] (cohort-batch-expression cohort-db batch opts)
                          :let [v (part g) sv (vec (sort vals))]
                          :when (some? v)]
                      {:hgnc-symbol g :rank (quant/insertion-index sv v) :cohort-n (count sv) :value v}))
                  (partition-all 100 (keys part)))))))

(defn cohort-zscores
  "z-score of each gene's sample value against the cohort distribution for
  that gene (mean / sample sd of log2(x + 1) when :log? default true):
  [{:hgnc-symbol :value :z :percentile :cohort-mean :cohort-n} ...] sorted by z desc.
  `genes` nil means all of the sample's genes (batched).

  The cohort sd is floored at :min-sd (default 0.25, in log2 units when
  :log?) so genes that are (near) constant in the cohort, e.g. almost always
  0, don't produce enormous z-scores from tiny denominators."
  ([db cohort-db sample-id genes] (cohort-zscores db cohort-db sample-id genes {}))
  ([db cohort-db sample-id genes {:keys [log? min-sd] :or {log? true min-sd 0.25} :as opts}]
   (let [part (into {} (if genes (geneset-expression db genes sample-id opts)
                           (sample->gene-expression db sample-id opts)))
         tf (if log? #(quant/log2 (inc (double %))) double)]
     (->> (partition-all 100 (keys part))
          (mapcat (fn [batch]
                    (for [[g vals] (cohort-batch-expression cohort-db batch opts)
                          :let [v (part g)]
                          :when (and (some? v) (> (count vals) 2))]
                      {:hgnc-symbol g :value v
                       :z (let [xs (map tf vals)]
                            (/ (- (double (tf v)) (quant/mean xs)) (max (double min-sd) (quant/sd xs))))
                       :percentile (quant/percentile-rank vals v)
                       :cohort-mean (quant/mean vals) :cohort-n (count vals)})))
          (sort-by (comp - :z))
          vec))))

(defn nearest-neighbors
  "Cosine distance from a target {gene value} map to every sample of a
  cohort database with expression measurement `:measurement` (default :tpm):
  [{:sample-id :distance} ...] ascending. Methodological caveat (from the
  source analysis): raw RNA-seq distance is sensitive to batch, vendor and
  pipeline effects."
  ([cohort-db target-gx] (nearest-neighbors cohort-db target-gx {}))
  ([cohort-db target-gx {:keys [measurement] :or {measurement :tpm}}]
   (let [db (pdb/as-db cohort-db)
         attr (mattr measurement)
         rows (pdb/q {:find '[?sid ?hgnc ?v]
                    :where [['?m attr '?v]
                            '[?m :measurement/sample ?s]
                            '[?s :sample/id ?sid]
                            '[?m :measurement/gene-product ?gp]
                            '[?gp :gene-product/gene ?g]
                            '[?g :gene/hgnc-symbol ?hgnc]]}
                   db)]
     (->> (group-by first rows)
          (map (fn [[sid rs]] {:sample-id sid
                               :distance (quant/cosine-distance target-gx (into {} (map (fn [[_ g v]] [g v])) rs))}))
          (sort-by :distance)
          vec))))

(defn examine-geneset
  "Participant samples vs one or more reference cohorts for a gene set, as a
  plotly spec (patternq.plot/vs-cohort; R: examine_geneset).

  cohort-dbs: {label db-name} or a seq of db names; genes: HGNC symbols or a
  gene set name (patternq.genesets); sample-ids: samples in `db`.
  opts: :measurement (:tpm), :cohort-measurement, :type (:violin | :box),
  :title. With :style (:box | :horizontal | :violin) the older
  patternq.plot/genex-vs-cohort rendering is used instead."
  [db cohort-dbs genes sample-ids {:keys [measurement cohort-measurement type title style]
                                   :or {measurement :tpm type :violin} :as opts}]
  (let [[genes title] (if (string? genes)
                        [((requiring-resolve 'patternq.genesets/geneset) genes) (or title genes)]
                        [genes title])
        cohorts (if (map? cohort-dbs) cohort-dbs (into (array-map) (map (fn [d] [(str d) d]) cohort-dbs)))]
    (if style
      (let [cohort-expr (vec (for [[_ c] cohorts]
                               (let [e (cohort-batch-expression c genes opts)]
                                 (mapv (fn [g] [g (get e g [])]) genes))))
            sample-expr (mapv (fn [sid] [sid (geneset-expression db genes sid opts)]) sample-ids)]
        (plot/genex-vs-cohort sample-expr cohort-expr
                              (merge {:cohort-names (vec (keys cohorts))
                                      :value-label (name (mattr measurement))
                                      :title (or title "Participant vs. cohort gene expression")}
                                     (dissoc opts :title))))
      (let [se (pqd/gene-expression db {:genes genes :samples sample-ids :measurement measurement})
            ce (vec (mapcat (fn [[label c]]
                              (map #(assoc % :cohort label)
                                   (pqd/gene-expression c {:genes genes :measurement (or cohort-measurement measurement)})))
                            cohorts))]
        (plot/vs-cohort se ce {:type type :title (or title "Samples vs cohort expression")
                               :xlab (str (name measurement) " (log scale)")})))))

(defn isoform-breakdown
  "Isoforms of a gene in a sample: [{:transcript-id :transcript-length
  :isoform-percent :effective-length} ...]."
  [db-or-name hgnc sample-id]
  (mapv (fn [[t l p e]] {:transcript-id t :transcript-length l :isoform-percent p :effective-length e})
        (pdb/q '[:find ?transcript ?tlen ?perc ?effective-len
               :in $ ?hgnc ?sample-id
               :where
               [?g :gene/hgnc-symbol ?hgnc]
               [?gp :gene-product/gene ?g]
               [?gp :gene-product/id ?transcript]
               [?gp :gene-product/transcript-length ?tlen]
               [?m :measurement/gene-product ?gp]
               [?m :measurement/isoform-percent ?perc]
               [?m :measurement/effective-transcript-length ?effective-len]
               [?m :measurement/sample ?s]
               [?s :sample/id ?sample-id]]
             (pdb/as-db db-or-name) hgnc sample-id)))

;; -- sample vs reference cohort, two-sample change (parity with R/Python
;;    analysis_expression: compare_to_cohort, top_by_zscore, compare_samples) --

(defn- sum-by-gene
  "gene-expression rows -> {gene summed value} (sums across measurement sets)."
  [rows]
  (reduce (fn [m {:keys [hgnc-symbol value]}] (update m hgnc-symbol (fnil + 0.0) (double value))) {} rows))

(defn sample-expression
  "{gene value} for one sample (summed across measurement sets).
  opts: :measurement (default :tpm), :measurement-set, :genes"
  ([db-or-name sample-id] (sample-expression db-or-name sample-id {}))
  ([db-or-name sample-id {:keys [measurement measurement-set genes] :or {measurement :tpm}}]
   (sum-by-gene (pqd/gene-expression db-or-name (cond-> {:samples [sample-id] :measurement measurement}
                                                  measurement-set (assoc :measurement-set measurement-set)
                                                  genes (assoc :genes genes))))))

(defn compare-to-cohort
  "Place a sample within a reference cohort's distribution, gene by gene:
  z-score and percentile of log2(1 + x) (with :log?, default true).

  opts:
    :measurement / :cohort-measurement  (default :tpm)
    :genes            restrict to these genes (default: all of the sample's)
    :fill-missing?    count cohort samples without a value for a gene as 0
                      (default true; the cohort size comes from :anchor-gene,
                      default \"GAPDH\")
    :sd-floor         lower bound on the cohort sd (default 0.25)
    :min-cohort       minimum cohort values per gene (default 5)
    :batch-size       genes per cohort query (default 2000)
  Returns maps sorted by |z| desc: :hgnc-symbol :value :z :percentile :cohort-n
  :cohort-observed :cohort-mean :cohort-sd :cohort-median (cohort stats on the
  transformed scale), with the comparison settings in metadata
  (:patternq/comparison)."
  [sample-id db cohort-db {:keys [measurement cohort-measurement genes log? fill-missing? anchor-gene
                                  sd-floor min-cohort batch-size]
                           :or {measurement :tpm log? true fill-missing? true anchor-gene "GAPDH"
                                sd-floor 0.25 min-cohort 5 batch-size 2000}}]
  (let [cohort-measurement (or cohort-measurement measurement)
        x (sample-expression db sample-id {:measurement measurement :genes genes})
        _ (when (empty? x) (throw (ex-info (str "No " (name measurement) " expression for sample " sample-id) {})))
        tr (if log? #(quant/log2 (+ 1.0 (max 0.0 (double %)))) double)
        n-all (when fill-missing?
                (let [n (count (distinct (map :sample-id (pqd/gene-expression cohort-db {:genes [anchor-gene]
                                                                                           :measurement cohort-measurement}))))]
                  (when (zero? n)
                    (throw (ex-info (str "No " (name cohort-measurement) " values for " anchor-gene " in cohort") {})))
                  n))
        rows (->> (partition-all batch-size (keys x))
                  (mapcat (fn [batch]
                            (let [cg (pqd/gene-expression cohort-db {:genes (vec batch) :measurement cohort-measurement})
                                  per (->> cg (group-by (juxt :sample-id :hgnc-symbol))
                                           (map (fn [[[_ g] rs]] [g (reduce + (map :value rs))]))
                                           (group-by first))]
                              (for [[g svs] per
                                    :let [vals (map (comp tr second) svs)
                                          observed (count vals)
                                          vals (if n-all (concat vals (repeat (max 0 (- (long n-all) observed)) (tr 0))) vals)
                                          v (x g)
                                          sdv (quant/sd vals)]
                                    :when (>= (count vals) (long min-cohort))]
                                {:hgnc-symbol g :value v
                                 :z (/ (- (double (tr v)) (quant/mean vals)) (max (double sd-floor) (if (Double/isNaN sdv) 0.0 sdv)))
                                 :percentile (quant/percentile-rank vals (tr v))
                                 :cohort-n (count vals) :cohort-observed observed
                                 :cohort-mean (quant/mean vals) :cohort-sd sdv :cohort-median (quant/median vals)}))))
                  (sort-by #(- (Math/abs (double (:z %)))))
                  vec)]
    (with-meta rows {:patternq/comparison {:sample sample-id :db db :cohort-db cohort-db :measurement measurement
                                           :cohort-measurement cohort-measurement :log? log? :cohort-size n-all}})))

(defn top-by-zscore
  "Top genes of a compare-to-cohort result. Keeps genes observed in at least
  :min-observed (default 0.5) of the cohort and expressed (sample value or
  cohort median >= :min-value, default 1).
  opts: :n (25), :direction (:both | :up | :down)"
  ([comparison] (top-by-zscore comparison {}))
  ([comparison {:keys [n direction min-value min-observed] :or {n 25 direction :both min-value 1 min-observed 0.5}}]
   (->> comparison
        (filter #(some? (:z %)))
        (filter #(>= (double (:cohort-observed %)) (* (double min-observed) (double (:cohort-n %)))))
        (filter #(or (>= (double (:value %)) (double min-value))
                     (>= (- (Math/pow 2.0 (double (:cohort-median %))) 1.0) (double min-value))))
        (filter #(case direction :up (pos? (double (:z %))) :down (neg? (double (:z %))) true))
        (sort-by (case direction :up #(- (double (:z %))) :down :z #(- (Math/abs (double (:z %))))))
        (take n)
        vec)))

(defn compare-samples
  "Expression change from sample a to sample b: [{:hgnc-symbol :value-a
  :value-b :lfc :avg-log10} ...] with lfc = log2((b + pc)/(a + pc)) and
  avg-log10 the mean of log10(1 + x); genes with avg-log10 < :min-avg (0.5)
  dropped; sorted by lfc desc. opts: :db-b, :measurement (:tpm), :min-avg,
  :pseudocount (1)."
  [sample-a sample-b db {:keys [db-b measurement min-avg pseudocount]
                         :or {measurement :tpm min-avg 0.5 pseudocount 1}}]
  (let [a (sample-expression db sample-a {:measurement measurement})
        b (sample-expression (or db-b db) sample-b {:measurement measurement})
        pc (double pseudocount)
        l10 #(Math/log10 (+ 1.0 (double %)))]
    (with-meta
      (->> (distinct (concat (keys a) (keys b)))
           (map (fn [g] (let [va (double (get a g 0.0)) vb (double (get b g 0.0))]
                          {:hgnc-symbol g :value-a va :value-b vb
                           :lfc (quant/log2 (/ (+ vb pc) (+ va pc)))
                           :avg-log10 (/ (+ (l10 va) (l10 vb)) 2.0)})))
           (filter #(>= (double (:avg-log10 %)) (double min-avg)))
           (sort-by #(- (double (:lfc %))))
           vec)
      {:patternq/comparison {:sample-a sample-a :sample-b sample-b :measurement measurement}})))
