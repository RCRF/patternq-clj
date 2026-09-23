(ns patternq.variants
  "Somatic variant analysis helpers, extracted from the unify-central variant
  forensics analysis (util.clj, template.clj, h37_cross_cohort.clj,
  unify-central.db.variant-queries). Participant-level fns take a subject
  id; cohort fns take database names/values explicitly."
  (:require [clojure.string :as str]
            [patternq.dataset :as pqd]
            [patternq.db :as pdb]
            [patternq.quant :as quant]
            [patternq.results :as res]))

(set! *warn-on-reflection* true)

(defn- exclude-clause [_] nil)

(defn- run-participant
  "Run a participant query whose second find element is the HGNC symbol,
  dropping `exclude-genes` (filtered client-side: the query service's parser
  has no set-membership predicate)."
  [db-or-name q participant-id exclude-genes]
  (let [ex (set exclude-genes)]
    (cond->> (pdb/q q (pdb/as-db db-or-name) participant-id)
      (seq ex) (remove #(ex (nth % 1))))))

(defn- in-clause [_] '[:in $ ?pid])

(defn participant-variants
  "Protein-altering variant measurements of a participant (variants with an
  HGVSp, excluding low impact/synonymous), with sample context:
  [{:variant-id :hgnc-symbol :HGVSp :sample-id :timepoint-id :vaf
    :metastasis :anatomic-site} ...]. Samples missing timepoint, metastasis or
  anatomic site still appear (values nil).

  opts: :exclude-genes (e.g. recurrent artifacts such as #{\"OR8U1\"})"
  ([db-or-name participant-id] (participant-variants db-or-name participant-id {}))
  ([db-or-name participant-id {:keys [exclude-genes]}]
   (let [q (vec (concat '[:find ?var-id ?hgnc ?hgvsp ?sid ?vaf ?s ?m] (in-clause exclude-genes)
                        '[:where
                          [?p :subject/id ?pid]
                          [?s :sample/subject ?p]
                          [?m :measurement/sample ?s]
                          [?m :measurement/variant ?v]
                          [?m :measurement/vaf ?vaf]
                          [?s :sample/id ?sid]
                          [?v :variant/id ?var-id]
                          [?v :variant/HGVSp ?hgvsp]
                          [?v :variant/gene ?g]
                          (not [?v :variant/impact :variant.impact/low])
                          [?g :gene/hgnc-symbol ?hgnc]]
                        (exclude-clause exclude-genes)))
         db (pdb/as-db db-or-name)
         rows (run-participant db q participant-id exclude-genes)
         smp (pdb/pull-map db '[:sample/metastasis :sample/freetext-anatomic-site
                                {:sample/timepoint [:timepoint/id]}]
                           (distinct (map #(nth % 5) rows)))]
     (vec (for [[vid hgnc hgvsp sid vaf s] rows
                :let [sm (smp s)]]
            {:variant-id vid :hgnc-symbol hgnc :HGVSp hgvsp :sample-id sid :vaf vaf
             :timepoint-id (get-in sm [:sample/timepoint :timepoint/id])
             :metastasis (:sample/metastasis sm)
             :anatomic-site (:sample/freetext-anatomic-site sm)})))))

(defn- impact-variants [db-or-name participant-id impact exclude-genes hgvsp?]
  (let [q (vec (concat [:find '?var-id '?hgnc '?sid '?vaf]
                       (when hgvsp? '[?hgvsp])
                       (in-clause exclude-genes)
                       '[:where
                         [?p :subject/id ?pid]
                         [?s :sample/subject ?p]
                         [?m :measurement/sample ?s]
                         [?m :measurement/variant ?v]
                         [?m :measurement/vaf ?vaf]
                         [?s :sample/id ?sid]
                         [?v :variant/id ?var-id]
                         [?v :variant/gene ?g]
                         [?g :gene/hgnc-symbol ?hgnc]]
                       [['?v :variant/impact impact]]
                       (when hgvsp? '[[?v :variant/HGVSp ?hgvsp]])
                       (exclude-clause exclude-genes)))]
    (mapv #(zipmap (cond-> [:variant-id :hgnc-symbol :sample-id :vaf] hgvsp? (conj :HGVSp)) %)
          (run-participant db-or-name q participant-id exclude-genes))))

(defn missense-variants
  "Alias of participant-variants (protein-altering, non-low-impact)."
  ([db-or-name participant-id] (participant-variants db-or-name participant-id {}))
  ([db-or-name participant-id opts] (participant-variants db-or-name participant-id opts)))

(defn modifier-variants
  "Modifier-impact variants of a participant (with SO consequences).
  opts: :exclude-genes"
  ([db-or-name participant-id] (modifier-variants db-or-name participant-id {}))
  ([db-or-name participant-id {:keys [exclude-genes]}]
   (let [db (pdb/as-db db-or-name)
         rows (impact-variants db participant-id :variant.impact/modifier exclude-genes false)
         conseq (update-vals (pdb/pull-by db '[{:variant/so-consequences [:so-sequence-feature/name]}]
                                          :variant/id (distinct (map :variant-id rows)))
                             #(res/join-many (map :so-sequence-feature/name (:variant/so-consequences %))))]
     (mapv #(assoc % :so-consequences (conseq (:variant-id %))) rows))))

(defn synonymous-variants
  "Low-impact (synonymous) variants of a participant. opts: :exclude-genes"
  ([db-or-name participant-id] (synonymous-variants db-or-name participant-id {}))
  ([db-or-name participant-id {:keys [exclude-genes]}]
   (impact-variants db-or-name participant-id :variant.impact/low exclude-genes true)))

(defn vaf-histogram-data
  "[[[start stop] count] ...] of a sample's VAFs in bins of `bin-size`
  (default 0.05)."
  ([db-or-name sample-id] (vaf-histogram-data db-or-name sample-id 0.05))
  ([db-or-name sample-id bin-size]
   (let [vals (pdb/q '[:find ?vaf (count ?var-id)
                     :in $ ?sample-id
                     :where
                     [?s :sample/id ?sample-id]
                     [?m :measurement/sample ?s]
                     [?m :measurement/variant ?v]
                     [?m :measurement/vaf ?vaf]
                     [?v :variant/id ?var-id]]
                   (pdb/as-db db-or-name) sample-id)]
     (quant/->histogram vals (vec (range 0.0 (+ 1.0 (double bin-size)) bin-size))))))

(defn trunk-candidates
  "Candidate clonal trunk variants: variants observed in at least
  `min-samples` (default 2) samples of a participant, ordered by sample count
  then summed VAF (desc). Input: participant-variants output.
  Returns [{:variant-id :hgnc-symbol :HGVSp :vafs :sample-count} ...]."
  ([variant-rows] (trunk-candidates variant-rows {}))
  ([variant-rows {:keys [min-samples] :or {min-samples 2}}]
   (->> (group-by :variant-id variant-rows)
        (map (fn [[vid rows]]
               {:variant-id vid :hgnc-symbol (:hgnc-symbol (first rows)) :HGVSp (:HGVSp (first rows))
                :vafs (mapv :vaf rows) :sample-count (count (distinct (map :sample-id rows)))}))
        (filter #(>= (long (:sample-count %)) (long min-samples)))
        (sort-by (juxt :sample-count #(reduce + (:vafs %))))
        reverse
        vec)))

(defn branch-variants
  "Non-trunk variants grouped by gene, with where they occur: [{:hgnc-symbol
  :HGVSp [..] :sites {timepoint site} :vafs {timepoint vaf}} ...] ordered by
  number of observations then summed VAF (desc)."
  [variant-rows trunk]
  (let [trunk-ids (set (map :variant-id trunk))]
    (->> variant-rows
         (remove (comp trunk-ids :variant-id))
         (group-by :hgnc-symbol)
         (map (fn [[g rows]]
                {:hgnc-symbol g
                 :HGVSp (vec (distinct (map :HGVSp rows)))
                 :sites (zipmap (map :timepoint-id rows) (map :anatomic-site rows))
                 :vafs (zipmap (map :timepoint-id rows) (map :vaf rows))}))
         (sort-by (juxt (comp count :vafs) #(reduce + (vals (:vafs %)))))
         reverse
         vec)))

(defn variants-by-impact
  "Variant measurements (all participants) at an impact level (\"high\",
  \"moderate\", \"low\", \"modifier\" or the ident), optionally restricted to
  timepoint ids: [{:sample-id :measurement-set :variant-id :hgnc-symbol
  :consequence :impact :vaf :alt :ref :timepoint-id :subject-id} ...]."
  ([db-or-name impact] (variants-by-impact db-or-name impact nil))
  ([db-or-name impact timepoint-ids]
   (let [impact (if (keyword? impact) impact (keyword "variant.impact" (name impact)))
         q (vec (concat '[:find ?sample-id ?ms-name ?var-id ?hgnc ?consequence ?impact-name ?vaf ?alt ?ref ?tp-id ?pid
                          :in $ ?impact]
                        (when (seq timepoint-ids) '[[?tp-id ...]])
                        '[:where
                          [?var :variant/impact ?impact]
                          [?impact :db/ident ?impact-name]
                          [?var :variant/id ?var-id]
                          [?var :variant/alt-allele ?alt]
                          [?var :variant/ref-allele ?ref]
                          [?var :variant/so-consequences ?soc]
                          [?soc :so-sequence-feature/name ?consequence]
                          [?var :variant/gene ?g]
                          [?g :gene/hgnc-symbol ?hgnc]
                          [?m :measurement/variant ?var]
                          [?m :measurement/vaf ?vaf]
                          [?ms :measurement-set/measurements ?m]
                          [?ms :measurement-set/name ?ms-name]
                          [?m :measurement/sample ?s]
                          [?s :sample/id ?sample-id]
                          [?s :sample/timepoint ?tp]
                          [?s :sample/subject ?p]
                          [?p :subject/id ?pid]
                          [?tp :timepoint/id ?tp-id]]))
         rows (apply pdb/q q (pdb/as-db db-or-name) impact (when (seq timepoint-ids) [(vec timepoint-ids)]))]
     (mapv #(update (zipmap [:sample-id :measurement-set :variant-id :hgnc-symbol :consequence :impact
                             :vaf :alt :ref :timepoint-id :subject-id] %)
                    :impact res/ident-name)
           rows))))

(defn patients-with-variants
  "Subject ids with at least one variant measurement."
  [db-or-name]
  (vec (sort (pdb/q '[:find [?pid ...]
                    :where [?m :measurement/variant] [?m :measurement/sample ?s]
                    [?s :sample/subject ?p] [?p :subject/id ?pid]]
                  (pdb/as-db db-or-name)))))

(defn variant-samples
  "Samples with variant measurements: [{:subject-id :sample-id :anatomic-site
  :timepoint-id :timepoint-order} ...]."
  [db-or-name]
  (let [db (pdb/as-db db-or-name)]
    (->> (pdb/q '[:find ?pid ?sid ?s
                :where [?m :measurement/variant] [?m :measurement/sample ?s]
                [?s :sample/id ?sid] [?s :sample/subject ?p] [?p :subject/id ?pid]]
              db)
         ((fn [rows]
            (let [pm (pdb/pull-map db '[:sample/freetext-anatomic-site
                                        {:sample/timepoint [:timepoint/id :timepoint/relative-order]}]
                                   (distinct (map #(nth % 2) rows)))]
              (map #(conj (vec %) (pm (nth % 2))) rows))))
         (map (fn [[pid sid s e]]
                (let []
                  {:subject-id pid :sample-id sid :anatomic-site (:sample/freetext-anatomic-site e)
                   :timepoint-id (get-in e [:sample/timepoint :timepoint/id])
                   :timepoint-order (get-in e [:sample/timepoint :timepoint/relative-order])})))
         (sort-by (juxt :subject-id #(or (:timepoint-order %) 0)))
         vec)))

(defn cohort-variants-by-gene
  "All variant observations of one gene across a (cohort) database:
  [{:subject-id :sample-id :anatomic-site :hgnc-symbol :HGVSp :variant-id} ...]
  (HGVSp \"non-coding\" when absent)."
  [db-or-name hgnc]
  (let [db (pdb/as-db db-or-name)]
    (->> (pdb/q '[:find ?pid ?sid ?s ?hgnc ?hgvsp ?var-id
                :in $ ?hgnc
                :where
                [?g :gene/hgnc-symbol ?hgnc]
                [?v :variant/gene ?g]
                [?v :variant/id ?var-id]
                (or-join [?v ?hgvsp]
                         [?v :variant/HGVSp ?hgvsp]
                         (and (not [?v :variant/HGVSp]) [(ground "non-coding") ?hgvsp]))
                [?m :measurement/variant ?v]
                [?m :measurement/sample ?s]
                [?s :sample/id ?sid]
                [?s :sample/subject ?p]
                [?p :subject/id ?pid]]
              db hgnc)
         ((fn [rows]
            (let [pm (pdb/pull-map db [:sample/freetext-anatomic-site] (distinct (map #(nth % 2) rows)))]
              (map #(conj (vec %) (pm (nth % 2))) rows))))
         (map (fn [[pid sid s g hgvsp vid e]]
                {:subject-id pid :sample-id sid :anatomic-site (:sample/freetext-anatomic-site e)
                 :hgnc-symbol g :HGVSp hgvsp :variant-id vid}))
         (sort-by (juxt :subject-id :sample-id))
         vec)))

(defn mutated-genes
  "[[hgnc variant-id] ...] for every variant entity in a database."
  [db-or-name]
  (vec (pdb/q '[:find ?hgnc ?var-id :where [?v :variant/id ?var-id] [?v :variant/gene ?g] [?g :gene/hgnc-symbol ?hgnc]]
            (pdb/as-db db-or-name))))

(defn variant-patient-count
  "Number of participants with a measured variant in a gene."
  [db-or-name hgnc]
  (or (pdb/q '[:find (count ?p) .
             :in $ ?hgnc
             :where [?g :gene/hgnc-symbol ?hgnc] [?v :variant/gene ?g] [?m :measurement/variant ?v]
             [?m :measurement/sample ?s] [?s :sample/subject ?p]]
           (pdb/as-db db-or-name) hgnc)
      0))

(defn variant-cohort-counts
  "Per gene: {:hgnc-symbol :participant-count :gene-size} in a cohort
  database, sorted by participant count (desc). Gene size helps normalize
  comparisons in high mutational burden cohorts."
  [db-or-name genes]
  (let [db (pdb/as-db db-or-name)
        size (fn [g] (when-let [[a b] (first (pdb/q '[:find ?start ?stop :in $ ?hgnc
                                                    :where [?g :gene/hgnc-symbol ?hgnc] [?g :gene/genomic-coordinates ?gc]
                                                    [?gc :genomic-coordinate/start ?start] [?gc :genomic-coordinate/end ?stop]]
                                                  db g))]
                       (- (long b) (long a))))]
    (->> genes
         (map (fn [g] {:hgnc-symbol g :participant-count (variant-patient-count db g) :gene-size (size g)}))
         (sort-by (comp - :participant-count))
         vec)))

(defn variant-genes-in-cohorts
  "For each cohort database name: {db-name {gene [[hgnc subject-id sample-id] ...]}}
  of observations of the given genes (cohorts without any are omitted)."
  [cohort-db-names genes]
  (into {}
        (for [n cohort-db-names
              :let [rows (pdb/q '[:find ?hgnc ?pid ?sid
                                :in $ [?hgnc ...]
                                :with ?m
                                :where [?g :gene/hgnc-symbol ?hgnc] [?v :variant/gene ?g] [?m :measurement/variant ?v]
                                [?m :measurement/sample ?s] [?s :sample/subject ?p] [?s :sample/id ?sid] [?p :subject/id ?pid]]
                              (pdb/db n) (vec genes))]
              :when (seq rows)]
          [n (group-by first rows)])))

(defn to-igv-coords
  "\"GRCh38:chr17:+:7675088:7675088/C/T\" -> \"chr17:7675088-7675088\"."
  [variant-id]
  (let [[_ chr _ start stop*] (str/split variant-id #":")]
    (format "%s:%s-%s" chr start (first (str/split stop* #"/")))))

(defn find-variants
  "Per-measurement variant rows for a participant with sample, timepoint and
  the copy number of CNV segments covering the gene in the same sample:
  [{:sample-id :timepoint-id :measurement-set :variant-id :hgnc-symbol :HGVSp
    :impact :classification :so-consequences :vaf :t-depth :cnv} ...].
  :cnv is \"cnv-id (absolute-cn)\" joined with \", \", or nil. opts: :samples"
  ([db-or-name participant-id] (find-variants db-or-name participant-id {}))
  ([db-or-name participant-id {:keys [samples]}]
   (let [db (pdb/as-db db-or-name)
         sids (or samples (pdb/q '[:find [?sid ...] :in $ ?pid
                                 :where [?p :subject/id ?pid] [?s :sample/subject ?p] [?s :sample/id ?sid]]
                               db participant-id))
         tps (into {} (map (juxt :sample-id :timepoint-id)) (pqd/samples db))
         vs (pqd/variants db {:samples sids})
         ann (update-vals (pdb/pull-by db '[{:variant/classification [:db/ident]}
                                            {:variant/so-consequences [:so-sequence-feature/name]}]
                                       :variant/id (distinct (map :variant-id vs)))
                          res/flatten-pull)
         genes (distinct (keep :hgnc-symbol vs))
         cnv (when (seq genes)
               (group-by (juxt :sample-id :hgnc-symbol)
                         (pqd/cnv-segments db {:genes genes :samples sids})))]
     (vec (for [v vs
                :let [a (ann (:variant-id v))
                      segs (get cnv [(:sample-id v) (:hgnc-symbol v)])]]
            (assoc v
              :timepoint-id (tps (:sample-id v))
              :classification (:variant-classification a)
              :so-consequences (res/join-many (:variant-so-consequences a))
              :cnv (when (seq segs)
                     (str/join ", " (map #(str (:cnv-id %) " (" (or (:absolute-cn %) (:segment-mean-lrr %)) ")")
                                         (sort-by :cnv-id segs))))))))))
