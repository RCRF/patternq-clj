(ns patternq.dataset
  "Canned queries over a dataset database: structure, samples & subjects,
  measurements, variants, gene expression, copy number.

  Every dataset is its own database, so queries start from samples,
  subjects, assays and measurement sets directly. Functions take a database
  value or a database name (plus an opts map), and return vectors of maps
  with provenance metadata (see patternq.db/provenance). Most have a *-query
  companion giving the query and its arguments as data.

  Mirrors the R/Python patternq catalog (see PARITY.md): R column
  `sample_id` is key :sample-id here."
  (:require [clojure.string :as str]
            [patternq.db :as pdb]
            [patternq.results :as res]))

(set! *warn-on-reflection* true)

(defn run
  "Run `query` with `args` against db-or-name; returns [result db db-name]."
  [db-or-name query args]
  (let [db (pdb/as-db db-or-name)]
    [(apply pdb/q query db args) db (when (string? db-or-name) db-or-name)]))

(defn- provenanced [rows db db-name]
  (pdb/with-provenance (vec rows) db db-name))

(defn pull-rows
  "Run a query whose find spec is a single pull expression; return flattened
  maps (enum refs resolved)."
  [db-or-name query args]
  (let [[r db db-name] (run db-or-name query args)]
    ;; enum refs resolved once for the whole result (one extra query at most)
    (provenanced (map res/flatten-pull (res/resolve-enum-refs db (mapv first r))) db db-name)))

(defn- measurement-attr [measurement]
  (keyword "measurement" (name measurement)))

;; -- samples & subjects --

(def samples-query
  '[:find (pull ?s [* {:sample/subject [:subject/id]}
                    {:sample/timepoint [:timepoint/id :timepoint/relative-order]}
                    {:sample/specimen [:db/ident]}
                    {:sample/type [:db/ident]}
                    {:sample/container [:db/ident]}
                    {:sample/study-day [:study-day/id]}
                    {:sample/gdc-anatomic-site [:gdc-anatomic-site/name]}])
    :where [?s :sample/id]])

(defn samples
  "One map per sample: :sample-id, :subject-id, :timepoint-id and the other
  sample attributes present."
  [db-or-name]
  (pull-rows db-or-name samples-query []))

(def subjects-query
  '[:find (pull ?s [* {:subject/sex [:db/ident]}
                    {:subject/race [:db/ident]}
                    {:subject/ethnicity [:db/ident]}
                    {:subject/smoker [:db/ident]}
                    {:subject/cause-of-death [:db/ident]}
                    {:subject/meddra-disease [:meddra-disease/preferred-name]}])
    :where [?s :subject/id]])

(defn subjects
  "One map per subject: :subject-id and demographic attributes present (enums
  as names, e.g. :subject-sex \"female\")."
  [db-or-name]
  (pull-rows db-or-name subjects-query []))

(defn timepoints
  "Timepoints ordered by :timepoint-relative-order: :timepoint-id,
  :timepoint-relative-order, :timepoint-type, :timepoint-offset, cycle/day."
  [db-or-name]
  (let [rows (pull-rows db-or-name '[:find (pull ?t [* {:timepoint/type [:db/ident]}])
                                     :where [?t :timepoint/id]] [])]
    (with-meta (vec (sort-by (fn [r] (or (:timepoint-relative-order r) Long/MAX_VALUE)) rows))
               (meta rows))))

;; -- dataset structure --

(def dataset-summary-query
  '[:find ?assay-name ?assay-technology ?measurement-set-name
    :where
    [?a :assay/name ?assay-name]
    [?a :assay/technology ?t]
    [?t :db/ident ?assay-technology]
    [?a :assay/measurement-sets ?ms]
    [?ms :measurement-set/name ?measurement-set-name]])

(defn dataset-summary
  "One map per measurement set: :assay-name, :assay-technology,
  :measurement-set-name."
  [db-or-name]
  (let [[r db db-name] (run db-or-name dataset-summary-query [])]
    (provenanced (->> r
                      (map (fn [[a t ms]] {:assay-name a
                                           :assay-technology (res/ident-name t)
                                           :measurement-set-name ms}))
                      (sort-by (juxt :assay-name :measurement-set-name)))
                 db db-name)))

(defn dataset-info
  "The dataset entity stored in the database: :dataset-name,
  :dataset-description, :dataset-doi, :dataset-url."
  [db-or-name]
  (first (pull-rows db-or-name '[:find (pull ?d [:dataset/name :dataset/description
                                                 :dataset/doi :dataset/url])
                                 :where [?d :dataset/name]] [])))

(defn schema-info
  "{:name :version} of the database's Unify schema."
  [db-or-name]
  (let [[r] (run db-or-name '[:find ?name ?version
                              :where
                              [?e :unify.schema/version ?version]
                              [?e :unify.schema/name ?name]] [])
        [n v] (first r)]
    {:name n :version v}))

(defn measurement-sets
  "Measurement sets with assay and size: :assay-name, :assay-technology,
  :measurement-set-name, :measurement-count."
  [db-or-name]
  (let [db (pdb/as-db db-or-name)
        counts (into {} (pdb/q '[:find ?name (count ?m)
                               :where
                               [?ms :measurement-set/name ?name]
                               [?ms :measurement-set/measurements ?m]]
                             db))]
    (provenanced (map #(assoc % :measurement-count (get counts (:measurement-set-name %) 0))
                      (dataset-summary db))
                 db (when (string? db-or-name) db-or-name))))

;; Measurement target references: the entity a measurement is "of", and the
;; clauses naming it. Keys of results are the var names without "?".
(def measurement-targets
  {:measurement/gene-product {:clauses '[[?tgp :gene-product/gene ?tg] [?tg :gene/hgnc-symbol ?hgnc-symbol]]
                              :ref '?tgp :var '?hgnc-symbol}
   :measurement/variant {:clauses '[[?tv :variant/id ?variant-id]] :ref '?tv :var '?variant-id}
   :measurement/cnv {:clauses '[[?tc :cnv/id ?cnv-id]] :ref '?tc :var '?cnv-id}
   :measurement/epitope {:clauses '[[?te :epitope/id ?epitope-id]] :ref '?te :var '?epitope-id}
   :measurement/cell-population {:clauses '[[?tcp :cell-population/name ?cell-population]]
                                 :ref '?tcp :var '?cell-population}
   :measurement/tcr {:clauses '[[?tt :tcr/id ?tcr-id]] :ref '?tt :var '?tcr-id}
   :measurement/otu {:clauses '[[?to :otu/id ?otu-id]] :ref '?to :var '?otu-id}
   :measurement/sgb {:clauses '[[?ts :sgb/metaphlan-id ?sgb-id]] :ref '?ts :var '?sgb-id}
   :measurement/pathway {:clauses '[[?tp :pathway/id ?pathway-id]] :ref '?tp :var '?pathway-id}
   :measurement/metabolite-feature {:clauses '[[?tmf :metabolite-feature/rt-mz-peak ?metabolite-feature]]
                                    :ref '?tmf :var '?metabolite-feature}
   :measurement/nanostring-signature {:clauses '[[?tns :nanostring-signature/name ?signature]]
                                      :ref '?tns :var '?signature}
   :measurement/atac-peak {:clauses '[[?tap :atac-peak/name ?atac-peak]] :ref '?tap :var '?atac-peak}
   :measurement/single-cell {:clauses '[[?tsc :single-cell/id ?single-cell-id]]
                             :ref '?tsc :var '?single-cell-id}})

(def enum-measurement-attrs #{:measurement/cnv-call :measurement/msi-status})

(defn measurement-types
  "Counts of every attribute across a measurement set's measurements: maps of
  :attribute (without namespace), :kind (\"value\" or \"target\"), :count.
  Sets can mix kinds of measurement, e.g. CyTOF populations with
  percent-of-parent alongside population x marker median-channel-value."
  [db-or-name measurement-set]
  (let [[r db db-name] (run db-or-name '[:find ?attr (count ?m)
                                         :in $ ?ms-name
                                         :where
                                         [?ms :measurement-set/name ?ms-name]
                                         [?ms :measurement-set/measurements ?m]
                                         [?m ?a]
                                         [?a :db/ident ?attr]]
                            [measurement-set])]
    (provenanced (->> r
                      (remove (comp #{:measurement/id :measurement/uid :measurement/sample} first))
                      (map (fn [[a n]] {:attribute (name a)
                                        :kind (if (contains? measurement-targets a) "target" "value")
                                        :count n}))
                      (sort-by (juxt :kind (comp - :count))))
                 db db-name)))

(defn measurement-set-attributes
  "Attributes of a random sample of `:n` (default 200) measurements of a set,
  optionally only those carrying `:measurement`. Cheap for big sets."
  ([db-or-name measurement-set] (measurement-set-attributes db-or-name measurement-set {}))
  ([db-or-name measurement-set {:keys [measurement n] :or {n 200}}]
   (let [db (pdb/as-db db-or-name)
         where (cond-> '[[?ms :measurement-set/name ?ms-name]
                         [?ms :measurement-set/measurements ?m]]
                 measurement (conj ['?m (measurement-attr measurement)]))
         eids (pdb/q {:find [(list 'sample n '?m) '.] :in '[$ ?ms-name] :where where} db measurement-set)]
     (->> (pdb/pull-many db '[*] (vec eids))
          (mapcat keys)
          (remove #{:db/id :measurement/id :measurement/uid})
          distinct
          sort
          vec))))

(defn measurements-query
  "Query and args for values of `measurement` (e.g. :percent-of-parent) with
  target columns for `targets` (target reference attributes)."
  [measurement {:keys [measurement-set samples targets]}]
  (let [attr (measurement-attr measurement)
        enum? (contains? enum-measurement-attrs attr)
        specs (map (fn [t] (or (measurement-targets t)
                               (throw (ex-info (str "Unknown measurement target " t) {:target t}))))
                   targets)
        where (concat [['?m attr (if enum? '?value-ref '?value)]
                       '[?m :measurement/sample ?s]
                       '[?s :sample/id ?sample-id]
                       '[?ms :measurement-set/measurements ?m]
                       '[?ms :measurement-set/name ?measurement-set]]
                      (when enum? '[[?value-ref :db/ident ?value]])
                      (mapcat (fn [t spec] (cons ['?m t (:ref spec)] (:clauses spec))) targets specs))
        ins (cond-> '[$]
              measurement-set (conj '?measurement-set)
              samples (conj '[?sample-id ...]))]
    {:query (vec (concat [:find '?sample-id '?measurement-set] (map :var specs) '[?value :with ?m :in]
                         ins [:where] where))
     :args (cond-> []
             measurement-set (conj measurement-set)
             samples (conj (vec samples)))
     :keys (vec (concat [:sample-id :measurement-set] (map (comp keyword #(subs % 1) str :var) specs) [:value]))}))

(defn measurements
  "Values of one measurement attribute (e.g. :tpm, :percent-of-parent,
  :olink-npx, :median-channel-value) with what each measurement is of.
  Targets are detected from a sample of the set's measurements carrying the
  attribute unless given.

  opts: :measurement-set (required unless :targets given), :samples, :targets
  Returns long maps: :sample-id, :measurement-set, target key(s), :value.
  Use patternq.results/->matrix for a wide matrix."
  ([db-or-name measurement] (measurements db-or-name measurement {}))
  ([db-or-name measurement {:keys [measurement-set targets] :as opts}]
   (let [db (pdb/as-db db-or-name)
         targets (or targets
                     (if measurement-set
                       (filterv (set (keys measurement-targets))
                                (measurement-set-attributes db measurement-set {:measurement measurement}))
                       (throw (ex-info "Give :measurement-set (or :targets) so targets can be detected" {}))))
         {:keys [query args keys]} (measurements-query measurement (assoc opts :targets targets))
         r (apply pdb/q query db args)]
     (provenanced (map (fn [row] (update (zipmap keys row) :value res/ident-name)) r)
                  db (when (string? db-or-name) db-or-name)))))

(defn sample-assays
  "Which samples were measured in which measurement sets: :subject-id,
  :sample-id, :timepoint-id, :assay-name, :assay-technology,
  :measurement-set-name."
  [db-or-name]
  (let [db (pdb/as-db db-or-name)
        smp (into {} (map (juxt :sample-id identity)) (samples db))]
    (provenanced
      (for [{:keys [measurement-set-name] :as ms} (dataset-summary db)
            sid (pdb/q '[:find [?sample-id ...]
                       :in $ ?ms-name
                       :where
                       [?ms :measurement-set/name ?ms-name]
                       [?ms :measurement-set/measurements ?m]
                       [?m :measurement/sample ?s]
                       [?s :sample/id ?sample-id]]
                     db measurement-set-name)]
        (merge (select-keys (smp sid) [:subject-id :timepoint-id])
               {:sample-id sid}
               ms))
      db (when (string? db-or-name) db-or-name))))

(def measurement-matrices-query
  '[:find ?assay-name ?measurement-set-name ?matrix-name ?measurement-type ?matrix-key
    :where
    [?a :assay/name ?assay-name]
    [?a :assay/measurement-sets ?ms]
    [?ms :measurement-set/name ?measurement-set-name]
    [?ms :measurement-set/measurement-matrices ?mm]
    [?mm :measurement-matrix/name ?matrix-name]
    [?mm :measurement-matrix/measurement-type ?mt]
    [?mt :db/ident ?measurement-type]
    [?mm :measurement-matrix/backing-file ?matrix-key]])

(defn measurement-matrices
  "File-backed measurement matrices: :assay-name, :measurement-set-name,
  :matrix-name, :measurement-type, :matrix-key. Download with
  patternq.http/measurement-matrix."
  [db-or-name]
  (let [[r db db-name] (run db-or-name measurement-matrices-query [])]
    (provenanced (map #(update (zipmap [:assay-name :measurement-set-name :matrix-name
                                        :measurement-type :matrix-key] %)
                               :measurement-type res/ident-name)
                      r)
                 db db-name)))

;; -- variants --

(defn variants-query
  "Query and args for variant measurements, optionally restricted to
  `samples`, `genes` (HGNC symbols) and a `measurement-set` name."
  [{:keys [samples genes measurement-set]}]
  (let [ins (cond-> '[$]
              samples (conj '[?sample-id ...])
              genes (conj '[?gene ...])
              measurement-set (conj '?measurement-set))
        args (cond-> []
               samples (conj (vec samples))
               genes (conj (vec genes))
               measurement-set (conj measurement-set))
        where (cond-> '[[?m :measurement/vaf ?vaf]
                        [?m :measurement/variant ?v]
                        [?m :measurement/sample ?s]
                        [?s :sample/id ?sample-id]
                        [?ms :measurement-set/measurements ?m]
                        [?ms :measurement-set/name ?measurement-set]]
                genes (into '[[?v :variant/gene ?g]
                              [?g :gene/hgnc-symbol ?gene]]))]
    {:query (into '[:find ?sample-id ?measurement-set ?vaf
                    (pull ?v [:variant/id :variant/HGVSp :variant/HGVSc
                              {:variant/gene [:gene/hgnc-symbol]}
                              {:variant/impact [:db/ident]}])
                    (pull ?m [:measurement/t-depth])
                    :in]
                  (concat ins [:where] where))
     :args args}))

(defn variants
  "Somatic variant measurements, one map per measurement: :sample-id,
  :measurement-set, :variant-id, :hgnc-symbol, :HGVSp, :HGVSc, :impact, :vaf,
  :t-depth (where present). HGVSp/HGVSc (cardinality many) are joined with
  \"; \".

  opts: :samples, :genes, :measurement-set"
  ([db-or-name] (variants db-or-name {}))
  ([db-or-name opts]
   (let [{:keys [query args]} (variants-query opts)
         [r db db-name] (run db-or-name query args)]
     (provenanced
       (map (fn [[sample-id ms vaf v m]]
              (let [fv (res/flatten-pull v)]
                (cond-> {:sample-id sample-id
                         :measurement-set ms
                         :variant-id (:variant-id fv)
                         :hgnc-symbol (:gene-hgnc-symbol fv)
                         :HGVSp (res/join-many (:variant-HGVSp fv))
                         :HGVSc (res/join-many (:variant-HGVSc fv))
                         :impact (:variant-impact fv)
                         :vaf vaf}
                  (:measurement/t-depth m) (assoc :t-depth (:measurement/t-depth m)))))
            r)
       db db-name))))

;; -- gene expression --

(defn gene-expression-query
  "Query and args for gene expression values of `measurement` (e.g. :tpm,
  \"rsem-normalized-count\"), optionally restricted to `genes`, `samples` and a
  `measurement-set` name."
  [{:keys [genes samples measurement measurement-set] :or {measurement :tpm}}]
  (let [attr (measurement-attr measurement)
        ins (cond-> '[$]
              genes (conj '[?hgnc-symbol ...])
              samples (conj '[?sample-id ...])
              measurement-set (conj '?measurement-set))
        args (cond-> []
               genes (conj (vec genes))
               samples (conj (vec samples))
               measurement-set (conj measurement-set))]
    {:query (into '[:find ?sample-id ?hgnc-symbol ?measurement-set ?value
                    :with ?m
                    :in]
                  (concat ins
                          [:where
                           '[?g :gene/hgnc-symbol ?hgnc-symbol]
                           '[?gp :gene-product/gene ?g]
                           '[?m :measurement/gene-product ?gp]
                           ['?m attr '?value]
                           '[?m :measurement/sample ?s]
                           '[?s :sample/id ?sample-id]
                           '[?ms :measurement-set/measurements ?m]
                           '[?ms :measurement-set/name ?measurement-set]]))
     :args args}))

(defn gene-expression
  "Gene expression, long format: maps of :sample-id, :hgnc-symbol,
  :measurement-set, :value.

  opts: :genes, :samples, :measurement (default :tpm), :measurement-set"
  ([db-or-name] (gene-expression db-or-name {}))
  ([db-or-name opts]
   (let [{:keys [query args]} (gene-expression-query opts)
         [r db db-name] (run db-or-name query args)]
     (provenanced (res/rows->maps [:sample-id :hgnc-symbol :measurement-set :value] r) db db-name))))

(defn isoforms
  "Isoform-level expression for a gene: :sample-id, :transcript-id,
  :transcript-length, :isoform-percent, :effective-length.

  opts: :samples"
  ([db-or-name gene] (isoforms db-or-name gene {}))
  ([db-or-name gene {:keys [samples]}]
   (let [q (vec (concat '[:find ?sample-id ?transcript-id ?transcript-length ?isoform-percent ?effective-length
                          :with ?m :in $ ?hgnc]
                        (when samples '[[?sample-id ...]])
                        '[:where
                          [?g :gene/hgnc-symbol ?hgnc]
                          [?gp :gene-product/gene ?g]
                          [?gp :gene-product/id ?transcript-id]
                          [?gp :gene-product/transcript-length ?transcript-length]
                          [?m :measurement/gene-product ?gp]
                          [?m :measurement/isoform-percent ?isoform-percent]
                          [?m :measurement/effective-transcript-length ?effective-length]
                          [?m :measurement/sample ?s]
                          [?s :sample/id ?sample-id]]))
         [r db db-name] (run db-or-name q (cond-> [gene] samples (conj (vec samples))))]
     (provenanced (res/rows->maps [:sample-id :transcript-id :transcript-length :isoform-percent
                                   :effective-length] r)
                  db db-name))))

;; -- copy number --
;;
;; CNV data is large (segments x samples, or genes x samples): the CNV
;; queries insist on a subset (genes, samples or subjects). Same rule in the
;; R and Python libraries.

(defn- require-cnv-subset! [{:keys [genes samples subjects]}]
  (when-not (or genes samples subjects)
    (throw (ex-info "CNV queries need a subset: give :genes, :samples or :subjects"
                    {:patternq/error :cnv-subset-required}))))

(defn- subset-clauses
  "[ins args where] for optional :samples / :subjects filters; `where` binds
  ?s (the sample entity)."
  [{:keys [samples subjects]}]
  [(cond-> [] samples (conj '[?sample-id ...]) subjects (conj '[?subject-id ...]))
   (cond-> [] samples (conj (vec samples)) subjects (conj (vec subjects)))
   (vec (concat (when subjects '[[?p :subject/id ?subject-id] [?s :sample/subject ?p]])
                (when samples '[[?s :sample/id ?sample-id]])))])

(defn cnv-segments-query
  [{:keys [genes] :as opts}]
  (require-cnv-subset! opts)
  (let [[ins args where] (subset-clauses opts)
        gene-where (when genes '[[?g :gene/hgnc-symbol ?hgnc-symbol] [?c :cnv/genes ?g]])]
    {:query (vec (concat '[:find ?sample-id ?measurement-set ?cnv-id]
                         (when genes '[?hgnc-symbol])
                         '[(pull ?c [{:cnv/genomic-coordinates [:genomic-coordinate/contig
                                                                :genomic-coordinate/start
                                                                :genomic-coordinate/end]}])
                           (pull ?m [:measurement/segment-mean-lrr :measurement/absolute-cn
                                     :measurement/a-allele-cn :measurement/b-allele-cn
                                     :measurement/loh])
                           :in $]
                         (when genes '[[?hgnc-symbol ...]])
                         ins
                         [:where]
                         ;; most selective subset first: genes, else samples/subjects
                         (if genes
                           (concat gene-where '[[?m :measurement/cnv ?c] [?m :measurement/sample ?s]] where)
                           (concat where '[[?m :measurement/sample ?s] [?m :measurement/cnv ?c]]))
                         '[[?c :cnv/id ?cnv-id]
                           [?s :sample/id ?sample-id]
                           [?ms :measurement-set/measurements ?m]
                           [?ms :measurement-set/name ?measurement-set]]))
     :args (vec (concat (when genes [(vec genes)]) args))}))

(defn cnv-segments
  "Copy number segments per sample (segment-level CNV measurements) with the
  genes each segment overlaps. Requires at least one of :genes, :samples,
  :subjects. Returns :sample-id, :measurement-set, :cnv-id, :hgnc-symbol (with
  :genes), :contig, :start, :end, :segment-mean-lrr, :absolute-cn (where
  present)."
  [db-or-name opts]
  (let [{:keys [query args]} (cnv-segments-query opts)
        genes? (boolean (:genes opts))
        [r db db-name] (run db-or-name query args)]
    (provenanced
      (map (fn [row]
             (let [[sid ms cnv-id & more] row
                   [hgnc c m] (if genes? more (cons nil more))
                   coords (get c :cnv/genomic-coordinates)]
               (cond-> (merge {:sample-id sid :measurement-set ms :cnv-id cnv-id
                               :contig (:genomic-coordinate/contig coords)
                               :start (:genomic-coordinate/start coords)
                               :end (:genomic-coordinate/end coords)}
                              (into {} (map (fn [[k v]] [(keyword (name k)) v])) m))
                 genes? (assoc :hgnc-symbol hgnc))))
           r)
      db db-name)))

(defn cnv-gene-calls
  "Gene-level copy number calls (e.g. GISTIC2 discrete calls
  :measurement/cnv-call-score, or enum :measurement/cnv-call). Requires at
  least one of :genes, :samples, :subjects. Returns :sample-id,
  :measurement-set, :hgnc-symbol, :value.

  opts: :genes, :samples, :subjects, :measurement (default :cnv-call-score)"
  [db-or-name {:keys [genes measurement] :or {measurement :cnv-call-score} :as opts}]
  (require-cnv-subset! opts)
  (let [attr (measurement-attr measurement)
        enum? (contains? enum-measurement-attrs attr)
        [ins args where] (subset-clauses opts)
        q (vec (concat '[:find ?sample-id ?measurement-set ?hgnc-symbol ?value :with ?m :in $]
                       (when genes '[[?hgnc-symbol ...]])
                       ins
                       [:where]
                       ;; most selective subset first: genes (gene-level CNV sets are
                       ;; genes x samples), else samples/subjects
                       (if genes
                         (concat '[[?g :gene/hgnc-symbol ?hgnc-symbol]
                                   [?gp :gene-product/gene ?g]
                                   [?m :measurement/gene-product ?gp]
                                   [?m :measurement/sample ?s]]
                                 where)
                         (concat where
                                 '[[?m :measurement/sample ?s]
                                   [?m :measurement/gene-product ?gp]
                                   [?gp :gene-product/gene ?g]
                                   [?g :gene/hgnc-symbol ?hgnc-symbol]]))
                       [['?m attr (if enum? '?value-ref '?value)]]
                       '[[?s :sample/id ?sample-id]
                         [?ms :measurement-set/measurements ?m]
                         [?ms :measurement-set/name ?measurement-set]]
                       (when enum? '[[?value-ref :db/ident ?value]])))
        [r db db-name] (run db-or-name q (vec (concat (when genes [(vec genes)]) args)))]
    (provenanced (map #(update (zipmap [:sample-id :measurement-set :hgnc-symbol :value] %)
                               :value res/ident-name)
                      r)
                 db db-name)))

(comment
  (str/join ", " (map :measurement-set-name (dataset-summary "H37001-2026-09-21a"))))

;; -- entities attached to measurement sets --

(defn- ms-entities [db-or-name ref-attr pattern measurement-set]
  (pull-rows db-or-name
             {:find [(list 'pull '?e pattern)]
              :in '[$ ?ms-name]
              :where [['?ms :measurement-set/name '?ms-name] ['?ms ref-attr '?e]]}
             [measurement-set]))

(defn cell-populations
  "Cell populations of a measurement set: :cell-population-name, cell type,
  positive/negative markers, parent."
  [db-or-name measurement-set]
  (ms-entities db-or-name :measurement-set/cell-populations
               '[* {:cell-population/cell-type [:cell-type/co-name]}
                 {:cell-population/positive-markers [:epitope/id]}
                 {:cell-population/negative-markers [:epitope/id]}
                 {:cell-population/parent [:cell-population/name]}]
               measurement-set))

(defn tcrs "TCRs of a measurement set." [db-or-name measurement-set]
  (ms-entities db-or-name :measurement-set/tcrs '[*] measurement-set))

(defn otus "OTUs of a measurement set." [db-or-name measurement-set]
  (ms-entities db-or-name :measurement-set/otus '[*] measurement-set))

(defn sgbs "SGBs (metagenomic species) of a measurement set." [db-or-name measurement-set]
  (ms-entities db-or-name :measurement-set/sgbs '[*] measurement-set))
