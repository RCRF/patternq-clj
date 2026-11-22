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

;; A query that binds every measurement of a set (count, sample, [?m ?a]) makes
;; the server hold the whole set in memory before any aggregate runs; on large
;; expression sets that runs the service out of heap. The functions below read
;; the set's measurement ids from the datoms index instead, in bounded pages,
;; and query only those ids.
(def ^:private page-size 10000) ; the server's datoms page cap; each page is also an audit write

(defn- measurement-set-eids [db measurement-set]
  (or (seq (pdb/q '[:find [?ms ...] :in $ ?ms-name :where [?ms :measurement-set/name ?ms-name]]
                  db measurement-set))
      (throw (ex-info (str "No measurement set named " (pr-str measurement-set))
                      {:measurement-set measurement-set :db (pdb/db-name db)}))))

(defn- set-page
  "Up to `limit` measurement ids of set `ms` from `pos`: [ids next-pos],
  next-pos nil at the end of the set. A position is [offset cursor]: the
  cursor when the service pages by cursor, else the offset. A cursor page
  answered without one (an instance from before cursor paging, during a
  deploy or a rollback) is read again by offset, never repeated."
  [db ms limit pos]
  (let [[offset after] (if (vector? pos) pos [pos nil])
        fetch (fn [opts] (pdb/datoms db :eavt [ms :measurement-set/measurements] (assoc opts :limit limit)))
        page (if after
               (let [p (fetch {:after after})]
                 (if (:patternq/cursor? (meta p)) p (fetch {:offset offset})))
               (fetch {:offset offset}))
        ids (mapv :v page)
        end (+ offset (count ids))]
    [ids (if (:patternq/cursor? (meta page))
           (when-let [nxt (:patternq/next (meta page))] [end nxt])
           (when (= (count ids) limit) [end nil]))]))

(defn- set-measurement-ids
  "Up to `limit` measurement ids of set `ms` in index order: [ids complete?],
  complete? when they reach the end of the set."
  [db ms limit]
  (loop [ids [] pos (identity 0)] ; boxed: pos becomes a cursor string
    (if (>= (count ids) limit)
      [ids (empty? (first (set-page db ms 1 pos)))]
      (let [[got pos] (set-page db ms (min page-size (- limit (count ids))) pos)
            ids (into ids got)]
        (if (nil? pos) [ids true] (recur ids pos))))))

(defn- count-set-measurements
  "[count exact?]: the set's measurement count when it is at most `cap`, else
  [cap false]. Counted in place when the service (or the peer) can;
  otherwise one page for a small set, and for a set filling the first page one
  probe past the cap, then one-datom probes binary-search the 10000-wide block
  it ends in and one page of that block finishes the count: about six round
  trips for any size, instead of a page per 10000."
  [db ms cap]
  (if-let [[n exact?] (pdb/datoms-count db :eavt [ms :measurement-set/measurements] cap)]
    (if (and exact? (<= n cap)) [n true] [cap false])
    (let [cap (min cap 100000) ; the deepest offset paging goes
          page (fn [offset limit] (count (pdb/datoms db :eavt [ms :measurement-set/measurements]
                                                     {:offset offset :limit limit})))
          step 10000
          n (page 0 step)]
      (cond
        (< n step) [n true]
        (pos? (page cap 1)) [cap false]
        :else (loop [lo 1 hi (long (Math/ceil (/ cap step)))]
                (if (> (- hi lo) 1)
                  (let [mid (quot (+ lo hi) 2)]
                    (if (pos? (page (* mid step) 1)) (recur mid hi) (recur lo mid)))
                  [(+ (* lo step) (page (* lo step) step)) true]))))))

(defn measurement-sets
  "Measurement sets with assay and size: :assay-name, :assay-technology,
  :measurement-set-name, :measurement-count, :measurement-count-exact.
  Counting never binds a whole set in one query: a set with more than
  :count-cap (default 1000000) measurements reports :count-cap with
  :measurement-count-exact false (read it as \"at least\"). The service
  counts in place; one without in-place counts is paged instead, capped at
  100000 (the deepest it pages by offset)."
  ([db-or-name] (measurement-sets db-or-name {}))
  ([db-or-name {:keys [count-cap] :or {count-cap 1000000}}]
   (let [db (pdb/as-db db-or-name)
         sets (vec (pdb/q '[:find ?ms ?name
                            :where
                            [?a :assay/measurement-sets ?ms]
                            [?ms :measurement-set/name ?name]]
                          db))
         ;; one or two round trips per set, and a dataset can have dozens of
         ;; sets: count them concurrently (back pressure still caps requests
         ;; in flight)
         found (mapv deref (mapv (fn [[ms]] (future (count-set-measurements db ms count-cap))) sets))
         counts (reduce (fn [m [[_ nm] [n exact?]]]
                          (let [[n0 e0] (get m nm [0 true])]
                            (assoc m nm [(+ n0 n) (and e0 exact?)])))
                        {}
                        (map vector sets found))]
     (provenanced (map (fn [{:keys [measurement-set-name] :as row}]
                         (let [[n exact?] (get counts measurement-set-name [0 true])]
                           (assoc row :measurement-count n :measurement-count-exact exact?)))
                       (dataset-summary db))
                  db (when (string? db-or-name) db-or-name)))))

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
  percent-of-parent alongside population x marker median-channel-value.

  Counts cover the set's first :limit (default 20000) measurements in index
  order; (:patternq/exact (meta result)) is false when the set is larger, and
  the counts are then a sample (kinds laid out after the first :limit can be
  missed)."
  ([db-or-name measurement-set] (measurement-types db-or-name measurement-set {}))
  ([db-or-name measurement-set {:keys [limit] :or {limit 20000}}]
   (let [db (pdb/as-db db-or-name)
         scans (for [ms (measurement-set-eids db measurement-set)] (set-measurement-ids db ms limit))
         exact? (every? second scans)
         counts (apply merge-with +
                       {}
                       (for [[ids] scans
                             chunk (partition-all page-size ids)]
                         (into {} (pdb/q '[:find ?attr (count ?m)
                                           :in $ [?m ...]
                                           :where [?m ?a] [?a :db/ident ?attr]]
                                         db (vec chunk)))))]
     (vary-meta
       (provenanced (->> counts
                         (remove (comp #{:measurement/id :measurement/uid :measurement/sample} key))
                         (map (fn [[a n]] {:attribute (name a)
                                           :kind (if (contains? measurement-targets a) "target" "value")
                                           :count n}))
                         (sort-by (juxt :kind (comp - :count))))
                    db (when (string? db-or-name) db-or-name))
       assoc :patternq/exact exact?))))

(defn measurement-set-attributes
  "Attributes of the set's first `:n` (default 200) measurements in index
  order, optionally only those carrying `:measurement` (looked for among the
  set's first 100000 measurements, the deepest the server pages by offset).
  Reads bounded pages of the datoms index, so it stays cheap for sets with
  millions of measurements."
  ([db-or-name measurement-set] (measurement-set-attributes db-or-name measurement-set {}))
  ([db-or-name measurement-set {:keys [measurement n] :or {n 200}}]
   (let [db (pdb/as-db db-or-name)
         attr (when measurement (measurement-attr measurement))
         carrying (fn [ids] (if (and attr (seq ids))
                              (vec (pdb/q {:find '[[?m ...]] :in '[$ [?m ...]] :where [['?m attr]]} db ids))
                              ids))
         eids (->> (measurement-set-eids db measurement-set)
                   (mapcat (fn [ms]
                             (loop [found [] pos (identity 0) pages 0] ; boxed: pos becomes a cursor string
                               (let [[ids pos] (set-page db ms (if attr page-size n) pos)
                                     found (into found (carrying ids))]
                                 ;; 100000 measurements at most, the service's offset cap
                                 (if (or (>= (count found) n) (nil? pos) (nil? attr) (>= (inc pages) 10))
                                   found
                                   (recur found pos (inc pages)))))))
                   (take n)
                   vec)]
     (->> (pdb/pull-many db '[*] eids)
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

(defn- sample-assays-by-walk
  "sample-assays for a database from before :measurement-set/includes-samples:
  walks each set's measurement -> sample refs."
  [db smp]
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
           ms)))

(def sample-assays-query
  "Samples x measurement sets from :measurement-set/includes-samples, the
  enrichment that lists each set's samples (schema 0.3.18), so nothing walks
  the measurements."
  '[:find ?sample-id ?assay-name ?assay-technology ?measurement-set-name
    :where
    [?a :assay/name ?assay-name]
    [?a :assay/technology ?t]
    [?t :db/ident ?assay-technology]
    [?a :assay/measurement-sets ?ms]
    [?ms :measurement-set/name ?measurement-set-name]
    [?ms :measurement-set/includes-samples ?s]
    [?s :sample/id ?sample-id]])

(defn- includes-samples? [db]
  (some? (pdb/q '[:find ?a . :in $ ?ident :where [?a :db/ident ?ident]]
                db :measurement-set/includes-samples)))

(defn sample-assays
  "Which samples were measured in which measurement sets: :subject-id,
  :sample-id, :timepoint-id, :assay-name, :assay-technology,
  :measurement-set-name. Reads :measurement-set/includes-samples; a database
  from before that enrichment falls back to walking each set's measurements,
  which can time out on large expression sets."
  [db-or-name]
  (let [db (pdb/as-db db-or-name)
        smp (into {} (map (juxt :sample-id identity)) (samples db))]
    (provenanced
      (if (includes-samples? db)
        (for [[sid a t ms] (pdb/q sample-assays-query db)]
          (merge (select-keys (smp sid) [:subject-id :timepoint-id])
                 {:sample-id sid
                  :assay-name a
                  :assay-technology (res/ident-name t)
                  :measurement-set-name ms}))
        (sample-assays-by-walk db smp))
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
        ;; Datomic joins in clause order: start from the genes or the samples
        ;; when given, not from every variant measurement in the database
        where (cond
                genes '[[?g :gene/hgnc-symbol ?gene]
                        [?v :variant/gene ?g]
                        [?m :measurement/variant ?v]
                        [?m :measurement/vaf ?vaf]
                        [?m :measurement/sample ?s]
                        [?s :sample/id ?sample-id]
                        [?ms :measurement-set/measurements ?m]
                        [?ms :measurement-set/name ?measurement-set]]
                samples '[[?s :sample/id ?sample-id]
                          [?m :measurement/sample ?s]
                          [?m :measurement/vaf ?vaf]
                          [?m :measurement/variant ?v]
                          [?ms :measurement-set/measurements ?m]
                          [?ms :measurement-set/name ?measurement-set]]
                :else '[[?m :measurement/vaf ?vaf]
                        [?m :measurement/variant ?v]
                        [?m :measurement/sample ?s]
                        [?s :sample/id ?sample-id]
                        [?ms :measurement-set/measurements ?m]
                        [?ms :measurement-set/name ?measurement-set]])]
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
                          [:where]
                          ;; samples without genes: start from the samples, not every gene
                          (if (and samples (not genes))
                            ['[?s :sample/id ?sample-id]
                             '[?m :measurement/sample ?s]
                             ['?m attr '?value]
                             '[?m :measurement/gene-product ?gp]
                             '[?gp :gene-product/gene ?g]
                             '[?g :gene/hgnc-symbol ?hgnc-symbol]]
                            ['[?g :gene/hgnc-symbol ?hgnc-symbol]
                             '[?gp :gene-product/gene ?g]
                             '[?m :measurement/gene-product ?gp]
                             ['?m attr '?value]
                             '[?m :measurement/sample ?s]
                             '[?s :sample/id ?sample-id]])
                          '[[?ms :measurement-set/measurements ?m]
                            [?ms :measurement-set/name ?measurement-set]]))
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
  ?s (the sample entity), samples first."
  [{:keys [samples subjects]}]
  [(cond-> [] samples (conj '[?sample-id ...]) subjects (conj '[?subject-id ...]))
   (cond-> [] samples (conj (vec samples)) subjects (conj (vec subjects)))
   (vec (concat (when samples '[[?s :sample/id ?sample-id]])
                (when subjects '[[?p :subject/id ?subject-id] [?s :sample/subject ?p]])))])

;; Starting from a few samples walks only their measurements. A subject can
;; own most of a dataset (a single-participant deep dive has one subject and
;; every sample), so subject-only CNV queries start from the measurement
;; attribute instead.

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
                         ;; most selective subset first: genes, else samples, else the
                         ;; CNV measurements themselves (subject-only)
                         (cond
                           genes (concat gene-where '[[?m :measurement/cnv ?c] [?m :measurement/sample ?s]] where)
                           (:samples opts) (concat where '[[?m :measurement/sample ?s] [?m :measurement/cnv ?c]])
                           :else (concat '[[?m :measurement/cnv ?c] [?m :measurement/sample ?s]] where))
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
                       (cond
                         genes (concat '[[?g :gene/hgnc-symbol ?hgnc-symbol]
                                         [?gp :gene-product/gene ?g]
                                         [?m :measurement/gene-product ?gp]
                                         [?m :measurement/sample ?s]]
                                       where)
                         (:samples opts) (concat where
                                                 '[[?m :measurement/sample ?s]
                                                   [?m :measurement/gene-product ?gp]
                                                   [?gp :gene-product/gene ?g]
                                                   [?g :gene/hgnc-symbol ?hgnc-symbol]])
                         :else (concat [['?m attr (if enum? '?value-ref '?value)]]
                                       '[[?m :measurement/sample ?s]]
                                       where
                                       '[[?m :measurement/gene-product ?gp]
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
