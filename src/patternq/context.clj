(ns patternq.context
  "Joining context (samples, subjects, outcomes, variants, CNVs) onto
  measurement rows, splitting by measurement set, and microbiome taxonomy
  helpers. Ported from wick (via the R library)."
  (:require [clojure.string :as str]
            [patternq.clinical :as clinical]
            [patternq.dataset :as pqd]
            [patternq.reference :as ref]
            [patternq.results :as res]))

(set! *warn-on-reflection* true)

(defn- join-by
  "Left join `rows` with `lookup` rows on key k; lookup values don't
  overwrite keys already present in the row."
  [rows k lookup-rows]
  (let [idx (into {} (map (juxt k identity)) lookup-rows)]
    (mapv (fn [r] (merge (get idx (get r k)) r)) rows)))

(defn add-subject-context
  "Join subject attributes (and optionally outcomes) on :subject-id.
  opts: :outcomes?"
  ([rows db-or-name] (add-subject-context rows db-or-name {}))
  ([rows db-or-name {:keys [outcomes?]}]
   (cond-> (join-by rows :subject-id (pqd/subjects db-or-name))
     outcomes? (join-by :subject-id (clinical/subject-outcomes db-or-name)))))

(defn add-sample-context
  "Join sample attributes on :sample-id, then (by default) subject attributes.
  opts: :subjects? (default true), :outcomes?"
  ([rows db-or-name] (add-sample-context rows db-or-name {}))
  ([rows db-or-name {:keys [subjects? outcomes?] :or {subjects? true}}]
   (cond-> (join-by rows :sample-id (pqd/samples db-or-name))
     subjects? (add-subject-context db-or-name {:outcomes? outcomes?}))))

(defn add-variant-context
  "Join variant annotations on :variant-id."
  [rows db-or-name]
  (join-by rows :variant-id (ref/variant-annotations db-or-name {:variant-ids (distinct (keep :variant-id rows))})))

(defn add-cnv-context
  "Join CNV entities (coordinates, genes) on :cnv-id."
  [rows db-or-name]
  (join-by rows :cnv-id (ref/cnvs db-or-name)))

(defn split-by-measurement-set
  "Group long rows by :measurement-set; with a `col` (target key or keys)
  each group becomes a matrix (patternq.results/->matrix)."
  ([rows] (group-by :measurement-set rows))
  ([rows col] (update-vals (group-by :measurement-set rows) #(res/->matrix % col))))

;; -- taxonomy (OTUs / SGBs) --

(def ^:private levels ["kingdom" "phylum" "class" "order" "family" "genus" "species"])

(defn- level-keys [taxa]
  (let [prefix (if (some #(contains? % :otu-kingdom) taxa) "otu-" "sgb-")]
    (filterv (fn [k] (some #(contains? % k) taxa)) (map #(keyword (str prefix %)) levels))))

(defn deduplicate-taxonomy
  "Make taxon names unique at each level: missing levels become `na-value`,
  and a name occurring under different parents is prefixed with its
  ancestors (joined by \"_\")."
  ([taxa] (deduplicate-taxonomy taxa "Unclassified"))
  ([taxa na-value]
   (let [lks (level-keys taxa)
         filled (mapv (fn [t] (reduce #(update %1 %2 (fn [v] (or v na-value))) t lks)) taxa)]
     (reduce
       (fn [out i]
         (let [cur (lks i)
               anc (subvec lks 0 (inc i))
               uniq (distinct (map #(select-keys % anc) filled))
               dup (set (map key (filter #(> (val %) 1) (frequencies (map cur uniq)))))]
           (mapv (fn [o f] (if (dup (get f cur)) (assoc o cur (str/join "_" (map f anc))) o)) out filled)))
       filled
       (range (count lks))))))

(defn aggregate-taxa
  "Sum measurements (rows with :sample-id, `id-key`, :value) to each taxonomic
  level. Returns {level-key [{:sample-id :taxon :value} ...]}, values as
  per-sample proportions when :normalize? (default true)."
  [rows taxa id-key {:keys [normalize? na-value] :or {normalize? true na-value "Unclassified"}}]
  (let [taxa (deduplicate-taxonomy taxa na-value)
        by-id (into {} (map (juxt id-key identity)) taxa)
        rows (keep #(when-let [t (by-id (get % id-key))] (merge t %)) (filter (comp some? :value) rows))]
    (into {}
          (for [lk (level-keys taxa)]
            (let [sums (->> rows (group-by (juxt :sample-id lk))
                            (map (fn [[[s t] rs]] {:sample-id s :taxon t :value (reduce + (map :value rs))})))
                  totals (update-vals (group-by :sample-id sums) #(reduce + (map :value %)))]
              [lk (vec (if normalize?
                         (map #(update % :value (fn [v] (/ (double v) (double (totals (:sample-id %)))))) sums)
                         sums))])))))
