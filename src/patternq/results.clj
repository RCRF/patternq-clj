(ns patternq.results
  "Shaping query results: pull flattening, enum idents, matrices.

  Canned queries return vectors of maps with unqualified kebab-case keys.
  These line up one-to-one with the snake_case column names of the R and
  Python libraries (:sample-id <-> sample_id)."
  (:require [clojure.string :as str]
            [clojure.walk :as walk]))

(set! *warn-on-reflection* true)

(defn ident-name
  ":variant.impact/high -> \"high\"; other values unchanged."
  [x]
  (if (keyword? x) (name x) x))

(defn attr-key
  ":subject/id -> :subject-id"
  [k]
  (keyword (str/replace (str (when (namespace k) (str (namespace k) "-")) (name k))
                        #"[./]" "-")))

(defn- ident-map? [v]
  (and (map? v) (= [:db/ident] (keys v))))

(defn- eid-map? [v]
  (and (map? v) (= [:db/id] (keys v))))

(defn resolve-enum-refs
  "Pulled refs without a nested pattern come back as {:db/id n}. Replace
  those that are enums (have a :db/ident) with {:db/ident ...}, so enum values
  read as names whichever pull pattern was used (e.g. `*`)."
  [db x]
  (let [eids (atom #{})
        _ (walk/postwalk (fn [v] (when (eid-map? v) (swap! eids conj (:db/id v))) v) x)
        idents (if (seq @eids)
                 (into {} ((requiring-resolve 'patternq.db/q)
                           '[:find ?e ?ident :in $ [?e ...] :where [?e :db/ident ?ident]]
                           db (vec @eids)))
                 {})]
    (walk/postwalk
      (fn [v]
        (if (eid-map? v)
          (if-let [ident (idents (:db/id v))] {:db/ident ident} v)
          v))
      x)))

(defn flatten-pull
  "Flatten one pulled entity into a map of unqualified keys.

  Scalar attributes keep their value under the attribute's key
  (:subject/id -> :subject-id). Enum refs ({:db/ident ...}) become the enum
  name under the referring attribute's key. Nested single refs are flattened
  recursively (colliding keys are prefixed with the referring attribute).
  Cardinality-many enums become vectors of names, cardinality-many refs to
  maps with a single attribute become vectors of that attribute's values,
  other cardinality-many refs become vectors of flattened maps. Unresolved
  refs ({:db/id n}) are dropped unless :exclude-ids? is false."
  ([m] (flatten-pull m {}))
  ([m {:keys [exclude-ids?] :or {exclude-ids? true} :as opts}]
   (reduce-kv
     (fn [out k v]
       (cond
         (and exclude-ids? (or (= :db/id k) (= "uid" (name k)))) out
         (nil? v) out
         (ident-map? v) (assoc out (attr-key k) (ident-name (:db/ident v)))
         (eid-map? v) (if exclude-ids? out (assoc out (attr-key k) (:db/id v)))
         (map? v) (reduce-kv (fn [o nk nv]
                               (assoc o (if (contains? o nk)
                                          (keyword (str (name (attr-key k)) "-" (name nk)))
                                          nk)
                                      nv))
                             out
                             (flatten-pull v opts))
         (and (coll? v) (every? ident-map? v)) (assoc out (attr-key k) (mapv (comp ident-name :db/ident) v))
         (and (coll? v) (every? map? v))
         (let [flat (mapv #(flatten-pull % opts) v)]
           (assoc out (attr-key k)
                  (if (every? #(= 1 (count %)) flat)
                    (mapv (comp val first) flat)
                    flat)))
         :else (assoc out (attr-key k) v)))
     {}
     m)))

(defn join-many
  "Collapse a cardinality-many value to a \"; \"-joined string."
  [v]
  (if (coll? v)
    (when (seq v) (str/join "; " (sort (map str v))))
    v))

(defn rows->maps
  "Zip relation rows with column keys."
  [ks rows]
  (mapv #(zipmap ks %) rows))

(defn rename-prefix
  "Strip `prefix` (a string like \"clinical-observation-\") from map keys."
  [m prefix]
  (into {} (map (fn [[k v]]
                  (let [n (name k)]
                    [(if (str/starts-with? n prefix) (keyword (subs n (count prefix))) k) v])))
        m))

;; -- matrices --
;;
;; A matrix is {:row-names [...] :col-names [...] :values [[row...] ...]}
;; with nil for missing cells (R: matrix; Python: DataFrame).

(defn ->matrix
  "Long rows (maps) -> matrix. `col` is a key or a vector of keys (joined with
  \"|\"). Duplicate cells are combined with `agg` (default mean).

  opts: :row (default :sample-id), :value (default :value), :agg"
  ([rows col] (->matrix rows col {}))
  ([rows col {:keys [row value agg] :or {row :sample-id value :value}}]
   (let [agg (or agg (fn [xs] (/ (double (reduce + xs)) (count xs))))
         col-fn (if (sequential? col)
                  (fn [r] (str/join "|" (map #(get r %) col)))
                  col)
         cells (->> rows
                    (group-by (juxt #(get % row) col-fn))
                    (reduce-kv (fn [m k rs] (assoc m k (agg (mapv #(get % value) rs)))) {}))
         rn (vec (sort (distinct (map first (keys cells)))))
         cn (vec (sort (distinct (map second (keys cells)))))]
     {:row-names rn
      :col-names cn
      :values (mapv (fn [r] (mapv (fn [c] (get cells [r c])) cn)) rn)})))

(defn ->long
  "Matrix -> long rows (maps), without nil cells."
  ([m] (->long m {}))
  ([{:keys [row-names col-names values]} {:keys [row col value]
                                          :or {row :sample-id col :target value :value}}]
   (vec (for [[r vs] (map vector row-names values)
              [c v] (map vector col-names vs)
              :when (some? v)]
          {row r col c value v}))))

(defn transpose-matrix
  [{:keys [row-names col-names values]}]
  {:row-names col-names
   :col-names row-names
   :values (if (seq values) (apply mapv vector values) [])})

(defn select-targets
  "Keep (:include) or drop (:exclude) matrix columns by name."
  [{:keys [col-names values] :as m} {:keys [include exclude]}]
  (let [keep? (fn [c] (and (or (nil? include) (contains? (set include) c))
                           (not (contains? (set exclude) c))))
        idx (keep-indexed (fn [i c] (when (keep? c) i)) col-names)]
    (assoc m
      :col-names (mapv #(nth col-names %) idx)
      :values (mapv (fn [row] (mapv #(nth row %) idx)) values))))
