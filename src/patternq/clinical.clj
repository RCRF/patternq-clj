(ns patternq.clinical
  "Clinical data: observation sets and observations, subject outcomes,
  adverse events, clinical interventions, and patient clinical timelines
  (queries + plotly traces, generalised from the unify-central analysis
  clinical timeline figures)."
  (:require [clojure.set]
            [patternq.dataset :as pqd]
            [patternq.db :as pdb]
            [patternq.plot :as plot]
            [patternq.results :as res]))

(set! *warn-on-reflection* true)

(defn clinical-observation-sets
  ":clinical-observation-set-name, :clinical-observation-set-description"
  [db-or-name]
  (pqd/pull-rows db-or-name '[:find (pull ?c [:clinical-observation-set/name
                                              :clinical-observation-set/description])
                              :where [?c :clinical-observation-set/name]]
                 []))

(defn clinical-observations
  "Clinical observations, either
  - :obs-type (e.g. :os, :pfs, :bor, :recist, :ldh): one map per observation
    with :subject-id, :timepoint-id, :study-day-id and a key named after the
    type (enums as names); or
  - :set-name: all attributes of the observations in that set, keys without
    the clinical-observation- prefix.

  opts: :obs-type | :set-name, :subjects"
  [db-or-name {:keys [obs-type set-name subjects]}]
  (cond
    obs-type
    (let [attr (keyword "clinical-observation" (name obs-type))
          q (vec (concat [:find '?subject-id
                          (list 'pull '?o [{:clinical-observation/timepoint [:timepoint/id]}
                                           {:clinical-observation/study-day [:study-day/id]}
                                           attr])
                          :in '$]
                         (when subjects '[[?subject-id ...]])
                         [:where ['?o attr]
                          '[?o :clinical-observation/subject ?p]
                          '[?p :subject/id ?subject-id]]))
          [r db db-name] (pqd/run db-or-name q (if subjects [(vec subjects)] []))
          k (res/attr-key attr)]
      (pdb/with-provenance
        ;; enum refs resolved once for the whole result
        (mapv (fn [sid o]
                (-> (res/flatten-pull o)
                    (clojure.set/rename-keys {k (keyword (name obs-type))})
                    (assoc :subject-id sid)))
              (map first r)
              (res/resolve-enum-refs db (mapv second r)))
        db db-name))

    set-name
    (let [rows (pqd/pull-rows
                 db-or-name
                 '[:find (pull ?o [* {:clinical-observation/subject [:subject/id]}
                                   {:clinical-observation/timepoint [:timepoint/id]}
                                   {:clinical-observation/study-day [:study-day/id]}
                                   {:clinical-observation/metastasis-gdc-anatomic-sites [:gdc-anatomic-site/name]}])
                   :in $ ?set-name
                   :where
                   [?cos :clinical-observation-set/name ?set-name]
                   [?cos :clinical-observation-set/clinical-observations ?o]]
                 [set-name])]
      (with-meta (mapv #(res/rename-prefix % "clinical-observation-") rows) (meta rows)))

    :else (throw (ex-info "Give :obs-type or :set-name" {}))))

(def ^:private recist-rank {"CR" 1 "PR" 2 "SD" 3 "PD" 4})

(defn subject-outcomes
  "One map per subject: :subject-id, :bor (from :clinical-observation/bor, or
  derived from RECIST observations when absent: CR > PR > SD > PD), :pfs,
  :pfs-event, :os, :os-event where present. Throws if a subject has more than
  one value of a single-valued outcome."
  [db-or-name]
  (let [db (pdb/as-db db-or-name)
        get1 (fn [t]
               (let [k (keyword (name t))
                     rows (filter #(some? (get % k)) (clinical-observations db {:obs-type t}))]
                 (when (not= (count rows) (count (distinct (map :subject-id rows))))
                   (throw (ex-info (str "More than one " (name t) " value for some subjects") {:outcome t})))
                 (into {} (map (juxt :subject-id #(get % k))) rows)))
        bor (let [b (get1 :bor)]
              (if (seq b)
                b
                (->> (clinical-observations db {:obs-type :recist})
                     (group-by :subject-id)
                     (into {} (map (fn [[sid rs]]
                                     (let [vals (keep (comp recist-rank :recist) rs)]
                                       [sid (if (seq vals)
                                              (key (first (filter #(= (apply min vals) (val %)) recist-rank)))
                                              "Unknown")])))))))
        outcomes {:bor bor :pfs (get1 :pfs) :pfs-event (get1 :pfs-event)
                  :os (get1 :os) :os-event (get1 :os-event)}]
    (pdb/with-provenance
      (mapv (fn [{:keys [subject-id]}]
              (reduce-kv (fn [m k lookup]
                           (if-some [v (get lookup subject-id)] (assoc m k v) m))
                         {:subject-id subject-id}
                         outcomes))
            (pqd/subjects db))
      db (when (string? db-or-name) db-or-name))))

(defn adverse-events
  "Adverse events (keys without the adverse-event- prefix).

  opts: :set-name (clinical observation set)"
  ([db-or-name] (adverse-events db-or-name {}))
  ([db-or-name {:keys [set-name]}]
   (let [pattern '[* {:adverse-event/subject [:subject/id]}
                   {:adverse-event/timepoint [:timepoint/id]}
                   {:adverse-event/meddra-adverse-event [:meddra-disease/preferred-name]}
                   {:adverse-event/ctcae-grade [:db/ident]}
                   {:adverse-event/ae-causality [:db/ident]}
                   {:adverse-event/study-day [:study-day/id]}]
         rows (if set-name
                (pqd/pull-rows db-or-name {:find [(list 'pull '?o pattern)]
                                           :in '[$ ?set-name]
                                           :where '[[?cos :clinical-observation-set/name ?set-name]
                                                    [?cos :clinical-observation-set/adverse-events ?o]]}
                               [set-name])
                (pqd/pull-rows db-or-name {:find [(list 'pull '?o pattern)]
                                           :where '[[?o :adverse-event/subject]]}
                               []))]
     (with-meta (mapv #(res/rename-prefix % "adverse-event-") rows) (meta rows)))))

(def intervention-pattern
  '[* {:clinical-intervention/subject [:subject/id]}
    {:clinical-intervention/timepoint [:timepoint/id :timepoint/relative-order :timepoint/offset]}
    {:clinical-intervention/treatment-regimen
     [:treatment-regimen/name
      {:treatment-regimen/drug-regimens [{:drug-regimen/drug [:drug/preferred-name]}
                                         :drug-regimen/freetext-drug]}]}
    {:clinical-intervention/biospecimen-derived-samples [:sample/id]}])

(defn clinical-interventions
  "Clinical interventions (treatments, surgeries, biopsies, ...), keys without
  the clinical-intervention- prefix; treatment regimen and drug names where
  present.

  opts: :subjects"
  ([db-or-name] (clinical-interventions db-or-name {}))
  ([db-or-name {:keys [subjects]}]
   (let [rows (pqd/pull-rows db-or-name
                             (vec (concat [:find (list 'pull '?ci intervention-pattern) :in '$]
                                          (when subjects '[[?subject-id ...]])
                                          '[:where [?ci :clinical-intervention/subject ?p]
                                            [?p :subject/id ?subject-id]]))
                             (if subjects [(vec subjects)] []))]
     (with-meta (mapv #(res/rename-prefix % "clinical-intervention-") rows) (meta rows)))))

;; -- clinical timeline --

(defn- timepoint-sort-key [e]
  (let [tp (or (:clinical-intervention/timepoint e) (:clinical-observation/timepoint e))
        id (:timepoint/id tp)]
    [(or (:timepoint/relative-order tp)
         (when (and id (re-matches #"-?\d+" id)) (Long/parseLong id))
         Long/MAX_VALUE)
     (str id)]))

(defn clinical-timeline
  "All clinical interventions and observations for one patient as raw pulled
  maps (namespaced keys, enums resolved), ordered by timepoint relative order
  (falling back to numeric timepoint ids)."
  [db-or-name patient-id]
  (let [db (pdb/as-db db-or-name)]
    (->> (pdb/q '[:find (pull ?ce [* {:clinical-observation/timepoint [:timepoint/id :timepoint/relative-order]
                                    :clinical-intervention/timepoint [:timepoint/id :timepoint/relative-order]
                                    :clinical-intervention/treatment-regimen
                                    [* {:treatment-regimen/drug-regimens [*]}]
                                    :clinical-intervention/biospecimen-derived-samples [:sample/id]}])
                :in $ ?patient-id
                :where
                [?p :subject/id ?patient-id]
                (or-join [?ce ?p]
                         [?ce :clinical-intervention/subject ?p]
                         [?ce :clinical-observation/subject ?p])]
              db patient-id)
         (mapv first)
         (res/resolve-enum-refs db)
         (sort-by timepoint-sort-key)
         vec)))

(defn observation-series
  "[[day value] ...] for one clinical-observation attribute (e.g.
  :clinical-observation/natera-signatera-mtm-per-ml) of a patient, with day
  taken from the observation's timepoint id (numeric day offsets) or
  study-day."
  [db-or-name patient-id attr]
  (->> (pdb/q {:find '[?tp-id ?v]
             :in '[$ ?pid]
             :where [['?p :subject/id '?pid]
                     ['?o :clinical-observation/subject '?p]
                     ['?o attr '?v]
                     ['?o :clinical-observation/timepoint '?tp]
                     ['?tp :timepoint/id '?tp-id]]}
            (pdb/as-db db-or-name) patient-id)
       (keep (fn [[tp v]] (when (re-matches #"-?\d+(\.\d+)?" tp) [(Double/parseDouble tp) v])))
       (sort-by first)
       vec))

(defn- axis-key [row] (if (= 1 row) :yaxis (keyword (str "yaxis" row))))

(defn- axis-style [th]
  {:showgrid true :gridcolor (:grid th) :zeroline false :showline true
   :linecolor (:grid th) :ticks "outside" :ticklen 4
   :tickfont {:family (:font th) :size 12 :color (:text_secondary th)}})

(defn series-trace
  "A numeric time series row (e.g. ctDNA, CEA, lab values). xy-data is
  [[day value] ...]. opts: :name, :log? (default true), :offset (added before
  log scaling, default 0.01), :color, :markers? . Returns {:traces :layout}."
  [xy-data row {:keys [name log? offset color markers?] :or {log? true offset 0.01}}]
  (let [th (plot/theme)]
    {:traces [(cond-> {:type "scatter"
                       :mode (if markers? "lines+markers" "lines")
                       :x (mapv first xy-data)
                       :y (mapv #(+ (double (second %)) (if log? offset 0.0)) xy-data)
                       :name name
                       :showlegend false
                       :xaxis "x"
                       :yaxis (str "y" row)
                       :line {:width 2}}
                color (assoc-in [:line :color] color))]
     :layout {(axis-key row) (merge (axis-style th)
                                    {:type (if log? "log" "linear")
                                     :title {:text (if (and log? (pos? offset)) (str offset " + " name) name)
                                             :font {:family (:font th) :size 12}}})}}))

(defn treatment-traces
  "Treatment spans as horizontal bars on one row. treatments: [{:label
  :day-range [start end]} ...]; colors follow the theme's categorical order."
  [treatments row]
  (let [th (plot/theme)]
    {:traces (vec (map-indexed
                    (fn [i {:keys [label day-range color]}]
                      {:type "scatter" :mode "lines"
                       :x [(first day-range) (second day-range)]
                       :y [label label]
                       :name label :showlegend false
                       :xaxis "x" :yaxis (str "y" row)
                       :line {:width 10 :color (or color (plot/series-color i))}})
                    treatments))
     :layout {(axis-key row) (merge (axis-style th) {:type "category" :showgrid false})}}))

(defn event-trace
  "Point events (samples, disease status, imaging, ...) on one row.
  events: [{:day :label :group} ...]; groups get theme colors in order."
  [events row {:keys [row-label symbol] :or {row-label "Events" symbol "diamond"}}]
  (let [th (plot/theme)
        groups (vec (distinct (map :group events)))]
    {:traces (vec (map-indexed
                    (fn [i g]
                      (let [es (filter #(= g (:group %)) events)]
                        {:type "scatter" :mode "markers+text"
                         :x (mapv :day es)
                         :y (vec (repeat (count es) row-label))
                         :text (mapv :label es)
                         :textposition "top right"
                         :xaxis "x" :yaxis (str "y" row)
                         :marker {:symbol symbol :size 12 :color (plot/series-color i)
                                  :line {:color (:surface th) :width 1.5}}
                         :name (str (or g row-label)) :showlegend false}))
                    groups))
     :layout {(axis-key row) (merge (axis-style th) {:type "category" :showgrid false})}}))

(defn timeline-figure
  "Assemble row components ({:traces :layout} from series-trace,
  treatment-traces, event-trace) into one figure with a shared day axis.
  opts: :title, :max-day, :height"
  [components {:keys [title max-day height]}]
  (let [th (plot/theme)
        n (count components)]
    {:data (vec (mapcat :traces components))
     :layout (merge (plot/base-layout title)
                    {:grid {:rows n :columns 1 :pattern "coupled"}
                     :height (or height (* 200 n))
                     :xaxis (merge (axis-style th)
                                   (cond-> {:title {:text "day"}}
                                     max-day (assoc :range [0 max-day])))}
                    (apply merge (map :layout components)))}))
