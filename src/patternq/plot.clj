(ns patternq.plot
  "plotly figure specs as plain data ({:data [...] :layout {...}}).

  Render with Kindly/Clay (kind/plotly), write to JSON for plotly.js, or
  embed in any report. No plotting dependency.

  Theme: resources/patternq/plotly-theme.json, shared with the R and Python
  libraries. The categorical order is colorblind-validated: colors follow
  entities in fixed order and are never cycled; series past 8 fold into
  \"Other\"."
  (:require [charred.api :as json]
            [clojure.java.io :as io]
            [patternq.quant :as quant]
            [patternq.results :as res]
            [patternq.survival :as surv]))

(set! *warn-on-reflection* true)

(def theme
  (memoize
    (fn []
      (json/read-json (slurp (io/resource "patternq/plotly-theme.json")) :key-fn keyword))))

(def other-color "#8a8983")

(defn series-color
  "Categorical color for series index i (0-based); neutral gray past 8."
  [i]
  (let [cs (:categorical (theme))]
    (if (< (long i) (count cs)) (nth cs i) other-color)))

(defn fold-other
  "Map values to at most `max-n` categories (most frequent kept, rest
  \"Other\"). Returns a fn value -> category."
  ([values] (fold-other values 8))
  ([values max-n]
   (let [freqs (sort-by (comp - val) (frequencies values))]
     (if (<= (count freqs) max-n)
       identity
       (let [keep (set (map key (take (dec max-n) freqs)))]
         (fn [v] (if (keep v) v "Other")))))))

(defn format-p
  "p-value formatted like R's format.pval(p, digits = 2)."
  [p]
  (let [p (double p)]
    (cond (< p 2.2e-16) "< 2.2e-16"
          (< p 1e-4) (format "%.1e" p)
          :else (let [bd (.round (java.math.BigDecimal. p) (java.math.MathContext. 2))]
                  (.toPlainString (.stripTrailingZeros bd))))))

(defn- alpha [^String hex a]
  (let [n (Long/parseLong (subs hex 1) 16)]
    (format "rgba(%d,%d,%d,%.2f)" (bit-shift-right n 16) (bit-and (bit-shift-right n 8) 255) (bit-and n 255) (double a))))

(defn sequential-scale []
  (let [s (:sequential (theme)) n (dec (count s))]
    (vec (map-indexed (fn [i c] [(/ (double i) n) c]) s))))

(defn diverging-scale []
  (let [{:keys [low mid high]} (:diverging (theme))]
    [[0 low] [0.5 mid] [1 high]]))

(defn- axis [m]
  (let [th (theme)]
    (merge {:gridcolor (:grid th) :zerolinecolor (:grid th) :linecolor (:grid th)
            :tickfont {:color (:text_secondary th)}}
           m)))

(defn base-layout
  "Theme layout with a title; merge figure-specific keys over it."
  [title]
  (let [th (theme)]
    {:title {:text title :x 0 :xanchor "left"}
     :font {:family (:font th) :color (:text_primary th)}
     :colorway (:categorical th)
     :paper_bgcolor (:surface th)
     :plot_bgcolor (:surface th)
     :hoverlabel {:font {:family (:font th)}}}))

(defn- layout [title m]
  (let [m (cond-> m
            (:xaxis m) (update :xaxis axis)
            (:yaxis m) (update :yaxis axis))]
    (merge (base-layout title) m)))

;; -- catalog plots (R/Python parity) --

(defn vaf-histogram
  "Overlaid VAF histograms per sample, from patternq.dataset/variants output.
  opts: :samples (subset), :bin-size (default 0.05), :title"
  ([variants] (vaf-histogram variants {}))
  ([variants {:keys [samples bin-size title] :or {bin-size 0.05 title "VAF histogram"}}]
   (let [by-sample (group-by :sample-id (if samples (filter (comp (set samples) :sample-id) variants) variants))
         ids (sort (keys by-sample))]
     {:data (vec (map-indexed
                   (fn [i sid]
                     {:type "histogram" :x (mapv :vaf (by-sample sid)) :name sid :opacity 0.7
                      :xbins {:start 0 :end 1 :size bin-size}
                      :marker {:color (series-color i) :line {:color (:surface (theme)) :width 1}}})
                   ids))
      :layout (layout title {:barmode "overlay" :showlegend (> (count ids) 1)
                             :xaxis {:title {:text "VAF"} :range [0 1]}
                             :yaxis {:title {:text "variants"}}})})))

(defn gene-expression
  "Grouped bars of expression per gene and sample, from
  patternq.dataset/gene-expression output. opts: :title, :ylab, :log?"
  ([expr] (gene-expression expr {}))
  ([expr {:keys [title ylab log?] :or {title "Gene expression" ylab "value"}}]
   (let [agg (->> expr
                  (group-by (juxt :sample-id :hgnc-symbol))
                  (map (fn [[[sid gene] rows]] {:sample-id sid :hgnc-symbol gene :value (reduce + (map :value rows))})))
         by-sample (group-by :sample-id agg)
         ids (sort (keys by-sample))]
     {:data (vec (map-indexed
                   (fn [i sid]
                     (let [rows (sort-by :hgnc-symbol (by-sample sid))]
                       {:type "bar" :x (mapv :hgnc-symbol rows) :y (mapv :value rows) :name sid
                        :marker {:color (series-color i)}}))
                   ids))
      :layout (layout title {:barmode "group" :bargap 0.2 :bargroupgap 0.05 :showlegend (> (count ids) 1)
                             :yaxis {:title {:text ylab} :type (if log? "log" "linear")}
                             :xaxis {:title {:text ""}}})})))

(defn sample-overview
  "Heatmap of sample counts per subject (rows) and measurement set
  (columns), from patternq.dataset/sample-assays output. Empty cells blank."
  ([sample-assays] (sample-overview sample-assays {}))
  ([sample-assays {:keys [title] :or {title "Samples per subject and measurement set"}}]
   (let [m (res/->matrix (map #(assoc % :n 1) (distinct (map #(select-keys % [:subject-id :sample-id :measurement-set-name])
                                                            sample-assays)))
                         :measurement-set-name {:row :subject-id :value :n :agg #(reduce + %)})]
     {:data [{:type "heatmap" :x (:col-names m) :y (:row-names m) :z (:values m) :zmin 0
              :colorscale (sequential-scale) :xgap 1 :ygap 1 :colorbar {:title {:text "samples"}}
              :hovertemplate "%{y}<br>%{x}<br>%{z} samples<extra></extra>"}]
      :layout (layout title {:xaxis {:title {:text ""} :tickangle -30 :automargin true}
                             :yaxis {:title {:text "subject"} :autorange "reversed"}})})))

(defn- box [x y name color]
  {:type "box" :x x :y y :name name :marker {:color color} :line {:color color}
   :boxpoints "all" :jitter 0.3 :pointpos 0})

(defn by-timepoint
  "Box plot of values per timepoint, optionally split by a grouping key.
  Rows as from patternq.context/add-sample-context.
  opts:
    :value      key holding the value (default :value; e.g. :change)
    :timepoints timepoint ids to show, in order (default: all, ordered by
                :timepoint-relative-order)
    :group      grouping key (e.g. :status-1y)
    :levels     group order
    :lines?     per-subject trajectories (colored by group, drawn under
                translucent boxes when grouped; needs :subject-id)
    :title :ylab"
  ([rows] (by-timepoint rows {}))
  ([rows {:keys [value timepoints group levels lines? title ylab] :or {value :value ylab "value"}}]
   (let [rows (->> rows
                   (map #(assoc % ::v (get % value)))
                   (filter #(and (some? (::v %)) (some? (:timepoint-id %))))
                   (filter #(or (nil? timepoints) (some #{(:timepoint-id %)} timepoints))))
         present (set (map :timepoint-id rows))
         order (if timepoints
                 (filterv present timepoints)
                 (->> rows (sort-by #(or (:timepoint-relative-order %) 0)) (map :timepoint-id) distinct vec))
         rank (zipmap order (range))
         trajectory (fn [rs color g]
                      (let [rs (sort-by (comp rank :timepoint-id) rs)]
                        (when (> (count rs) 1)
                          (cond-> {:type "scatter" :mode "lines" :x (mapv :timepoint-id rs) :y (mapv ::v rs)
                                   :line {:color color :width 1} :showlegend false :hoverinfo "text"
                                   :text (str (:subject-id (first rs)) (when g (str "<br>" g)))}
                            g (assoc :legendgroup g)))))
         traces
         (if group
           (let [fold (fold-other (keep #(some-> (get % group) str) rows))
                 rows (keep #(when-some [g (get % group)] (assoc % ::g (fold (str g)))) rows)
                 present-g (set (map ::g rows))
                 lv (if levels (filterv present-g levels) (vec (sort present-g)))]
             (vec (concat
                    (when lines?
                      (for [[i g] (map-indexed vector lv)
                            [_ rs] (group-by :subject-id (filter #(= g (::g %)) rows))
                            :let [t (trajectory rs (alpha (series-color i) 0.35) g)]
                            :when t]
                        t))
                    (for [[i g] (map-indexed vector lv)
                          :let [rs (filter #(= g (::g %)) rows) c (series-color i)]]
                      {:type "box" :x (mapv :timepoint-id rs) :y (mapv ::v rs) :name g :legendgroup g
                       :marker {:color c :size 4} :line {:color c} :fillcolor (alpha c 0.15)
                       :boxpoints (if lines? false "all") :jitter 0.3 :pointpos 0}))))
           (vec (concat
                  (when lines?
                    (keep (fn [[_ rs]] (trajectory rs "rgba(82,81,78,0.25)" nil)) (group-by :subject-id rows)))
                  [(assoc (box (mapv :timepoint-id rows) (mapv ::v rows) ylab (series-color 0)) :showlegend false)])))]
     {:data traces
      :layout (layout title (cond-> {:xaxis {:title {:text "timepoint"} :type "category"
                                             :categoryorder "array" :categoryarray order}
                                     :yaxis {:title {:text ylab}}}
                              group (assoc :boxmode (if lines? "overlay" "group"))))})))

(defn by-group
  "Box (or violin) plot of :value per group key. opts: :violin?, :title, :ylab"
  ([rows group] (by-group rows group {}))
  ([rows group {:keys [violin? title ylab] :or {ylab "value"}}]
   (let [rows (filter #(and (some? (:value %)) (some? (get % group))) rows)
         fold (fold-other (map #(str (get % group)) rows))
         gs (sort (distinct (map #(fold (str (get % group))) rows)))]
     {:data (vec (map-indexed
                   (fn [i g]
                     (let [ys (mapv :value (filter #(= g (fold (str (get % group)))) rows))
                           c (series-color i)]
                       (if violin?
                         {:type "violin" :x (vec (repeat (count ys) g)) :y ys :name g
                          :line {:color c :width 1.5} :fillcolor (alpha c 0.25)
                          :box {:visible true :fillcolor (:surface (theme)) :line {:color c} :width 0.15}
                          :meanline {:visible false} :points "all" :jitter 0.4 :pointpos 0
                          :marker {:color c :size 5 :opacity 0.8}}
                         (box (vec (repeat (count ys) g)) ys g c))))
                   gs))
      :layout (layout title {:showlegend false :xaxis {:title {:text (name group)}} :yaxis {:title {:text ylab}}})})))

(defn survival
  "Kaplan-Meier curves (hand-rolled, patternq.quant/kaplan-meier) with
  censoring ticks and, for two or more groups, a log-rank p annotation
  (patternq.survival/logrank-test). Rows as from
  patternq.clinical/subject-outcomes or patternq.survival/survival-by-median.
  opts: :time (default :os), :event (default :os-event), :group, :levels
  (group order), :pvalue? (default true), :title, :xlab"
  ([rows] (survival rows {}))
  ([rows {:keys [time event group levels pvalue? title xlab]
          :or {time :os event :os-event title "Survival" pvalue? true}}]
   (let [rows (filter #(and (some? (get % time)) (some? (get % event))) rows)
         fold (if group (fold-other (keep #(some-> (get % group) str) rows)) identity)
         rows (if group (keep #(when-some [g (get % group)] (assoc % ::g (fold (str g)))) rows)
                  (map #(assoc % ::g "all") rows))
         present (set (map ::g rows))
         gs (if levels (filterv present levels) (vec (sort present)))
         lr (when (and pvalue? (> (count gs) 1))
              (surv/logrank-test rows {:time time :event event :group ::g}))]
     {:data (vec (apply concat
                        (map-indexed
                          (fn [i g]
                            (let [rs (filter #(= g (::g %)) rows)
                                  km (quant/kaplan-meier (map (fn [r] {:time (double (get r time))
                                                                       :event (let [e (get r event)] (if (number? e) (pos? (double e)) (boolean e)))})
                                                              rs))
                                  nm (format "%s (n=%d)" g (count rs))
                                  c (series-color i)
                                  cens (filter #(pos? (long (:n-censor %))) km)]
                              (cond-> [{:type "scatter" :mode "lines" :name nm
                                        :x (into [0] (map :time km)) :y (into [1.0] (map :surv km))
                                        :line {:shape "hv" :color c :width 2}
                                        :hovertemplate (str nm "<br>t=%{x:.1f}<br>S=%{y:.2f}<extra></extra>")}]
                                (seq cens) (conj {:type "scatter" :mode "markers" :showlegend false :hoverinfo "skip"
                                                  :x (mapv :time cens) :y (mapv :surv cens)
                                                  :marker {:symbol "line-ns-open" :size 9 :color c}}))))
                          gs)))
      :layout (layout title (cond-> {:showlegend (> (count gs) 1)
                                     :xaxis {:title {:text (or xlab (name time))} :rangemode "tozero"}
                                     :yaxis {:title {:text "survival probability"} :range [0 1.02]}}
                              (:p lr) (assoc :annotations
                                             [{:xref "paper" :yref "paper" :x 0.02 :y 0.04 :xanchor "left" :showarrow false
                                               :text (str "log-rank p = " (format-p (:p lr)))
                                               :font {:color (:text_secondary (theme))}}])))})))

(defn heatmap
  "Clustered heatmap of a matrix ({:row-names :col-names :values}, see
  patternq.results/->matrix). opts: :scale (:none | :row | :column; scaled
  values use the diverging scale), :cluster-rows? :cluster-cols? (default
  true), :col-groups {col-name group} (drawn as an annotation strip above
  the heatmap, legend as colored labels), :title, :zlab"
  ([m] (heatmap m {}))
  ([{:keys [row-names col-names values]} {:keys [scale cluster-rows? cluster-cols? col-groups title zlab]
                                          :or {scale :none cluster-rows? true cluster-cols? true}}]
   (when (or (empty? row-names) (empty? col-names))
     (throw (ex-info "heatmap: empty matrix" {})))
   (let [values (case scale
                  :row (quant/scale-rows values)
                  :column (apply mapv vector (quant/scale-rows (apply mapv vector values)))
                  values)
         ro (if cluster-rows? (quant/hclust-order values) (range (count row-names)))
         co (if cluster-cols? (quant/hclust-order (apply mapv vector values)) (range (count col-names)))
         div? (not= scale :none)
         z (mapv (fn [r] (mapv #(get-in values [r %]) co)) ro)
         xs (mapv col-names co)
         lim (when div? (apply max 1e-9 (map #(Math/abs (double %)) (remove nil? (flatten z)))))
         main (cond-> {:type "heatmap" :x xs :y (mapv row-names ro) :z z
                       :colorscale (if div? (diverging-scale) (sequential-scale))
                       :colorbar {:title {:text (or zlab (if div? "z-score" "value"))}}
                       :hovertemplate "%{y}<br>%{x}<br>%{z:.3g}<extra></extra>"}
                div? (assoc :zmin (- (double lim)) :zmax lim))]
     (if-not col-groups
       {:data [main]
        :layout (layout title {:xaxis {:title {:text ""} :tickangle -45 :type "category" :automargin true}
                               :yaxis {:title {:text ""} :type "category" :autorange "reversed" :automargin true}})}
       (let [grp (mapv #(get col-groups %) xs)
             lv (vec (sort (distinct (remove nil? grp))))
             idx (zipmap lv (range 1 (inc (count lv))))
             n (max 1 (count lv))
             strip {:type "heatmap" :x xs :y ["group"] :z [(mapv idx grp)] :text [(mapv #(or % "") grp)]
                    :showscale false :xgap 1 :zmin 1 :zmax n :xaxis "x" :yaxis "y"
                    :colorscale (vec (map-indexed (fn [i _] [(if (= 1 n) 0 (/ (double i) (dec n))) (series-color i)]) lv))
                    :hovertemplate "%{x}<br>%{text}<extra></extra>"}]
         {:data [strip (assoc main :xaxis "x" :yaxis "y2")]
          :layout (layout title {:margin {:t (+ 60 (* 18 (count lv)))}
                                 :xaxis {:title {:text ""} :tickangle -45 :type "category" :anchor "y2"}
                                 :yaxis (axis {:domain [0.955 1.0] :showticklabels false :showgrid false})
                                 :yaxis2 (axis {:domain [0 0.95] :type "category" :autorange "reversed" :automargin true})
                                 :annotations (vec (map-indexed
                                                     (fn [i g] {:xref "paper" :yref "paper" :x 1
                                                                :y (+ 1.02 (* 0.045 (- (count lv) i 1)))
                                                                :xanchor "right" :yanchor "bottom" :showarrow false
                                                                :text (str "\u25A0 " g)
                                                                :font {:color (series-color i) :size 12}})
                                                     lv))})})))))

(def ^:private impact-levels ["modifier" "low" "moderate" "high"])

(defn mutation-landscape
  "Genes (most frequently mutated first) x samples; cells show the most
  severe variant impact, or simply mutated when impact is not annotated.
  From patternq.dataset/variants output. opts: :n-genes (25), :genes, :title"
  ([variants] (mutation-landscape variants {}))
  ([variants {:keys [n-genes genes title] :or {n-genes 25 title "Mutation landscape"}}]
   (let [vs (filter :hgnc-symbol variants)
         impact? (some :impact vs)
         sev (fn [v] (if impact? (inc (max 0 (.indexOf ^java.util.List impact-levels (or (:impact v) "modifier")))) 1))
         freq (into {} (map (fn [[g rs]] [g (count (distinct (map :sample-id rs)))])) (group-by :hgnc-symbol vs))
         genes (or genes (take n-genes (map key (sort-by (comp - val) freq))))
         m (res/->matrix (map #(assoc % :severity (sev %)) (filter (comp (set genes) :hgnc-symbol) vs))
                         :sample-id {:row :hgnc-symbol :value :severity :agg #(apply max %)})
         gidx (into {} (map-indexed (fn [i g] [g i]) (:row-names m)))
         genes (filter gidx genes)
         rows (mapv #(get-in m [:values (gidx %)]) genes)
         cols (range (count (:col-names m)))
         so (sort-by (fn [c] (mapv #(if (nil? (get-in rows [% c])) 1 0) (range (count rows)))) cols)
         z (mapv (fn [row] (mapv #(get row %) so)) rows)
         s (:sequential (theme))
         labels (mapv #(format "%s (%d)" % (freq %)) genes)]
     {:data [(if impact?
               {:type "heatmap" :x (mapv (:col-names m) so) :y labels :z z :zmin 1 :zmax 4 :xgap 1 :ygap 1
                :text (mapv (fn [row] (mapv #(when % (impact-levels (dec (long %)))) row)) z)
                :colorscale [[0 (s 1)] [0.33 (s 2)] [0.34 (s 3)] [0.66 (s 4)] [0.67 (s 5)] [1 (s 6)]]
                :colorbar {:title {:text "impact"} :tickvals [1 2 3 4] :ticktext impact-levels}
                :hovertemplate "%{y}<br>%{x}<br>%{text}<extra></extra>"}
               {:type "heatmap" :x (mapv (:col-names m) so) :y labels :z z :xgap 1 :ygap 1
                :colorscale [[0 (s 4)] [1 (s 4)]] :showscale false
                :hovertemplate "%{y}<br>%{x}<br>mutated<extra></extra>"})]
      :layout (layout title {:xaxis {:title {:text (format "samples (%d)" (count so))}
                                     :showticklabels (<= (count so) 40) :type "category"}
                             :yaxis {:title {:text ""} :autorange "reversed" :type "category"
                                     :tickmode "linear" :dtick 1 :automargin true}})})))

;; -- cohort comparison plots (from the variant forensics reports) --

(defn- floor-values [xs floor] (mapv #(max (double %) (double floor)) xs))

(defn genex-vs-cohort
  "Participant sample expression vs reference cohort distributions per gene.

  sample-expr: [[sample-id [[gene value] ...]] ...] (one entry per sample)
  cohort-expr: [[[gene [values...]] ...] ...] (one entry per cohort) or a
               single cohort [[gene [values...]] ...]
  opts: :style (:box vertical boxes | :horizontal boxes | :violin horizontal
        half-violins; default :box), :cohort-names, :log? (default true),
        :floor (value floor before log scaling; 0.01, violins 1),
        :title, :value-label"
  ([sample-expr cohort-expr] (genex-vs-cohort sample-expr cohort-expr {}))
  ([sample-expr cohort-expr {:keys [style cohort-names log? floor title value-label]
                             :or {style :box log? true title "Participant vs. cohort gene expression"
                                  value-label "expression"}}]
   (let [cohorts (if (string? (ffirst cohort-expr)) [cohort-expr] cohort-expr)
         floor (or floor (if (= style :violin) 1 0.01))
         fl #(if log? (floor-values % floor) (vec %))
         horizontal? (not= style :box)
         n-genes (count (first cohorts))
         cohort-traces
         (for [[ci cohort] (map-indexed vector cohorts)
               [gi [gene vals]] (map-indexed vector cohort)
               :let [c (series-color (+ ci (count sample-expr)))
                     cname (get (vec cohort-names) ci (str "cohort " (inc ci)))
                     xs (fl vals)
                     common {:name cname :legendgroup cname :showlegend (zero? gi)}]]
           (merge common
                  (case style
                    :violin {:type "violin" :orientation "h" :side "negative" :scalemode "width" :width 1.6
                             :x xs :y (vec (repeat (count xs) gene))
                             :line {:color c :width 1} :fillcolor (alpha c 0.3) :points false
                             :box {:visible true} :meanline {:visible true}}
                    :horizontal {:type "box" :orientation "h" :x xs :y (vec (repeat (count xs) gene))
                                 :marker {:color c} :line {:color c} :boxpoints false}
                    {:type "box" :y xs :x (vec (repeat (count xs) gene))
                     :marker {:color c} :line {:color c} :boxpoints false})))
         sample-traces
         (map-indexed
           (fn [i [sid sexpr]]
             (let [genes (mapv first sexpr) vals (fl (map second sexpr))]
               (merge {:type "scatter" :mode "markers" :name sid
                       :marker {:size 13 :color (series-color i) :symbol "diamond"
                                :line {:width 1 :color (:text_primary (theme))}}}
                      (if horizontal? {:x vals :y genes} {:x genes :y vals}))))
           sample-expr)
         val-axis {:type (if log? "log" "linear")
                   :title {:text (if log? (str value-label " (log scale)") value-label)}}]
     {:data (vec (concat cohort-traces sample-traces))
      :layout (layout title (merge {:showlegend true}
                                   (if horizontal?
                                     {:xaxis val-axis :yaxis {:title {:text "gene"} :type "category"}
                                      :height (+ 180 (* 40 n-genes))}
                                     (cond-> {:yaxis val-axis :xaxis {:type "category"} :height 600}
                                       (> (count cohorts) 1) (assoc :boxmode "group")))
                                   (when (= style :violin) {:violinmode "overlay"})))})))

(defn samples-vs-histogram
  "Cohort distribution of one gene's expression as a histogram, with each
  participant sample's value as a dashed vertical line.
  cohort-values: [numbers]; sample-values: [[sample-id value] ...]"
  [gene cohort-values sample-values {:keys [log? title] :or {log? true}}]
  {:data [{:type "histogram" :x (vec cohort-values) :name (str gene " cohort distribution")
           :marker {:color (:grid (theme)) :line {:color (:text_secondary (theme)) :width 0.5}}}]
   :layout (layout (or title (str "Expression of " gene))
                   {:shapes (vec (map-indexed (fn [i [_ v]]
                                                {:type "line" :x0 v :x1 v :y0 0 :y1 1 :yref "paper"
                                                 :line {:color (series-color i) :width 2 :dash "dash"}})
                                              sample-values))
                    :annotations (vec (map-indexed (fn [i [sid v]]
                                                     {:x v :y 1 :yref "paper" :text sid :showarrow false
                                                      :yanchor "bottom" :font {:color (series-color i)}})
                                                   sample-values))
                    :xaxis {:title {:text "expression"} :type (if log? "log" "linear")}
                    :yaxis {:title {:text "samples"}}})})

(defn ma-plot
  "MA plot from patternq.quant/ma-values output [[gene M A] ...]: M
  (log2 fold change) vs A (mean log2 expression). Genes with |M| >=
  :threshold (default 1) are colored up/down; the top :n-labels (default 8,
  split between up and down) by |M| x A among genes with A >= :min-a
  (default 2, i.e. reasonably expressed) are labeled."
  ([ma] (ma-plot ma {}))
  ([ma {:keys [threshold n-labels min-a title x-label]
        :or {threshold 1.0 n-labels 8 min-a 2.0 title "MA plot" x-label "A: mean log2(expression + 1)"}}]
   (let [th (theme)
         {:keys [high low]} (:diverging th)
         cls (fn [[_ m _]] (cond (>= (double m) (double threshold)) :up
                                 (<= (double m) (- (double threshold))) :down
                                 :else :ns))
         groups (group-by cls ma)
         tr (fn [k nm color]
              (let [pts (groups k)]
                {:type "scattergl" :mode "markers" :name (format "%s (%d)" nm (count pts))
                 :x (mapv #(nth % 2) pts) :y (mapv second pts) :text (mapv first pts)
                 :marker {:color color :size (if (= k :ns) 4 6) :opacity (if (= k :ns) 0.4 0.8)}
                 :hovertemplate "%{text}<br>M=%{y:.2f}<br>A=%{x:.2f}<extra></extra>"}))
         top (fn [k n] (->> (groups k) (filter #(>= (double (nth % 2)) (double min-a)))
                            (sort-by #(- (* (Math/abs (double (second %))) (double (nth % 2))))) (take n)))
         labeled (concat (top :up (quot (long n-labels) 2)) (top :down (- (long n-labels) (quot (long n-labels) 2))))]
     {:data [(tr :ns "not changed" (:text_secondary th)) (tr :up "up" high) (tr :down "down" low)]
      :layout (layout title {:xaxis {:title {:text x-label}}
                             :yaxis {:title {:text "M: log2 fold change"} :zeroline true}
                             :shapes [{:type "line" :xref "paper" :x0 0 :x1 1 :y0 threshold :y1 threshold
                                       :line {:color (:grid th) :dash "dot"}}
                                      {:type "line" :xref "paper" :x0 0 :x1 1 :y0 (- (double threshold)) :y1 (- (double threshold))
                                       :line {:color (:grid th) :dash "dot"}}]
                             :annotations (mapv (fn [[g m a]] {:x a :y m :text g :showarrow true :arrowsize 0.5
                                                               :ax 15 :ay -15 :font {:size 10}})
                                                labeled)})})))

;; -- expression comparison plots (parity with R plot_zscores, plot_vs_cohort,
;;    plot_ma, plot_fold_change) --

(defn zscores
  "Horizontal bars of the top genes of a patternq.expression/compare-to-cohort
  result (via top-by-zscore rules: :n, :min-value, :min-observed), colored
  up/down on the diverging scale."
  ([comparison] (zscores comparison {}))
  ([comparison {:keys [n min-value min-observed title] :or {n 30 min-value 1 min-observed 0.5}}]
   (let [top ((requiring-resolve 'patternq.expression/top-by-zscore)
              comparison {:n n :min-value min-value :min-observed min-observed})
         top (sort-by :z top)
         {:keys [high low]} (:diverging (theme))
         info (:patternq/comparison (meta comparison))]
     {:data [{:type "bar" :orientation "h" :x (mapv :z top) :y (mapv :hgnc-symbol top)
              :marker {:color (mapv #(if (>= (double (:z %)) 0) high low) top)}
              :text (mapv #(format "value %.3g · cohort median %.3g · %.1f pct" (double (:value %))
                                   (- (Math/pow 2.0 (double (:cohort-median %))) 1.0) (double (:percentile %)))
                          top)
              :textposition "none"
              :hovertemplate "%{y}<br>z = %{x:.2f}<br>%{text}<extra></extra>"}]
      :layout (layout (or title (when info (format "%s vs %s: top genes by z-score" (:sample info) (:cohort-db info))))
                      {:showlegend false :bargap 0.25 :height (+ 140 (* 20 (count top)))
                       :xaxis {:title {:text "z-score vs cohort (log2(1+x))"} :zeroline true}
                       :yaxis {:title {:text ""} :type "category" :categoryorder "array"
                               :categoryarray (mapv :hgnc-symbol top)}})})))

(def ^:private marker-symbols ["diamond" "square" "circle" "triangle-up" "x" "star"])

(defn vs-cohort
  "Sample values (markers) against cohort distributions (violins or boxes), one
  row per gene. On a log scale values are drawn as log10 on a linear axis, so
  violin densities are estimated on the log scale, with power-of-ten tick
  labels.
  sample-expr: rows of :sample-id :hgnc-symbol :value
  cohort-expr: rows of :sample-id :hgnc-symbol :value and :cohort (label)
  opts: :type (:violin | :box), :log? (true), :floor (0.01), :title, :xlab"
  ([sample-expr cohort-expr] (vs-cohort sample-expr cohort-expr {}))
  ([sample-expr cohort-expr {:keys [type log? floor title xlab]
                             :or {type :violin log? true floor 0.01 title "Samples vs cohort expression" xlab "expression"}}]
   (let [th (theme)
         sum-rows (fn [rows ks] (->> rows (group-by (apply juxt ks))
                                     (map (fn [[k rs]] (assoc (zipmap ks k) :value (reduce + (map :value rs)))))))
         ce (sum-rows (map #(update % :cohort (fnil identity "cohort")) cohort-expr) [:cohort :sample-id :hgnc-symbol])
         se (sum-rows sample-expr [:sample-id :hgnc-symbol])
         genes (vec (distinct (concat (map :hgnc-symbol se) (map :hgnc-symbol ce))))
         fl (fn [v] (if log? (Math/log10 (max (double v) (double floor))) (double v)))
         cohorts (sort (distinct (map :cohort ce)))
         ctraces (map-indexed
                   (fn [i c]
                     (let [rs (filter #(= c (:cohort %)) ce) col (series-color i)]
                       (if (= type :violin)
                         {:type "violin" :orientation "h" :x (mapv (comp fl :value) rs) :y (mapv :hgnc-symbol rs)
                          :name c :legendgroup c :line {:color col :width 1} :fillcolor (alpha col 0.25)
                          :points false :spanmode "hard" :scalemode "width" :width 0.8
                          :box {:visible true :fillcolor (:surface th) :line {:color col} :width 0.2}
                          :meanline {:visible false} :hoverinfo "y+name"}
                         {:type "box" :orientation "h" :x (mapv (comp fl :value) rs) :y (mapv :hgnc-symbol rs)
                          :name c :legendgroup c :marker {:color col :size 3} :line {:color col}
                          :fillcolor (alpha col 0.2) :boxpoints "outliers"})))
                   cohorts)
         sids (sort (distinct (map :sample-id se)))
         straces (map-indexed
                   (fn [j sid]
                     (let [rs (filter #(= sid (:sample-id %)) se)]
                       {:type "scatter" :mode "markers" :x (mapv (comp fl :value) rs) :y (mapv :hgnc-symbol rs)
                        :name sid :customdata (mapv :value rs)
                        :marker {:symbol (marker-symbols (mod j (count marker-symbols))) :size 11
                                 :color (:text_primary th) :line {:color (:surface th) :width 1.5}}
                        :hovertemplate (str sid "<br>%{y}: %{customdata:.3g}<extra></extra>")}))
                   sids)
         all-x (filter #(Double/isFinite %) (map (comp fl :value) (concat ce se)))
         xaxis (cond-> {:title {:text xlab}}
                 (and log? (seq all-x))
                 (merge (let [ticks (range (long (Math/floor (apply min all-x))) (inc (long (Math/ceil (apply max all-x)))))]
                          {:tickvals (vec ticks)
                           :ticktext (mapv #(let [v (Math/pow 10.0 %)]
                                              (if (>= (long %) 3) (format "%,d" (long v))
                                                  (.toPlainString (.stripTrailingZeros (java.math.BigDecimal. (str v))))))
                                           ticks)})))]
     {:data (vec (concat ctraces straces))
      :layout (layout title {:height (+ 160 (* 42 (count genes))) :boxmode "group" :xaxis xaxis
                             :yaxis {:title {:text ""} :type "category" :categoryorder "array"
                                     :categoryarray (vec (reverse genes)) :autorange true}})})))

(defn ma
  "MA plot of a patternq.expression/compare-samples result: log2 fold change
  vs average log10(1 + x); genes beyond :lfc-threshold (2.5) colored
  up/down, the top :label-top (8) each way with avg >= :label-min-avg (1.5)
  labeled (thinned so labels don't pile up), plus any :highlight genes."
  ([change] (ma change {}))
  ([change {:keys [lfc-threshold highlight label-top label-min-avg title]
            :or {lfc-threshold 2.5 label-top 8 label-min-avg 1.5}}]
   (let [th (theme)
         {:keys [high low]} (:diverging th)
         t (double lfc-threshold)
         info (:patternq/comparison (meta change))
         cls (fn [r] (cond (>= (double (:lfc r)) t) :up (<= (double (:lfc r)) (- t)) :down :else :unchanged))
         spec {:unchanged [(alpha other-color 0.35) "within threshold"]
               :down [low (str "down (lfc ≤ -" lfc-threshold ")")]
               :up [high (str "up (lfc ≥ " lfc-threshold ")")]}
         traces (for [k [:unchanged :down :up]
                      :let [rs (filter #(= k (cls %)) change) [c nm] (spec k)]
                      :when (seq rs)]
                  {:type "scattergl" :mode "markers" :x (mapv :avg-log10 rs) :y (mapv :lfc rs)
                   :text (mapv :hgnc-symbol rs) :name nm
                   :marker {:color c :size (if (= k :unchanged) 5 7)}
                   :customdata (mapv (juxt :value-a :value-b) rs)
                   :hovertemplate "%{text}<br>lfc %{y:.2f}<br>a %{customdata[0]:.3g} → b %{customdata[1]:.3g}<extra></extra>"})
         cand (filter #(>= (double (:avg-log10 %)) (double label-min-avg)) change)
         lab (->> (concat (take label-top (sort-by #(- (double (:lfc %))) cand))
                          (take label-top (sort-by :lfc cand)))
                  (filter #(>= (Math/abs (double (:lfc %))) t)))
         lab (distinct (concat (filter #(some #{(:hgnc-symbol %)} highlight) change) lab))
         xs (map :avg-log10 change) ys (map :lfc change)
         xr (- (double (apply max 1 xs)) (double (apply min 0 xs)))
         yr (- (double (apply max 1 ys)) (double (apply min -1 ys)))
         kept (reduce (fn [kept r]
                        (if (some #(and (< (Math/abs (- (double (:avg-log10 %)) (double (:avg-log10 r)))) (* 0.06 xr))
                                        (< (Math/abs (- (double (:lfc %)) (double (:lfc r)))) (* 0.05 yr)))
                                  kept)
                          kept (conj kept r)))
                      [] lab)
         guide (fn [y] {:type "line" :xref "paper" :x0 0 :x1 1 :y0 y :y1 y
                        :line {:color (:text_secondary th) :width 1 :dash (if (zero? (double y)) "solid" "dot")}})]
     {:data (vec traces)
      :layout (layout (or title (when (:sample-a info) (format "Expression change: %s → %s" (:sample-a info) (:sample-b info))))
                      {:hovermode "closest"
                       :annotations (mapv (fn [r] {:x (:avg-log10 r) :y (:lfc r) :text (:hgnc-symbol r)
                                                   :showarrow true :arrowhead 0 :arrowwidth 1
                                                   :arrowcolor (:text_secondary th) :ax 18
                                                   :ay (if (pos? (double (:lfc r))) -16 16)
                                                   :font {:size 11 :color (:text_primary th)}})
                                          kept)
                       :shapes [(guide 0) (guide t) (guide (- t))]
                       :xaxis {:title {:text "average expression, log10(1 + x)"}}
                       :yaxis {:title {:text "log2 fold change"}}})})))

(defn fold-change
  "Horizontal bars of log2 fold change for :genes (default the 15 most down
  and 15 most up) of a compare-samples result."
  ([change] (fold-change change {}))
  ([change {:keys [genes title] :or {title "Change in gene expression"}}]
   (let [s (sort-by :lfc change)
         genes (set (or genes (concat (map :hgnc-symbol (take 15 s)) (map :hgnc-symbol (take-last 15 s)))))
         s (filter #(genes (:hgnc-symbol %)) s)
         {:keys [high low]} (:diverging (theme))]
     {:data [{:type "bar" :orientation "h" :x (mapv :lfc s) :y (mapv :hgnc-symbol s)
              :marker {:color (mapv #(if (>= (double (:lfc %)) 0) high low) s)}
              :hovertemplate "%{y}: %{x:.2f}<extra></extra>"}]
      :layout (layout title {:showlegend false :height (+ 140 (* 22 (count s)))
                             :xaxis {:title {:text "log2 fold change"}}
                             :yaxis {:title {:text ""} :type "category" :categoryorder "array"
                                     :categoryarray (mapv :hgnc-symbol s)}})})))
