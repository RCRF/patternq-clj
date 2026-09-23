^{:kindly/hide-code true
  :clay {:title "PRINCE: clinical and immune correlates with patternq"
         :quarto {:title "PRINCE: clinical and immune correlates with patternq"
                  :subtitle "prince-2022 on the Pattern Data Commons · Clojure"
                  :format {:html {:theme :cosmo :toc true :embed-resources true :code-overflow :wrap}}}}}
(ns prince-demo
  (:require [patternq.clinical :as clinical]
            [patternq.context :as ctx]
            [patternq.dataset :as pqd]
            [patternq.http :as http]
            [patternq.plot :as plot]
            [patternq.results :as res]
            [patternq.survival :as surv]
            [scicloj.kindly.v4.kind :as kind]))

;; This notebook re-creates, on the Pattern Data Commons copy of the PRINCE
;; trial, the kinds of analyses in Padron et al., *Nat Med* 2022 (PRINCE:
;; sotigalimab and/or nivolumab with gemcitabine/nab-paclitaxel in first-line
;; metastatic pancreatic cancer): overall survival, circulating immune
;; populations over treatment, survival stratified at the median of a
;; baseline biomarker, marker expression on those populations, serum
;; proteins, and tumor gene expression signatures.
;;
;; It is a demonstration of the `patternq` workflow, not a re-analysis.
;;
;; > **Treatment arm is not in the `prince-2022` database.** The paper
;; > stratifies most analyses by arm (nivo/chemo, sotiga/chemo,
;; > sotiga/nivo/chemo). Arm is in the source data (`PICI0002_ph2_clinical.csv`,
;; > `Arm` / `Actual Arm`) but was not imported, so here every analysis pools
;; > patients across arms and groups by survival status at 1 year. P-values
;; > are unadjusted and exploratory.
;;
;; The same notebook exists in R (knitr), Python (marimo) and Julia in
;; `prince_demo.qmd` in [patternq-r](https://github.com/RCRF/patternq-r), `prince_demo.qmd` in [patternq](https://github.com/RCRF/patternq) and
;; `prince_demo.qmd` in [PatternQ.jl](https://github.com/RCRF/PatternQ.jl).

^{:kindly/hide-code true :kindly/kind :kind/hidden}
(defn- fmt [v]
  (if (and (float? v) (not (Double/isNaN v)))
    (let [a (Math/abs (double v))]
      (cond (zero? a) "0" (or (>= a 1e5) (< a 1e-3)) (format "%.3g" v) :else (format "%.4g" v)))
    v))

^{:kindly/hide-code true :kindly/kind :kind/hidden}
(defn table
  ([rows] (table rows 10))
  ([rows n] (kind/table {:row-maps (mapv #(update-vals % fmt) (take n rows))})))

^{:kindly/hide-code true :kindly/kind :kind/hidden}
(defn crosstab
  "Counts of rows by (row-key, col-key) as a table."
  [rows row-key col-key]
  (let [label #(or (some-> % str) "NA")
        cols (sort (distinct (map #(label (get % col-key)) rows)))
        freq (frequencies (map (juxt #(label (get % row-key)) #(label (get % col-key))) rows))]
    (kind/table {:column-names (into [(name row-key)] cols)
                 :row-vectors (vec (for [r (sort (distinct (map #(label (get % row-key)) rows)))]
                                     (into [r] (map #(get freq [r %] 0) cols))))})))

;; ## The dataset
;;
;; Every dataset is its own database. `http/resolve-db` maps the stable
;; dataset name to the current database.

(def db (http/resolve-db "prince-2022"))

db

(table (pqd/dataset-summary db))

;; Which subjects were profiled by which assays:

(def sa (pqd/sample-assays db))

(kind/plotly (plot/sample-overview sa))

;; Samples by timepoint and specimen (blood draws at C1D1, C1D15, C2D1, C4D1;
;; tumor biopsies at screening):

(def smp (pqd/samples db))

(crosstab smp :timepoint-id :sample-specimen)

;; ## Outcomes
;;
;; Best overall response, progression-free and overall survival (months) per
;; subject, and survival status at 1 year (the paper's primary endpoint is
;; 1-year OS). Subjects censored before 1 year have unknown status.

(def outcomes
  (mapv #(assoc % :status-1y (surv/survival-status % {:time :os :event :os-event :at 12}))
        (clinical/subject-outcomes db)))

(def status-levels ["alive at 1 year" "died within 1 year"])

(frequencies (map :status-1y outcomes))

(frequencies (map :bor outcomes))

;; Overall survival (Fig. 2a pools arms here) and by best overall response:

(kind/plotly (plot/survival outcomes {:title "Overall survival" :xlab "months"}))

(kind/plotly (plot/survival outcomes {:group :bor :title "Overall survival by best overall response" :xlab "months"}))

;; Progression-free survival (Extended Data Fig. 2a):

(kind/plotly (plot/survival outcomes {:time :pfs :event :pfs-event :title "Progression-free survival" :xlab "months"}))

;; ## Circulating immune populations
;;
;; Frequencies of gated populations from the X50 (flow) and CyTOF panels,
;; with sample context (subject, timepoint) and 1-year survival status joined
;; on.

(def status (into {} (map (juxt :subject-id :status-1y)) outcomes))

^:kind/hidden
(defn with-status [rows]
  (->> (ctx/add-sample-context rows db {:subjects? false})
       (filter #(contains? status (:subject-id %)))
       (mapv #(assoc % :status-1y (status (:subject-id %))))))

(def x50 (with-status (pqd/measurements db :percent-of-parent {:measurement-set "PICI X50 Immune Profiling"})))

(def cytof (with-status (concat (pqd/measurements db :percent-of-parent {:measurement-set "PICI CyTOF Immune Profiling"})
                                (pqd/measurements db :percent-of-leukocytes {:measurement-set "PICI CyTOF Immune Profiling"}))))

(def blood-timepoints ["C1D1" "C1D15" "C2D1" "C4D1"])

(sort (distinct (map :cell-population x50)))

;; Baseline (C1D1) value per subject for one population:

^:kind/hidden
(defn baseline
  ([rows population] (baseline rows population "C1D1"))
  ([rows population timepoint]
   (->> rows
        (filter #(and (= population (:cell-population %)) (= timepoint (:timepoint-id %))))
        (group-by :subject-id)
        (into {} (map (fn [[s rs]] [s (/ (reduce + (map :value rs)) (count rs))]))))))

;; ### CD38+ effector memory CD8 T cells (Fig. 3a, c)
;;
;; OS stratified by the pre-treatment frequency, above vs below the median:

(def cd38-pop "CD38+ effector memory CD8+ T cells (% effector memory CD8 T cells)")

(def cd38-km (surv/survival-by-median (baseline x50 cd38-pop) outcomes))

(:p (surv/logrank cd38-km))

(kind/plotly (plot/survival cd38-km {:group :group :levels ["low" "high"] :xlab "months"
                                     :title "OS by baseline CD38+ EM CD8 T cells (median split)"}))

;; Frequencies before and on treatment, grouped by 1-year survival, with each
;; patient's trajectory:

(kind/plotly (plot/by-timepoint (filter #(= cd38-pop (:cell-population %)) x50)
                                {:group :status-1y :lines? true :timepoints blood-timepoints :levels status-levels
                                 :title "CD38+ EM CD8 T cells over treatment" :ylab "% of EM CD8 T cells"}))

;; ### PD-1+CD39+ EM1 CD4 T cells (Fig. 3d) and PD-1+Tbet+ non-naive CD4 T cells (Fig. 4f, h)

(kind/fragment
  (for [pop ["PD-1+CD39+ effector memory 1 CD4+ T cells (% of effector memory 1 CD4 T cells)"
             "PD-1+Tbet+ non-Naive CD4+ T cells (% of CD4 not naive T cells)"]]
    (kind/plotly (plot/survival (surv/survival-by-median (baseline x50 pop) outcomes)
                                {:group :group :levels ["low" "high"] :xlab "months"
                                 :title (str "OS by baseline " (first (clojure.string/split pop #" \(")))}))))

;; ### Cross-presenting dendritic cells (Fig. 4b)

(kind/plotly (plot/survival (surv/survival-by-median (baseline cytof "CD1C+ CD141+ DC (% of leukocytes)") outcomes)
                            {:group :group :levels ["low" "high"] :xlab "months"
                             :title "OS by baseline CD1c+ CD141+ DCs (CyTOF)"}))

;; ### Proliferating T cells over treatment (Extended Data Fig. 3a)
;;
;; Change relative to C1D1 (log2 fold change of the frequency):

(def ki67 (surv/change-from-baseline (filter #(= "Ki67+ Non-Naive CD8 T Cells" (:cell-population %)) x50)
                                     {:baseline "C1D1" :method :log2-ratio :pseudocount 0.001}))

(kind/plotly (plot/by-timepoint ki67 {:value :change :group :status-1y :lines? true :timepoints blood-timepoints
                                      :levels status-levels
                                      :title "Ki-67+ non-naive CD8 T cells, change from C1D1"
                                      :ylab "log2 fold change vs C1D1"}))

;; ### Screening all baseline populations
;;
;; Every X50 and CyTOF population at C1D1, median split, log-rank p (the
;; paper's Fig. 5 summarizes such associations per arm):

^:kind/hidden
(defn screen [rows panel]
  (for [pop (distinct (map :cell-population rows))
        :let [km (surv/survival-by-median (baseline rows pop) outcomes)
              {:keys [p observed expected]} (surv/logrank km)]
        :when p]
    {:panel panel :cell-population pop :n (count km) :p p
     :higher-is (if (< (/ (double (observed "high")) (double (expected "high")))
                       (/ (double (observed "low")) (double (expected "low"))))
                  "longer survival" "shorter survival")}))

(def assoc-screen (sort-by :p (concat (screen x50 "X50") (screen cytof "CyTOF"))))

(table assoc-screen 12)

;; ## Marker expression on populations (Extended Data Figs. 7b, 8b)
;;
;; Median channel values of markers on CD38+ non-naive CD8 T cells,
;; pre-treatment, per patient, z-scored per marker, patients annotated by
;; 1-year survival:

(def mcv (with-status (pqd/measurements db :median-channel-value {:measurement-set "PICI X50 Immune Profiling"})))

(def cd8 (filter #(and (= "CD38+ non-naive CD8+ T cells (% non-naive CD8 T cells)" (:cell-population %))
                       (= "C1D1" (:timepoint-id %)))
                 mcv))

(kind/plotly (plot/heatmap (res/->matrix cd8 :subject-id {:row :epitope-id})
                           {:scale :row :col-groups status
                            :title "Markers on CD38+ non-naive CD8 T cells (C1D1, z-scored MFI)"}))

;; CCR7+ CD11b+ CD27- B cells (CyTOF) on treatment (C1D15):

(def mcv-b (with-status (pqd/measurements db :median-channel-value {:measurement-set "PICI CyTOF Immune Profiling"})))

(kind/plotly (plot/heatmap (res/->matrix (filter #(= "C1D15" (:timepoint-id %)) mcv-b) :subject-id {:row :epitope-id})
                           {:scale :row :col-groups status
                            :title "Markers on CCR7+ CD11b+ CD27- B cells (C1D15, z-scored MFI)"}))

;; ## Serum proteins (Extended Data Fig. 3c–f)
;;
;; Olink NPX is on a log2 scale, so change from C1D1 is a difference:

(def olink (with-status (pqd/measurements db :olink-npx {:measurement-set "PICI Olink Proteomics"})))

(def olink-change (surv/change-from-baseline olink {:baseline "C1D1" :method :difference}))

;; Proteins as `[label Olink-id]` (Olink abbreviates CXCL10 as CXL10):

(def proteins [["IFNG" "IFNG"] ["PD-1" "PDCD1"] ["CXCL9" "CXCL9"] ["CXCL10" "CXL10"]])

(kind/fragment
  (for [[label id] proteins]
    (kind/plotly (plot/by-timepoint (filter #(= id (:epitope-id %)) olink-change)
                                    {:value :change :group :status-1y :timepoints blood-timepoints
                                     :levels status-levels :title (str label " change from C1D1")
                                     :ylab "log2 NPX change"}))))

;; ## Tumor: gene expression signatures and multiplex IF
;;
;; Pre-treatment tumor RNA-seq signature scores, clustered, patients
;; annotated by 1-year survival (Extended Data Figs. 4a, 6a). The RNA-seq data
;; comes from the cBioPortal / iAtlas release, which labels the pre-treatment
;; timepoint `0` (the PICI modalities use `screening` / `C1D1`):

(def sig (with-status (pqd/measurements db :nanostring-signature-score
                                        {:measurement-set "cBioPortal Gene Expression Signatures"})))

(def sig-pre (filter #(contains? #{"screening" "0"} (:timepoint-id %)) sig))

(kind/plotly (plot/heatmap (res/->matrix sig-pre :subject-id {:row :signature})
                           {:scale :row :col-groups status
                            :title "Tumor gene expression signatures (pre-treatment, z-scored)"}))

;; OS by the Th1 signature (Extended Data Fig. 6b):

^:kind/hidden
(defn subject-means [rows]
  (->> rows (group-by :subject-id)
       (into {} (map (fn [[s rs]] [s (/ (reduce + (map :value rs)) (count rs))])))))

(kind/plotly (plot/survival (surv/survival-by-median (subject-means (filter #(= "Bindea_Th1_Cells" (:signature %)) sig-pre))
                                                     outcomes)
                            {:group :group :levels ["low" "high"] :xlab "months"
                             :title "OS by tumor Th1 signature (median split)"}))

;; OS by intratumoral iNOS+ macrophages from multiplex IF (Extended Data
;; Fig. 4c):

(def vectra (with-status (pqd/measurements db :percent-of-total-cells {:measurement-set "PICI Vectra Multiplex IF"})))

(def inos (filter #(and (= "iNOS+ macrophages (% of total cells)" (:cell-population %))
                        (= "screening" (:timepoint-id %)))
                  vectra))

(kind/plotly (plot/survival (surv/survival-by-median (subject-means inos) outcomes)
                            {:group :group :levels ["low" "high"] :xlab "months"
                             :title "OS by tumor iNOS+ macrophages (median split)"}))

;; ## Tumor mutations

(kind/plotly (plot/mutation-landscape (pqd/variants db) {:n-genes 20 :title "PRINCE mutation landscape (top 20 genes)"}))

;; ## Provenance
;;
;; Every canned query result records the database and the database basis t
;; it was read at, in metadata:

(:patternq/provenance (meta (pqd/measurements db :percent-of-parent {:measurement-set "PICI X50 Immune Profiling"})))
