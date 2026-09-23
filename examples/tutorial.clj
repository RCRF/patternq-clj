^{:kindly/hide-code true
  :clay {:title "patternq tutorial"
         :quarto {:title "patternq tutorial"
                  :subtitle "Querying and analyzing the Pattern Data Commons · Clojure"
                  :format {:html {:theme :cosmo :toc true :embed-resources true :code-overflow :wrap}}}}}
(ns tutorial
  (:require [patternq.clinical :as clinical]
            [patternq.context :as ctx]
            [patternq.dataset :as pqd]
            [patternq.db :as pdb]
            [patternq.expression :as expr]
            [patternq.http :as http]
            [patternq.plot :as plot]
            [patternq.results :as res]
            [patternq.survival :as surv]
            [scicloj.kindly.v4.kind :as kind]))

;; The basic patternq workflow, step by step. The same steps, with the same
;; datasets and outputs, are in `tutorial.qmd` in [patternq-r](https://github.com/RCRF/patternq-r) (R, knitr),
;; `tutorial.qmd` in [patternq](https://github.com/RCRF/patternq) (marimo) and `tutorial.qmd` in [PatternQ.jl](https://github.com/RCRF/PatternQ.jl) (Julia).
;;
;; | step | what |
;; |---|---|
;; | 1 | connect and find a dataset |
;; | 2 | explore a dataset: assays, measurement sets, samples, subjects |
;; | 3 | write a Datalog query |
;; | 4 | canned queries: variants and gene expression |
;; | 5 | any measurement, reshaped and joined with sample context |
;; | 6 | clinical outcomes and survival |
;; | 7 | a sample against a reference cohort |
;; | 8 | expression change between two samples |
;;
;; By default the Clojure library queries the Pattern Data Commons query
;; service over HTTP, like the R and Python libraries: it needs only
;; `PATTERNQ_API_KEY`. (Internal users can switch to the Datomic peer; see the
;; library README.)

^{:kindly/hide-code true :kindly/kind :kind/hidden}
(defn- fmt [v]
  (if (and (float? v) (not (Double/isNaN v)))
    (let [a (Math/abs (double v))]
      (cond (zero? a) "0" (or (>= a 1e5) (< a 1e-3)) (format "%.3g" v) :else (format "%.4g" v)))
    v))

^{:kindly/hide-code true :kindly/kind :kind/hidden}
(defn table
  "Show rows (maps) as a table: first n rows, optionally only `cols`, floats
  rounded for display."
  ([rows] (table rows 10))
  ([rows n] (table rows n nil))
  ([rows n cols]
   (kind/table {:row-maps (mapv (fn [r] (update-vals (if cols (select-keys r cols) r) fmt)) (take n rows))})))

;; ## 1. Connect and find a dataset
;;
;; Set `PATTERNQ_API_KEY` (and optionally `PATTERNQ_ENDPOINT`) in the
;; environment (or call `http/set-token!` / `http/set-endpoint!`).

(http/endpoint)

(def datasets (http/list-datasets))

(table (map #(select-keys % [:dataset :db :patient-count :sample-count]) datasets) 6)

;; Every dataset is its own database. Dataset names are stable; database
;; names change when a dataset is re-imported, so resolve them:

(def h37001 (http/resolve-db "H37001"))
(def uvm (http/resolve-db "tcga-uvm"))
[h37001 uvm]

;; ## 2. Explore a dataset
;;
;; Functions take a database name (or a database value) as their first
;; argument.

(table (pqd/dataset-summary h37001))

(table (pqd/measurement-types h37001 "RSEM RNASeq"))

(table (pqd/samples h37001))

(table (pqd/subjects h37001))

;; ## 3. Write a Datalog query
;;
;; In Clojure, queries are plain Datomic Datalog written as data, run with
;; `patternq.db/q` against a database handle (the query is converted to the
;; query service's JSON form behind the scenes):

(def db (pdb/db h37001))

(table (map #(zipmap [:sample-id :site] %)
            (pdb/q '[:find ?sample-id ?site
                   :where
                   [?s :sample/id ?sample-id]
                   [?s :sample/freetext-anatomic-site ?site]]
                 db)))

;; Inputs bind with `:in`; a vector binds a collection:

(def vaf-rows
  (pdb/q '[:find ?sample-id ?hgnc ?vaf
         :in $ [?hgnc ...]
         :where
         [?g :gene/hgnc-symbol ?hgnc]
         [?v :variant/gene ?g]
         [?m :measurement/variant ?v]
         [?m :measurement/vaf ?vaf]
         [?m :measurement/sample ?s]
         [?s :sample/id ?sample-id]]
       db ["BAP1" "GNAQ" "SF3B1"]))

(table (map #(zipmap [:sample-id :hgnc :vaf] %) vaf-rows))

;; Results of the canned queries record where they came from (database and
;; basis t), in metadata:

(:patternq/provenance (meta (pqd/samples h37001)))

;; ## 4. Variants and gene expression

(def v (pqd/variants h37001))

(count v)

(table (->> v (filter #(= "high" (:impact %))) (map #(select-keys % [:sample-id :hgnc-symbol :HGVSp :vaf]))) 6)

(kind/plotly (plot/vaf-histogram v {:samples ["H37001-001" "H37001-006"]}))

(def gx (pqd/gene-expression h37001 {:genes ["BAP1" "GNAQ" "PRAME" "PMEL" "MLANA"] :measurement :tpm}))

(table gx 6)

(kind/plotly (plot/gene-expression gx {:title "Melanocytic genes (TPM)" :ylab "TPM" :log? true}))

;; Copy number queries always take a subset (genes, samples or subjects):

(table (pqd/cnv-segments h37001 {:genes ["BAP1" "GNAQ"]}) 20
       [:sample-id :hgnc-symbol :contig :start :end :segment-mean-lrr :absolute-cn])

;; ## 5. Any measurement, reshaped and joined with context
;;
;; `measurements` works for any measurement attribute and detects what each
;; measurement is of (gene, cell population, epitope, ...). Here: CyTOF cell
;; population frequencies in PRINCE.

(def prince (http/resolve-db "prince-2022"))

(table (pqd/measurement-types prince "PICI CyTOF Immune Profiling"))

(def cy (pqd/measurements prince :percent-of-parent {:measurement-set "PICI CyTOF Immune Profiling"}))

(table cy 6)

;; Wide (samples x populations):

(def m (res/->matrix cy :cell-population))

[(count (:row-names m)) (count (:col-names m))]

;; With sample and subject context:

(def cy-ctx (ctx/add-sample-context cy prince))

(table (map #(select-keys % [:sample-id :subject-id :timepoint-id :cell-population :value]) cy-ctx) 6)

(kind/plotly
  (plot/by-timepoint (filter #(= "HLA-DR+ Non-Naive CD8 T Cells" (:cell-population %)) cy-ctx)
                     {:timepoints ["C1D1" "C1D15" "C2D1" "C4D1"] :lines? true
                      :title "HLA-DR+ non-naive CD8 T cells" :ylab "% of parent"}))

;; ## 6. Clinical outcomes and survival

(def outcomes (clinical/subject-outcomes prince))

(table outcomes 6)

(kind/plotly (plot/survival outcomes {:group :bor :xlab "months"
                                      :title "Overall survival by best overall response"}))

;; Survival split at the median of a baseline biomarker, with a log-rank
;; test:

(def base
  (->> cy-ctx
       (filter #(and (= "HLA-DR+ Non-Naive CD8 T Cells" (:cell-population %)) (= "C1D1" (:timepoint-id %))))
       (group-by :subject-id)
       (into {} (map (fn [[s rs]] [s (/ (reduce + (map :value rs)) (count rs))])))))

(def km (surv/survival-by-median base outcomes))

(:p (surv/logrank km))

(kind/plotly (plot/survival km {:group :group :levels ["low" "high"] :xlab "months"
                                :title "OS by baseline HLA-DR+ non-naive CD8 T cells"}))

;; ## 7. A sample against a reference cohort
;;
;; Place a patient sample within a reference cohort's distribution, gene by
;; gene (z-score and percentile on log2(1 + TPM)). Only compare comparable
;; measurements, and check `:cohort-observed`.

(def genes ["BAP1" "GNAQ" "PRAME" "MITF" "PMEL" "TYR" "MLANA" "CD8A" "GZMB"])

(def cmp (expr/compare-to-cohort "H37001-003" h37001 uvm {:genes genes}))

(table (map #(select-keys % [:hgnc-symbol :value :z :percentile :cohort-observed]) cmp) 20)

(kind/plotly (expr/examine-geneset h37001 {"TCGA-UVM" uvm} genes ["H37001-003" "H37001-001"]
                                   {:title "H37001 vs TCGA-UVM"}))

;; ## 8. Expression change between two samples

(def change (expr/compare-samples "H37001-003" "H37001-001" h37001 {}))

(table change 6)

(kind/plotly (plot/ma change))
