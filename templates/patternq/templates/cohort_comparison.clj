^{:kindly/hide-code true
  :clay {:quarto {:title "Participant vs. Reference Cohort Expression"
                  :format {:html {:embed-resources true :theme :cosmo}}}}}
(ns patternq.templates.cohort-comparison
  "Participant gene expression vs a reference cohort (e.g. TCGA): top genes
  by z-score, gene set comparisons, and a two-sample MA plot."
  (:require [patternq.expression :as expr]
            [patternq.genesets :as gs]
            [patternq.http :as http]
            [patternq.plot :as plot]
            [patternq.quant :as q]
            [scicloj.kindly.v4.kind :as kind]))

;; # Participant vs. Reference Cohort
(def db-name (http/resolve-db "H37001"))
(def cohort-db (http/resolve-db "tcga-uvm"))
(def sample-ids ["H37001-001" "H37001-004"])

;; Pick an expression measurement both databases share:
(def genex (expr/resolve-genex-attrs db-name cohort-db))
genex
(def measurement (:participant-attr genex))

;; ## Top genes by z-score vs the cohort
(def zs (expr/cohort-zscores db-name cohort-db (first sample-ids) nil {:measurement measurement}))

^:kindly/hide-code
(kind/table {:row-maps (vec (take 25 zs))})

;; ## Gene sets
^:kindly/hide-code
(kind/plotly (expr/examine-geneset db-name [cohort-db] (gs/geneset "melanoma-proliferative") sample-ids
                                   {:measurement measurement :style :violin :cohort-names ["reference"]}))

^:kindly/hide-code
(kind/plotly (expr/examine-geneset db-name [cohort-db] gs/glutamine-vulnerable sample-ids
                                   {:measurement measurement :title "Glutaminase inhibitor vulnerability"}))

;; ## Two-sample difference (MA plot)
(def ma (q/ma-values (expr/sample->gene-expression db-name (first sample-ids) {:measurement measurement})
                     (expr/sample->gene-expression db-name (second sample-ids) {:measurement measurement})))

^:kindly/hide-code
(kind/plotly (plot/ma-plot ma {:title (str (second sample-ids) " vs " (first sample-ids))}))
