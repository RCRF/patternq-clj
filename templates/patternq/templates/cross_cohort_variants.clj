^{:kindly/hide-code true
  :clay {:quarto {:title "Cross-Cohort Variant Checks"
                  :format {:html {:embed-resources true :theme :cosmo}}}}}
(ns patternq.templates.cross-cohort-variants
  "Looking across a multi-participant database: which samples have variant
  data, observations of commonly mutated genes (TERT, TP53), and recurrent
  identical variants that may be artifacts."
  (:require [patternq.http :as http]
            [patternq.variants :as pv]
            [scicloj.kindly.v4.kind :as kind]))

(def db-name (http/resolve-db "tcga-uvm"))

;; ## Samples with variant data
^:kindly/hide-code
(kind/table {:row-maps (pv/variant-samples db-name)})

;; ## TERT and TP53
^:kindly/hide-code
(kind/table {:row-maps (pv/cohort-variants-by-gene db-name "TERT")})

^:kindly/hide-code
(kind/table {:row-maps (pv/cohort-variants-by-gene db-name "TP53")})

;; ## Genes mutated in the most participants
^:kindly/hide-code
(kind/table {:row-maps (pv/variant-cohort-counts db-name ["GNAQ" "GNA11" "BAP1" "SF3B1" "EIF1AX" "TP53"])})
