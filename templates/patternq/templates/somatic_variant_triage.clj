^{:kindly/hide-code true
  :clay {:quarto {:title "Somatic Variant Triage"
                  :format {:html {:embed-resources true :theme :cosmo}}}}}
(ns patternq.templates.somatic-variant-triage
  "Generic participant somatic variant triage (from the unify-central variant
  forensics template): samples & timepoints, protein-altering variants,
  candidate clonal trunk variants, branch variants by site. Set `db-name`
  and `participant-id` below, then render with Clay:

    clojure -M:examples -e \"(require '[scicloj.clay.v2.api :as clay]) (clay/make! {:source-path \\\"templates/patternq/templates/somatic_variant_triage.clj\\\"})\""
  (:require [patternq.dataset :as pqd]
            [patternq.http :as http]
            [patternq.plot :as plot]
            [patternq.reference :as ref]
            [patternq.variants :as pv]
            [scicloj.kindly.v4.kind :as kind]))

;; # Somatic Variant Triage
;;
;; A starting point for examining a participant's variants to identify
;; plausible clonal trunk variants and primary drivers. Focused on coding
;; (protein-altering) mutations.
;;
;; _Set the dataset and participant here:_
(def db-name (http/resolve-db "H37001"))
(def participant-id "H37001")

;; ## Samples and timepoints
^:kindly/hide-code
(kind/table
  {:row-maps (->> (pqd/samples db-name)
                  (filter #(= participant-id (:subject-id %)))
                  (sort-by :timepoint-relative-order)
                  (mapv #(select-keys % [:sample-id :sample-freetext-anatomic-site :timepoint-id
                                         :timepoint-relative-order :sample-metastasis])))})

;; ## Protein-altering variants
;;
;; Including HGVSp restricts to variants that alter proteins; low impact
;; (synonymous) variants are excluded. Recurrent artifacts can be excluded
;; with `:exclude-genes`.
(def all-variants (pv/participant-variants db-name participant-id {:exclude-genes #{"OR8U1"}}))

^:kindly/hide-code
(kind/hiccup [:p (format "Variant measurements: %d" (count all-variants))])

^:kindly/hide-code
(kind/plotly (plot/vaf-histogram all-variants))

;; ## Candidate trunk variants
;;
;; Variants present in several samples/timepoints are candidate clonal
;; trunk mutations; sorted by sample count, then summed VAF.
(def trunk (pv/trunk-candidates all-variants {:min-samples 2}))

^:kindly/hide-code
(kind/table {:row-maps (mapv #(update % :vafs (fn [vs] (mapv (fn [v] (format "%.2f" (double v))) vs))) trunk)})

;; ### UniProt links for trunk candidates
^:kindly/hide-code
(kind/table
  {:row-maps (vec (for [{:keys [hgnc-symbol HGVSp]} trunk
                        :let [[_ uniprot acc] (first (ref/hgnc->uniprot db-name hgnc-symbol))]]
                    {:hgnc hgnc-symbol :HGVSp HGVSp
                     :uniprot (if acc (kind/hiccup [:a {:href (ref/uniprot-link acc)} uniprot]) "")}))})

;; ## Branch variants
;;
;; Variants not in the trunk, grouped by gene, with the sites and
;; timepoints they were seen at.
^:kindly/hide-code
(kind/table {:row-maps (mapv #(update % :HGVSp (partial clojure.string/join "; "))
                             (pv/branch-variants all-variants trunk))})
