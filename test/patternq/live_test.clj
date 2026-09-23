(ns patternq.live-test
  "Live tests against dev dataset databases through the peer (needs AWS
  credentials with read access to the dev storage, and PATTERNQ_API_KEY for
  the HTTP listing API). Test datasets per SCOPE.md 5.1."
  (:require [clojure.test :refer [deftest is testing]]
            [patternq.clinical :as clinical]
            [patternq.context :as ctx]
            [patternq.dataset :as pqd]
            [patternq.db :as pdb]
            [patternq.expression :as expr]
            [patternq.genesets :as gs]
            [patternq.http :as http]
            [patternq.plot :as plot]
            [patternq.reference :as ref]
            [patternq.results :as res]
            [patternq.variants :as pv]))

(def db-name
  (memoize (fn [dataset] (http/resolve-db dataset))))

(deftest ^:live h37001-sanity
  (let [db (db-name "H37001")]
    (testing "samples & subjects"
      (is (= 6 (count (pqd/samples db))))
      (is (= ["H37001"] (map :subject-id (pqd/subjects db))))
      (is (= db (:db (:patternq/provenance (meta (pqd/samples db)))))))
    (testing "variants"
      (let [vs (pqd/variants db)]
        (is (< 100 (count vs)))
        (is (some #(= "BAP1" (:hgnc-symbol %)) vs))
        (is (= 1 (count (:data (plot/vaf-histogram vs {:samples ["H37001-001"]})))))))
    (testing "gene expression"
      (let [gx (pqd/gene-expression db {:genes ["BAP1" "GNAQ"] :measurement :rsem-normalized-count})]
        (is (= #{"BAP1" "GNAQ"} (set (map :hgnc-symbol gx))))))
    (testing "participant variant analysis"
      (let [pvs (pv/participant-variants db "H37001" {:exclude-genes #{"OR8U1"}})
            trunk (pv/trunk-candidates pvs {:min-samples 2})]
        (is (seq pvs))
        (is (not-any? #(= "OR8U1" (:hgnc-symbol %)) pvs))
        (is (some #(= "BAP1" (:hgnc-symbol %)) trunk))
        (is (every? #(>= (:sample-count %) 2) trunk))
        (is (seq (pv/find-variants db "H37001")))))))

(deftest ^:live tcga-uvm-sanity
  (let [db (db-name "tcga-uvm")]
    (is (= 80 (count (pqd/subjects db))))
    (is (= #{"WES" "RNA-seq"} (set (map :assay-technology (pqd/dataset-summary db)))))
    (is (< 100 (count (pqd/gene-expression db {:genes ["BAP1" "PRAME"]}))))
    (is (= 80 (count (pqd/samples db))))
    (is (pos? (pv/variant-patient-count db "GNAQ")))))

(deftest ^:live prince-catalog
  (let [db (db-name "prince-2022")]
    (testing "measurement types and generic measurements"
      (let [mt (pqd/measurement-types db "PICI CyTOF Immune Profiling")
            n-pp (:count (first (filter #(= "percent-of-parent" (:attribute %)) mt)))
            pp (pqd/measurements db :percent-of-parent {:measurement-set "PICI CyTOF Immune Profiling"})
            mcv (pqd/measurements db :median-channel-value {:measurement-set "PICI CyTOF Immune Profiling"})]
        (is (every? (set (map :attribute mt)) ["percent-of-parent" "median-channel-value" "cell-population" "epitope"]))
        (is (= n-pp (count pp)))
        (is (= #{:sample-id :measurement-set :cell-population :value} (set (keys (first pp)))))
        (is (contains? (first mcv) :epitope-id))
        (is (seq (:col-names (res/->matrix (pqd/measurements db :olink-npx {:measurement-set "PICI Olink Proteomics"})
                                            :epitope-id))))))
    (testing "clinical"
      (let [oc (clinical/subject-outcomes db)]
        (is (= (count (pqd/subjects db)) (count oc)))
        (is (some :bor oc))
        (is (some :os oc))
        (is (seq (:data (plot/survival oc {:group :bor}))))))
    (testing "timepoints ordered"
      (let [tp (pqd/timepoints db)]
        (is (= (map :timepoint-relative-order tp) (sort (map :timepoint-relative-order tp))))))
    (testing "context + plots"
      (let [rows (ctx/add-sample-context (take 20 (pqd/measurements db :percent-of-parent {:measurement-set "PICI CyTOF Immune Profiling"}))
                                         db {:outcomes? true})]
        (is (every? :subject-id rows))
        (is (seq (:data (plot/by-timepoint rows {:group :bor}))))))
    (testing "reference"
      (is (every? #(= "KRAS" (:hgnc-symbol %)) (ref/variant-annotations db {:genes ["KRAS"]})))
      (is (= {"p53" "TP53" "HER2" "ERBB2"} (ref/map-gene-symbols db ["p53" "HER2"])))
      (is (= 2 (count (ref/gene-coordinates db {:genes ["KRAS" "TP53"]})))))
    (testing "sample assays + overview"
      (let [sa (pqd/sample-assays db)]
        (is (seq sa))
        (is (= "heatmap" (get-in (plot/sample-overview sa) [:data 0 :type])))))
    (is (seq (clinical/clinical-observations db {:obs-type :os})))
    (is (= "PICI Survival Outcomes (OS/PFS/BOR)" (:clinical-observation-set-name (first (clinical/clinical-observation-sets db)))))
    (is (= 16 (count (pqd/cell-populations db "PICI CyTOF Immune Profiling"))))))

(deftest ^:live h37004-cnv-and-more
  (let [db (db-name "H37004")]
    (is (thrown-with-msg? clojure.lang.ExceptionInfo #"subset" (pqd/cnv-segments db {})))
    (is (thrown-with-msg? clojure.lang.ExceptionInfo #"subset" (pqd/cnv-gene-calls db {})))
    (let [cs (pqd/cnv-segments db {:genes ["TP53"]})
          cs2 (pqd/cnv-segments db {:subjects ["H37004"]})]
      (is (seq cs))
      (is (every? #(= "TP53" (:hgnc-symbol %)) cs))
      (is (every? #(contains? % :contig) cs))
      (is (> (count cs2) (count cs))))
    (is (seq (pqd/isoforms db "TP53")))
    (is (seq (clinical/clinical-interventions db)))
    (is (= 2 (count (pqd/measurement-matrices db))))
    (is (seq (pqd/sgbs db "Taxa Abundance")))
    (is (seq (clinical/clinical-timeline db "H37004")))))

(deftest ^:live painter-spotty-cnv
  (let [db (db-name "painter-2025-angiosarc")
        calls (pqd/cnv-gene-calls db {:genes ["MYC" "CDKN2A"]})]
    (is (every? #(contains? #{-2 -1 0 1 2} (:value %)) calls))
    (is (< (count (distinct (map :sample-id calls))) (count (pqd/samples db))))))

(deftest ^:live cohort-comparison
  (let [db (db-name "H37001")
        cohort (db-name "tcga-uvm")
        {:keys [shared? participant-attr]} (expr/resolve-genex-attrs db cohort)]
    (is shared?)
    (let [zs (expr/cohort-zscores db cohort "H37001-001" ["BAP1" "GNAQ" "PRAME"] {:measurement participant-attr})]
      (is (= 3 (count zs)))
      (is (every? #(number? (:z %)) zs)))
    (is (= 3 (count (expr/pathway-gene-percentiles db cohort "H37001-001" ["BAP1" "GNAQ" "PRAME"] {:measurement participant-attr}))))
    (let [spec (expr/examine-geneset db [cohort] (take 5 (gs/geneset "melanoma-proliferative")) ["H37001-001"]
                                     {:measurement participant-attr :style :violin})]
      (is (some #(= "violin" (:type %)) (:data spec))))
    (is (= #{db cohort} (set (map :db (pdb/across-dbs pqd/dataset-summary [db cohort])))))))

(deftest ^:live prince-survival-matches-r
  ;; reference p-values from the R library (patternq::survival_by_median on prince-2022)
  (let [db (db-name "prince-2022")
        outcomes (clinical/subject-outcomes db)
        x50 (ctx/add-sample-context (pqd/measurements db :percent-of-parent {:measurement-set "PICI X50 Immune Profiling"})
                                    db {:subjects? false})
        base (fn [pop] (->> x50
                            (filter #(and (= pop (:cell-population %)) (= "C1D1" (:timepoint-id %))))
                            (group-by :subject-id)
                            (into {} (map (fn [[s rs]] [s (/ (reduce + (map :value rs)) (count rs))])))))
        p (fn [pop] (:p ((requiring-resolve 'patternq.survival/logrank)
                          ((requiring-resolve 'patternq.survival/survival-by-median) (base pop) outcomes))))]
    (is (< (Math/abs (- 0.01151863 (p "PD-1+Tbet+ non-Naive CD4+ T cells (% of CD4 not naive T cells)"))) 1e-6))
    (is (< (Math/abs (- 0.7055291 (p "CD38+ effector memory CD8+ T cells (% effector memory CD8 T cells)"))) 1e-6))))
