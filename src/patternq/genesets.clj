(ns patternq.genesets
  "Gene sets used in the variant forensics reports, as data. Curated sets
  are inline; MSigDB hallmark sets and the structured antibody-target
  reference are resources under patternq/genesets/."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.string :as str]))

(defn- resource-list [fname]
  (->> (str/split (str/trim (slurp (io/resource (str "patternq/genesets/" fname)))) #"[,\s]+")
       (remove str/blank?)
       vec))

(def hallmark-sets
  {"hallmark-apoptosis" "hallmark_apoptosis.txt"
   "hallmark-dna-repair" "hallmark_dna_repair.txt"
   "hallmark-hypoxia" "hallmark_hypoxia.txt"
   "hallmark-inflammatory" "hallmark_inflammatory.txt"
   "antibody-therapy-targets" "antibody_therapy_targets.txt"})

(defn structured-reference
  "Antibody-target groups: {\"Antibody Targets (...)\" [genes]}."
  []
  (edn/read-string (slurp (io/resource "patternq/genesets/structured-reference.edn"))))

(def breast-cancer
  ["PTEN" "BRCA1" "ERBB2" "TP53" "PALB2" "CHEK2" "STK11" "BRCA2"
   "BRIP1" "ATM" "MRE11" "RAD50" "NBN" "XRCC3" "CYP19A1"])

;; germline multi-cancer panel
(def germline-multi-cancer
  ["AIP" "ALK" "APC" "ATM" "AXIN2" "BAP1" "BARD1" "BLM" "BMPR1A" "BRCA1" "BRCA2" "BRIP1" "CDC73" "CDH1" "CDK4" "CDKN1B" "CDKN2A" "CHEK2" "CTNNA1" "DICER1" "EGFR" "EPCAM" "FH" "FLCN" "GREM1" "HOXB13" "KIT" "LZTR1" "MAX" "MBD4" "MEN1" "MET" "MITF" "MLH1" "MSH2" "MSH3" "MSH6" "MUTYH" "NF1" "NF2" "NTHL1" "PALB2" "PDGFRA" "PMS2" "POLD1" "POLE" "POT1" "PRKAR1A" "PTCH1" "PTEN" "RAD51C" "RAD51D" "RB1" "RET" "SDHA" "SDHAF2" "SDHB" "SDHC" "SDHD" "SMAD4" "SMARCA4" "SMARCB1" "SMARCE1" "STK11" "SUFU" "TMEM127" "TP53" "TSC1" "TSC2" "VHL"])

;; melanoma phenotype markers (loosely: high-variance genes in TCGA-SKCM)
(def melanoma-proliferative ["MITF" "TYPR1" "DCT" "MLANA" "SOX10" "PAX3" "PRAME" "ZEB2"])
(def melanoma-mesenchymal-invasive ["AXL" "MMP2" "MMP16" "VCAN" "LUM" "FSCN1" "COL1A1" "COL1A2" "COL3A1" "FN1" "THBS1" "ZEB1"])
(def melanoma-neural-crest ["SOX9" "NES" "MCAM"])

(def glutamine-vulnerable ["GLS" "MYC" "JUN" "HIF1A" "NFE2L2" "VIM" "SLC1A5" "SLC7A11"])
(def dna-repair-agents ["WEE1" "ATR" "POLQ" "PRKDC" "CHEK1"])

(defn housekeeping
  "HSIAO_HOUSEKEEPING_GENES (MSigDB): control set of genes expressed in most cells."
  [] (resource-list "housekeeping.txt"))

(defn adult-kidney
  "LAKE_ADULT_KIDNEY_C8_DECENDING_THIN_LIMB (MSigDB): normal tissue control."
  [] (resource-list "adult_kidney.txt"))

(defn lee-ncc-up "Lee neural crest cell up-regulated genes." [] (resource-list "lee_ncc_up.txt"))

(def ^:private named
  {"breast-cancer" (constantly breast-cancer)
   "germline-multi-cancer" (constantly germline-multi-cancer)
   "melanoma-proliferative" (constantly melanoma-proliferative)
   "melanoma-mesenchymal-invasive" (constantly melanoma-mesenchymal-invasive)
   "melanoma-neural-crest" (constantly melanoma-neural-crest)
   "glutamine-vulnerable" (constantly glutamine-vulnerable)
   "dna-repair-agents" (constantly dna-repair-agents)
   "housekeeping" housekeeping
   "adult-kidney" adult-kidney
   "lee-ncc-up" lee-ncc-up})

(defn geneset-names [] (sort (concat (keys named) (keys hallmark-sets))))

(defn geneset
  "Genes of a named set (see geneset-names)."
  [set-name]
  (if-let [f (named set-name)]
    (f)
    (if-let [fname (hallmark-sets set-name)]
      (resource-list fname)
      (throw (ex-info (str "Unknown gene set " set-name) {:geneset set-name :known (geneset-names)})))))
