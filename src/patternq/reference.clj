(ns patternq.reference
  "Reference data queries. Reference entities (genes, gene products,
  proteins, epitopes, cell types, ...) are included in each dataset
  database, so these take a db like every other query."
  (:require [clojure.set]
            [clojure.string :as str]
            [patternq.dataset :as pqd]
            [patternq.db :as pdb]
            [patternq.results :as res]))

(set! *warn-on-reflection* true)

(defn- names-of [db-or-name attr]
  (let [[r db db-name] (pqd/run db-or-name {:find '[[?v ...]] :where [['_ attr '?v]]} [])]
    (pdb/with-provenance (vec (sort r)) db db-name)))

(defn gene-symbols "HGNC symbols." [db-or-name] (names-of db-or-name :gene/hgnc-symbol))
(defn gdc-anatomic-sites [db-or-name] (names-of db-or-name :gdc-anatomic-site/name))
(defn epitopes [db-or-name] (names-of db-or-name :epitope/id))
(defn cell-types [db-or-name] (names-of db-or-name :cell-type/co-name))
(defn meddra-diseases [db-or-name] (names-of db-or-name :meddra-disease/preferred-name))
(defn drugs [db-or-name] (names-of db-or-name :drug/preferred-name))

(def genes-query
  '[:find (pull ?g [:gene/hgnc-symbol :gene/hgnc-id :gene/hgnc-name :gene/ensembl-id
                    :gene/previous-hgnc-symbols :gene/alias-hgnc-symbols
                    {:gene/hgnc-locus-group [:db/ident]}])
    :where [?g :gene/hgnc-symbol]])

(defn genes
  "Genes: :gene-hgnc-symbol, :gene-hgnc-name, ids, :gene-previous-hgnc-symbols
  and :gene-alias-hgnc-symbols (vectors)."
  [db-or-name]
  (pqd/pull-rows db-or-name genes-query []))

(def gene-products-query
  '[:find ?gene-product-id ?hgnc-symbol
    :where
    [?gp :gene-product/id ?gene-product-id]
    [?gp :gene-product/gene ?g]
    [?g :gene/hgnc-symbol ?hgnc-symbol]])

(defn gene-products
  ":gene-product-id, :hgnc-symbol"
  [db-or-name]
  (let [[r db db-name] (pqd/run db-or-name gene-products-query [])]
    (pdb/with-provenance (res/rows->maps [:gene-product-id :hgnc-symbol] r) db db-name)))

(defn gene-coordinates
  ":hgnc-symbol, :assembly, :contig, :strand, :start, :end.

  opts: :genes"
  ([db-or-name] (gene-coordinates db-or-name {}))
  ([db-or-name {:keys [genes]}]
   (let [q (vec (concat '[:find ?hgnc-symbol ?assembly ?contig ?strand ?start ?end :in $]
                        (when genes '[[?hgnc-symbol ...]])
                        '[:where
                          [?g :gene/hgnc-symbol ?hgnc-symbol]
                          [?g :gene/genomic-coordinates ?gc]
                          [?gc :genomic-coordinate/assembly ?a]
                          [?a :db/ident ?assembly]
                          [?gc :genomic-coordinate/contig ?contig]
                          [?gc :genomic-coordinate/strand ?strand]
                          [?gc :genomic-coordinate/start ?start]
                          [?gc :genomic-coordinate/end ?end]]))
         [r db db-name] (pqd/run db-or-name q (if genes [(vec genes)] []))]
     (pdb/with-provenance
       (mapv #(update (zipmap [:hgnc-symbol :assembly :contig :strand :start :end] %) :assembly res/ident-name) r)
       db db-name))))

(defn variant-annotations
  "Variant reference entities: :variant-id, :hgnc-symbol, :HGVSp, :HGVSc,
  :impact, :classification, :type, :so-consequences, alleles, ids.
  Cardinality-many strings are joined with \"; \".

  opts: :variant-ids, :genes"
  ([db-or-name] (variant-annotations db-or-name {}))
  ([db-or-name {:keys [variant-ids genes]}]
   (let [q (vec (concat '[:find (pull ?v [:variant/id :variant/HGVSp :variant/HGVSc :variant/ref-allele
                                           :variant/alt-allele :variant/coordinate-string :variant/dbSNP
                                           :variant/max-af :variant/external-ids
                                           {:variant/gene [:gene/hgnc-symbol]}
                                           {:variant/impact [:db/ident]}
                                           {:variant/classification [:db/ident]}
                                           {:variant/type [:db/ident]}
                                           {:variant/so-consequences [:so-sequence-feature/name]}])
                          :in $]
                        (when variant-ids '[[?variant-id ...]])
                        (when genes '[[?gene ...]])
                        '[:where [?v :variant/id ?variant-id]]
                        (when genes '[[?g :gene/hgnc-symbol ?gene] [?v :variant/gene ?g]])))
         args (cond-> [] variant-ids (conj (vec variant-ids)) genes (conj (vec genes)))
         rows (pqd/pull-rows db-or-name q args)]
     (with-meta
       (mapv (fn [r]
               (-> r
                   (res/rename-prefix "variant-")
                   (clojure.set/rename-keys {:id :variant-id :gene-hgnc-symbol :hgnc-symbol})
                   (update-vals (fn [v] (if (and (coll? v) (every? string? v)) (res/join-many v) v)))))
             rows)
       (meta rows)))))

(defn cnvs
  "CNV reference entities: :cnv-id, :genomic-coordinate-contig/-start/-end,
  :cnv-genes (vector of HGNC symbols)."
  [db-or-name]
  (pqd/pull-rows db-or-name
                 '[:find (pull ?c [:cnv/id
                                   {:cnv/genomic-coordinates [:genomic-coordinate/contig
                                                              :genomic-coordinate/start
                                                              :genomic-coordinate/end]}
                                   {:cnv/genes [:gene/hgnc-symbol]}])
                   :where [?c :cnv/id]]
                 []))

(defn proteins
  ":protein-preferred-name, :protein-uniprot-name, :protein-uniprot-accessions,
  :gene-hgnc-symbol"
  [db-or-name]
  (pqd/pull-rows db-or-name
                 '[:find (pull ?p [:protein/preferred-name :protein/uniprot-name :protein/uniprot-accessions
                                   {:protein/gene [:gene/hgnc-symbol]}])
                   :where [?p :protein/preferred-name]]
                 []))

(defn map-gene-symbols
  "Map gene symbols (current, previous or alias; case-insensitive) to current
  HGNC symbols. Returns a map input-symbol -> HGNC symbol (nil if unmapped).
  Current symbols take precedence over previous/alias ones.

  `all-genes` is output of `genes`; fetched from `db-or-name` when nil."
  ([db-or-name symbols] (map-gene-symbols db-or-name symbols nil))
  ([db-or-name symbols all-genes]
   (let [all-genes (or all-genes (genes db-or-name))
         up str/upper-case
         current (into {} (map (fn [g] [(up (:gene-hgnc-symbol g)) (:gene-hgnc-symbol g)])) all-genes)
         others (into {} (for [g all-genes
                               k [:gene-previous-hgnc-symbols :gene-alias-hgnc-symbols]
                               s (get g k)]
                           [(up s) (:gene-hgnc-symbol g)]))
         lookup (merge others current)]
     (into {} (map (fn [s] [s (get lookup (up s))])) symbols))))

(defn resolve-aliases
  "Direct Datalog alias resolution (from unify-central analysis): pairs of
  [input-symbol hgnc-symbol] for inputs matching a previous/alias symbol or
  name, or an Ensembl id."
  [db-or-name symbols]
  (pdb/q '[:find ?gene ?hgnc
         :in $ [?gene ...]
         :where
         (or-join [?g ?gene]
                  [?g :gene/previous-hgnc-symbols ?gene]
                  [?g :gene/previous-hgnc-names ?gene]
                  [?g :gene/ensembl-id ?gene]
                  [?g :gene/hgnc-prev-symbols ?gene]
                  [?g :gene/alias-hgnc-symbols ?gene]
                  [?g :gene/alias-hgnc-names ?gene])
         [?g :gene/hgnc-symbol ?hgnc]]
       (pdb/as-db db-or-name) (vec symbols)))

(defn remap-gene-names
  "Map for renaming gene names to HGNC symbols: valid symbols map to
  themselves, aliases to their HGNC symbol (unmappable names are absent)."
  [db-or-name gene-names]
  (let [db (pdb/as-db db-or-name)
        confirmed (pdb/q '[:find [?hgnc ...] :in $ [?hgnc ...] :where [_ :gene/hgnc-symbol ?hgnc]]
                       db (vec gene-names))]
    (merge (zipmap confirmed confirmed)
           (into {} (resolve-aliases db gene-names)))))

(defn hgnc->uniprot
  "[hgnc uniprot-name accession] tuples for a gene's proteins."
  [db-or-name hgnc]
  (pdb/q '[:find ?hgnc ?uniprot ?acc
         :in $ ?hgnc
         :where
         [?g :gene/hgnc-symbol ?hgnc]
         [?p :protein/gene ?g]
         [?p :protein/uniprot-accessions ?acc]
         [?p :protein/uniprot-name ?uniprot]]
       (pdb/as-db db-or-name) hgnc))

(defn uniprot-link [uniprot-accession]
  (str "https://www.uniprot.org/uniprotkb/"
       (java.net.URLEncoder/encode (str uniprot-accession) "UTF-8") "/entry"))

(defn gene-size
  "Genomic span (bp, end - start) of a gene."
  [db-or-name hgnc]
  (when-let [[start stop] (first (pdb/q '[:find ?start ?stop
                                        :in $ ?hgnc
                                        :where
                                        [?g :gene/hgnc-symbol ?hgnc]
                                        [?g :gene/genomic-coordinates ?gc]
                                        [?gc :genomic-coordinate/start ?start]
                                        [?gc :genomic-coordinate/end ?stop]]
                                      (pdb/as-db db-or-name) hgnc))]
    (- (long stop) (long start))))
