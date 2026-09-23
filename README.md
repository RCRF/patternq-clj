# patternq (Clojure)

Read-only query and analysis library for the Pattern Data Commons. By
default it queries the Pattern Data Commons query service over HTTP (like the
R, Python and Julia libraries); internal users can switch to the Datomic
peer. The patternq family: [patternq](https://github.com/RCRF/patternq) (Python), [patternq-r](https://github.com/RCRF/patternq-r) (R), [patternq-clj](https://github.com/RCRF/patternq-clj) (Clojure) and [PatternQ.jl](https://github.com/RCRF/PatternQ.jl) (Julia) share one function catalog, the same result columns and the same plots.

- Every dataset is its own Datomic database; functions take a database
  **name** (or a handle from `patternq.db/db`) as their first argument.
- Strictly read-only: no transactions, no database creation/deletion, never
  the admin database (`test/patternq/read_only_test.clj` enforces this).
- Results are vectors of maps with unqualified kebab-case keys (`:sample-id`
  here is `sample_id` in R/Python), with provenance in metadata:
  `(:patternq/provenance (meta rows))` -> `{:db :basis-t :timestamp}`.
- Plots are plain plotly spec maps (`{:data [...] :layout {...}}`), no
  plotting dependency. Render with Kindly/Clay (`kind/plotly`) or plotly.js.

## Setup

deps.edn coordinate (git dependency; use the latest commit sha):

```clojure
{:deps {io.github.RCRF/patternq-clj {:git/url "https://github.com/RCRF/patternq-clj"
                                      :git/sha "<commit sha>"}}}
```

Set **`PATTERNQ_API_KEY`** (from the user settings page of the Pattern Data
Commons dashboard) and optionally **`PATTERNQ_ENDPOINT`** (default
`https://data-commons.rcrf-dev.org`), or call `patternq.http/set-token!` /
`set-endpoint!`. That is all the default HTTP transport needs: no Datomic
jars and no cloud credentials.

### Transports

All queries go through `patternq.db/q` (Datomic Datalog written as Clojure
data, the same shape as `datomic.api/q`), so every function works over both
transports:

| transport | select | needs |
|---|---|---|
| `:http` (default) | nothing, or `PATTERNQ_TRANSPORT=http` | `PATTERNQ_API_KEY` |
| `:peer` | `PATTERNQ_TRANSPORT=peer` or `(patternq.db/set-transport! :peer)` | the `:peer` deps alias (Datomic peer), `PATTERNQ_DATOMIC_URI` (the storage base URI; database names are appended) or `patternq.db/set-base-uri!`, and read access to that storage |

Over HTTP the query is converted client-side to the JSON form the query
service parses (the service accepts no EDN), results come back through the
service's S3 result cache by default (bind `patternq.http/*cache*` to false
for inline, uncached results, `*refresh-cache*` to recompute, `*timeout-ms*`
for the query timeout), and are converted back so they look exactly like peer
results (idents as keywords, find specs `[?x ...]`, `?x .` and `[?a ?b]`
reshaped). Requesting the peer without the peer on the classpath, or without
`PATTERNQ_DATOMIC_URI`, is an error rather than a silent fallback.

To use the peer, add the alias to your command (or copy its deps):

```sh
PATTERNQ_TRANSPORT=peer PATTERNQ_DATOMIC_URI=<storage base uri> clojure -M:peer ...
```

Differences between transports: the query service restricts expression
functions to a whitelist (comparisons, arithmetic, `str`, `re-find`,
`ground`, `get-else`, `missing?`, ...) and cannot take string literals that
start with `:` or `?` (pass those as arguments); `:db.type/float` values come
back as doubles over HTTP (floats from the peer); relations are vectors over
HTTP (sets from the peer); and the basis t in provenance is that of the
latest query through a handle.

JDK 21+ prints netty warnings with the peer unless you add
`--enable-native-access=ALL-UNNAMED` (the `:peer` and `:test` aliases do).

## Quick start

```clojure
(require '[patternq.http :as http]
         '[patternq.dataset :as pqd]
         '[patternq.clinical :as clinical]
         '[patternq.plot :as plot])

(def db (http/resolve-db "prince-2022"))      ; dataset name -> current db name
(pqd/dataset-summary db)                       ; assays and measurement sets
(pqd/measurement-types db "PICI CyTOF Immune Profiling")
(pqd/measurements db :percent-of-parent {:measurement-set "PICI CyTOF Immune Profiling"})
(plot/survival (clinical/subject-outcomes db) {:group :bor})

;; CNV queries require a subset (genes, samples or subjects)
(pqd/cnv-segments (http/resolve-db "H37004") {:genes ["TP53"]})
```

## Namespaces

| namespace | contents |
|---|---|
| `patternq.db` | transport selection, `db` handles, `q` (the single query entry point), `pull`, `pull-many`, `pull-by`, `pull-map`, `across-dbs`, provenance |
| `patternq.http` | HTTP transport (`q`, `query->wire`, cache options), `list-datasets`, `resolve-db`, `measurement-matrix`, `set-token!`, `set-endpoint!` |
| `patternq.dataset` | samples, subjects, timepoints, dataset-summary/info, schema-info, measurement-sets/-types/-set-attributes, `measurements` (targets auto-detected), sample-assays, measurement-matrices, variants, gene-expression, isoforms, cnv-segments, cnv-gene-calls, cell-populations/tcrs/otus/sgbs; `*-query` companions |
| `patternq.clinical` | clinical observation sets/observations, subject-outcomes (BOR/PFS/OS), adverse-events, clinical-interventions, clinical-timeline + timeline traces (series, treatments, events) |
| `patternq.reference` | gene-symbols, genes, gene-products, gene-coordinates, variant-annotations, cnvs, proteins, epitopes, cell-types, meddra-diseases, drugs, anatomic sites, map-gene-symbols, resolve-aliases, remap-gene-names, hgnc->uniprot, gene-size |
| `patternq.context` | add-sample/subject/variant/cnv-context, split-by-measurement-set, deduplicate-taxonomy, aggregate-taxa |
| `patternq.results` | flatten-pull, resolve-enum-refs, ->matrix / ->long / select-targets, join-many |
| `patternq.plot` | theme; vaf-histogram, gene-expression, sample-overview, by-timepoint (`:value`, `:timepoints`, `:levels`, trajectories colored by group), by-group, survival (log-rank p annotation, `:levels`), heatmap (clustered, `:col-groups` annotation strip), mutation-landscape, vs-cohort, zscores, ma, fold-change, genex-vs-cohort (older box / horizontal / violin style), samples-vs-histogram, ma-plot |
| `patternq.survival` | logrank-test (Mantel-Haenszel, k groups; = R `logrank_test` = `survival::survdiff`), chisq-upper, median-split, survival-status (landmark), survival-by-median (log-rank in metadata, `logrank`), change-from-baseline (`:log2-ratio` / `:difference` / `:ratio`) |
| `patternq.quant` | mean, variance, sd, median, quantile, percentile-rank, insertion-index, z-scores, log-fold-change, ma-values, l2/cosine distance, top-varying, ssgsea-score, ->histogram, kaplan-meier, hclust-order, scale-rows |
| `patternq.expression` | compare-to-cohort, top-by-zscore, compare-samples, sample-expression (R/Python parity), geneset-expression, sample->gene-expression, cohort (batch) expression, available/resolve-genex-attrs, pathway-gene-percentiles, rank-expression, cohort-zscores, nearest-neighbors, examine-geneset, isoform-breakdown |
| `patternq.variants` | participant-variants (missense), modifier/synonymous variants, vaf-histogram-data, trunk-candidates, branch-variants, variants-by-impact, find-variants (with CNV), variant-samples, cohort-variants-by-gene, mutated-genes, variant-patient-count, variant-cohort-counts, variant-genes-in-cohorts, to-igv-coords |
| `patternq.genesets` | named gene sets (hallmark, housekeeping, germline, melanoma phenotypes, antibody targets, ...) |

## Examples and report templates

`examples/` has the same notebooks as every language in the family, section
for section, as Clay namespaces:

- `tutorial.clj`: the basic workflow (connect, explore, Datalog, variants and
  expression, measurements and context, survival, cohort comparison, MA plot)
- `prince_demo.clj`: PRINCE trial clinical and immune correlates

Render with Clay (over HTTP by default; needs `quarto` on the PATH and
`PATTERNQ_API_KEY`):

```sh
clojure -M:examples -e "(require '[scicloj.clay.v2.api :as clay]) (clay/make! {:source-path \"examples/tutorial.clj\" :format [:quarto :html]})"
```

`templates/patternq/templates/` holds the earlier generic report templates
ported from the unify-central variant forensics analyses
(`somatic_variant_triage`, `cohort_comparison`, `cross_cohort_variants`).
They are Clojure-only starting points for patient reports, not cross-language
examples; render one with
`clojure -M:examples -e "(require '[scicloj.clay.v2.api :as clay]) (clay/make! {:source-path \"templates/patternq/templates/somatic_variant_triage.clj\" :format [:quarto :html]})"`.

## Tests

```sh
clojure -X:test                     # unit + live tests over HTTP (live tests need PATTERNQ_API_KEY)
PATTERNQ_TRANSPORT=peer clojure -X:test:peer   # the same over the peer (needs PATTERNQ_DATOMIC_URI)
clojure -X:test :nses '[patternq.quant-test patternq.results-plot-test patternq.read-only-test]'   # offline only
```

## Contributing

patternq is one library in four languages: [patternq](https://github.com/RCRF/patternq)
(Python), [patternq-r](https://github.com/RCRF/patternq-r) (R),
[patternq-clj](https://github.com/RCRF/patternq-clj) (Clojure) and
[PatternQ.jl](https://github.com/RCRF/PatternQ.jl) (Julia). All four are generated
from a common source and published here, so this repository does not accept pull
requests.

Please report bugs and feature requests as
[issues](https://github.com/RCRF/patternq-clj/issues). Code is welcome in an issue: a
minimal example (the call or query, the dataset, what you expected and what you
got), or a proposed change as a snippet, is the most useful way to suggest one.
