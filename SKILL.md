---
name: patternq-clj
description: Query the Pattern Data Commons (PDC, RCRF datasets) from Clojure with patternq-clj. Use when asked about Pattern Data Commons datasets, subjects, samples, measurements, variants, copy number, gene expression, or clinical outcomes, or when Clojure code requires patternq namespaces.
---

# Pattern Data Commons from Clojure (patternq-clj)

patternq queries the Pattern Data Commons query service over HTTPS. You act on
behalf of the user whose API key is configured. The key can see only that
user's datasets, the data is de-identified, and every call is read-only.

## Setup

deps.edn (git dependency; use the latest commit sha):

```clojure
{:deps {io.github.RCRF/patternq-clj {:git/url "https://github.com/RCRF/patternq-clj"
                                      :git/sha "<commit sha>"}}}
```

```sh
export PATTERNQ_ENDPOINT=https://data-commons.rcrf-dev.org   # the default
export PATTERNQ_API_KEY=...   # the user's key, from the dashboard's user settings
```

The user supplies the key. Read it from the environment and never print it,
log it, write it into code or notebooks, or commit it. The HTTP transport is
the default and needs nothing else; ignore the `:peer` transport.

```clojure
(require '[patternq.http :as http]
         '[patternq.db :as pdb]
         '[patternq.dataset :as pqd]
         '[patternq.clinical :as clinical]
         '[patternq.context :as ctx])
```

## Every dataset is its own database

Every function takes a database name as its first argument. Resolve the
dataset name to its current database first. Database names change when a
dataset is re-imported, so don't hard-code them.

```clojure
(http/list-datasets)                     ; maps with :dataset, :db, counts, :assays, :tags
(def db (http/resolve-db "tcga-uvm"))    ; dataset name -> current database name
```

To cover several datasets, use `(pdb/across-dbs f db-names)`, which runs
`f` on each database and adds `:db` to every row.

## The calls you need most

| task | call |
|---|---|
| what's in a dataset | `(pqd/dataset-summary db)`, `(pqd/measurement-sets db)`, `(pqd/measurement-types db set-name)` |
| subjects (demographics, disease) | `(pqd/subjects db)` |
| samples (subject, timepoint, type) | `(pqd/samples db)`, `(pqd/timepoints db)`, `(pqd/sample-assays db)` |
| measurement values | `(pqd/measurements db :tpm {:measurement-set set-name :samples [...]})` |
| somatic variants | `(pqd/variants db {:genes [...] :samples [...]})` |
| gene expression | `(pqd/gene-expression db {:genes [...]})` |
| copy number | `(pqd/cnv-gene-calls db {:genes [...]})`, `(pqd/cnv-segments db {:genes [...]})` (genes, samples, or subjects required) |
| clinical | `(clinical/subject-outcomes db)` (BOR, PFS, OS), `(clinical/clinical-observation-sets db)`, `(clinical/clinical-observations db {:obs-type :os})` |
| sample or subject context | `(ctx/add-sample-context rows db)`, `(ctx/add-subject-context rows db)` |
| large matrices | `(pqd/measurement-matrices db)`, `(http/measurement-matrix db matrix-key)` |

Results are vectors of maps with unqualified kebab-case keys, and enum
values come back as keywords. `(:patternq/provenance (meta rows))` gives the
database, basis t, and query time. Cite these when you report numbers.

## Raw queries are the fallback

Use the catalog first. When nothing in it fits, write Datomic Datalog as
Clojure data and run it with `patternq.db/q`, the same shape as
`datomic.api/q`:

```clojure
(pdb/q '[:find ?id ?sex
         :where [?s :subject/id ?id]
                [?s :subject/sex ?x]
                [?x :db/ident ?sex]]
       db)
```

Every catalog function has a `*-query` companion that returns its query, a
good starting point to adapt. `(pqd/schema-info db)` gives the schema
version. Bind an attribute or entity in every clause: the server refuses full
scans (`query-too-broad`). Over HTTP, expression functions are limited to a
whitelist, and string literals that start with `:` or `?` must be passed as
arguments.

## Rate limits

Each API key gets about 600 requests a minute and 4 queries in flight.
patternq handles HTTP 429 and 503 itself: it honors Retry-After, retries up
to 5 times, and caps concurrent queries at 4 across threads and futures.
Don't add your own retry loop. An `ex-info` with `:patternq/throttled true`
means the limit outlasted the retries, so slow down. A server-side query
timeout isn't retried; narrow or page the query instead (the server caps
timeouts at 120 s). Prefer a few broad catalog calls over many tiny ones.

## More for agents

- Agent guide: https://data-commons.rcrf-dev.org/docs/agents (also
  https://data-commons.rcrf-dev.org/llms.txt).
- Without a client library, the commons also serves MCP at
  https://data-commons.rcrf-dev.org/mcp (Streamable HTTP, `Authorization:
  Bearer <API key>`), with `list_datasets`, `describe_dataset`, and `query`
  (EDN Datalog, at most 1,000 rows, no matrices). For larger results,
  matrices, and the analysis functions, use this library.
