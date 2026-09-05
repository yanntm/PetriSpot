# `test/mcc/` — reading a cluster campaign

The campaign itself is described in `BENCH.md` at the repository root: one
`oarsub` per model instance and examination, each writing an
`OAR.<jobid>.stdout` log into a directory named after the examination
(`L/`, `OS/`, `QL/`, `RD/`, `SM/`, `UB/`, ...). This folder holds what turns
those logs into tables we can analyse.

* `mcclogs2csv.py` — the collector. Takes examination directories, writes
  `runs.csv` (one row per log) and `verdicts.csv` (one row per formula).
* `csv/` — the tables of the campaigns we have collected, one folder per
  campaign, each with a `README.md` reporting what that campaign showed.

## What the collector extracts

A log carries both sides of the comparison, so no oracle file is needed: the
`Control values :` block is the oracle, the `FORMULA` / `STATE_SPACE` lines are
the tool's answers. The rules of `run_test.pl` are reproduced — the last
verdict a formula receives wins, an oracle `?` constrains nothing, and an
oracle `inf` accepts `+inf`.

`runs.csv` columns: the log path, model, examination, the host that ran it,
the ITS-Tools build version, the timeout, the extra arguments of the call, the
start stamp, the teamcity test counts, the suite duration and the time to the
first verdict; then what PetriSpot did — `petrispot` invocations, of which
`petrispot fail` could not even be spawned and `petrispot timeout` exceeded
their `PFLOWS`/`PSEMIFLOWS` budget, and the walker census `walk calls`,
`walk asked`, `walk solved`, `walk ms`, `walk killed`; then `failure`, the
first fatal signature the log carries (`no_input`, `overlarge_marking`,
`eclipse_fatal`, `out_of_memory`, `its_abort`); then the verdict census: `oracle` values expected, of which `known` are not
`?`, `answered` by the tool, and among those `ok`, `wrong`, `missed` (the
oracle knows, we do not), `bonus` (the oracle does not know, we do) and `extra`
(a verdict for a formula the oracle never mentions).

`verdicts.csv` has one row per formula with the same `status` vocabulary, plus
the `TECHNIQUES` the tool claimed for that answer — the column that says
whether a random walk, an SMT solver or the decision diagrams did the work.

## Use

```
python3 mcclogs2csv.py L OS QL RD SM UB -o csv/<campaign>/
```

Logs are not kept in this repository: they are rsynced from the cluster into
`/data/ythierry/MCC26run/`, the CSVs are what we commit.
