# `test/mcc/` — reading a cluster campaign

The campaign itself is described in `BENCH.md` at the repository root: one
`oarsub` per model instance and examination, each writing an
`OAR.<jobid>.stdout` log into a directory named after the examination
(`L/`, `OS/`, `QL/`, `RD/`, `SM/`, `UB/`, ...). This folder holds what turns
those logs into tables we can analyse.

* `mcclogs2csv.py` — the collector. Takes examination directories, writes
  `runs.csv` (one row per log) and `verdicts.csv` (one row per formula).
* `totallogs2csv.py` — the collector of the total examinations
  (`TOTAL_QUERIES.md`): one row per log with the atoms answered, the
  witnessed/proved split, the engine that closed them and the time at which a
  quarter, half, three quarters and all were closed; `--oracles DIR` also
  writes the answers as vector oracles in the `pnmcc-models-2026` format.
* `log2oracle.py` — one harness log of a total examination into its oracle
  file, the run's answers as the vector and `?` where it said nothing; the
  same writer as the collector's `--oracles`, for a single run.
* `toolsupport.py` — joins a `verdicts.csv` with the contest's
  `raw-result-analysis.csv` (in `pnmcc-models-20xx/website/`) and writes
  `support.csv`: for every value, which tools produced the oracle verdict
  (`who`, `backing`) and which produced a different one (`dissent`). It is what
  separates a value nobody has from a value we alone are missing.
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
`eclipse_fatal`, `out_of_memory`, `its_abort`); then the verdict census:
`oracle` values expected, of which `known` are not `?`, `answered` by the tool,
and among those `ok`, `wrong`, `missed` (the
oracle knows, we do not), `none` (neither the oracle nor we have it), `bonus`
(the oracle does not know, we do) and `extra` (a verdict for a formula the
oracle never mentions).

`missed` and `none` are both silences, and only the first is a gap of ours:
`missed + none` is what an ideal tool would have answered, `missed` is what the
field already answers and we do not.

`verdicts.csv` has one row per formula, with the oracle value, the answer and
the same `status` vocabulary.

## Use

```
python3 mcclogs2csv.py L OS QL RD SM UB -o csv/<campaign>/
python3 totallogs2csv.py QLA SMA UBA -o csv/<campaign>/ --oracles csv/<campaign>/oracles/
python3 toolsupport.py ~/git/pnmcc-models-2026/website/raw-result-analysis.csv \
        csv/<campaign>/verdicts.csv -o csv/<campaign>/support.csv
```

Logs are not kept in this repository: they are rsynced from the cluster into
`/data/ythierry/MCC26run/`, the CSVs are what we commit.
