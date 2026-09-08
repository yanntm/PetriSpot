# libHSC on the contest: first experiments

The spec the campaign reports against; results go to `libHSC_in_MCC.md`
(dated sections) and the collected tables to `Petri/test/mcc/csv/`. Plan:
`HSC_PLAN.md`; tool: libHSC `tools/README.md`; harness: `docs/CLUSTER.md`.

## Readiness

| piece | state |
|---|---|
| `MCC-drivers/hsc/` on `hsc-pn`, six examinations, `install.sh` from `HSC-Linux` | committed, not pushed; deploy clone must pull |
| `run_test.pl` pass in the deploy tree (`libHSC_in_MCC.md`, "The driver") | to do before a submit |
| memory: no bound in `hsc`; `ulimit -v` in the driver, a kill is no answer | acceptable |
| ITS-Tools `-hscBenchReduce` | committed, not pushed; needs the ITS-Tools CI product |

## E1. StateSpace baseline, whole P/T corpus

`BK_TOOL=hsc`, examination StateSpace, confinement 300 s, 4 cores (the
portfolio: NUPN, FORCE, Louvain, Louvain+FORCE). Baseline per instance from
the contest `raw-result-analysis.csv` (`~/git/pnmcc-models-2026/website/`):
ITS-Tools (SDD over GAL), GreatSPN (Meddly), plus the user's older Marcie
(IDD) data. Collect: answered, exact (all four values against the oracle),
time, winning configuration, and the failure class from the driver's
stderr tails: timeout, `overflow_error` (leaf bound), memory kill.
Expected: wins on safe process-like nets, losses on high markings,
correlated with the maximum marking (E4).

## E2. Dev set: PNET + properties pairs

About thirty instances by regime: safe NUPN (Raft, Philosophers,
ShieldRVt), high marking (FMS, Kanban, SwimmingPool, Angiogenesis,
Diffusion2D, GPPP), structurally large (DES, Vasy2003, Peterson, Sudoku).
Each as the raw PNET (`petri64 --exportNet`) and the ITS-Tools reduced PNET
(kept by `its-tools ... -hscBenchReduce -Dhsc.debug=1`), with RC and RF
forms (`--printProps=sexpr-index`, or the runner's files) and the oracle
lines. Committed under libHSC `examples/mcc/` with a manifest, run by
`tests/pn_samples.sh`. The bench for `select` cost, shape choice, and the
semiflow shapes of `HSC_PLAN.md` section 8; a minute on a laptop.

## E3. The value of the pipeline (RC, RF)

Three columns on the same instances: the `hsc` driver on the raw net;
ITS-Tools `-hscBenchReduce` (reductions, then hsc-pn); full ITS-Tools
`-its`. First gap = structural reductions; second gap = what walks and SMT
settle before any engine. Column one is cluster ready; two and three need
the ITS-Tools product with the new plugins.

## E4. Leaf theory diagnosis

For every E1 timeout: maximum marking, places, transitions, which
configuration won on the family's smaller instances. Losses clustering at
high markings point to roadmap R3 (bounded-integer leaf; IDD-style
intervals the candidate); losses elsewhere point to shape.
