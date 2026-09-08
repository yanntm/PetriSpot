# libHSC as an MCC competitor — a first exploration (2026-09-07)

A sideline of the campaign session, to be picked up in a dedicated one. What
`hsc-mcc` can answer today, how it does on contest models, what a driver looks
like, and what the next session should build. Read-only on libHSC: nothing was
changed there; the driver prototype lives in `MCC-drivers/hsc/`.

## What exists (libHSC 6fc4342, `build/tools/`)

* `hsc-mcc model.pnml -mcc StateSpace|OneSafe [--decompose]`: PNML P/T import
  (PetriSpot's loader, vendored), the NUPN unit tree as the shape or a Louvain
  decomposition with `--decompose` (a flat shape when neither), saturation,
  MCC-formatted lines. StateSpace prints the STATES value only. Deadlock was in
  the sample oracle but dropped from the pipeline (noted as debt in
  `tools/check_samples.cmake`). No timeout or memory flag.
* `nupn2hsc model.pnml`: the same import as `.hsc` text, ending with
  `(reach R saturate) (count R) (nodes R) (bill)`; strip those and append
  directives and queries to drive the engine by hand (`hsc file.hsc -e FORM`).
* Surface queries (`doc/hsc_manual.md` §8): `reach`, `select` (atoms over
  leaves, boolean forms, comparisons between leaves), `count` (a double),
  `max-value`, `states` (MCC format), `get-witness`, `xreach` (explicit),
  `cegar`. Directives: `reorder-force`, `decompose-louvain`,
  `simplify-constants`, `hotbit`, `flatten`.
* Binaries are dynamically linked (5 shared libraries); no CI publishes them.

## Measurements (one core, 45-60 s cap, oracle = MCC StateSpace)

Default `hsc-mcc` on 37 P/T models: 31 exact, mostly under 100 ms
(Raft-PT-04 3·10^12 states in 36 ms, ShieldRVt-010A 2·10^12 in 33 ms, DES-00a
2.4·10^10 in 8.8 s, Vasy2003 10^21 in 3 s, Kanban-00050 10^16 in 14 s); 1
off by one (Vasy2003: 9794739147610899087360 for …361, the count is a
`double`); 5 timeouts: Philosophers-000020, ParamProductionCell-0, Sudoku-AN04,
FMS-00500, Angiogenesis-10.

The directives change everything, in both directions (`/data/ythierry/hsc-exp/`):

| model | NUPN shape (default) | FORCE | Louvain | Louvain+FORCE |
|---|---|---|---|---|
| Philosophers-000020 | timeout | 35 ms | 29 ms | 31 ms |
| ParamProductionCell-0 | timeout | 96 ms | 189 ms | 50 ms |
| Kanban-00050 | 14.7 s | 10.7 s | 0.20 s | 0.20 s |
| DES-00a | 11.3 s | 15.2 s | 4.2 s | 0.10 s |
| Peterson-3 | 0.47 s | 0.21 s | 7.3 s | 1.4 s |
| SwimmingPool-01 | 0.15 s | 0.12 s | 5.5 s | 5.4 s |
| SafeBus-03 | 60 ms | 45 ms | 192 ms | 187 ms |
| Sudoku-AN04, Angiogenesis-10 | timeout | timeout | timeout | timeout |

(`--decompose` of `hsc-mcc` behaves like Louvain without FORCE.) FORCE alone
never lost more than a factor 1.4 and often won; Louvain wins big on the
process-like nets and loses big on SwimmingPool and Peterson. No configuration
dominates: the driver runs them as a portfolio. Vasy2003: 3 s by default,
0.1 s under Louvain (with or without FORCE), 6.5 s under FORCE alone, the
count off by one in every case (the `double`). FMS-00500 (500 tokens per
place, an `int_set` leaf per value) times out under every configuration: the
high-token nets are the bounded-integer leaf theory of roadmap R3.

## What it could compete in

* **StateSpace**: STATES now (exact below 2^53). The other three values are
  expressible with existing queries, driver side: `max-value` is
  MAX_TOKEN_IN_PLACE; TRANSITIONS is the sum over transitions t of
  `count(select(R, guard_t))` (edges = enabled pairs); MAX_TOKEN_PER_MARKING is
  a maximum of a sum, which needs the `lia` layer to accept a sum in a query.
* **OneSafe**: now; an overflow past the leaf bound answers FALSE, correctly.
* **ReachabilityDeadlock**: `select(R, (not (or guard_1 … guard_n)))`
  non-empty; the dropped pipeline is that formula.
* **RC, RF**: EF phi is `select(R, phi)` non-empty, AG the negation. PetriSpot's
  `Petri/src/parse/mcc` turns the contest XML into the expression tree and
  `parse/sexpr` prints s-expressions: a printer target for the surface's atom
  syntax, provided sums of leaves in comparisons are accepted (guards take
  `(< a c)`; a sum is the same `lia` path, to confirm).
* **QuasiLiveness**: every `select(R, guard_t)` non-empty. **StableMarking**:
  per-leaf domains. **UpperBounds** on a place: `max-value`.
* Not reachable soon: Liveness, CTL, LTL, COL.

## Correctness and safety for a sweep

* The `double` count is wrong above 2^53 (Vasy): a big-integer `count` in
  `src/surface_query.cc` before any scored StateSpace run.
* Memory: the tables are unbounded (roadmap R2 open); under 16 GB the failure
  is a crash, never a wrong answer; `ulimit -v` in the driver makes it clean.
* Non-safe nets: the leaf domain is `[0, max initial + 1)`; a place growing past
  it raises `overflow_error`: OneSafe reads it as FALSE, StateSpace aborts with
  no output. Kanban-50, Diffusion2D, SwimmingPool, GPPP answered exactly, so the
  bound as set covers them.
* Unsupported examinations are refused (`DO_NOT_COMPETE`).

## The driver (`MCC-drivers/hsc/`, prototype committed)

`BenchKit_head.sh`: `nupn2hsc` once, then four `hsc` runs in parallel on the
same `.hsc` (NUPN shape, FORCE, Louvain, Louvain+FORCE), the first with an
answer wins and the others are killed; `ulimit -v` from the memory
confinement, `timeout` from the time confinement, `CANNOT_COMPUTE` when none
answers. OneSafe is decided from `max-value`. `SupportedExamination.txt`:
StateSpace PT, OneSafe PT. `install.sh` copies binaries from a local build
(`HSC_BUILD`). Checked through `run_test.pl` locally
(`/data/ythierry/MCC26deploy/smoke-hsc-*.log`): Philosophers-000020 STATES
exact in 0.06 s (winner Louvain+FORCE), Kanban-00050 in 0.23 s, SwimmingPool-01
in 0.17 s (winner FORCE), Raft-02 and Angiogenesis-05 OneSafe right, Sudoku-AN04
`CANNOT_COMPUTE` at the 55 s confinement. `run_test.pl` marks StateSpace
runs failed because they carry one of the four expected values; the value
itself is right. The deploy clone `/data/ythierry/MCC26deploy/MCC-drivers`
holds a copy of `hsc/` but cannot `git pull` (untracked local oracle files
collide with ones upstream now tracks); move them aside, then pull.

## Next session

1. libHSC: exact `count` (big integer or `__int128` plus a string), and static
   linking of the tools (a `-static` option like PetriSpot's) so a cluster copy
   is one file each; a CI job publishing them (Inv-Linux style) later.
2. Driver: TRANSITIONS and MAX_TOKEN_IN_PLACE from existing queries; then
   ReachabilityDeadlock and the RC/RF atoms through PetriSpot's parser and a
   surface printer.
3. Sweep: `BK_TOOL=hsc ./run_oar.sh "oracle/*-SS.out"` on the cluster after
   the Airplane warmup, collected with `mcclogs2csv.py`; the contest's
   StateSpace table gives the comparison with ITS-Tools and TEDD.

## 2026-09-08, evening: superseded

Items 1 and 2 of "Next session" are done: exact `count` (GMP), static
binaries from the libHSC CI (`HSC-Linux`), `hsc-pn` answering RC, RF, RD, UB,
StateSpace (four values) and OneSafe on PNML + MCC XML or PNET +
s-expressions, the driver `MCC-drivers/hsc/` rewritten on it, and ITS-Tools
plugins with `-hsc` / `-hscBench`. Design and decisions: `HSC_PLAN.md`;
the campaign spec (item 3 and beyond): `HSC_EXPERIMENTS.md`; current state
and next actions: libHSC `handoff_mcc.md`. New measurements go in dated
sections below this one.

## 2026-09-08, campaign 20260908-hsc: StateSpace, first partial read

`BK_TOOL=hsc`, StateSpace, 300 s, 4 cores, `tall%`, 1391 instances submitted
(binaries from the libHSC CI branch `HSC-Linux`; driver: four configurations
in parallel, the memory confinement split among them). Read at 163
collected logs while the campaign drains; collected into
`/data/ythierry/MCC26run/20260908-hsc/` (`csv/REPORT.md`), on the pages as
the set `libHSC 20260908 SS`.

| | |
|---|---:|
| instances read | 163 of 1391 |
| oracle values known | 652 |
| answered | 312 |
| ok | 309 |
| wrong | 3 |
| runs at the 300 s wall | 41 |
| median run | 60 s |
| longest run | 300.6 s (the walltime of 10 min was never binding) |

**The three wrong values are all TRANSITIONS, on BART-COL-002, -005, -010,
and they are not the engine's.** A coloured instance reaches us unfolded:
the harness runs ITS-Tools on it because we declare P/T only, and that
unfolder fuses symmetric bindings without reporting their multiplicity. The
unfolded net has the right reachable states (STATES, MAX_TOKEN_IN_PLACE and
MAX_TOKEN_PER_MARKING all matched the oracle on those three) but fewer arcs
than the coloured semantics, so an arc count over it is an undercount — a
constant 167/202 of the oracle across the three instances. The driver now
leaves TRANSITIONS unanswered on a coloured input rather than reporting a
number we know to be biased. Answering it properly needs the unfolder to
report a per-transition binding multiplicity `m(t)`, after which
TRANSITIONS is `Σ_t m(t) · |{s ∈ R : s enables t}|`: a feature request on
the ITS-Tools side, noted in `HSC_PLAN.md` section 4.

Also found: `hsc-pn` dumped core on a coloured PNML (the vendored reader
throws a string literal, which no handler caught). Fixed in libHSC: such a
net is now a clean error, `Net is not a P/T net-> Colors are not supported
currently.`, exit 1.

Two collector notes for whoever reads the tables: 58 of the 163 logs are
counted as "truncated (no trailer)" although no run exceeded its walltime —
`collect.sh` excludes `*.stderr`, where the `time -p` trailer lives; and 3
runs carry an `eclipse_fatal` signature, which is the unfolder failing on a
coloured instance, not our tool.

The full read (all 1391 instances, the failure classes, the winning
configuration per family, the comparison against ITS-Tools, TEDD and the
2025 gold) belongs to the next session: `HSC_EXPERIMENTS.md` E1, with the
ITS-Tools 300 s baseline of the same examination beside it.
