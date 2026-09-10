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

The directives change everything, in both directions (`/data/ythierry/MCC26logs/_local/hsc-exp/`):

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

### The full read (1391 instances, campaign 20260908-hsc)

| set | answered | ok | wrong | at the wall | median s |
|---|---:|---:|---:|---:|---:|
| **libHSC, 300 s** | 3612 | 3589 | 23 | 668 | 225.6 |
| ITS-Tools 2026 (contest run) | 3485 | 3485 | 0 | 202 | 69.2 |
| TEDD 2026 | 5524 | 5524 | 0 | 0 | 4.9 |
| 2025 gold | 5512 | 5512 | 0 | 0 | 5.6 |
| TY 2026 | 1252 | 1252 | 0 | 314 | 332.1 |
| petrivet 2026 | 460 | 431 | 29 | 659 | 2143.4 |

5572 oracle values are known; we answered 3612 of them on our first attempt,
slightly more than the contest's ITS-Tools run (whose StateSpace answers omit
TRANSITIONS on many instances) and well short of TEDD, which answers 99 %.
Values right by kind: STATES 1066, MAX_TOKEN_IN_PLACE 1066,
MAX_TOKEN_PER_MARKING 779, TRANSITIONS 681. 1067 of 1391 runs produced at
least one value; 668 hit the 300 s wall. The gap to TEDD is time and memory,
not correctness — but note the contest sets ran at the contest's much longer
confinement, so the honest comparison is the ITS-Tools 300 s baseline
submitted beside this one (`SS.itstools`).

**All 23 wrong values are understood and fixed.**

* 22 TRANSITIONS on coloured instances (BART, DrinkVendingMachine,
  GlobalResAllocation, PhilosophersDyn, UtilityControlRoom): the unfolding
  path fuses bindings, as analysed in `HSC_PLAN.md` sections 10 to 13. The
  driver now leaves TRANSITIONS unanswered on a coloured input.
* 1 STATES and 1 MAX_TOKEN_IN_PLACE on GPPP-PT-C0010N1000000000, whose
  markings pass four billion: the PNML parse assigned a `long` into the
  build's integer type and silently shortened, so we answered STATES 2
  against an oracle of 176894515156. Fixed at the parse (`castExact` in
  `Arithmetic.hpp`, `PTNetHandler` throws): a value beyond the build's width
  is refused, and `petri32` refuses it too where it used to truncate.
  `petri64` is unaffected.

**Two lessons for the next campaign.** The winning configuration was printed
on stderr, which `collect.sh` drops, so this run cannot say which shape and
order answered; the driver now prints the attribution on stdout. And the
per-configuration memory split is untested under real pressure: with four
configurations each holding a quarter of the confinement, a memory-bound
instance may lose all four.

### TRANSITIONS is exact before the reductions (2026-09-08, evening)

Decisive measurement on BART-COL-002. Handed the *unreduced* unfolding that
ITS-Tools produces (10865 places, 646 transitions, written as a PNET by
`HscRunner` with `-Dhsc.debug=1`), `hsc-pn --states` answers

```
STATE_SPACE STATES 17424          (oracle 17424)
STATE_SPACE TRANSITIONS 53328     (oracle 53328)
```

so the arc count is exact when we receive the net before the arc-destroying
rules, and the `TMULT` block is all ones there, honestly. The campaign's
44088 came from the harness path, where the reducer had already removed 302
of those transitions. Consequences:

* The harness path stays a defect for this value: the driver leaves
  TRANSITIONS unanswered on a coloured input, and `MCC-drivers` needs no
  further change.
* `hsc-pn` now reports TRANSITIONS only with evidence that the net's arcs are
  those of the net it came from: always on the PNML path, and on the PNET
  path only when a `TMULT` block says the producer accounted for what it
  dropped. Without evidence it prints the other three values and says why.
* The way to answer it on coloured models is therefore ITS-Tools calling us
  as a StateSpace back end on the unreduced unfolding.

Note on a neighbouring metric: ITS-Tools prints `STATE_SPACE
UNIQUE_TRANSITIONS`, which is not this value — it counts two transitions
with the same effect and different read arcs once.

### The StateSpace back end answers all four on a coloured model (2026-09-08)

`its-tools -pnfolder BART-COL-002 -examination StateSpace -hscBench -timeout 240`
now calls `HscRunner.runStateSpace` on the net as the StateSpace branch builds
it, before the rules that remove transitions:

| value | ours | oracle |
|---|---:|---:|
| STATES | 17424 | 17424 |
| MAX_TOKEN_IN_PLACE | 1 | 1 |
| MAX_TOKEN_PER_MARKING | 274 | 274 |
| TRANSITIONS | 53328 | 53328 |

All four right, in 120 s of a 240 s budget. This is what removes the arc-count
penalty from a diagram run: the examination, or that one value, can be
delegated. The cost is the unreduced net — the diagram path answers the other
three in 13 s on the reduced one — so the production shape is a reduced run
for the three cheap values beside an unreduced run for the arcs, which a
portfolio can do in parallel.

The metric deserves its reputation: everywhere else the contest reads a Petri
net as an unlabelled Kripke structure (no transition labels in the temporal
logics, no "t fired"), and only here does it ask for labelled edges, which is
why the value is fragile under every transformation that touches the
transition set.

### Equal budget against ITS-Tools (300 s, same corpus and hardware)

The baseline campaign `20260908-its` (`SS.itstools`, product 202609081701,
the native image) beside ours, right answers by value:

| value | libHSC | ITS-Tools |
|---|---:|---:|
| STATES | 1066 | 982 |
| MAX_TOKEN_IN_PLACE | 1066 | 982 |
| MAX_TOKEN_PER_MARKING | 779 | 982 |
| TRANSITIONS | 681 | 0 |
| **total right** | **3592** | **2946** |
| wrong | 24 (all fixed since) | 0 |
| runs at the wall | 668 | 367 |
| wall hours | 60.5 | 45.9 |

Two readings. We answer **84 more state counts** than the SDD engine at the
same budget, and the arc count is ours alone (ITS-Tools prints
`UNIQUE_TRANSITIONS`, a different metric). And we lose
MAX_TOKEN_PER_MARKING 779 to 982, which is our own doing: we binary-search
it with repeated `select`s over a sum spanning every place, a crossing atom,
where the diagram engine reads it off in one pass. The fix is a query that
folds the maximum of a sum over the diagram bottom-up, the same shape as
`cardinal_as`, turning log(bound x P) symbolic searches into one traversal.
That is the cheapest win visible in this table.

### The counting record on reduced nets (2026-09-08)

Through `its-tools -examination StateSpace -hscBenchReduce`, all four values
on nets the reductions have worked over:

| model | reduction | ours | oracle |
|---|---|---|---|
| AutonomousCar-PT-01a | 9 duplicate transitions fused, 2 no-effect kept | 227, 1, 6, 654 | same |
| Dekker-PT-010 | 2 no-effect kept | 6144, 1, 20, 171530 | same |
| AirplaneLD-PT-0010 | 32 constant places removed, holding 32 tokens | 43463, 1, 38, 183664 | same |

AirplaneLD is the one that shows `PDROP` working: the exchanged net carries
`TMULT` and `PDROP: 32 removed constant places holding 32 tokens`, and
without that block the per-marking maximum would read 6 instead of 38.
Dekker shows the no-effect rule: its arcs are self-loops nothing else
carries, so the rule now keeps those transitions while a record exists.
AutonomousCar shows the fusion: the survivor of each duplicate inherits its
weight.

The redundant composition rule is skipped while a record is tracked rather
than allowed to destroy the arc count — it prefers longer firing paths and
buys little. A plain run with no HSC flags applies every rule as before and
creates no record.
