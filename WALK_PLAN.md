# PetriSpot explicit walk engine: survey, architecture, heuristics, plan

Status: design document, revision 3. Iteration 1 is implemented (commits up
to the relaxed-plan strategy, 2026-09-03): folders and CMake, property AST and
MCC parser, sparse walk engine, random / best-first / structural / relaxed
strategies, `--props --query` driver. Section 8 records the state.

Goal: extend PetriSpot with an explicit, (almost) memoryless, multi-threaded
walk engine whose job is to find counter-examples to reachability properties
quickly on large P/T nets. It complements the Java ITS-Tools tool chain (which
does the reductions, the SMT and symbolic work and today's walks); it does not
replace it. Being exhaustive is not a goal. Being very good at finding a
witness when one exists is.

Design axioms fixed by the user:

- Sparse everywhere. Transitions touch a handful of places out of tens of
  thousands; `SparseArray` / `MatrixCol` are the substrate, no dense
  per-place structures in the hot path.
- Threaded by design. Hot data is thread-local; what threads share goes
  through an explicit Knowledge Base, whose content is deliberately flexible.
- Properties are parsed into one AST; several concrete syntaxes are a side
  goal we will reach (MCC XML, PetriVizu text, s-expressions).
- Models and benchmarks live outside the repo (cluster, MCC corpora,
  ITS-Tools dumps); the repo holds a few small examples.
- Parikh hints and other ITS-Tools inputs come later.
- The source tree gets reorganised into folders.

---

## 1. What exists today

### 1.1 master (the invariant tool)

| Component | File | State / notes |
|---|---|---|
| Sparse vector | `src/SparseArray.h` | Template on integral `T`. Sorted keys, merge-based `sumProd`, `greaterOrEqual` with galloping binary search, `manhattanDistance`, `countContainsPos`, `hash`/`==`. Solid, the substrate. |
| Column matrix | `src/MatrixCol.h` | Vector of `SparseArray` columns, `transpose`, `sumProd`. flowPT / flowTP / incidence. |
| Net | `src/SparsePetriNet.h` | `marks` (dense `vector<T>`, initial only), `flowPT`, `flowTP`, names. |
| PNML parser | `src/PTNetHandler.h`, `src/PTNetLoader.h` | expat SAX, P/T only, deferred arc patching. Recently optimised. Reused as is. |
| Walker | `src/Walker.h` | 2020 deadlock prototype, ported from the Java `RandomExplorer` of that era. Design reference, not code (1.3). |
| Invariants | `Invariant*.h`, `RowSigns*.h`, `MixedSignsUniqueTable.h`, `Heuristic.h` | Mature. Out of scope here except as a later heuristic input (semi-flows give bounds and cycles). |
| IO | `SparseMatrixIO.h` (KERS), `MatrixExporter.h`, `PNMLExport.h`, `FlowPrinter.h` | Exporters. |
| Timeout pattern | `InvariantMiddle::computePInvariants` | `std::thread` + promise/future `wait_for`. |
| Build | autotools, C++23, static, only dependency expat. `petri32/64/128` via `-DVAL`. `subdir-objects` already enabled, so folders cost nothing in `Makefile.am`. |

Working tree has an uncommitted `--check` feature (`Petri.cpp`,
`InvariantHelpers.h`), unrelated; commit it first.

### 1.2 origin/er/link-spot (Jan 2021, 5 years stale)

Merge base `dbe0ee9`, before the template refactoring. Adds:

| Component | Reusable? | Comment |
|---|---|---|
| `expr/Op.h` | yes | Operator enum: MCC atoms (CARD, ENABLED, BOUND, DEAD), comparisons, boolean ops, CTL and LTL operators. |
| `expr/Expression.h`, `BinOp.h`, `NaryOp.h`, `VarRef.h`, `Constant.h`, `BoolConstant.h` | yes, after porting | Virtual `eval(const SparseIntArray&)`, `print`, child navigation. Needs `T` template, evaluation against the walker's marking, and `is-fireable`. The Java simplifier is present only as a comment block: it is the spec for ours. |
| `expr/Property.h` | yes | name + body. Needs a type (reach / invariant / deadlock / bounds / ctl / ltl / atom). |
| `expr/parse/ExprHandler.h`, `ExprLoader.*` | yes, after porting | expat SAX parser for the MCC property XML: tokens-count, integer-le, boolean ops, CTL and LTL path operators. `is-fireable` throws: we need it. `ext_hash_map` (libDDD) to be replaced by std containers. |
| `expr/AtomicPropManager.h` | later | Unique atomic propositions for LTL. |
| `Petricube.h` | later, as design | Spot `kripkecube` adapter over a marking and an enabled list. Good template for a future LTL back-end; pulls in Spot and libDDD. Not now. |
| `Petri.cpp` (branch) | no | Spot CNDFS driver. |

Verdict: cherry-pick the `expr/` tree and the MCC parser by hand (`git show`),
modernise, leave Spot on the shelf. The walk engine's successor iterator must
stay adaptable to a `kripkecube`-style interface later.

### 1.3 Lessons from the current `Walker.h`

Keep:

- `combFlow`: one effect column per transition (post minus pre).
- `behaviorMap`: transitions grouped by identical effect; random choice over
  distinct behaviours, duplicates never fired twice from a state.
- Incremental `updateEnabled`: only consumers of places that grew can become
  newly enabled.
- "Repeat last": iterate a transition while enabled if its effect strictly
  decreases some place (no divergence).
- Dead-end handling: reset to the initial marking, count resets.
- The "fewest successors" rule for deadlock, which is a distance heuristic in
  disguise (distance to deadlock = number of enabled transitions).
- `dropUnavailable(enabled, parikh)`: vestige of the Parikh-guided walk.

Fix:

- `updateEnabled` allocates and zeroes two arrays of size |T| per step, so the
  incremental update is O(|T|) anyway; `computeEnabled` on reset is
  O(|T| x |pre|) with a `greaterOrEqual` per transition.
- One `sumProd` allocation per step for the new marking; raw `new int[]`
  lists and `memcpy` per candidate in the lookahead.
- Deadlock hard-wired; no property, no trace, no parseable output, one
  strategy, one thread.

### 1.4 Neighbouring code worth knowing

- libHSC (`~/git/libHSC`) vendors our `SparseArray`, `MatrixCol`,
  `SparsePetriNet` and PNML parser unchanged (see its
  `include/hsc/petri/README.md`), and has a 84-line s-expression reader
  (`include/hsc/surface/sexpr.hh`: a `datum` is an atom or a list with a
  source line; `parse`, `write`). Its manual section 6 defines the
  SMT-flavoured expression syntax: `(and BEXP*) (or BEXP*) (not BEXP)
  (imply a b) (CMP EXPR EXPR) (OP EXPR EXPR)`.
- PetriVizu (`~/git/PetriVizu/public/syntax.md`) defines the infix syntax:
  `property "Name" [reach|atom|ctl|ltl|bounds] : EF P1 + P2 <= 5;`, `@Name`
  references to earlier properties, quoted identifiers when they clash with
  keywords.

---

## 2. Requirements recap

- Reachability first: given a property set, answer as many as possible with
  a witness (EF phi true, AG phi false), stay silent or say UNKNOWN on the
  rest. Deadlock is a target like the others.
- Mostly memoryless: current marking, enabled set, scratch, static per-net
  tables, a bounded Knowledge Base. Never a growing set of visited states.
- Per-step cost proportional to the arcs actually touched, independent of
  |T| and |P|.
- Guided: pluggable heuristics, restarts, portfolio of strategies.
- Threaded from the first version that matters (Phase 4), designed for it
  from Phase 2.
- Subprocess of ITS-Tools: deterministic CLI, parseable output, witness
  trace on request.
- 32/64/128-bit variants kept, static link, no new dependencies (expat only).

---

## 3. Proposed architecture

### 3.1 Source tree

```
Petri/src/
  core/        SparseArray.h SparseBoolArray.h MatrixCol.h SparsePetriNet.h
               Arithmetic.hpp Rational.h InvariantHelpers.h   (moved, unchanged)
  parse/       PTNetHandler.h PTNetLoader.h                  (moved)
               mcc/   PropertyHandler.h PropertyLoader.h     (ported from link-spot)
               sexpr/ Sexpr.h SexprProperties.h              (new, small)
               vizu/  PropertyLexer.h PropertyParser.h       (new, later phase)
  expr/        Op.h Expression.h BinOp.h NaryOp.h VarRef.h Constant.h
               BoolConstant.h Property.h PropertySet.h Simplify.h Printer.h
  invariants/  InvariantCalculator.h InvariantMiddle.h InvariantsTrivial.h
               RowSigns.h RowSignDomination.h MixedSignsUniqueTable.h
               Heuristic.h                                   (moved, unchanged)
  io/          SparseMatrixIO.h MatrixExporter.h PNMLExport.h FlowPrinter.h
  walk/        WalkNet.h EnabledSet.h Marking.h Target.h Distance.h
               Strategy.h strategies/*.h Walker.h Trace.h
               KnowledgeBase.h Checkpoint.h Portfolio.h
  Petri.cpp    CLI (thin; option parsing may move to cli/ if it grows)
  kersconv.cpp
```

Each folder gets a `README.md`; `walk/` and `expr/` also get an
`algorithm.md` written before the code, per the repo conventions. Headers are
included as `"core/SparseArray.h"` with `-I$(srcdir)`. `Makefile.am` lists
sources per folder; `subdir-objects` is already on.

### 3.2 expr: one AST, several syntaxes

The AST is the link-spot tree, modernised:

- `Expression<T>` (or `Expression` evaluating in `T` through a marking
  accessor; decision below) with `Op`, children, `eval(const Marking&)`,
  `print(Printer&)`.
- Leaves: place reference, transition reference (for `is-fireable`),
  integer constant, boolean constant, deadlock, property reference `@Name`
  (resolved at load time into the referenced body).
- `Property`: name, type (atom, reach, bounds, ctl, ltl), body.
- `PropertySet`: ordered list with name index.
- `Simplify`: constant folding, flattening of n-ary and/or, double negation,
  negation pushing, `EF !phi` versus `AG phi` normalisation. Transcribed from
  the Java block commented in link-spot's `Expression.h`.
- `Printer`: one printer per syntax (MCC XML, PetriVizu, s-expr). Printing
  then re-parsing in every syntax is the cross-validation of the parsers.

Parsers, in order of arrival:

1. **MCC XML** (`parse/mcc/`). The MCC corpora and the ITS-Tools dumps are
   PNML plus this XML, so it is mandatory and it is what the benchmark
   pipeline feeds us. Ported from link-spot, plus `is-fireable`, plus the
   property type inferred from the file name and the top operator. Its
   defects (no true/false leaves, verbose) are absorbed by the AST and the
   simplifier.
2. **s-expressions** (`parse/sexpr/`). A reader in the style of libHSC's
   `datum` (about 100 lines) and a datum-to-AST pass (about 150 lines).
   Syntax: `(property Name reach (EF (<= (+ P1 P2) 5)))`, `(atom Name ...)`,
   `(and ...) (or ...) (not ...) (imply ...)`, `(fireable t1 t2)`,
   `(deadlock)`, `(ref Name)` or `@Name`. Trivial to parse and to emit, so it
   is the natural format for tests, for hand-written properties in the
   examples folder, and for program-to-program exchange with ITS-Tools if we
   ever prefer text over the XML (KERS already covers the matrix side).
3. **PetriVizu infix** (`parse/vizu/`). Human-friendly, already specified,
   and it is the syntax the user reads and writes. Needs a lexer and a
   precedence-climbing parser with the keyword-quoting rule (about 400
   lines). Scheduled after the walker works, unless it turns out we want to
   hand-write many properties earlier.

Recommendation: MCC XML and s-expr in Phase 1, PetriVizu in Phase 6 or when
needed. All three print; the CLI gets `--props=<file>` with the syntax chosen
by extension (`.xml`, `.sexpr`/`.lisp`, `.prop`) and `--printProps=<syntax>`
to convert.

### 3.3 WalkNet: the compiled net, read-only, shared by threads

Built once from `SparsePetriNet<T>`:

- `pre[t]`, `post[t]`: the existing sparse columns.
- `effect[t]`: post minus pre (the `combFlow`); empty-effect transitions
  flagged (only matter for `is-fireable` and deadlock semantics).
- `effectClass[t]`: identical-effect class id (the `behaviorMap`) and a map
  from an effect to its exact negation, if any (anti-oscillation).
- `consumers[p]`: transitions with an arc from `p`, with weight, **sorted by
  weight** (transpose of `flowPT`, one sort per row).
- `producers[p]`: transitions with a positive effect on `p`.
- `maxConsumeWeight[p]`, `nbPre[t]`.
- Per-property structural tables (4.3) computed by `Target`, not here.

Memory: three to four sparse copies of the arc set. Tens of MB for a net with
half a million transitions. Fine.

### 3.4 Marking: sparse, thread-local, updated in place

The marking is a `SparseArray<T>` (sorted keys, values). Two update schemes,
both sparse, both allocation-free after warm-up; we keep both behind one
interface and measure:

- **In-place merge.** For each `(p, d)` in `effect[t]`: locate `p` in the
  marking (binary search, or a merge walk since both are sorted), add `d`,
  remove the key if it hits zero, insert it if it was absent. Insert and
  remove shift the tail of the key and value arrays (`memmove`); the tail is
  short on average since effects are small and keys are spread.
- **Double buffer.** `sumProd(1, m, 1, effect[t])` written into a
  preallocated second buffer of the same capacity, then swap. Zero
  allocation, cost O(|m| + |effect|), branch-friendly, no shifting.

The in-place merge is asymptotically better when |m| is large and effects
are tiny, which is the common case on big nets; the double buffer is
simpler and probably faster on small nets. The interface exposes what the
enabled-set update needs: for each touched place, its old and new value.

`unfire(t)` is `fire` with the negated effect, giving cheap one-step
lookahead without a copy. Initial marking and checkpoints are also
`SparseArray<T>`, so copying a state is one `operator=`.

### 3.5 EnabledSet: incremental maintenance with unsatisfied-arc counters

- `unsat[t]`: number of pre-arcs `(p, w)` of `t` with `m[p] < w`. Enabled iff
  `unsat[t] == 0`.
- `enabled`: packed list of enabled transition ids plus `posInEnabled[t]`
  for O(1) swap-remove; random choice O(1).
- After `fire(t)`, for each touched place `p` going from `a` to `b`:
  - decrease (`b < a`): consumers of `p` with weight `w` in `(b, a]` get
    `unsat++`; those reaching 1 leave `enabled`.
  - increase (`b > a`): consumers with `w` in `(a, b]` get `unsat--`; those
    reaching 0 join `enabled`.
  - `consumers[p]` is sorted by weight, so the affected range is a binary
    search and only the transitions that actually flip are visited. When
    `a` and `b` are both at least `maxConsumeWeight[p]` (or both below the
    smallest weight), nothing is visited.

Per-step cost: `O(sum over touched p of log(fanout(p)) + number of
transitions that flip)`. Independent of |T| and |P|.

`unsat`, `posInEnabled` and `enabled` are per-thread arrays of size |T|
(a few MB on the largest nets). They are static working storage, not state
memory; `reset()` copies them from a precomputed initial snapshot (O(|T|)
memcpy, amortised by the restart schedule). If |T| ever makes this hurt, a
`SparseBoolArray` snapshot of the initially enabled set plus lazy `unsat`
recomputation is the fallback; not expected to be needed.

`is-fireable(t)` atoms are then O(1) reads of `unsat[t]`.

### 3.6 Target: compiled goal predicate and distance

A reachability property normalises to "find a marking where `phi` holds"
(from `EF phi`, or `AG psi` with `phi = !psi`) or to deadlock (`enabled`
empty, with the MCC convention about empty-effect transitions decided
explicitly). `Target` compiles `phi` into:

- an evaluator over the marking and `unsat` (booleans, comparisons of linear
  sums of places and constants, fireability);
- a distance `dist(m) >= 0`, zero iff `phi` holds (4.2), evaluated
  incrementally: `atomsOfPlace[p]` (sparse: only places that appear in
  atoms) and `atomsOfTransition[t]` for fireability atoms, so that firing `t`
  re-evaluates only atoms touching `effect[t]`;
- goal places and goal transitions with a direction (need more / need less
  tokens) for the structural heuristics (4.3).

Several properties are handled together: every open target is checked at
each step (cheap, indexed by touched places). A strategy is told which
property it currently aims at; each has its own distance.

### 3.7 Strategy: transition choice policy

```
struct Strategy {
  virtual int choose(WalkContext& ctx) = 0;  // transition id, or -1 for "reset"
  virtual void onFired(int t) {}
  virtual void onReset() {}
};
```

`WalkContext`: marking, enabled set, current target, RNG, recent history,
WalkNet tables, read access to the Knowledge Base. Concrete strategies in
section 4; `Composite` mixes them (epsilon-greedy, phase scheduling).

### 3.8 Walker: the loop (one per thread)

```
reset(from initial marking or from a KB checkpoint)
loop while budget and not stopped:
  if an open target holds: publish witness to KB (trace replayed for verification)
  if enabled is empty: deadlock target? publish : reset
  t = strategy.choose(); if t < 0: reset
  fire(t); enabled.update(touched places); trace.push(t); strategy.onFired(t)
  every K steps: KB.pollStop(); maybe publish a checkpoint (3.9)
  restart schedule (Luby or geometric) -> reset
```

`Trace`: transition ids since the run's origin. A run started from a KB
checkpoint prepends the checkpoint's trace at output time (3.9), so a run's
own trace stays bounded by the restart schedule.

### 3.9 Knowledge Base: what threads share

Everything hot is thread-local. The Knowledge Base (KB) is the one shared
object, coarse-grained (a mutex per section, atomics for flags), touched only
every K steps or at reset. Its content is deliberately open; the initial
sections:

- **Verdicts**: per property, an atomic solved flag and the winning witness
  (marking + trace). Threads stop aiming at solved properties.
- **Checkpoints**: a bounded pool of promising states per property: marking,
  distance score, the trace from the initial marking, the strategy that
  produced it. Insertion when a thread's best distance in the current run
  beats the pool's worst; eviction keeps diversity (distinct markings by
  hash, distinct distance levels). Restarts draw from the pool with some
  probability instead of the initial marking: population-style search with
  a few hundred states of memory at most. Traces are stored as a parent
  checkpoint id plus a suffix, forming a checkpoint tree, so shared prefixes
  are stored once.
- **Structural knowledge** computed once and shared: goal tables (4.3),
  effect classes, and later invariant-derived facts (place bounds from
  semi-flows, transitions provably dead).
- **Statistics**: per-property best distance seen, restarts, steps;
  per-transition fire counts and "useful" counts (fired on a witness or on a
  checkpoint's path), which heuristics may read as priors. Cheap, bounded.
- **Hints** (later): Parikh vectors and other ITS-Tools inputs live here so
  every strategy can read them the same way.

Whether checkpoints actually help is an empirical question; the design keeps
them optional (`--kb-checkpoints=0` disables the pool).

### 3.10 Portfolio and timeout

N threads, each with its own `Marking`, `EnabledSet`, RNG seed, `Strategy`
and restart schedule, sharing the WalkNet, the compiled targets and the KB.
Global stop through an atomic checked every K steps; wall clock via the
existing promise/future idiom in the driver. Default assignment of strategies
to threads in 4.7.

### 3.11 CLI and output

```
--props=<file>            properties (MCC XML, .sexpr, later PetriVizu .prop)
--findDeadlock            kept; becomes a deadlock target
--walkSteps=<n>           per-run budget before restart (default e.g. 1M)
--walkThreads=<n>
--strategy=<random|bestfirst|structural|portfolio>
--trace=<file|->          print witnesses as transition sequences
--seed=<n>
--printProps=<mcc|sexpr|vizu>   convert and exit
-t <seconds>              wall clock, as today
```

Output: one MCC line per solved property, `FORMULA <name> TRUE|FALSE
TECHNIQUES EXPLICIT RANDOM_WALK ...`, nothing for unsolved ones unless
`--printUnknown`. Statistics after `-q` gating, as today. Witnesses are
replayed on a fresh marking before printing.

---

## 4. Heuristics to guide the search

A toolbox of cheap, composable signals. Each is O(touched arcs) per step or
precomputed once per property.

### 4.1 Baseline: random walk with restarts, made cheap

- Uniform choice over enabled effect classes, not transitions.
- Saturation: a transition whose effect strictly decreases some place and
  stays enabled can be fired `k` times in one update (`countContainsPos`).
- Restart schedule (Luby or geometric), fresh seed each time, optional
  restart from a KB checkpoint.
- Cycle signal under fixed memory: a Zobrist-style 64-bit marking hash
  maintained incrementally (XOR of `H(p, value)` for touched places) fed to
  a small fixed-size table or Bloom filter. A hit is a signal to perturb or
  restart, never a pruning decision. This is the only per-thread "visited"
  memory and it is constant-size.

### 4.2 Marking distance, best-first (TAPAAL style)

Jensen, Nielsen, Oestergaard, Srba, "TAPAAL and reachability analysis of P/T
nets", ToPNoC 2016; link-spot's `Expression.h` already points at it.

- `p >= k`: `max(0, k - m[p])`; `p <= k`: `max(0, m[p] - k)`; `p == k`:
  `|m[p] - k|`; `p != k`: `1 if m[p] == k else 0`; linear sums likewise.
- `AND`: sum of children; `OR`: min; negation pushed to atoms.
- `fireable(t)`: sum over pre-arcs of `max(0, w - m[p])`; `!fireable(t)`:
  min over pre-arcs of `max(0, m[p] - w + 1)`.
- Deadlock: `|enabled|`, possibly weighted by how hard each enabled
  transition is to disable.

Use:

- Candidate `t` scored by `dist(m + effect[t])`, computed incrementally on
  the atoms touching `effect[t]`. Minimum wins, ties random, epsilon-greedy
  or softmax over `-dist`.
- Candidate sampling: score a random sample of K (16 to 64) enabled
  transitions plus the enabled goal transitions (4.3), never the whole
  enabled list. Bounded cost and useful noise.
- Plateau escape: no improvement for N steps switches to random for a while
  or restarts; a small tabu list of recent transitions.

### 4.3 Structural distance: firings away in the net graph

Marking distance is blind when tokens must first travel along a chain, a
lock must be released, or a counter must reach a threshold through an
unrelated cycle. A static per-property signal:

- Goal places with a direction, goal transitions: producers of "need more"
  places, consumers of "need less" places, the transition itself for
  fireability.
- Backward BFS in the bipartite graph from goal transitions: `hops[t]` =
  length of the shortest chain `t -> post place -> consumer -> ... -> goal
  transition`. One pass, O(arcs), once per property, one small integer per
  transition. Unreachable transitions get infinity and are deprioritised,
  never forbidden.
- Strategy: prefer enabled transitions with smallest `hops`, tie-break with
  marking distance, epsilon randomness; or a weighted sum, tuned empirically.
- Stall-time re-rooting ("unlocking"): when no goal transition is enabled,
  take the disabled ones with smallest `hops`, read their unsatisfied
  pre-places (`unsat` plus a scan of `pre[t]`), and temporarily promote the
  producers of those places. A one-level dynamic re-rooting of the BFS,
  recomputed only on stall.

### 4.4 Knowledge-Base driven heuristics

- Checkpoint restarts (3.9): "go with the winners". A thread restarting from
  a checkpoint inherits a state already close to the goal under some
  distance; different threads restarting from the same checkpoint with
  different strategies and seeds is a cheap beam.
- Transition priors: `useful[t]` counts from witnesses and checkpoint paths
  of other properties on the same net bias the sampling weights.
- Cross-property reuse: a witness for property A is a state; it enters the
  checkpoint pools of the other properties if its distance to them is good.

### 4.5 Anti-oscillation and novelty

- Do not immediately undo: penalise the effect class that is the exact
  negation of the last fired one; tabu window of the last L effect classes.
- Place-level novelty with constant memory: per run, sparse min and max seen
  for goal-related places; a candidate that pushes such a place outside its
  seen range gets a bonus.
- The marking-hash filter (4.1) as a "we are looping" signal.

### 4.6 Deadlock specifics

Deadlock is a target with `dist = |enabled|` and structural preferences:
consume from places with a single consumer; prefer transitions whose firing
disables the most other transitions (computable from `consumers[p]` for the
places it decreases). The old forward walk and the "backward" (reversed net)
walk both become strategies; whether the reversed walk is worth keeping is
an empirical question.

### 4.7 Portfolio

Default with 4 threads:

1. random with restarts (fast, unbiased, feeds checkpoints);
2. best-first on marking distance, epsilon 0.1, K = 32;
3. structural hops with marking-distance tie-break;
4. best-first with softmax, restarting from KB checkpoints with high
   probability.

Each thread rotates its focus over the open properties at every restart.
First verified witness wins for a property.

### 4.8 Later: hints from ITS-Tools

Parikh vectors from the SMT over-approximation as a soft or hard restriction
of the choice (the vestigial `dropUnavailable`), T-semiflows to recognise
useless cycles, place bounds from P-semiflows to prune. All enter through
the KB hints section, none in the first deliverables.

---

## 5. Validation and benchmarking

- **Development**: the handful of nets in `Petri/examples/` plus a few small
  MCC instances with their property files and known verdicts, checked into
  `Petri/test/` (kilobytes, not megabytes). Every phase below names its
  check.
- **Benchmarking**: on the cluster, over the MCC corpora
  (`~/git/pnmcc-models-20xx`) and over ITS-Tools dumps (nets already
  reduced, agglomerated, with the properties ITS-Tools could not settle by
  other means). Those dumps are the real target: a property ITS-Tools
  already solves is not interesting. Verdict oracles are the MCC reference
  results and ITS-Tools' own answers.
- The repo holds the campaign scripts (`Petri/test/bench/`), the result
  tables and the analysis, never the models.
- Metrics: solved count per strategy and time-to-witness; steps per second
  per net size class; sensitivity to seeds. Disagreements with the oracle
  are bugs until proven otherwise.

---

## 6. Plan of attack

Each phase ends with something that runs and can be measured. Line counts
are rough effort calibration.

### Phase 0: housekeeping and tree reorganisation

- Commit the pending `--check`.
- `.gitignore`: binaries, objects, `res.txt`, the 245 MB PNML in
  `examples/` (or move it out of the tree).
- `git mv` into `core/ parse/ invariants/ io/`, fix includes, `Makefile.am`,
  build the three binaries, run the existing invariant examples as
  regression. One commit for the move.
- `README.md` per folder; `CLAUDE.md` at the root (done, revision 1).

### Phase 1: properties in (about 900 lines, half ported)

- `expr/`: AST, `Property`, `PropertySet`, `Simplify`, printers for the
  three syntaxes. `expr/algorithm.md` first (AST shape, normal forms,
  simplifier rules).
- `parse/mcc/`: port of link-spot's `ExprHandler`, `is-fireable`, property
  types.
- `parse/sexpr/`: reader and datum-to-AST.
- CLI: `--props` and `--printProps`. Check: parse the MCC files of a few
  models, print as s-expr, re-parse, compare ASTs; parse hand-written s-expr
  properties for the `examples/` nets.

### Phase 2: walk core, single thread (about 1000 lines)

- `walk/algorithm.md` first: WalkNet, marking update schemes, enabled
  maintenance, invariants of the data structures.
- `WalkNet`, `Marking` (both update schemes), `EnabledSet`, `Trace`.
- Rewire `--findDeadlock` onto the new core with the existing random
  strategy, as regression: same verdicts, steps per second at least as good.
  Delete the old `Walker.h`.
- Micro-benchmark: steps per second on `Airplane.pnml`, a mid-size MCC net,
  and the half-million-transition FamilyReunion net (kept out of the tree),
  for both marking update schemes. This validates the "cost independent of
  |T|" claim before any heuristic work.

### Phase 3: reachability end to end, single thread (about 600 lines)

- `Target` compilation and incremental evaluation; deadlock as a target.
- Random strategy with restarts and the hash-based cycle signal.
- MCC output, witness replay and printing, exit codes.
- Check: verdicts against known results on a batch of small MCC models where
  random walks are known to succeed.

### Phase 4: threads and Knowledge Base (about 600 lines)

- `KnowledgeBase` (verdicts, checkpoints, statistics), `Checkpoint` tree,
  `Portfolio`, global stop and timeout.
- Check: N threads of the random strategy give the same verdicts as one, no
  data races under `-fsanitize=thread` on the small tests, throughput scales.

### Phase 5: heuristics (about 800 lines, iterative)

- Marking distance with sampling and epsilon-greedy; structural hops and
  stall re-rooting; anti-oscillation and novelty; checkpoint restarts.
- Evaluation harness for the cluster (`Petri/test/bench/`): solved counts and
  time-to-witness per strategy versus the random baseline, on MCC subsets
  and ITS-Tools dumps. Most of the phase is tuning, not code.

### Phase 6: syntax completion and integration (about 500 lines)

- PetriVizu infix parser.
- Output protocol frozen and documented; Java side hook in ITS-Tools (other
  repo). Superseded by `INTEROP.md`, which designs the exchange formats and
  plans the replacement of `RandomExplorer` by PetriSpot.

### Later, explicitly out of scope now

- Parikh and invariant hints through the KB.
- LTL through Spot using the link-spot `kripkecube` design over the walker's
  successor iterator; bounds; CTL.
- Any exhaustive or stateful exploration.

---

## 7. Decisions taken (2026-09-03)

1. Marking update: sparse, in place; the double-buffer variant only if
   profiling asks for it.
2. Enabled maintenance by unsatisfied-arc counters over the place-to-transition
   matrix. Delta updates only; no transition-to-transition structure, ever.
3. Syntax: MCC XML first, AST independent of MCC. Bricks: true, false, and,
   or, not, and linear atoms `sum(coeff * place) cmp constant` (sums are
   required for MCC cardinality). `is-fireable(t)` is desugared by the parser
   into the conjunction of `place >= weight` over the pre-arcs of `t`, so it
   never reaches the AST. No deadlock atom in the AST; deadlock is a separate
   target.
4. Template on `T` kept end to end, effort on the 64-bit build. No template
   metaprogramming; the AST is plain data, evaluation is a small template
   function over the marking.
5. No Knowledge Base in the first deliverables; only a stop flag, a solved
   flag and a witness slot. One property per run for now (`--query=<n>`).
6. Threads after the single-thread engine works.
7. `--findDeadlock` moves to the MCC `FORMULA` output line.
8. Folder layout done (commit d5b034d); the build is CMake only (autotools
   removed on 2026-09-04, CI and `buildPetriSpot.sh` converted).

## 8. State after iteration 1 (2026-09-03)

Iteration 2 (2026-09-04) is the tool-to-tool interface of `INTEROP.md`: CLI
on CLI11 (`cli/`), s-expression properties (`parse/sexpr/`), the PNET binary
net (`io/PNETIO.h`), the result protocol, and multi-target walking
(`walk/TargetSet.h`: one walk checks every open property, a sweep round
precedes the focused rounds). The notes below describe iteration 1.

* `petri64 -i model.pnml --props=ReachabilityFireability.xml --query=3
  --strategy=relaxed --stall=300` answers the challenge query in about 100
  steps. On the whole challenge file, 9 of the 11 queries that have a witness
  are solved in under half a second each (queries 12 and 15 remain); the 5
  others have no witness. Every verdict produced on the three development
  models matches the MCC 2026 consensus.
* What made the difference is not raw speed but the heuristic: marking
  distance (TAPAAL style) and hop distance both stall on this net, whose goal
  is a coordination of 9 processes gated by a control token. The planning
  relaxed plan (h_add with helpful transitions, as in directed unfolding)
  sees through it. See `Petri/src/walk/algorithm.md`.
* Speed: 12k to 40k steps/ms on small nets with the random strategy. On the
  challenge net the enabled-set delta visits about 39k consumer arcs per
  step (places with 6561 consumers each flip between 0 and 1) for 2 real
  status changes, giving 80 steps/ms; the relaxed plan costs 0.5 to 3 ms per
  step there.
* Portfolio (2026-09-04): `--threads=N --strategies=name[:eps[:stall]],...`
  runs N walkers with thread-local state, round-robin strategies, first
  verified witness wins. On Angiogenesis fireability a 4-thread portfolio
  solves 15/16 in 1 s each (the union of what each strategy solved alone);
  on the hand-written "bridge full" query of BridgeAndVehicles-PT-V50P50N50
  (`Petri/test/props/`) best-first finds the witness in 75 ms where the
  random walk fails. ThreadSanitizer clean.
* `--share=N` is the first Knowledge Base piece: a bounded shared pool of
  best-of-run markings that restarts may draw from, with eviction of entries
  that fail to lead to improvements. Not decisive so far on the two open
  challenge queries (the relaxed plan runs too few restarts there to feed
  it), to be evaluated on the cluster.
* MCC harness (2026-09-04): `~/git/MCC-drivers/petrispot/` is the driver
  (PT reachability only, `BK_TOOL=petrispot`, or `petrispotxred` to run the
  ITS-Tools reducer first). Harness and oracles now follow
  pnmcc-models-2026 (formula ids verbatim). On the challenge model with a
  300 s confinement: alone 9/16 (all witnesses except 12 and 15); with the
  reducer 15/16, ITS-Tools settling the 5 queries without witness in its 30 s
  slice and the walk taking 10 residual ones on the reduced net. The residual
  properties come as `p <= 1 && -p <= -1`, which exposed a normal-form gap
  (negative leading coefficients) now fixed. `--totalTime` schedules open
  properties in rounds of growing budget. Logs in `Petri/test/logs/`.
* Open: the last challenge query (12); h_add cost on wide nets
  (bounded/backward variants).

Development models (outside the repo, `bench/models/`, git-ignored):
AirplaneLD-PT-0010, Angiogenesis-PT-05, and the challenge
ErlangenMainframeV1-PT-bP09C09 ReachabilityFireability (333 places, 59403
transitions, 16 properties, 209 fireability atoms; ITS-Tools solves 4/16,
TAPAAL all). The win condition is time-to-counterexample on such queries.

## 9. Hub places (2026-09-06): where the per-step cost stops being sparse

Measured on the QuasiLivenessAll campaign (`Petri/test/probes/qla_coverage.py`
over `Petri/test/mcc/csv/2026-09-06/total-runs.csv`, 1680 P/T instances, one
target per transition, 1800 s): completion is 1.0 at the median in every
transition-count decile up to 27 934 transitions; in the top decile it is
0.376 and 86 of 168 runs hit the wall. Among the 249 runs at the wall, the 110
with at least 20 000 transitions have a median completion of 0.294, the others
0.935. The top decile is also where the arcs per place jump: median 413
against 14 in the decile below, and 4 000 to 24 000 on ErlangenMainframe. The
two phenomena at the wall are distinct: small nets with a handful of dead
transitions that need a proof, not a walk (DoubleExponent 196/198, LamportFastMutEx
529/536, SieveSingleMsgMbox 745/749), and hub-dense nets where the walk itself
crawls.

ResIsolation-PT-N10P4, 445 places, 147 855 transitions, 2.06 M place-to-transition
arcs: twenty places are each consumed by 98 305 transitions. Walker alone,
4 threads, random sweep, 20 s:

| targets | steps/ms per thread | target checks per step | arc visits per step |
| ---: | ---: | ---: | ---: |
| 1 000 | 27 | 208 to 263 | 65 530 |
| 10 000 | 6 to 8 | 3 100 to 3 700 | 65 500 |
| 95 000 | 0.4 | 41 000 to 43 000 | 65 000 |

Two costs, both from the hubs, both violating the "proportional to the arcs
touched" promise of 3.5 and 3.6:

1. **Enabled set maintenance.** A hub place moving between 0 and 1 crosses the
   weight-1 threshold of its 98 305 consumers; `onPlaceChanged` visits each
   one. The weight-sorted consumer list bounds the visit to the weights
   crossed, which is everything when every consumer has weight 1. This is the
   65 000 arc visits per step, present with 1 000 targets as with 95 000, and
   it caps the walk at about 27 steps/ms where a sparse net does 20 000. Read
   arcs are not the issue: `Marking::apply` iterates the effect, so a
   transition that only tests a place never touches it.
2. **Target checks.** `fireable(t)` is desugared onto the pre-places of `t`
   (decision 7.3), so the place-to-targets index puts the target of every
   consumer of a hub on that hub's list; a step that moves the hub re-evaluates
   them all, 876 000 candidate ids at the median firing on this net. Indexing
   `fireable` targets by transition was considered and rejected: the whole
   chain reasons on state predicates, and the walker must stay that way.

Directions to weigh, none decided:

* **Hub metric.** Consumers of a place against the tokens it can hold: a
  place with 98 305 weight-1 consumers and a 0/1 marking is a hub; the arcs
  per place of the parsed net (413 at the median of the top decile) flags the
  nets where any of this matters. Only such nets pay for a mitigation.
* **Blinders per thread.** Each thread walks a sub-net where a random
  fraction of a hub's consumers is switched off (removed from its consumer
  list and never chosen). Every state reached is reachable, so verdicts stay
  sound; threads take different fractions, the shared pool (3.9) carries
  markings across them, so the exploration is no more incomplete than the
  restart schedule already makes it. Cost scales with the fraction kept.
* **Lazy hubs.** Do not update a hub's consumers on a move; mark the hub
  stale and let the strategy verify enabledness at choice time (one `pre`
  scan), drawing also among the stale hub's consumers so that newly enabled
  ones are not lost. Sound (a fired transition is verified), biased in
  distribution, no fan-out visit per step.
* **Delta target checks.** Per-target unsatisfied-atom counters, the 3.5
  scheme applied to the atoms of a conjunction: a place move flips only the
  atoms whose threshold it crosses, a target is reached when its counter hits
  zero. This keeps state predicates, turns cost 2 into the same order as
  cost 1, and needs a fallback to full evaluation for targets whose normal
  form is not a conjunction of atoms or whose atoms carry several places.
* **Done (2026-09-06): per-walker up/down index.** Each walker owns its
  target lists, split per place by the direction of change that can make a
  target hold, and drops solved ids as it meets them. On ResIsolation-PT-N10P4,
  4 threads, 20 s random sweep:

  | targets | index | steps/ms | checks/step | claimed in 20 s |
  | ---: | --- | ---: | ---: | ---: |
  | 95 000 | shared, by place | 0.4 | 42 000 | 129 |
  | 95 000 | own, up/down | 0.9 | 21 000 | 129 |
  | 95 000 | own, up/down, split in 4 | 1.5 to 3 | 5 100 | 92 |
  | 10 000 | shared, by place | 7 | 3 500 | 124 |
  | 10 000 | own, up/down | 9 to 15 | 1 600 | 125 |
  | 10 000 | own, up/down, split in 4 | 17 to 21 | 450 | 88 |

  The direction split halves the checks and doubles the step rate at equal
  claims. Splitting the set between threads (`--partition`) divides the
  checks again and buys steps, but a thread walks past the finds another
  owns and the claims drop by 30 %: off by default, kept as a knob. The arc
  visits, 65 000 per step, are untouched by either.
* **Effort share.** On the cluster this run got two walker calls totalling
  62 s of 1800; the rest was flattening and decision diagrams. Whatever the
  walker's speed, the portfolio starves it here (`PORTFOLIO.md`).

## 10. Awareness (design, 2026-09-06): firing statistics, rare events, shared restart points

The number behind this section: on ResIsolation-PT-N10P4 with 1 000
`fireable` targets, 4 threads, a 20 s random sweep, each thread fires 87 to 96
distinct transitions in about 550 000 steps, out of 147 855; the union claims
about 125 targets whatever the size of the target set. The walk is not slow
at exploring, it is not exploring: uniform choice among the enabled
transitions keeps it cycling in a region of ninety transitions. The net says
why: one token in one of its 445 places initially, one transition enabled;
399 transitions have a single input place and 147 456 have fourteen, so the
uniform walk lives among the cheap ones and almost never assembles the
fourteen-place conjunction a big one needs. The runs also report 0 resets: with `runLength` at a million steps and 27 steps/ms, a thread
never restarts inside the 30 to 60 s the cluster gives a call, so the shared
pool never receives or serves a state there. Whatever the target cost, this
is where the twenty minutes of a long run would go.

The direction: a little memory in the walker, bounded and mostly counters, so
that threads know what has been done and steer toward what has not.

### Firing statistics

* Per thread, `fired[t]`: firings of `t` in this walk (an array of |T|
  32-bit counters, 600 KB on the largest nets; `distinctFired` in the report
  already comes from it). Incremented in `fire`, O(1).
* Shared, `firedAll[t]`: the merge of every thread's counts, updated at run
  end and at every reset (O(|T|) per merge, amortised by the restart
  schedule; on hub nets where resets are rare, also every K steps). Plain
  atomics, no lock; a stale read is harmless.
* Rarity of `t` = `1 / (1 + firedAll[t])`; a never-fired enabled transition
  is a rare event.

### Rare-event choice

A strategy (`RandomStrategy` variant, or the epsilon share of the guided
ones) that prefers rarely fired enabled transitions. The enabled set can hold
tens of thousands of transitions on hub nets, so no scan: draw `k` candidates
uniformly from the enabled list (`--sample`, k about 8) and fire the least
fired of them; with probability epsilon fire the first. O(k) per step. This
is the novelty pressure of 4.5 with the memory that 4.5 lacked: it pushes the
walk out of the ninety-transition cycle toward transitions never fired, which
on a `fireable` sweep are the open targets themselves.

### Which states to share

The `SharedPool` (3.9) exists with a bounded capacity, a score and a draw at
restart. Publish a state when the walk has reason to believe it is a branching
point it could not explore:

* **a rare event just happened**: a transition with `firedAll[t] == 0` (or
  below a rank) was fired; the marking after it is new ground;
* **it enables much**: the enabled count is high relative to the walk's
  running average, or many of the enabled transitions are rarely fired (the
  sample of k gives an estimate for free);
* **the heuristic tied**: a guided strategy found several candidates at the
  best score and chose at random; the state is where guidance ran out.

The score is the count of rarely fired enabled transitions at the state
(estimated from the sample); the pool keeps the best-scored few. Memory stays
bounded by the pool capacity; there is no state hashing and no growth.

### Restarts that matter

Two changes to the restart schedule, both cheap:

* **a time-based run length**: a run ends after about a second of wall clock
  as well as after `runLength` steps, so that a 30 s call on a hub net
  restarts a few dozen times and the pool is actually used;
* **novelty stall**: a run that fired no new transition (per thread) in the
  last K steps restarts, drawing from the pool; the epsilon of guided
  strategies applies to the rare-event choice rather than to uniform random.

### Threads and related goals

Splitting a sweep by id modulo the thread count (section 9) loses the finds a
thread walks past. If a split is wanted, it should follow the support: targets
whose atoms share places (for `fireable`, transitions sharing pre-places) go
to the same thread, so that the region a thread explores is where its targets
live, and a guided strategy can aim at "the nearest of my open targets" (the
minimum over the own targets of the unsatisfied-atom count, kept incrementally
by the delta counters of section 9). Grouping by the hub places a target
needs, or by a hash of its place set, is enough for a first version.

### Measure

`distinct transitions fired` in the thread report, and claims, on the
ResIsolation 1 000-target sweep as the yardstick (today 90 and 125), then on
ErlangenMainframe and StigmergyElection from the QLA campaign. The goal is
thousands of distinct transitions in 20 s and a claims curve that keeps rising
across the minutes, not the first 20 s.
