# libHSC in the ITS-Tools chain: plan

Design document, revision 2 (2026-09-08, end of the first session).
Companion of `libHSC_in_MCC.md` (measurements), `INTEROP.md` (the
tool-to-tool protocol) and `HSC_EXPERIMENTS.md` (the campaign spec). The
**current state and next actions live in libHSC's `handoff_mcc.md`**; this
file carries the design and the decisions. Phases 0 to 4 of section 5 are
done; what remains is the campaign and the shape work of sections 7 to 9.

Goal: make libHSC a subprocess of ITS-Tools the way PetriSpot is, on the same
protocol, so that the symbolic engine replaces its-reach on the Petri net
examinations of the MCC, and test its limits there. GAL is not a target: the
s-expression formats supersede it, and cutting the Petri to GAL to
decomposition path written in Java is part of the point.

---

## 1. Where we stand

* libHSC (`~/git/libHSC`): `hsc` runs `.hsc` models (surface s-expressions,
  saturation, `select` with linear atoms, sums and weighted sums included,
  `count` with an `exact` GMP variant, `max-value`, exact `states`);
  `hsc-pn` (section 2, done) answers Petri net properties on both input
  pairs; `nupn2hsc` and `hsc-mcc` remain until the driver no longer needs
  them. Shape choices: the NUPN unit tree, FORCE reordering, Louvain
  decomposition; no configuration dominates (`libHSC_in_MCC.md`), the
  driver runs a portfolio.
* CI (done): `build_hsc.sh`, workflow on `ubuntu-24.04` (an `macos-15` one
  exists, not maintained this session), static binaries `hsc`, `hsc-pn`,
  `hsc-mcc`, `nupn2hsc`, `dve2hsc`, `fsp2hsc` deployed to the branch
  `HSC-Linux` of `yanntm/libHSC`; the test suite is the gate; GMP from the
  runner's packages. `libHSC/doc/ci.md`.
* Vendoring (done): PetriSpot's `core/`, `parse/` (MCC XML with CTL,
  s-expressions), `expr/`, `io/` as byte-exact copies under
  `include/hsc/petri/`, `vendor.sh` re-copies and checks; the edits they
  needed are upstream here (section 9).
* PetriSpot: the protocol of `INTEROP.md`: PNET binary net, s-expression
  properties over indices (`reach`, `invariant`, `deadlock`, `bound`, `ctl`
  forms), `FORMULA` lines on stdout; the MCC XML property parser
  (`parse/mcc/`, CTL included) and the s-expression reader
  (`parse/sexpr/`), printers for both syntaxes (`expr/SexprPrinter.h`).
* ITS-Tools (sections 3 and 4 done, pushed to `lip6/ITSTools`): the shared
  formats live in `interop/fr.lip6.move.gal.interop` (`KERSFormatIO`,
  `PNETFormatIO`, `SexprPropertyPrinter`), used by the PetriSpot runner and
  by `hsc/fr.lip6.move.hsc.runner` (`HscRunner`); `hsc/fr.lip6.hsc.binaries`
  downloads `hsc-pn` at Maven build time. In the MCC application `-hsc`
  starts `HscSolverRunner` beside the decision diagrams, `-hscBench` /
  `-hscBenchReduce` run it alone right after the model is read (verified
  against the oracle on Raft-PT-02 RC, RD, UB and Angiogenesis-PT-05 RC).
  its-reach stays on `ITSRunner` over GAL for CTL and LTL.

---

## 2. One tool: `hsc-pn` (done; the reference is libHSC `tools/README.md`)

One binary that eats either the contest inputs or the tool-to-tool inputs,
then does the same thing: build the model, ask the questions, print the
answers as they fall. The MCC examination protocol (examination names, the
model folder, the four StateSpace values) is a wrapper's job:
`MCC-drivers/hsc/BenchKit_head.sh` composes `hsc-pn` calls; it must not shape
the tool beyond expressivity.

```
hsc-pn (-i model.pnml | --net model.pnet) [--props FILE] [options]
```

| option | meaning |
|---|---|
| `-i FILE` | PNML P/T net; the NUPN unit tree read from `toolspecific` when present |
| `--net FILE` | PNET binary net (`INTEROP.md` section 3), places `p<i>`, transitions `t<i>` |
| `--props FILE` | MCC XML (`.xml`) or s-expression forms (anything else); `--propsSyntax=auto|mcc|sexpr` overrides |
| `--shape nupn|flat|louvain` | the hierarchy: the unit tree (flat when absent), a flat spine, Louvain clustering |
| `--force` | FORCE reordering after the shape |
| `--bound N` | leaf domain `[0, N)` (default the max initial marking plus one, at least 2) |
| `--states` | the four `STATE_SPACE` lines of the MCC examination |
| `--max-tokens` | the `MAX_TOKEN_IN_PLACE` line alone (OneSafe) |
| `--deadlock NAME` | a deadlock query without a property file |
| `--totalTime S` | wall-clock budget; unanswered properties are reported `UNKNOWN` at exit when `--printUnknown` |
| `--export-hsc FILE` | write the model (and the queries) as `.hsc`, the debugging path; what `nupn2hsc` does today |
| `-q` | quiet |

Output is the line protocol of `INTEROP.md` section 5, flushed per line:
`FORMULA <name> TRUE|FALSE TECHNIQUES DECISION_DIAGRAMS SATURATION ...`,
`FORMULA <name> <k> TECHNIQUES ...` for a bound (exact, the diagram holds
the whole set, so no `BOUND` lines), `STATE_SPACE STATES <n> TECHNIQUES ...`
for `--states`, `UNKNOWN <name>` for what the budget did not close. Exit 0
when the run ends, non-zero for errors (bad file, unresolved reference,
unknown option). A symbolic answer is a proof either way: `reach` answers
`TRUE` or `FALSE`, unlike the walker whose `FALSE` is rare; the Java reader
already honours the polarity (`INTEROP.md` section 5).

### 2.1 Properties to surface queries

Every property form maps to `select` atoms over the reachable set `R`, the
atom syntax of the surface (manual section 8), which already accepts what
the grammar of `INTEROP.md` 4.2 produces: `(and|or|not ...)`,
`(CMP EXPR EXPR)` with sums and `(* k p)`, integers, place references.

| form | query | verdict |
|---|---|---|
| `(reach N B)` | `select S R B` | TRUE iff S non-empty |
| `(invariant N B)` | `select S R (not B)` | FALSE iff S non-empty |
| `(deadlock N)` | `select S R (not (or g_1 ... g_n))`, g_t the guard of t | TRUE iff S non-empty |
| `(bound N E [K])` | `max-value` of E over R (a sum needs the lia layer's max of a sum: to confirm, else enumerate the leaf) | `FORMULA N <max>` |
| `(ctl N F)` | parsed, kept in the tree | `UNKNOWN`, until a CTL engine exists |

`fireable` is desugared from the pre-arcs as PetriSpot does. Place
references: `p<i>` is an index on both input paths; names resolve against the
net on the PNML path. The leaf names of the emitted model are the place ids,
so the translation is a printer from the property tree to surface atoms.

The MCC XML path goes through PetriSpot's parser into the same property tree
(`expr/Property.h`) and the same printer, so the two inputs meet before any
translation; `is-fireable`, `tokens-count`, `place-bound`, and the CTL
operators come for free. A round trip check exists on the PetriSpot side
(`--printProps=sexpr` from XML, re-parse, compare) and holds here too.

### 2.2 What libHSC must gain

1. **Vendored PetriSpot code** beside the loader already in
   `include/hsc/petri/` (its README lists the vendored files and the local
   edits, which stay minimal): `io/SparseMatrixIO.h`, `io/PNETIO.h` (KERS
   and PNET), `expr/Expression.h`, `Property.h`, `CtlFormula.h`,
   `CtlSimplify.h`, `Simplify.h`, `SexprPrinter.h`, `Hint.h`,
   `parse/mcc/PropertyHandler.h`, `PropertyLoader.h`,
   `parse/sexpr/Sexpr.h`, `PropertyReader.h`, `HintReader.h`,
   `parse/NetResolver.h`, `parse/PropertyFile.h`, `core/Log.h`: about 3000
   lines, header-only. Kept in PetriSpot's subfolders so their mutual
   includes stay verbatim, with `include/hsc/petri` on the include path of
   `petri_import` and the one rewrite `"core/X.h"` to `"X.h"` (the loader's
   files are flat there). `Sexpr.h` is PetriSpot's copy of our own reader;
   it stays as the property tree's reader, our `datum` as the surface's.
2. **Exact `count`.** `diagram_engine::cardinal` accumulates a `double`
   (`src/diagram.cc`); the MCC wants the integer. An exact variant beside it
   with a small unsigned big-integer (add, multiply, decimal print) and the
   same memo; `(count)` prints it, the `double` stays for the meters.
   Vasy2003 (10^21 states, off by one today) is the check.
3. **The property printer** `expr` tree to surface atoms, and the query
   driver that binds `R`, runs each query, prints the protocol lines.
   Deadlock re-enters through it (the pipeline dropped from `hsc-mcc`).
4. **`hsc-pn`** in `tools/`, on CLI11 like `hsc`; `hsc-mcc` and `nupn2hsc`
   are retired once the driver uses it (`--export-hsc` covers `nupn2hsc`).
5. Tests: the `examples/mcc` corpus gets property files (the contest's
   `ReachabilityCardinality.xml` and `ReachabilityFireability.xml` of the
   same instances) and their oracle values in `expected.csv`; the check
   runs `hsc-pn` on PNML plus XML and on PNET plus s-expressions produced by
   `petri64 --exportNet --printProps=sexpr-index`, both against the oracle.

### 2.3 Bounds and shape on the PNET path

On the PNET path there is no unit tree: `--shape louvain` or `flat` plus
`--force`. ITS-Tools could pass its own partition later (it has
`fr.lip6.move.gal.louvain` and the NUPN units of the original model,
lost after reductions); a `--shape FILE` taking a unit list is cheap to add
when measurements say it matters.

The leaf domain: the marking bound of the reduced net is unknown to the
Java side too. `hsc` raises `overflow_error` when a place grows past the
bound; `hsc-pn` catches it and reports the affected properties `UNKNOWN`
(never a wrong answer), and the driver may retry with a larger `--bound`.
A bounded-integer leaf theory (roadmap R3 of libHSC) is the real fix.

---

## 3. ITS-Tools: plugins

* **`interop/fr.lip6.move.gal.interop`** (new): the formats shared by every
  native companion: `KERSFormatIO`, `PNETFormatIO`, `SexprPropertyPrinter`
  moved out of `fr.lip6.move.petrispot.runner`, which then requires it.
  Depends on `fr.lip6.move.gal.structural` only.
* **`hsc/fr.lip6.hsc.binaries`** (new): `pom.xml` with the ant `<get>` of
  `hsc-pn` from `https://github.com/yanntm/libHSC/raw/HSC-Linux/hsc-pn` (and
  `HSC-OSX`), `BinaryToolsPlugin` copied from the PetriSpot one with the
  bundle name changed, `bin/` in `build.properties`.
* **`hsc/fr.lip6.move.hsc.runner`** (new): `HscRunner`, shaped like
  `PetriSpotWalker`: writes the PNET and the forms, runs the binary through
  `fr.lip6.move.gal.process.Runner` with a budget, reads the stream in a
  thread, fills verdicts, honours the polarity; `-Dhsc.bin=<path>` names a
  binary outside OSGi; returns null when the binary cannot be used so the
  caller falls back.
* Registration: the three modules in `fr.lip6.move.gal.parent/pom.xml`, the
  two plugins in `pnmcc/fr.lip6.move.gal.feature.pnmcc/feature.xml`, the
  runner in the `Require-Bundle` of `fr.lip6.move.gal.application.pnmcc`.

Check: `mvn -o install` in `fr.lip6.move.gal.parent` (`docs/BUILD.md`), the
product tarball carries `plugins/fr.lip6.hsc.binaries_*/bin/hsc-pn` with the
hash of the `HSC-Linux` binary.

---

## 4. The drop-in: `-hsc` beside `-its`

`Application` gets a `-hsc` flag. Where `-its` starts `ITSRunner` on the GAL
export for the reachability family (ReachabilityCardinality, Fireability,
Deadlock, UpperBounds, OneSafe, QuasiLiveness, StableMarking, StateSpace),
`-hsc` starts an `HscRunner` implementing `IRunner` on the PNET export of
the same reduced net, in the same portfolio slot: structural reductions and
SMT run before as they do today, the verdicts land in `DoneProperties` as
they stream, the runner is interrupted when everything is answered. Both
flags together run both engines side by side, which is how the comparison
is made. CTL and LTL stay on its-ctl and its-ltl.

Feature request found by the first campaign
(`libHSC_in_MCC.md`, 2026-09-08): the unfolder fuses transitions with
identical pre and post vectors without reporting how many bindings each
stands for, so an arc count over the unfolded net undercounts the coloured
semantics. Fusing is right — it is behaviourally neutral and it is what
keeps the net tractable for every engine — so the fix is a per-transition
multiplicity `m(t)` travelling with the net, not a de-fused net.

The contest value is a multigraph count, arcs labelled by transitions: two
fused bindings enabled in one state lead to the *same* successor, yet the
oracle counts them separately (our fused count is 167/202 of it), so
`TRANSITIONS = Σ_t m(t) · |{s ∈ R : s enables t}|` is exact with the
multiplicities and unanswerable without them. `m(t)` is well defined at
unfolding and does not survive transition-level structural reductions
(agglomeration and removal change the graph), which is consistent with
StateSpace running unreduced; on P/T instances every `m(t)` is 1.

Cheapest path, the harness already orchestrating the unfolding: the
unfolder writes a multiplicity file beside the unfolded model, the driver
passes it to `hsc-pn` (a `--mult FILE` flag), the tool weights its sum. A
PNET field or a PNML tool-specific annotation is the tidier home, and would
serve PetriSpot too if it ever counts arcs.

Measurement: the `itstools` MCC driver (`~/git/MCC-drivers/itstools/`) with
`-hsc` in place of `-its` on the reachability examinations of the 2026
corpus, collected with the campaign scripts of `Petri/test/mcc/`, read
against the `-its` baseline on the same product: answered, wrong, time.

---

## 5. Plan of attack (phases 0 to 4 done; the campaign is `HSC_EXPERIMENTS.md`)

Each phase ends with something that runs and is checked.

### Phase 0 (done): CI
Binaries on `HSC-Linux`; `MCC-drivers/hsc/install.sh` downloads them
instead of copying a local build.

### Phase 1 (done): vendoring and exact count (libHSC)
The files of 2.2 item 1 under `include/hsc/petri/`, README updated; exact
`count`. Check: the suite; Vasy2003 STATES equals the oracle.

### Phase 2 (done): `hsc-pn` (libHSC)
Items 3 and 4 of 2.2, the `examples/mcc` property files and the double
check of 2.2 item 5; `--states`, `--max-tokens`, deadlock. Check: every
oracle value of `examples/mcc` on both input paths; the MCC driver rewritten
on `hsc-pn` passes `run_test.pl` as the prototype did (`libHSC_in_MCC.md`).

### Phase 3 (done): plugins (ITS-Tools)
Section 3. Check: the local product runs `hsc-pn` on Airplane through
`HscRunner` from a small main, PNET byte-identical to `hsc-pn`'s own export.

### Phase 4 (done, `-hsc` and `-hscBench`): ITS-Tools
Section 4. Check: `its-tools -pnfolder . -examination ReachabilityCardinality
-hsc` answers Airplane and Angiogenesis like `-its`; then the harness
campaign.

---

## 6. Decisions taken

* Branch `HSC-Linux` (OSX and Windows later). C++23 kept, the build
  machine is constrained, not the client (static link). GMP accepted as a
  library dependency for exact counts; the double stays the default.
* One tool `hsc-pn` for both input pairs; MCC protocol is the wrapper's.
* Shared Java formats in a new `interop` plugin, used by both companions.
* No GAL frontend.
* CTL parsing vendored now, the engine later.

---

## 7. The composite machinery, and where it went

Notes from the first end-to-end runs (`-hscBench`), to steer the campaign.

**Bypassed.** Past the reductions and the property bookkeeping, the
`-hscBench` path writes a PNET and s-expressions; the GAL export,
`GALRewriter.flatten`, array-to-variable rewriting, `SumRewriter`, the Java
Louvain plugin and its veto, `CompositeBuilder` (variables into types), the
label synchronisation that makes one transition a set of matching local
transitions, the order file, and its-reach's GAL parser and
`-reachable-file` syntax are not on it. libHSC clusters from the flow
matrices itself; the hierarchy is a shape, not a type system.

**Matching is free for nets.** An event is compiled to a product term along
the shape, each separable piece on its leaf; the tree is the
synchronisation. A P/T transition is a conjunction of `place >= w` guards
and per-place increments, hence entirely separable (libHSC petri README):
no labels to invent, no nesting to keep consistent. This is a property of
nets, not of GAL: a guard over two variables or a sum is a crossing piece
and goes through the case engine.

**"Or of alts" moves into the algebra.** ITS-Tools splits a disjunctive
guard into alternative synchronisations because a composite transition is
a conjunction of local pieces. For nets the disjunctions are in the
properties, and a `select` of an `or` is a union of selections on the
fixpoint: done once on the result, not multiplied before it. Where a
disjunction must be an event (an `is-fireable` operand under a CTL
operator, later) the surface has `alt`, and the same choice returns.

**Not free: sums across components.** ITS-Tools keeps comparison supports
inside one component by feeding them to the Louvain graph (veto when it
cannot). libHSC builds the graph from the net alone and lets the case
engine resolve crossing atoms whatever the shape: more general, and the
cost is visible on Angiogenesis-PT-05 (a `select` on a sum about 0.5 s on
the 42M-state flat diagram, the all-places sum of MAX_TOKEN_PER_MARKING over
15 s). Two options, to measure: property supports as optional hyperedges of
the decomposition (a property-aware shape in the portfolio, ITS-Tools'
idea), and a better case engine for linear atoms (calculus work).

**Still owed to Java.** The structural reductions; the walk and SMT
portfolio that settles most contest properties before any engine runs
(why plain `-hsc` never fired on Raft); the hierarchy source on the PNET
path, where the NUPN units are lost and we recluster. A `--shape FILE`
taking ITS-Tools' partition would compare the two decompositions on equal
footing.

## 8. Flows, fusion, renaming: notes for the shape work

**Flows as shape information.** Ciardo's argument (a place determined by a
P-invariant sits right after its determiners, or every level in between
duplicates the pending value) survives shapes and sharpens: a semiflow's
support as a sub-sort is a small diagram of token distributions,
referenced from above by one arc; the determined place is redundancy
confined to that sort. Nested supports give nested sorts. GreatSPN's
invariant orders are the linear shadow of this. Three ways in, least
intrusive first: supports as weighted hyperedges of the Louvain graph
(the bounded-contribution mechanism exists); minimal supports as unit
seeds with Louvain arbitrating overlaps, in the portfolio; a FORCE-like
objective on invariant spans within a sort. The solver is PetriSpot's
`invariants/`, the one folder the vendoring left out: one line in
`vendor.sh`, in-process, time-capped. Do not eliminate the determined
place: its guards would become sums, the crossing atoms we pay for.

**Fusion of transitions.** ITS-Tools fuses transitions with one local
effect into one label by hand, so a synchronisation is a product of
per-component alternatives. In the calculus this is the normal form of a
sum of terms: `sum_at` folds `node(A,id) ⊕ node(A',id)` into
`node(A⊕A',id)` recursively, one summand per level mirroring the shape,
and the static saturation pass partitions events at every cut into F, L
and G, recursively. Locality and hierarchy are automatic. The mixed pair
`node(A,id) ⊕ node(id,B)` and crossing terms stay flat: the same boundary
the composite builder hits.

**Renaming.** Variables are positional, codes position-relative:
isomorphic subcomponents anywhere in the tree intern to the same terms and
diagrams. GAL names globally, so a composite shares nothing across
isomorphic components. N philosophers cost one component's nodes.

**Sums across components.** Supported exactly whatever the shape (SDD could
not), paid for today. Knobs: property-aware clustering (keeps frequent sums
local, at the price of a net-unnatural shape), and a case engine that
treats a sum constraint as a weighted-count sub-shape, the semiflow idea
from the other side. Measure before choosing.

## 9. Who depends on whom

By content the arrow is libHSC → PetriSpot: libHSC consumes the net layer
(`core/`, `parse/`, `expr/`, `io/`), PetriSpot consumes one small reader.
"PetriSpot depends on libHSC" is a product decision (one binary that also
saturates), which is fusion under another name and would create a cycle
through the vendored files. Fusion is not wanted: symbolic and explicit
have different performance disciplines, tests and failure modes; the
shared philosophy is in the formats, which are small.

A third project would carry PetriSpot minus `walk/`, `ctl/`, `lp/`,
`cli/`: the net, the formats, the parsers, `invariants/` (the shape work
wants semiflows), about 4000 header-only lines on expat alone. The Java
`fr.lip6.move.gal.interop` plugin is its other side; `KERS.md` and
`INTEROP.md` are its docs. It does not need its own repository or CI now.

Duplication is disciplined: one direction, byte-identical, scripted with
an identity check, edits upstream first. Guards to add: record the
PetriSpot commit of the snapshot in `vendor.sh`, and have the libHSC CI
fetch PetriSpot at that commit and run the check, so a patched copy fails
the build; script the reverse copy (`parse/sexpr/Sexpr.h` from libHSC's
`surface/sexpr.hh`) the same way.

When `invariants/` is needed inside libHSC, switch to a real dependency
without a new repository: PetriSpot exposes `Petri/src` as an INTERFACE
library target, libHSC fetches PetriSpot at a pinned commit at configure
time as it does sparsehash. Same arrow, no copies, no new CI.

## 10. Multiplicities: transitions now, places as weighted counting

**Transitions.** The `m(t)` of section 4 travels as a one-column KERS: PNET
already ends with the initial marking in that shape, so a fourth optional
block (a header flags bit saying it is there, old readers unaffected) or a
sibling `.kers` file does it. Sparse by nature: on P/T instances every
coefficient is 1 and the block is empty.

**Where the count is lost, measured on BART-COL-002.** The harness runs
`its-tools -pnfolder . -examination StateSpace --reduce-single STATESPACE`;
the log reports "Unfolded HLPN to a Petri net with 10865 places and 646
transitions" and no reduction rule, so the arcs went in the unfolder's
**unique tables**, which never keep a duplicate binding. Three collapse
points, all exact to instrument where they happen:
`fuseEqualParameters` (two parameters constrained `x = y` are iterated as
one: unconditional and load bearing, since the alternative enumerates the
full product to keep only its diagonal — multiplicity-neutral, because the
bindings it never generates have an unsatisfiable guard: they are not
fused transitions, they do not exist, so this rule must never carry a
coefficient), the `i<j` canonicalisation of symmetric parameters (orbit
size is a closed-form multinomial over the repeated values; not applied
under `ReductionType.STATESPACE`), and the unique table (on a hit, add the
incoming weight to the surviving entry — one line, given a table that
reports the hit).

**The counter composes part by part.** Unfolding binds one parameter at a
time, so a partially bound `HLTrans` carries a weight starting at 1: a step
that skips values on a false guard contributes nothing, a step that
collapses an orbit multiplies by its size, and `transformHLtoPT`
accumulates into the surviving entry instead of overwriting, so a unique
table hit adds. Product under composition, sum over alternatives: a
counting semiring, no global bookkeeping. 64-bit, and stale on overflow
rather than wrap (the final sum is GMP anyway). Places need none of this:
the same log shows 10865 places for 646 transitions, so places unfold in
full and place weights are purely a reduction-side notion (free-SCC
agglomeration), not a decolouring one.

The cost is plumbing, not arithmetic: `dropTransitions` (13 call sites in
`StructuralReduction`) renumbers, so permuting the weight vector there
covers most rules cheaply, while rules that *create* transitions get the
stale bit rather than an invented composed weight. Carried on
`SparsePetriNet` as the user's design: default all ones, stored as
weight-1 so 0 means 1 and an unmodified net is an empty sparse vector (no
transition can weigh 0); a marker naming what the weights are relative to
(the full unfolding), because chained transformations otherwise make them
meaningless; separate staleness bits for place and transition weights,
which different rules preserve. Check that costs nothing: a debug flag
disabling the unique table, then assert Σ m(t) equals the transition count
of the un-interned unfolding on a few COL models.

**Places.** A merged place standing for K places over which tokens travel
freely (a free SCC) is *not* a scalar factor: a state with m tokens there
represents C(m+K-1, K-1) real states, so the count is a sum over reachable
merged states of a product of per-place binomials. That is the shape our
counting already has — `cardinal_as` folds bottom-up and asks each leaf for
the size of its value set, so a per-leaf *weight function of the value*
gives weighted counting with the recursion unchanged, exact under the GMP
instantiation (about 30 lines; a flat cardinality cannot express it).

**Why it matters more than the coefficient.** StateSpace runs unreduced
today because reductions change the count. With per-leaf weights we accept
the reducer's net and still report exact STATES, which is often the
difference between answering and timing out. The four values differ:
MAX_TOKEN_PER_MARKING is untouched, MAX_TOKEN_IN_PLACE survives when the
free SCC lets tokens gather in one place, STATES needs the weights,
TRANSITIONS needs per-state arc weights and is the hardest.

**Boundary.** Exact when the fibre over a reduced state factorises into
independent per-place counts (free SCC, agglomerated place). A reduction
tying several places by an invariant gives a general polytope: Berthomieu's
setting, lattice-point counting, no closed form. So this is a cheap, well
delimited subset — closed form, composable by product, refusable when the
reducer reports a fibre it cannot factor. It needs the ITS-Tools reducer to
emit the weights, composed through chained reductions, and the surface to
carry them (a leaf-weight declaration, so a `.hsc` file stays
self-describing).

## 11. The multiplicity contract (spec, before code)

Evidence from one run (BART-COL-002, `--reduce-single STATESPACE`): the arc
count is lost at three sites, not one.

```
Unfolded HLPN to a Petri net with 10865 places and 646 transitions
Reduce places removed 10281 places and 212 transitions
Drop transitions (Redundant composition of simpler transitions.) removed 70 transitions
Reduce places removed 28 places and 20 transitions
```

Hence the counters belong to the `SparsePetriNet` contract, not to the
unfolder: every modifier inherits the obligation.

**Fields.** `tmult` per transition, `pmult` per place, each with a validity
flag. Values stored as weight minus one, so an all-ones vector is empty and
sparse (no element can weigh 0). Semantics: how many elements of the
*baseline* net this one stands for, the baseline being the net when the
counters were created (the full unfolding).

**Presence is the switch.** The record is attached at the top by the step
that knows the values are wanted — the unfolder under
`ReductionType.STATESPACE`, the examination that counts objects rather than
deciding a property — and every later step is written as "if the record is
there, maintain it, else skip this code". No flag travels through the rules,
nothing is allocated for a run that does not track (the field is null), and
a consumer weights its counts when the blocks are present.

**Java side, implemented (2026-09-08).** `SparsePetriNet` carries optional
named matrices (`getExtra`/`putExtra`/`getExtras`/`clearExtras`, null until
something is attached, copied by the copy constructor);
`SparseHLPetriNet.unfold` attaches an empty `TMULT` under STATESPACE, which
declares tracking; `StructuralReduction` copies the record in, hands it back
through `SparsePetriNet.readFrom`, and **drops it in `dropTransitions`** with
a logged reason, since a dropped transition's arcs cannot be attributed to a
survivor; `PNETFormatIO` writes whatever the net carries as named blocks.
So today a reduction that removes transitions leaves the value unanswered
(honest), and the ghost contributors below are what will make it survive.

**Obligation, three options.** Maintain, invalidate, or — for any rule
nobody has audited — invalidate by default. That default is what makes the
feature safe to land incrementally: an unaudited rule costs an unanswered
value, never a wrong one.

* Maintainable exactly: fusing identical transitions (additive); dropping a
  transition proven never enabled (no-op, it contributes no arcs);
  renumbering (permutation).
* Provably not maintainable for arc counting: "redundant composition of
  simpler transitions" (the removed transition's arcs come from the states
  where *it* was enabled, which is no survivor's enabling set), the
  agglomerations, transition splitting. These clear the flag.
* Places: constant-place removal is arc-neutral but changes STATES;
  free-SCC agglomeration is the factorisable fibre (binomial weight,
  section 10); anything invariant-tied clears the flag (Berthomieu).

**What STATESPACE actually applies** (`StructuralReduction.reduce`, the
`rt == STATESPACE` branch, "pretty basic stuff only"): constant-place
removal, no-effect transition removal (`flowPT == flowTP`), redundant
compositions. Two gates are the point: a constant place *holding tokens* is
skipped, and `ensureUnique` (duplicate fusion) is not run at all. Those
refusals exist because the consumer cannot interpret the result otherwise —
the same post-interpretation of simplified results that ITS-Tools does by
hand for the token values, historically a source of bugs.

**Ghost contributors make TRANSITIONS survive reductions.** For arcs a
transition need not survive: what is needed is its enabling predicate and
its multiplicity. A transition dropped as no-effect or as a redundant
composition can be kept as a ghost contributor (its pre-arc vector plus a
weight) and evaluated against the reachable set at the end, so its arcs are
counted exactly though it is gone from the net. Likewise a removed constant
place holding `c` tokens is an additive `c` in token sums and a candidate
`c` in the per-place maximum. This *lifts* the reducer's gates instead of
working around them: StateSpace can then be answered on a more reduced net,
not a less reduced one. (It supersedes the earlier conclusion here that
TRANSITIONS forces an unreduced run.)

**One vocabulary for all four values.** Each place carries three functions
of its value: multiplicative count weight, additive token contribution, max
contribution. Each transition carries a multiplicity. The net carries ghost
contributors and constants. Then STATES is the weighted count,
MAX_TOKEN_PER_MARKING the max of a weighted sum, MAX_TOKEN_IN_PLACE the max
over per-place functions, TRANSITIONS the weighted sum over surviving and
ghost contributors — four folds over one diagram, with the producer stating
what each object stands for instead of the consumer guessing. A closed
vocabulary of kinds (identity, constant, affine, binomial-in-K) keeps
serialisation to a kind tag plus one or two integers per object, hence a
named KERS block per kind.

**Format.** PNET keeps its three blocks; a flags bit says named blocks
follow. Each: 8-byte ASCII name, uint32 byte length (an unknown block is
skipped without parsing), then a KERS block. `TMULT`, `PMULT` to start.
Absence is the staleness signal, so no validity bit in the format: a
producer that cannot maintain omits the block.

**Done (2026-09-08).** PNET carries optional named blocks (8-byte name,
uint32 length, KERS payload; `KERS.md` and `INTEROP.md` §3), read only when
asked for, unknown names skipped, absence meaning "not available".
`hsc-pn` honours `TMULT` in its arc count and leaves transitions that cannot
change the marking out of the fixpoint, reading their guards from the net
instead (`tools/README.md`). Checked by `tests/pn_samples.sh`: weights of 2
on Raft-PT-02 double TRANSITIONS (55824 to 111648) and leave STATES at 7381;
`tests/pnet_block.py` writes a block without needing a producer.

**Next, and it wants the same treatment:** the NUPN unit tree as blocks, so
the PNET path stops losing the hierarchy (today it reclusters with Louvain).
A places-by-units incidence matrix plus a one-column parent index per unit
are two ordinary KERS payloads — one format for the net, its counting
record and its shape.

**Prototype, in this order.** (1) Weights in the unfolder only, on an
unreduced unfold, written as a standalone one-column KERS. (2) `hsc-pn
--mult FILE` weights its TRANSITIONS sum. (3) Test: BART-COL-002 (oracle
TRANSITIONS 53328, the other three values unchanged), then -005 and -010;
inside ITS-Tools assert Σ m(t) equals the number of bindings generated.
(4) Only then promote the vector to the net contract with the audit above,
and the file to a PNET named block.

## 12. Placing the counting record: what it costs where

**Producer side, least intrusive placement.** The reducer keeps deleting as
it does today and *additionally* appends what it deleted to a counting
record on the net. No rule has to honour a new kind of object, no index
moves, the verification path is unchanged; the cost is bookkeeping in the
three STATESPACE rules plus an accumulator, and default-invalidate leaves
every unaudited rule alone. Flagging transitions "counting-only" inside the
working net is the alternative and is worse: every consumer would then have
to respect the flag.

**KERS stays a dumb matrix codec — no flags there.** The framing belongs to
PNET: a header flags bit says named blocks follow; each block is an 8-byte
name, a uint32 length, then an ordinary KERS payload. Unknown names are
skipped by length; a producer emits only what it can maintain, so absence is
the staleness signal.

**One file, not side files.** The record is not only per-object scalars: the
dropped transitions carry pre-arc vectors, so those blocks reference place
indices and are meaningful against one exact net. Side files express that by
convention and fail silently when the pairing goes stale. Blocks make the
coupling structural. Side files stay acceptable for the prototype, where one
command produces net and record together.

**Blocks (all ordinary KERS payloads).**

| name | shape | content |
|---|---|---|
| `TMULT` | 1 column, T rows | weight − 1 per surviving transition (0 means 1) |
| `PCOEF` | 1 column, P rows | K − 1 per surviving place (the free-SCC binomial) |
| `PDROP` | 1 sparse column | constant markings of dropped token-holding places: their additive token contribution and their candidates for the per-place max |
| `GHOSTPT` | P × G matrix | pre-arcs of transitions dropped but still contributing arcs |
| `GHOSTMULT` | 1 column, G rows | their weights |

**Loading versus using.** Read every block you understand, always, and make
the weights part of the semantics of counting queries rather than an option:
a query that silently ignores a present weight returns a wrong number. A
caller wanting unweighted counts asks for raw explicitly — default-correct
with an opt-out.

**Staging.** Three of the four StateSpace values need no engine change:
TRANSITIONS is arithmetic over enabled-state counts (surviving transitions
and ghosts alike), the two token values arithmetic over `max-value` and the
dropped constants. Only weighted STATES touches the calculus (a per-leaf
weight in the counting fold), and it is also the only piece needing the
weights to reach the surface language instead of staying out of band in
`hsc-pn`.

## 13. PNET named blocks: framing and impact radius (measured on the code)

PNET is days old and only we produce and consume it, so there is no
versioning question: the feature is simply added, no version bump, no
header flag, no per-block "required" bit. Extensibility comes from the
names — an unknown block is skipped by its length, and anything needed
later is a new name rather than a re-framing.

**What the code already gives.** Both sides are stream-composable:
`SparseMatrixIO::read` reads one KERS block from a stream, the C++ PNET
reader (`io/PNETIO.h`) reads exactly three and stops, and on the Java side
PNET is **write-only** with a stream-based `KERSFormatIO`. So blocks
append with no change to KERS and no change to the existing header
validation (`version == 1`, `flags == 0`) — we never set the flags byte.

**Framing.** After the three mandatory blocks, zero or more of: 8 bytes
zero-padded ASCII name, uint32 payload length, then the ordinary KERS
payload. Append-only, each name at most once, canonical order so bytes stay
reproducible; read until the stream ends; an unknown name is skipped by its
length (worth a line on stderr in verbose mode). Row counts are validated
per name against P or T with the `expect` helper already there, which
catches a stale pairing at once.

**Impact radius, two files:**

| where | change | size |
|---|---|---|
| PetriSpot `io/PNETIO.h` (vendored into libHSC by `vendor.sh`) | write blocks; read them into a small extras struct; one new overload keeping both existing entry points | ~60 lines |
| ITS-Tools `PNETFormatIO.java` | optional record argument, same framing, lengths via a byte buffer | ~30 lines |
| `SparseMatrixIO.h`, `KERSFormatIO.java`, every existing call site | untouched | 0 |

`kersconv --decode-net` can dump the blocks in ~15 lines. Docs: the framing
table into `KERS.md` / `INTEROP.md` §3.

**Residual risk, noted not engineered.** A deployed *old* binary reading a
net that carries weights ignores them and prints an unweighted count. Not
worth machinery: the only value affected is TRANSITIONS on coloured
instances, which we leave unanswered today, and the deploy chain is ours
end to end.
