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
