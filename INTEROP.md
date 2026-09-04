# ITS-Tools and PetriSpot: tool-to-tool reachability

Status: design document, revision 2 (2026-09-04). Revision 1 was discussed
with the user; section 10 records the decisions. Phases 0 to 6a (section 9)
are implemented: the PetriSpot side here, the Java side in ITS-Tools
(`petrispot/fr.lip6.move.petrispot.runner`). Sections 3 to 5 are the
reference for the formats until they move next to `KERS.md`.

Goal: let ITS-Tools (Java) delegate its explicit walks to PetriSpot with the
same ease it already delegates invariant computation, so that
`fr.lip6.move.gal.structural.RandomExplorer` and its call sites in
`ReachabilitySolver`, `DeadlockSolver`, `AtomicReducer(SR)` and later
`UpperBoundsSolver` become calls to the PetriSpot binary. The exchange must be
cheap enough to be invoked ten or twenty times per ITS-Tools run on nets with
10^5 transitions: binary net, index-based properties, streamed results.

---

## 1. Where we stand

### 1.1 PetriSpot

* Invariants over KERS matrices (`--loadKERS`, `--basisKERS`, `KERS.md`).
  This is the existing tool-to-tool path and its conventions (binary,
  little-endian, sparse columns, indices instead of names) are the model for
  everything below.
* Walk engine (`Petri/src/walk/`, `WALK_PLAN.md`): PNML in, MCC property XML
  in, `FORMULA` lines out. One target per walk; the driver schedules open
  properties in rounds (`--totalTime`). Strategies random, bestfirst,
  structural, relaxed; N threads racing on one target; a shared restart pool.
* MCC harness: `~/git/MCC-drivers/petrispot/BenchKit_head.sh` drives the
  binary on MCC inputs (PNML plus property XML) and checks verdicts against
  the oracles; `petrispotxred` runs the ITS-Tools reducer first.

### 1.2 ITS-Tools

* `petrispot/fr.lip6.move.petrispot.runner`: `PetriSpotRunner` writes the
  incidence matrix as KERS (`KERSFormatIO`, both directions), runs the binary
  found by `BinaryToolsPlugin` through `fr.lip6.move.gal.process.Runner`, reads
  the basis back. Temp files, timeout, exit code handling and a one-entry
  cache are there. The walk runner is a second entry point in this plugin,
  not a new one.
* `Runner.runTool` either inherits stdout or redirects it to a file and waits
  for exit; `ITSRunner` shows the other pattern, a thread reading
  `process.getInputStream()` line by line. Streaming is what we want for
  verdicts (section 5).
* Export to MCC formats exists (`MCCExporter`: `StructuralToPNML` plus
  `PropertiesToPNML`), but it is the slow path: XML, names, and integer
  constants of the properties turned into extra places to stay
  MCC-compliant. Good for debugging and for the `petrispotxred` harness, not
  for the hot loop.

### 1.3 What the call sites need

| Call site | Today | Needs from the interface |
|---|---|---|
| `ReachabilitySolver.randomCheckReachability` | 4 threads random over all open predicates (steps, 30 s); then best-first on each predicate (5 s each, capped at 50 misses); then probabilistic / exhaustive Bloom walk | many predicates in one request, all checked at every step; a per-predicate heuristic phase; budgets in steps and seconds; verdict per predicate |
| `ReachabilitySolver.tryReplayParikh` | `runGuidedReachabilityDetection` per distinct Parikh vector, 100 x |parikh| steps, 1 to 30 s | Parikh vector and partial order per predicate as hints (later, section 6.3) |
| `DeadlockSolver` | random, then Parikh-guided deadlock walks | deadlock target; Parikh hint later |
| `AtomicReducer`, `AtomicReducerSR` | random 100k steps, 30 s over up to hundreds of candidate atoms | many predicates in one request |
| `UpperBoundsSolver` | max mode: track the largest value of each expression; `CoverWalker` | bound targets returning the max seen (later, section 6.4) |

Verdicts on the Java side are `int[]`, one per predicate: 1 if a witness was
seen (or the max value in bound mode). Predicates are already normalised by
Java to "find a state where this holds" (`EF phi` gives `phi`, `AG psi` gives
`not psi`), `is-fireable` is already desugared, and transitions of the walked
net map to the SMT representatives through `repr`; that mapping stays in Java,
anything sent to PetriSpot is expressed over the transmitted net.

---

## 2. The exchange in one picture

```
ITS-Tools                                   PetriSpot (one process per request)
  net (reduced, in memory)  --PNET file-->   --net=req.pnet
  predicates                --sexpr file-->  --props=req.sexpr
  budget                    --CLI flags-->   --threads=4 --totalTime=30 --walkSteps=1000000 -q
  verdicts                  <--stdout----    FORMULA lines, flushed as they fall
```

One request is one process. Loading a binary net of a few million arcs takes
tens of milliseconds, well under the walk budgets (seconds), and the net
changes between ITS-Tools iterations anyway (reductions), so a persistent
server process is a non-goal. Files rather than pipes: they are debuggable
(`DEBUG` keeps them, as `PetriSpotRunner` does today) and match the existing
runner.

The exchange is tool to tool: nothing in it is meant to be read by a person.
On the PNET path places are `p0, p1, ...` and transitions `t0, t1, ...` by
construction, and that is the traceability; whoever wants readable traces
uses the PNML path. Properties carry a name because it costs nothing and
ITS-Tools has one (`prop0, prop1, ...` or the MCC id).

---

## 3. Net: the PNET container

A KERS matrix carries one matrix and the incidence matrix loses read arcs, so
the net needs `flowPT`, `flowTP` and the initial marking. PNET is a header
followed by three KERS blocks, each exactly as specified in `KERS.md` (16-byte
header, non-empty columns in ascending order, `0xFFFFFFFF` terminator), so
both sides reuse their KERS reader and writer unchanged. It is a separate
format with its own magic: a KERS reader stays a matrix reader.

All integers little-endian.

| Offset | Size | Type | Content |
|---|---|---|---|
| 0 | 4 | char[4] | magic `P` `N` `E` `T` |
| 4 | 1 | uint8 | version `1` |
| 5 | 1 | uint8 | flags, reserved, must be zero |
| 6 | 4 | uint32 | number of places `P` |
| 10 | 4 | uint32 | number of transitions `T` |
| 14 | 2 | | padding, zero |
| 16 | | KERS | `flowPT`, `P` rows x `T` columns (column `t` = pre-arcs of `t`) |
| | | KERS | `flowTP`, `P` rows x `T` columns (post-arcs) |
| | | KERS | initial marking, `P` rows x 1 column |

Places and transitions are identified by their index in this file, on both
sides, for properties and traces. The loaded net names them `p<i>` and
`t<i>`. Values are int64 in the file as in KERS; `petri32` warns and
truncates as it does for KERS. The dimensions of the three blocks must agree
with the header; a mismatch is an error.

Size: about `12 x arcs + 8 x T + 16 x 3 + 16` bytes; a net with 500k
transitions and 2M arcs is about 50 MB and loads in tens of milliseconds
(the PNML of the same net is several hundred MB and parses in seconds).

PetriSpot side: `--net=<file.pnet>` as an alternative to `-i model.pnml`;
`--exportNet=<file.pnet>` from any loaded net, which is also the round-trip
test (PNML to PNET to walk gives the same verdicts as PNML to walk).
`kersconv --decode-net` prints a PNET as text for inspection.

Java side: `PNETFormatIO.write(ISparsePetriNet, Path)`: header plus three
`KERSFormatIO.write(matrix, stream)` calls, the marking as a one-column
`IntMatrixCol` built from `getMarks()`. About 40 lines.

---

## 4. Properties: s-expressions over indices

### 4.1 Why not the MCC XML

The XML path exists on both sides but costs what we want to avoid: an XML
writer that mangles constants into places, expat, name resolution, and files
whose size is dominated by markup. Property files are kilobytes; what matters
is a parser that is a few hundred lines, allocation-light and total, and a
printer that is trivial to write in Java. S-expressions are that, they are the
syntax `WALK_PLAN.md` 3.2 already planned for `parse/sexpr/`, and libHSC uses
the same shape (`(and ...) (or ...) (not ...) (<= e e) (+ e e)`) and the same
reader (`datum`: atom or list, with a source line), copied here so the two
projects stay compatible without depending on each other.

### 4.2 Grammar

One top-level form per property; `;` starts a comment to end of line;
whitespace is free. A NAME is an atom or a `"quoted string"`.

```
FORM  ::= (reach NAME BEXP)          ; find a marking satisfying BEXP; verdict TRUE  (EF)
        | (invariant NAME BEXP)      ; find a marking violating BEXP;  verdict FALSE (AG)
        | (deadlock NAME)            ; find a dead marking;             verdict TRUE
        | (bound NAME EXPR [INT])    ; maximise EXPR; the optional INT is a known upper bound
BEXP  ::= true | false
        | (and BEXP+) | (or BEXP+) | (not BEXP)
        | (CMP EXPR EXPR)            ; CMP in == != <= >= < >
        | (fireable TREF+)           ; some listed transition is enabled
EXPR  ::= INT | PREF
        | (+ EXPR+) | (- EXPR EXPR) | (* INT EXPR)
PREF  ::= p<digits> | NAME           ; place by index, or by name
TREF  ::= t<digits> | NAME           ; transition by index, or by name
```

`p<digits>` and `t<digits>` are indices into the net; any other atom is a
name resolved against the net (a place really named like `p12` in a PNML net
must be quoted). ITS-Tools writes indices only. A `(CMP EXPR EXPR)` is
linearised into one `LinearAtom` (`sum coeff x place cmp constant`); an
expression that is not linear (product of two places) is a parse error, which
is fine: MCC has none and neither does ITS-Tools after `toPredicates`.
`fireable` is desugared from `flowPT` into the disjunction of conjunctions of
`place >= weight`, as the MCC parser does today; the AST does not change.

A `bound` target is pure exploration: the walk reports the largest value of
`EXPR` it saw. The optional integer is a hint the caller knows to be an upper
bound (structural, from invariants or the state equation; it may
over-approximate but is never below the truth): reaching it ends the target
with a `FORMULA` verdict, since the value is then exact. Without it the target
stays open for the whole budget and only a `BOUND` line is reported. Where the
hint comes from is the caller's business; PetriSpot does not compute it.

Example request written by ITS-Tools for three predicates, a deadlock and two
bounds, one with a structural hint:

```
(reach prop0 (and (>= (+ p12 (* 2 p40)) 5) (not (fireable t7 t9))))
(reach prop1 (<= (- p3 p4) 0))
(invariant prop2 (<= p1 1))
(deadlock prop3)
(bound prop4 (+ p5 p6) 3)
(bound prop5 p9)
```

### 4.3 Printer and converter

`--printProps=sexpr` prints a property set in this syntax; combined with the
MCC parser it converts MCC XML to s-expressions, which is how the parser is
tested: parse an MCC file, print, re-parse, compare ASTs, over the property
files of a few MCC models. `--printProps=infix` keeps the current
human-readable dump.

Java side: `SexprPropertyPrinter implements ExprVisitor<Void>` after the
model of `CExpressionPrinter`, about 80 lines, emitting indices.

### 4.4 Hints: `--hints=<file>`

A hint does not change the question, it enables more strategies. Hints come
as a separate input in the same s-expression syntax, each form naming the
property it helps (a name absent from `--props` is a warning):

```
HINT  ::= (parikh NAME (TREF INT)+)    ; a Parikh vector: transition, firing count
```

The Parikh vector comes from the SMT state equation over *representatives*:
transitions with identical effect are one variable there. PetriSpot groups
transitions by identical effect itself, so a count on any member of a class
is the class count, and firing any member consumes it. ITS-Tools sends the
vector as the solver returned it, keyed by the representative transitions of
the walked net, and never sends its `repr` mapping.

The known upper bound of a `bound` target is not a hint in this sense: it is
part of the form (4.2), since it decides when the target is answered.

---

## 5. Results: a line protocol on stdout

Every record is one line, starts with a keyword, and is flushed when written.
Anything else on stdout is a log and may be ignored. Nothing is written twice
for the same property. Property names are printed verbatim and must not
contain whitespace (MCC ids and ITS-Tools names never do).

```
FORMULA <name> TRUE|FALSE TECHNIQUES <words>     verdict (unchanged MCC line)
FORMULA <name> <k> TECHNIQUES <words>            bound target whose hint k was reached: the bound is k
BOUND <name> <max>                               every bound target, at exit: the largest value seen
WITNESS <name> <k> t3 t9 t3 ...                  only with --trace: the k transitions fired from the initial marking
UNKNOWN <name>                                   only with --printUnknown: at exit, one per property without verdict
```

A bound target without a hint never gets a `FORMULA` line nor an `UNKNOWN`
line; its answer is the `BOUND` line.

Witnesses are not part of the fast path: without `--trace` the walker never
records anything, and the MCC does not ask for traces. With `--trace` the
trace is recorded, replay-verified and printed, indices on the PNET path and
names on the PNML path. The Java side trusts the verdicts; the guarantee that
PetriSpot does not lie is its own testing, not a replay on the Java side.

Exit code 0 when the run ends normally, whether or not anything was solved;
non-zero only for errors (bad file, unknown option, unresolved reference).

Java side reads stdout in a thread (the `ITSRunner` pattern) and records
every `FORMULA` into `DoneProperties` as it arrives. Verdicts already
recorded survive a kill on timeout, which today's redirect-to-file-then-parse
would not guarantee.

---

## 6. What changes in PetriSpot

### 6.0 CLI on CLI11, `cli/` folder

`Petri.cpp` parses options by hand and mixes parsing, the invariant driver and
the walk driver. Before adding `--net`, `--exportNet`, `--printUnknown` and a
syntax switch for `--props`, move to CLI11:

* Single header (v2.7.2), vendored under `Petri/third_party/CLI11/` with its
  BSD-3 licence file; no build or link dependency, no network at configure
  time (FetchContent would break the static CI build). Only the driver
  includes it.
* Both `--opt=value` and `--opt value` accepted, typed values, generated
  `--help`, validation of enumerations (`--strategy`), grouping of options
  (invariants, walk, I/O). Existing flags keep their spelling so
  `PetriSpotRunner` and the MCC driver do not change.
* `cli/Options.h`: one struct with every option and the CLI11 setup;
  `cli/InvariantDriver.h` and `cli/WalkDriver.h` take the current bodies of
  `Petri.cpp`, which becomes a thin `main`. `cli/README.md` documents the
  layout; `README.md` keeps the user table.

### 6.1 Loaders and printers

* `io/PNETIO.h`: PNET read and write, built on `SparseMatrixIO` stream
  variants (today it takes a file name only).
* `parse/sexpr/Sexpr.h` (the libHSC `datum` reader, about 100 lines) and
  `parse/sexpr/PropertyReader.h` (datum to `Property`, about 200 lines).
  `--props` picks the syntax by extension (`.xml` MCC, anything else
  s-expressions); `--propsSyntax=mcc|sexpr` overrides.
* `expr/SexprPrinter.h`: s-expression printer over indices or names.

### 6.2 Multi-target walking (implemented)

Before this work a `Target` was one goal and the portfolio raced threads on it. The
ITS-Tools use is dozens to hundreds of predicates on one net (one-safeness,
quasi-liveness and the like produce one query per net object), so the
engine must be able to check many targets in one walk. On the other hand
the directed strategies are built to solve one given property, and in some
cases focusing on a single property is what works. Which policy to apply
when is a matter of profiling; the architecture stays flexible:

* `TargetSet`: the open targets, each with an atomic `solved` flag (the
  "verdicts" section of the Knowledge Base of `WALK_PLAN.md` 3.9). A step
  checks every open target whose atoms touch a place the fired transition
  changed: a sparse `place -> targets` index built once, so the check stays
  proportional to the arcs touched and to the targets that depend on them,
  never to the number of targets.
* A walker aims its heuristic at one target (its focus) but every open
  target reached on the way is claimed; the focus rotates at restarts over
  the still-open targets.
* The driver schedule is a policy on top: a sweep round of random
  multi-target walks on all threads (`--sweepTime`, what ITS-Tools does by
  hand today), then heuristic rounds per open target with the growing
  per-property budget of `--totalTime`, in milliseconds so that many
  properties under a short budget still each get a walk. `FORMULA` lines
  stream as targets fall. Focus rotation inside a run is not done: a walker
  keeps its focus, the driver picks the next one. Options select the policy;
  profiling on the MCC harness decides the defaults.

### 6.3 Parikh strategy (implemented)

The hint-driven walk of `RandomExplorer.runGuidedReachabilityDetection`, as
the `parikh` strategy: the remaining count per effect class starts at the
hint; a step chooses uniformly among the enabled transitions whose class
still has count and decrements it; when no such transition is enabled the
run restarts. (The partial order the SMT side can produce is not used: it is
practically disabled there, too expensive.) The restriction is relaxed with probability
growing with the number of restarts (one per mille per restart, the Java
decay), so a slightly wrong vector still guides. A focus with a hint runs the
`--hintStrategies` pool (default `parikh,parikh,relaxed,bestfirst`) instead
of `--strategies`; a `parikh` spec on a target without hint is a random walk.

### 6.4 Bound targets (implemented)

`PropertyKind::Bound`: a linear form to maximise and an optional hint. In the
walk a bound is a `Target` with the form and a limit (the hint, or none); it
lives in the `TargetSet` like the others, with a shared atomic running maximum
per target updated by the walkers when they improve on it (each walker keeps
its own last published value to stay off the atomic). The place-to-targets
index covers bound forms, so the value is re-read only when a place of the
form changed. Reaching the limit claims the target as any other goal.

Heuristics are inherited from reachability: with a hint K the focus goal is
the atom `form >= K`, and best-first, structural and relaxed strategies work
unchanged. Without a hint the goal moves: `BoundDistance` aims at the running
maximum plus one, so best-first climbs; structural and relaxed strategies
need a fixed atom and fall back to that best-first. The shared restart pool
naturally keeps the highest markings as restart points.

Not done: witnessing unboundedness (a pumping segment along the path, what
`CoverWalker` does in Java); the Java side keeps its structural methods and
`CoverWalker`, PetriSpot reports the best value seen.

### 6.5 Documentation

`KERS.md` gets a PNET section, the s-expression grammar and the result
protocol live in this document until frozen, then move next to `KERS.md`;
`parse/sexpr/README.md`, `cli/README.md`, `walk/algorithm.md` (multi-target
step) and `README.md` follow the code.

---

## 7. What changes in ITS-Tools

All in `petrispot/fr.lip6.move.petrispot.runner`, beside `PetriSpotRunner`:

* `PNETFormatIO` (section 3) and `SexprPropertyPrinter` (section 4.3).
* `PetriSpotWalker`: the replacement API, deliberately shaped like
  `RandomExplorer` so call sites change one line:
  `int[] runReachability(ISparsePetriNet net, List<Expression> predicates,
  long steps, int seconds)`, `boolean runDeadlock(...)`, later hints and
  `int[] runBounds(...)`. It writes the two files, runs the binary with
  `--threads`, `--walkSteps`, `--totalTime`, `-q`, reads the stream, fills
  the verdict array, and falls back to `RandomExplorer` if the binary is
  missing or fails.
* Call sites: `ReachabilitySolver.randomCheckReachability` (one request
  replaces the random plus best-first phases; the probabilistic and
  exhaustive Bloom walks stay on `RandomExplorer`), `DeadlockSolver`,
  `AtomicReducer(SR)`; then, with hints, `tryReplayParikh`; last
  `UpperBoundsSolver`. Technique words in `DoneProperties` come from the
  `FORMULA` line.

This work happens in the ITS-Tools repository, after the PetriSpot side is
validated on the MCC harness; it is listed here so the plan is whole.

---

## 8. What stays in Java, and other honest limits

* The probabilistic (Bloom filter) and exhaustive walks that can prove
  absence of a witness on small state spaces: stateful exploration, outside
  the memoryless engine by design. ITS-Tools keeps them; nothing prevents a
  bounded explicit search in PetriSpot later, it is simply not this plan.
* `CoverWalker` (coverability with omega) for bounds.
* Structural reductions, SMT, the `repr` mapping and everything that
  produces hints.
* PetriSpot never proves unreachability; every verdict it emits is a
  witness. The Java side's `AG` true / `EF` false verdicts keep coming from
  its own methods.

---

## 9. Plan of attack

Each phase ends with something that runs and is checked; line counts are
rough. Phases 0 to 4 are done (2026-09-04); 5 and 6 follow in this order.

### Phase 0 (done): CLI11 and `cli/`

Vendor CLI11, `cli/Options.h`, split the drivers out of `Petri.cpp`, same
flags. Check: the invariant and walk outputs on `Petri/examples/` are
identical before and after (fixed seed); the MCC driver and `PetriSpotRunner`
run unchanged.

### Phase 1 (done): s-expression properties

Reader, property reader, printer, `--props` by extension, `--printProps=`
syntax. Check: MCC XML to sexpr to AST equals MCC XML to AST on the property
files of the development models; hand-written `.sexpr` files for
`Petri/examples/` nets.

### Phase 2 (done): PNET

`io/PNETIO.h`, `--net`, `--exportNet`, `kersconv --decode-net`; `KERS.md`
section. Check: PNML to PNET to walk reproduces the PNML walk verdicts and
step counts with a fixed seed on the development models.

### Phase 3 (done): result protocol

`WITNESS` line, `--printUnknown`, flush discipline, exit codes. Check: the
MCC driver still passes the harness.

### Phase 4 (done): multi-target walking

`TargetSet`, the place-to-targets index, the sweep round, millisecond
budgets. Check: on Angiogenesis the set of verdicts of a fixed-seed random run
equals the per-property baseline; a one-second sweep answers 12 of the 16
fireability properties at 8700 steps/ms with 3 goal evaluations per step.
ThreadSanitizer clean with 4 threads.

### Phase 5 (done): Java runner and call sites (ITS-Tools repository)

`PNETFormatIO`, `SexprPropertyPrinter` and `PetriSpotWalker` in the existing
wrapper plugin `fr.lip6.move.petrispot.runner`; call sites
`ReachabilitySolver.randomCheckReachability` (one request replaces the random
and best-first phases), `AtomicReducer`, `AtomicReducerSR`, `DeadlockSolver`
(random walks). `PetriSpotWalker.USE_PETRISPOT` is the compile-time switch
back to the Java walker; every entry point returns null when the binary cannot
be used and the caller falls back. `-Dpetrispot.bin=<path>` names the binary
outside OSGi. Checked: the PNET written from Java is byte-identical to
PetriSpot's own export of the same PNML (Airplane, Angiogenesis), and a
request with one `p >= 1` predicate per place plus an unreachable one solves
exactly the reachable ones. Still to do: run the `itstools` MCC driver on the
reachability examinations and compare with today's results. Found on the
way: the structural plugin's SAX `PTNetHandler` added post arcs with place
and transition swapped (fixed; production loads nets through the `nupn`
reader, which was correct).

### Phase 6a (done): bounds

`bound` targets (6.4) on both sides: `PetriSpotWalker.runBounds` sends the
bound expressions with the structural bound as hint when ITS-Tools knows one,
and `UpperBoundsSolver.randomCheckReachability` uses it in place of the
random and best-first phases. Checked through the MCC harness on the
reinstalled release: Angiogenesis and Airplane UpperBounds all match the
oracle, the walker reporting the 16 bounds in 20 to 35 ms.

### Phase 6b: hints

`--hints`, the Parikh strategy (6.3), then their Java call sites
(`tryReplayParikh`, `DeadlockSolver` guided walks).

---

## 10. Decisions (2026-09-04)

1. PNET is its own format with its own magic; KERS is not extended, its
   reuse as a block is the proof it serves its purpose.
2. No name table, no readability goal in the exchange: PNET nets are named
   `p<i>` / `t<i>`, pretty traces go through the PNML path. Property names
   are kept because they cost nothing.
3. Hints (Parikh vectors, known bounds) are a side input, not part of the
   query, and come later: integrate the random walks first, then the
   directed heuristics, then hints.
4. Traces and witnesses are optional and off the fast path; no replay or
   checking on the Java side, trust plus thorough testing of PetriSpot.
5. Multi-target versus single-focus is a policy question to be settled by
   profiling; the architecture supports both.
6. The s-expression reader is the libHSC `datum` reader, copied (no
   dependency between the projects), so the syntaxes stay compatible.
