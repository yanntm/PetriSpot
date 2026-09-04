# ITS-Tools and PetriSpot: tool-to-tool reachability

Status: design document, revision 1 (2026-09-04). Nothing in sections 3 to 5
is implemented yet; section 1 records what exists. The user peruses, comments,
and code follows validation.

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
* `Petri.cpp` (661 lines): hand-rolled `substr` option parsing, invariants and
  walk driver in one file.
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
  verdicts (1.3).
* Export to MCC formats exists (`MCCExporter`: `StructuralToPNML` plus
  `PropertiesToPNML`), but it is the slow path: XML, names, and integer
  constants of the properties turned into extra places to stay
  MCC-compliant. Good for debugging and for the `petrispotxred` harness, not
  for the hot loop.

### 1.3 What the call sites need

| Call site | Today | Needs from the interface |
|---|---|---|
| `ReachabilitySolver.randomCheckReachability` | 4 threads random over all open predicates (steps, 30 s); then best-first on each predicate (5 s each, capped at 50 misses); then probabilistic / exhaustive Bloom walk | many predicates in one request, all checked at every step; a per-predicate heuristic phase; budgets in steps and seconds; verdict per predicate |
| `ReachabilitySolver.tryReplayParikh` | `runGuidedReachabilityDetection` per distinct Parikh vector, 100 x |parikh| steps, 1 to 30 s | Parikh vector and partial order per predicate as hints |
| `DeadlockSolver` | random, then Parikh-guided deadlock walks | deadlock target, Parikh hint |
| `AtomicReducer`, `AtomicReducerSR` | random 100k steps, 30 s over up to hundreds of candidate atoms | many predicates in one request |
| `UpperBoundsSolver` | max mode: track the largest value of each expression; `CoverWalker` | bound targets returning the max seen; a known upper bound to stop at |

Verdicts on the Java side are `int[]`, one per predicate: 1 if a witness was
seen (or the max value in bound mode). Predicates are already normalised by
Java to "find a state where this holds" (`EF phi` gives `phi`, `AG psi` gives
`not psi`), `is-fireable` is already desugared, and transitions of the walked
net map to the SMT representatives through `repr`; that mapping stays in Java,
hints are expressed over the transmitted net.

---

## 2. The exchange in one picture

```
ITS-Tools                                   PetriSpot (one process per request)
  net (reduced, in memory)  --PNET file-->   --net=req.pnet
  predicates + hints        --sexpr file-->  --props=req.sexpr
  budget                    --CLI flags-->   --threads=4 --totalTime=30 --walkSteps=1000000 -q
  verdicts                  <--stdout----    FORMULA / WITNESS / STATE / BOUND lines, flushed
```

One request is one process. Loading a binary net of a few million arcs takes
tens of milliseconds, well under the walk budgets (seconds), and the net
changes between ITS-Tools iterations anyway (reductions), so a persistent
server process is a non-goal. Files rather than pipes: they are debuggable
(`DEBUG` keeps them, as `PetriSpotRunner` does today) and match the existing
runner; `-` for stdin can come later if ever wanted.

---

## 3. Net: the PNET container

A KERS matrix carries one matrix and the incidence matrix loses read arcs, so
the net needs `flowPT`, `flowTP` and the initial marking. PNET is a header
followed by three KERS blocks, each exactly as specified in `KERS.md` (16-byte
header, non-empty columns in ascending order, `0xFFFFFFFF` terminator), so
both sides reuse their KERS reader and writer unchanged.

All integers little-endian.

| Offset | Size | Type | Content |
|---|---|---|---|
| 0 | 4 | char[4] | magic `P` `N` `E` `T` |
| 4 | 1 | uint8 | version `1` |
| 5 | 1 | uint8 | flags (below) |
| 6 | 4 | uint32 | number of places `P` |
| 10 | 4 | uint32 | number of transitions `T` |
| 14 | 2 | | padding, zero |
| 16 | | KERS | `flowPT`, `P` rows x `T` columns (column `t` = pre-arcs of `t`) |
| | | KERS | `flowTP`, `P` rows x `T` columns (post-arcs) |
| | | KERS | initial marking, `P` rows x 1 column |
| | | names | only if flag `NAMES`: `T` then `P` strings, each uint32 length + UTF-8 |

Flags: bit 0 `NAMES` (name table present), bit 1 `SAFE` (the sender knows the
net is 1-safe; a hint, never trusted for correctness). Other bits reserved,
must be zero.

Conventions: places and transitions are identified by their index in this
file, on both sides, for properties, hints, traces and witness markings. Names
are optional and only serve human-readable output (`--trace` with names);
ITS-Tools never sends them. Values are int64 in the file as in KERS; `petri32`
warns and truncates as it does for KERS.

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

## 4. Properties and hints: s-expressions over indices

### 4.1 Why not the MCC XML

The XML path exists on both sides but costs what we want to avoid: an XML
writer that mangles constants into places, expat, name resolution, and files
whose size is dominated by markup. Property files are kilobytes; what matters
is a parser that is a few hundred lines, allocation-light and total, and a
printer that is trivial to write in Java. S-expressions are that, they are the
syntax `WALK_PLAN.md` 3.2 already planned for `parse/sexpr/`, and libHSC uses
the same shape (`(and ...) (or ...) (not ...) (<= e e) (+ e e)`), so the
family stays coherent.

### 4.2 Grammar

One top-level form per property or hint; `;` starts a comment to end of line;
whitespace is free. Names are atoms or `"quoted strings"`.

```
FORM  ::= (reach NAME BEXP)          ; find a marking satisfying BEXP; verdict TRUE  (EF)
        | (invariant NAME BEXP)      ; find a marking violating BEXP;  verdict FALSE (AG)
        | (deadlock NAME)            ; find a dead marking;             verdict TRUE
        | (bound NAME EXPR [INT])    ; report the largest EXPR seen; stop at INT if given
        | (parikh NAME (TREF INT [INT])+)   ; hint: transition, count, optional rank
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
name resolved against the net, and a real place named like `p12` must be
quoted. ITS-Tools writes indices only. A `(CMP EXPR EXPR)` is linearised into
one `LinearAtom` (`sum coeff x place cmp constant`); an expression that is not
linear (product of two places) is a parse error, which is fine: MCC has none
and neither does ITS-Tools after `toPredicates`. `fireable` is desugared from
`flowPT` into the disjunction of conjunctions of `place >= weight`, as the MCC
parser does today; the AST does not change.

`bound` and `parikh` are new property and hint kinds; `PropertyKind` gains
`Bound`, `Property` gains an optional hint. A `parikh` form names the property
it belongs to and may appear before or after it in the file; a hint for an
unknown property is a warning, not an error.

Example request written by ITS-Tools:

```
(reach p3 (and (>= (+ p12 (* 2 p40)) 5) (not (fireable t7 t9))))
(invariant p5 (<= p1 1))
(deadlock dl)
(parikh p3 (t7 2 0) (t12 1 1) (t3 1 1))
```

### 4.3 Printer and converter

`printProps=sexpr` prints a property set in this syntax; combined with the MCC
parser it converts MCC XML to s-expressions, which is how the parser is
tested: parse an MCC file, print, re-parse, compare ASTs, over the property
files of a few MCC models. `printProps=infix` keeps the current
human-readable dump.

Java side: `SexprPropertyPrinter implements ExprVisitor<Void>` after the
model of `CExpressionPrinter`, about 80 lines, emitting indices.

---

## 5. Results: a line protocol on stdout

Every record is one line, starts with a keyword, and is flushed when written.
Anything else on stdout is a log and may be ignored. Nothing is written twice
for the same property.

```
FORMULA <name> TRUE|FALSE TECHNIQUES <words>     verdict (unchanged MCC line)
WITNESS <name> <k> t3 t9 t3 ...                  the k fired transitions, from the initial marking
STATE <name> p1:3 p7:1                           the witness marking, sparse, places in ascending index
BOUND <name> <max seen>                          bound targets, at exit (FORMULA instead when the known bound is hit)
UNKNOWN <name>                                   at exit, one per property without verdict, with --printUnknown
```

`WITNESS` needs `--trace` (recorded and replay-verified before printing, as
today); `STATE` is printed whenever a witness exists. With `--names` both use
the net's names instead of indices (PNML input, or PNET with a name table).
Exit code 0 when the run ends normally, whether or not anything was solved;
non-zero only for errors (bad file, unknown option, unresolved reference).

Java side reads stdout in a thread (the `ITSRunner` pattern), records every
`FORMULA` into `DoneProperties` as it arrives, and keeps `WITNESS` for an
optional replay on its own net when `DEBUG` is on. Verdicts already recorded
survive a kill on timeout, which today's redirect-to-file-then-parse would
not guarantee. `Runner` gets a variant that hands the stream to a line
consumer, or the runner uses `ProcessBuilder` directly as `ITSRunner` does.

---

## 6. What changes in PetriSpot

### 6.0 CLI on CLI11, `cli/` folder

`Petri.cpp` parses options by hand and mixes parsing, the invariant driver and
the walk driver. Before adding `--net`, `--exportNet`, `--printUnknown`,
`--names` and a syntax switch for `--props`, move to CLI11:

* Single header, vendored under `Petri/third_party/CLI11/` with its
  BSD-3 licence file; no build or link dependency, no network at configure
  time (FetchContent would break the static CI build). Only `Petri.cpp`
  includes it, so compile time is unaffected elsewhere.
* Both `--opt=value` and `--opt value` accepted, typed values, generated
  `--help`, validation of enumerations (`--strategy`), grouping of options
  (invariants, walk, I/O, formats). Existing flags keep their spelling so
  `PetriSpotRunner` and the MCC driver do not change.
* `cli/Options.h`: one struct with every option and the CLI11 setup;
  `cli/InvariantDriver.h` and `cli/WalkDriver.h` take the current bodies of
  `Petri.cpp`, which becomes a thin `main`. `cli/README.md` documents the
  options; `README.md` keeps the user table.

### 6.1 Loaders and printers

* `io/PNETIO.h`: PNET read and write, built on `SparseMatrixIO` stream
  variants (today it takes a file name only).
* `parse/sexpr/Sexpr.h` (reader, about 100 lines, the `datum` style of
  libHSC) and `parse/sexpr/PropertyReader.h` (datum to `Property`, about 200
  lines). `--props` picks the syntax by extension (`.xml` MCC, `.sexpr`
  s-expressions), `--propsSyntax=` overrides.
* `expr/Printer.h`: s-expression printer over indices or names.

### 6.2 Multi-target walking

Today a `Target` is one goal and the portfolio races threads on it. The
ITS-Tools use is dozens to hundreds of predicates on one net, so:

* `TargetSet`: the open targets, each with an atomic `solved` flag and a
  witness slot (the "verdicts" section of the Knowledge Base of
  `WALK_PLAN.md` 3.9, finally). A step checks every open target whose atoms
  touch a place the fired transition changed: a sparse `place -> targets`
  index built once, so the check stays proportional to the arcs touched and
  to the targets that depend on them, never to the number of targets.
* A walker aims its heuristic at one target (its focus) but every open
  target reached on the way is claimed; the focus rotates at restarts over
  the still-open targets, as `WALK_PLAN.md` 4.7 says.
* The driver schedule becomes what ITS-Tools does by hand today: a first
  round of cheap random multi-target walks on all threads, then heuristic
  rounds per open target with the growing per-property budget of
  `--totalTime`. `FORMULA` lines stream as targets fall.
* Deadlock is one target among the others; `bound` targets never solve
  unless a known bound is given, they record a maximum instead.

### 6.3 Parikh strategy

The hint-driven walk of `RandomExplorer.runGuidedReachabilityDetection`, as a
`Strategy`: choose among enabled transitions with a positive remaining count
in the hint, decrement on firing, prefer the smallest rank when a partial
order is given, cycle through rank-first / random / last / first modes across
restarts, and relax the restriction with probability growing with the number
of restarts (the Java decay), so a slightly wrong Parikh vector still guides.
Later the same hint can weight the relaxed-plan or best-first choice rather
than restrict it; first reproduce what works. Step budget defaults to
100 x |parikh| when the property has a hint.

### 6.4 Bound targets

`Target` kind `Bound`: a linear expression, the maximum seen so far, an
optional known upper bound that ends the target. Best-first for bounds is
greedy on the increase of the expression, with the usual epsilon and stall.
Last phase; the `UpperBoundsSolver` also relies on `CoverWalker`, which stays
in Java.

### 6.5 Documentation

`KERS.md` gets a PNET section (or PNET gets its own short file), the
s-expression grammar and the result protocol live in this document until
frozen, then move next to `KERS.md`; `parse/sexpr/README.md`,
`cli/README.md`, `walk/algorithm.md` (multi-target step, Parikh strategy)
and `README.md` follow the code.

---

## 7. What changes in ITS-Tools

All in `petrispot/fr.lip6.move.petrispot.runner`, beside `PetriSpotRunner`:

* `PNETFormatIO` (section 3) and `SexprPropertyPrinter` (section 4.3).
* `PetriSpotWalker`: the replacement API, deliberately shaped like
  `RandomExplorer` so call sites change one line:
  `int[] runReachability(ISparsePetriNet net, List<Expression> predicates,
  List<SparseIntArray> parikh, List<SparseIntArray> order, long steps, int
  seconds)`, `boolean runDeadlock(...)`, later `int[] runBounds(...)`. It
  writes the two files, runs the binary with `--threads`, `--walkSteps`,
  `--totalTime`, `-q`, reads the stream, fills the verdict array, and falls
  back to `RandomExplorer` if the binary is missing or fails.
* Call sites: `ReachabilitySolver.randomCheckReachability` (one request
  replaces the random plus best-first phases; the probabilistic and
  exhaustive Bloom walks stay on `RandomExplorer`), `tryReplayParikh` (one
  request with all hints instead of one walk per Parikh vector),
  `DeadlockSolver`, `AtomicReducer(SR)`, and last `UpperBoundsSolver`.
  Technique words in `DoneProperties` come from the `FORMULA` line.

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
rough.

### Phase 0: CLI11 and `cli/` (about 400 lines moved, 150 new)

Vendor CLI11, `cli/Options.h`, split the drivers out of `Petri.cpp`, same
flags. Check: `Petri/test/` invariant and walk scripts give byte-identical
`FORMULA` and `Computed ...` lines; the MCC driver and `PetriSpotRunner` run
unchanged.

### Phase 1: s-expression properties (about 350 lines)

Reader, property reader, printer, `--props` by extension, `--printProps=`
syntax. Check: MCC XML to sexpr to AST equals MCC XML to AST on the property
files of the development models; hand-written `.sexpr` files for
`Petri/examples/` nets replace or complement `Petri/test/props/`.

### Phase 2: PNET (about 250 lines, plus 40 in Java)

`io/PNETIO.h`, `--net`, `--exportNet`, `kersconv --decode-net`; `KERS.md`
section. Check: PNML to PNET to walk reproduces the PNML walk verdicts and
step counts with a fixed seed on the development models; a PNET written by
the Java writer (a small `main` in the runner plugin) loads and walks.

### Phase 3: result protocol (about 100 lines)

`STATE`, `WITNESS` over indices, `--names`, `--printUnknown`, flush
discipline, exit codes. Check: the MCC driver still passes the harness; a
script in `Petri/test/` parses the stream.

### Phase 4: multi-target walking (about 500 lines)

`TargetSet`, the place-to-targets index, focus rotation, the two-stage
schedule. Check: on the development models the union of verdicts equals
today's per-property runs at equal budget, and a request with the 200-odd
fireability atoms of the challenge model as separate `reach` forms answers in
one walk what today takes 200 walks. ThreadSanitizer clean.

### Phase 5: Parikh strategy (about 300 lines)

`(parikh ...)` hints, `ParikhStrategy`, budget defaults. Check: Parikh vectors
dumped from ITS-Tools on a few models (SMT solutions, saved as `.sexpr` in
`Petri/test/props/`) are realised at least as often as by the Java walk.

### Phase 6: Java runner and call sites (ITS-Tools repository)

`PNETFormatIO`, `SexprPropertyPrinter`, `PetriSpotWalker`, then the call
sites one at a time, `ReachabilitySolver` first. Check: the `itstools` MCC
driver on the reachability examinations solves at least what it solves today,
with the walk time per property logged on both sides.

### Phase 7: bounds (about 250 lines)

`bound` targets, greedy strategy, `BOUND` lines, `UpperBoundsSolver` call
site.

---

## 10. Questions for the user

1. PNET as a separate magic and file, or as a KERS version 2 with a "net"
   flag? Separate is proposed: a KERS reader stays a matrix reader.
2. Do we want the name table at all, or PNML for anything human? Proposed:
   optional flag, never sent by ITS-Tools, cheap to have for debugging.
3. `bound` and `parikh` in the same file as the properties (proposed) or a
   separate `--hints=` file?
4. Should the Java side replay every witness on its own net as a check, or
   only under `DEBUG`? Proposed: `DEBUG` only, once the boundary is trusted.
5. Order of phases 4 and 5: multi-target first serves `AtomicReducer` and the
   random phase; Parikh first serves the SMT loop, which is where ITS-Tools
   gains most today. Proposed as written, both before the Java work.
