# KERS — Kernel/Elimination Result Sparse format

KERS is a compact binary format for sparse integer matrices.
It is used by PetriSpot both for **input** (incidence matrices) and **output** (invariant bases),
enabling efficient program-to-program transfer without going through PNML or ASCII representations.

## Design goals

- Single format for input and output: the caller reads the result the same way it wrote the input.
- Suitable for large models: 10⁵ × 10⁵ matrices with millions of non-zeros load in milliseconds.
- No external dependencies: plain binary, no compression, no serialisation library.
- Self-describing dimensions: a reader can allocate correctly before parsing any entries.

## File structure

All multi-byte integers are **little-endian**.

### Header (16 bytes)

| Offset | Size | Type   | Description                        |
|--------|------|--------|------------------------------------|
| 0      | 4    | char[4]| Magic: `K` `E` `R` `S` (0x4B 0x45 0x52 0x53) |
| 4      | 1    | uint8  | Version: `1`                       |
| 5      | 1    | uint8  | Flags: reserved, must be `0`       |
| 6      | 4    | uint32 | `nrows` — number of rows           |
| 10     | 4    | uint32 | `ncols` — number of columns        |
| 14     | 2    | —      | Padding (zero)                     |

### Body — column entries

Only non-empty columns are written, in **ascending column order**.

For each non-empty column:

| Size | Type   | Description                              |
|------|--------|------------------------------------------|
| 4    | uint32 | Column index (0-based)                   |
| 4    | uint32 | `nnz` — number of non-zero entries in this column |
| 4×nnz | uint32 | Row indices (0-based, contiguous, sorted ascending) |
| 8×nnz | int64  | Values (signed, contiguous, matching row order) |

Row indices are written as a contiguous block first, then all values as a contiguous block.
This separates structure from data and allows a reader to load the index array and value array
directly without interleaving.

Row indices within a column are **sorted in ascending order**.
This allows readers to use `append()` (amortised O(1)) rather than `put()` (O(log n)).

### Terminator

After the last column, a 4-byte sentinel marks end of data:

| Size | Type   | Value        |
|------|--------|--------------|
| 4    | uint32 | `0xFFFFFFFF` |

## Semantic conventions

### Input matrix (incidence matrix)

When used as input to PetriSpot via `--loadKERS`:

- `nrows` = number of **variables** (places for P-flows, transitions for T-flows)
- `ncols` = number of **constraints** (transitions for P-flows, places for T-flows)
- The matrix is the **incidence matrix** C = flowTP − flowPT
- For T-flows/T-semiflows, pass the same C; PetriSpot transposes internally

### Output matrix (invariant basis)

When written by PetriSpot via `--basisKERS`:

- `nrows` = number of variables (same as input)
- `ncols` = number of basis vectors (invariants)
- Each column is one basis vector: a sparse integer vector in the null space of C

No constant terms are stored (constants depend on the initial marking and are not part of the basis).

## PNET: a net as three KERS blocks

PNET is the tool-to-tool net format of the reachability side (`INTEROP.md`
section 3): a 16-byte header followed by three KERS blocks, `flowPT` (places x
transitions, pre-arcs), `flowTP` (post-arcs) and the initial marking (places x
1), each exactly as above. Places and transitions are identified by index; the
loaded net names them `p<i>` and `t<i>`.

| Offset | Size | Type | Description |
|--------|------|------|-------------|
| 0 | 4 | char[4] | Magic: `P` `N` `E` `T` |
| 4 | 1 | uint8 | Version: `1` |
| 5 | 1 | uint8 | Flags: reserved, must be `0` |
| 6 | 4 | uint32 | Number of places |
| 10 | 4 | uint32 | Number of transitions |
| 14 | 2 | | Padding (zero) |
| 16 | | KERS | `flowPT`, then `flowTP`, then the marking |

After those three blocks a net may carry zero or more optional **named
blocks**. Each is a 12-byte framing followed by an ordinary KERS payload:

| Offset | Size | Type | Description |
|--------|------|------|-------------|
| 0 | 8 | char[8] | Name, ASCII, zero-padded on the right (e.g. `TMULT\0\0\0`) |
| 8 | 4 | uint32 | Payload length in bytes |
| 12 | | KERS | The payload, a matrix as specified above |

Rules of the container:

* Blocks follow the three mandatory ones and nothing follows them, so a
  reader that wants none stops after the marking.
* Each name appears at most once. Order is the enum order of the producer and
  carries no meaning, but staying canonical keeps written bytes reproducible.
* An unknown name is skipped by its length, without parsing its payload. New
  information is therefore a new name, and the framing never changes: PNET
  has no version or flag for this.
* **A producer emits only what it can vouch for.** The absence of a block
  means "not available", never a default. This is the whole discipline: a
  transformation that cannot maintain a block removes it, and a consumer that
  needs it then leaves its value unanswered rather than reporting a wrong
  one.

### Why any of this exists

A net rarely is what a consumer wants to count. It may be the unfolding of a
coloured net, or the result of structural reductions, and then its objects
stand for several objects of the net the question was really about, or hide
objects that were removed. These blocks carry exactly that: what each object
of *this* net stands for in the net it came from, and what was dropped along
the way. See `HSC_PLAN.md` sections 10 to 13 in this repository for the
reasoning and the maintenance rules.

### `TMULT` — transition multiplicities

**Shape.** One column, as many rows as this net has transitions. The value at
row `t` is `m(t) - 1`, so a row with no entry means `m(t) = 1`. Values are
non-negative.

**Semantics.** `m(t)` is how many transitions of the producer's baseline net
this transition `t` stands for. Fusing two transitions with identical pre and
post vectors is the case that produces `m(t) > 1`.

**Use.** A consumer counting arcs of the reachability graph computes
`Σ_t m(t) · |{ s reachable : s enables t }|`. That is the MCC `TRANSITIONS`
value, which counts arcs labelled by transitions.

**The presence rule, which is the subtle part.** Inside a present block, a
missing row means a multiplicity of one. A *missing block* means something
else entirely: that the net carries no evidence its arcs are the arcs of the
net it came from, and a consumer must then not report an arc count at all.
The two must not be conflated: the block is both the weights and the
producer's statement that it accounted for everything it removed.

**Maintenance.** Fusing duplicates adds the dropped weights to the survivor,
following the chain when three or more transitions are identical. Removing a
transition that could fire, without a survivor standing for it, cannot be
accounted for and removes the block. Removing a transition that can never
fire contributes no arc and only re-indexes the column.

### `PDROP` — what removed constant places held

**Shape.** One column. It is a list, not an indexed vector: the rows are
`0..k-1` for the `k` places recorded so far, and the value of a row is the
marking that place held. A place holding no token contributes nothing and is
not recorded. Later removals append.

**Semantics.** Each value is a number of tokens present in *every* reachable
marking of the net that remains, held by a place that no longer exists.

**Use.** `MAX_TOKEN_PER_MARKING` is the maximum over reachable markings of
the sum over the surviving places, plus the sum of this block. Each value is
also a candidate for `MAX_TOKEN_IN_PLACE`, which is the maximum of the
engine's answer and the largest value here. The block does not affect
`STATES`, since a constant place has one value, nor arc counts.

**Maintenance.** Appending composes, so chained removals need nothing else. A
removal that leaves a place whose marking is not constant is not this rule
and must not be recorded here.

### `PCOEF` — place multiplicities (declared, not yet produced)

**Shape.** One column, as many rows as this net has places; value `K - 1` at
row `p`, so a missing row means `K = 1`.

**Semantics.** Place `p` stands for `K` places of the baseline net over which
tokens travel freely, a free strongly connected component collapsed to one
place holding the component's total.

**Use.** A marking of `m` there represents `C(m+K-1, K-1)` markings of the
baseline net, so a consumer counting states folds that weight into its count
rather than multiplying at the end: the correction depends on the marking.
Token totals are unchanged, since the component's total is what the place
holds.

### `GHOSTPT` and `GMULT` — transitions removed but still counted (declared, not yet produced)

**Shape.** `GHOSTPT` is a matrix of as many rows as this net has places and
one column per ghost: the pre-arcs of a transition that was removed. `GMULT`
is one column of as many rows as `GHOSTPT` has columns, holding `m - 1` for
each.

**Semantics.** A transition may be removable from the transition relation
while its arcs still belong to the graph being counted, because its firing
changed nothing or was achievable otherwise. What a consumer needs to count
it is not the transition but its enabling condition, which is what these
blocks carry.

**Use.** The arc count adds `Σ_g m(g) · |{ s reachable : s enables g }|`. A
ghost never takes part in the fixpoint.

The block dimensions must agree with the header. Size is about
`12 x arcs + 8 x transitions` bytes: the 59k-transition, 1.26M-arc
ErlangenMainframe net is 16 MB against 66 MB of PNML and loads in 40 ms.

```
petri64 -i model.pnml --exportNet=model.pnet          # PNML -> PNET
petri64 --net=model.pnet --props=req.sexpr --threads=4 # walk from PNET
petri64 --net=model.pnet --Pflows                      # invariants from PNET
kersconv --decode-net model.pnet                       # inspect as text
```

## Usage with PetriSpot

```
# Export incidence matrix of a Petri net in KERS format
petri64 -i model.pnml --exportAsKERS=model.kers

# Compute P-semiflows from a KERS matrix, export result as KERS
petri64 --loadKERS=model.kers --Psemiflows --basisKERS=basis.kers

# Equivalent: compute directly from PNML, export basis
petri64 -i model.pnml --Psemiflows --basisKERS=basis.kers
```

## kersconv — ASCII conversion utility

`kersconv` is a small companion tool for inspecting and constructing KERS files without PetriSpot.
It is built alongside the `petri*` binaries.

### Usage

```
kersconv --decode <input.kers> [output.txt]      # KERS binary → ASCII (stdout if no output file)
kersconv --encode <input.txt>  <output.kers>     # ASCII → KERS binary
kersconv --decode-net <input.pnet> [output.txt]  # PNET → ASCII (counts, flowPT, flowTP, marking)
```

### ASCII format

The text representation is the same as `--exportAsMatrix`:

```
<nrows> <ncols>
<rowIdx> <colIdx>:<val> <colIdx>:<val> ...
...
```

- First line: matrix dimensions.
- Each subsequent line encodes one **non-empty row**: the row index followed by space-separated `colIdx:value` pairs.
- Empty rows are omitted.
- Values are signed integers.

### Example

```sh
# Inspect a basis file
kersconv --decode basis.kers

# Roundtrip: decode then re-encode
kersconv --decode basis.kers basis.txt
kersconv --encode basis.txt basis2.kers
```

### Consistency with PNML-based output

When comparing `--loadKERS` output to direct PNML computation, the invariant **coefficients and
count are identical**. The following differences are expected and by design:

| Property | PNML path | KERS path |
|---|---|---|
| Variable names | Real place/transition names from PNML | `x0`, `x1`, … (no name table in binary format) |
| RHS constant | Dot product with initial marking (e.g., `= 1`) | Always `0` (no marking stored in KERS) |
| Log lines | Includes parsing and model info | Only computation lines |
| "Computed" line | `Computed N P flows in T ms.` | `Computed N P flows in T ms.` (same) |

For pure **flows** (as opposed to semiflows), `= 0` on the RHS is mathematically correct regardless.
For **semiflows**, the PNML path annotates each invariant with the conserved quantity value
(the invariant dotted with the initial marking), which is unavailable in the KERS path.

## Size estimate

For a matrix with `N` non-zero entries:
- Header: 16 bytes
- Per non-empty column header: 8 bytes
- Per entry: 12 bytes (4 row + 8 value)
- Terminator: 4 bytes

Total ≈ 16 + 8×ncols_nonempty + 12×N bytes.

For N = 10⁶ entries: approximately **12 MB**.

## Overflow behaviour

PetriSpot ships three binaries: `petri32` (int), `petri64` (long), `petri128` (long long).
All three read and write values as `int64` in the file.
When reading into a narrower type (e.g. `petri32`), values outside `[-2³¹, 2³¹-1]` trigger
a warning on stderr and are truncated. Use `petri64` or `petri128` for large coefficients.
