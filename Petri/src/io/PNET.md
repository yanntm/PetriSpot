# PNET: a binary P/T net, and what its producer knows about it

PNET is the tool-to-tool container for a Petri net: a small header, three
KERS blocks (`KERS.md`), then optional named blocks carrying what the
producer knows about the net's provenance. It reuses KERS as its payload
codec and is otherwise a format of its own; a KERS reader stays a matrix
reader.

Producers and consumers today: PetriSpot (`Petri/src/io/PNETIO.h`, both
directions, `--net` and `--exportNet`, `kersconv --decode-net`), ITS-Tools
(`interop/fr.lip6.move.gal.interop/PNETFormatIO`, write), libHSC
(`tools/hsc-pn --net`, read). The exchange it serves is `INTEROP.md`; the
reasoning behind the optional blocks is `HSC_PLAN.md` sections 10 to 13.

All integers are little-endian. Values inside a KERS payload are int64 as
KERS specifies; a 32-bit build refuses a value that does not fit rather than
truncating it.

## Header and the three mandatory blocks
A 16-byte header followed by three KERS blocks: `flowPT` (places x
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

## Named blocks

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
