# Reconstructing transition counts over a reduced state diagram

Status: in-memory prototype implemented; no serialization. Scope: ordinary weighted P/T nets
under the certified STATESPACE reduction schedule, initially free-SCC fusion,
constant-place removal, duplicate transitions and dead-transition removal.

## What is counted

Count pairs (reachable marking, original enabled transition). Different
transitions remain distinct even when they have the same effect. Choosing
different individual tokens does not make several occurrences of one transition.
Self-loops count too.

A reduced marking represents a family of original markings. For a free
component of K original places whose total is M, that family consists of every
nonnegative distribution of M tokens over the K places. The free unit transfers
can reach every such distribution without changing any outside coordinate.
Families of disjoint free components therefore combine independently at a fixed
reduced marking. Their totals can still be correlated in the reduced diagram.

This family property is the certificate needed by both counting methods below.
A more general reduction cannot use this record unless its reconstruction has
the same property or supplies another explicitly supported reconstruction.

## One shared reconstruction record

Retain a baseline view before the first transformation that loses enabling
information:

* A sparse preset for each baseline transition. The PNML-seeded prototype
  retains each transition separately, with implicit multiplicity one. Retain
  internal transitions and no-effect transitions even when the working net
  subsequently deletes them. Identical presets may share a counting query,
  adding multiplicities, regardless of their effects.
* A partition of baseline places into reconstruction groups. Each group has
  its number K of baseline places and a total supplied either by a current
  reduced place or by a fixed constant. An unfused place is a group with K=1.
* Stable group identities and a baseline-place-to-group map, separate from
  the compacted indices of the working net.

This is a sparse enabling view, not a runnable PNET: postsets and an initial
marking are unnecessary to count enabled transitions over an already computed
reachable set. A full source PNET remains a convenient development oracle and
would support reconstructing events for other operations later.

Do not make independent transition copies in each SCC record. A transition
can consume from several components and ordinary places; it must be counted
once with all its requirements combined. Keep one shared transition table and
let group references express its component dependencies.

## Maintaining the record

Working-net edits do not destroy the baseline transition table.

Free-component fusion unions the absorbed groups, adds their sizes, and points
the resulting group at the surviving place. This must flatten repeated fusions:
K counts baseline places, not the number of current representatives. Compaction
changes only the current-place reference. Preserve the original sparse presets;
derive demands against the current partition when counting.

If a representative becomes constant, retain the group and replace its current
place reference with the fixed total. This includes ordinary constant places
and the components currently recorded by PCONST. No correspondence is inferred
from the order of PDROP or PCONST entries.

Duplicate-transition fusion in the working net does not merge away baseline
presets with different reconstruction weights. Initially, keeping all baseline
transitions avoids this complication. A transition proven dead need not be
removed from the counting view: its enabled count must be zero. This also makes
the first implementation easier to compare against the original net.

Unknown transformations invalidate the transition certificate. A PNET arriving
with PCOEF/PCONST but no reconstruction view cannot recover lost transitions;
absence of the new record must not be treated as an identity record. Existing
TMULT may seed multiplicities only while its current certificate is valid.

## Analytic evaluation

For a transition t and group g, let d(t,g) be the sum of its baseline preset
weights over places in that group. Every place has its own lower bound from
the preset. Subtracting those lower bounds gives a bijection to distributions
of M-d(t,g) tokens over K places. Thus its enabling factor is

    W(t,g,M) = binomial(M - d(t,g) + K - 1, K - 1), if M >= d(t,g),
               0 otherwise.

This formula includes untouched groups (d=0), ordinary places (K=1, a Boolean
guard), and internal unit transfers (d=1 in their component). For each reduced
reachable marking q, multiply the factors of all groups, including constant
groups. Sum over q, then over baseline transitions with their multiplicities.
Use exact integers throughout.

The DD supplies the correlations between totals and outside guards: multiply
head and child counts on each existing arc and sum over arcs. Do not replace
the DD by independent marginal counts. Leaf evaluation sums W over that leaf's
values. Constant-group factors are computed once per distinct counting query.

Memoize by canonical DD node, shape/coordinate context, and the interned sparse
group-demand query. A per-query cache can make the query identity implicit.
Share identical group-demand queries across transitions, adding multiplicities.
Poll the common deadline; publish the metric only after the complete exact sum.

## Relationship to diagram expansion

An explicit reconstruction replaces each component total M with the diagram
of its K-place distributions, restores constant groups, and applies the
baseline transition guards before ordinary counting. The analytic evaluator
computes exactly the cardinality of those filtered expansions without building
them. Neither approach needs to explore reachability again.

The first comparison should use AutonomousCar-PT-01b: the reduced diagram plus
the in-memory reconstruction record must reproduce all four unreduced values,
especially 521442 transitions. Compare per-transition enabled counts as well
as the total to locate errors. Then check another real model with successive
fusion and a constant component. Serialization and general transformation tracing are outside this feature.

## Prototype boundaries and interfaces

Only a directly loaded PNML seeds the reconstruction record, before native
StateSpace reduction. PNET inputs retain their existing counting behavior; no
lost ancestry is invented from old blocks. The original presets and names are
immutable shared data. A disjoint-set forest tracks original-place groups; a
separate current-slot map follows compaction. Fixed group totals survive even
when no DD coordinate survives. Repeated reductions carry this typed optional
record alongside Counting, never through its PNET block codec.

The consumer compiles original presets into sparse group demands and caches
identical resulting queries. Fixed groups contribute exact scalar factors. Live
groups are evaluated through `(count-enabled R (PLACE DEMAND) ...)`, using
existing leaf-weight declarations for their component sizes. That operation
uses a fresh memo table keyed by shape, coordinate offset and diagram handle;
its demands and weights are fixed for its lifetime. It includes ordinary
places, so no additional selector construction is needed.

Count every group exactly once: the consumer must not multiply the existing
PCONST state-count factor again. An empty reduced net has one empty residual
marking and still requires counting its fixed groups and original transitions.
An incomplete reachable DD cannot supply a StateSpace transition metric.
Interruption publishes no partial metric; deadline checks cover query assembly,
cache hits, DD recursion and combinatorial evaluation.

For a source-level audit, an opt-in diagnostic reports each original transition
name and its exact enabled count, for reduced and unreduced execution alike.
Compare those vectors, not only their sums, on real inputs. Normal verbose
logging reports the aggregate query count without printing every transition.
