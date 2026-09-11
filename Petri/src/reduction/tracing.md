# Optional reduction tracing and visual study

Design only. Preserve ITS-Tools' high-debug capability as an explicitly
selectable service: local net exports, rule labels, and highlighted objects.
A later renderer produces a multipage PDF that can be stepped through as an
animation of reduction. This is valuable for studying rules and explaining
their applications, independently of solver witnesses or Java's LTL image.

## Capture has no work when disabled

Tracing defaults to off. Select the traced or untraced coordinator once at
entry, using a small compile-time trace policy (`NoTrace` / enabled capture).
The same rule code runs in both; `if constexpr` removes observation hooks in
the untraced instantiation. Keep this policy limited to instrumentation rather
than generating a specialization per goal or trace option.

No trace-only allocations, event histories, net copies, label formatting,
highlight-set construction, neighborhood traversal, transpose construction,
per-match clocks or renderer calls occur on the disabled path. Passing eagerly
built arguments to a no-op logger does not meet this requirement. Existing
mutation bookkeeping, names and minimal coordinator progress remain necessary
regardless of tracing; extra diagnostic data does not become required state.

Trace levels are independent of rule eligibility and reduction goals:

* Off: no event capture.
* Applications: streamed compact rule-application descriptions and names.
* Local views: selected before/after neighborhoods with semantic highlights.
* Diagnostic: selected candidate/guard explanations and optional expensive
  consistency checks. Rejected matches are not logged by default.

A sink may decline an event using only cheap data already available (rule ID,
step number, selected names/IDs). Consult filters and remaining limits before
creating its payload. Diagnostic checks are separately requested: asking for
a PDF does not enable full-net invariant checks after every application.

## What a rule exposes

Each rule supplies its stable name/variant, focus objects and highlight roles
from its match. Common edit helpers know which objects/arcs are removed,
created, redirected or retained. Together these produce the view without
duplicating the rule algorithm in a visualization module.

For a selected application, capture the bounded local before-view after its
guards and arithmetic checks, before mutation. Commit the edit, capture the
after-view, then emit the pair with its goal, phase/pass and application number.
Include relevant marking and arc weights, protected places, rule-specific
conditions and composed names. Distinguish removed, added, modified and context
objects with labels/styles as well as color. Failed edits do not appear as
successful applications; discard an uncommitted before-view on failure.

Batch rules may emit one labeled batch event to avoid forcing a slower rewrite
strategy. State the granularity in the event. Compaction is an index operation,
not a mathematical reduction step; keep visual identities stable across it.
Names remain human traceability; capture-only stable keys disambiguate repeated
names and index reuse, without adding tracking to the untraced workspace.

Do not keep references into mutable columns after a hook returns. Consume views
synchronously or own only their bounded captured contents. No full-net snapshot
per application, and no retention of retired slots for rendering. The stream is
a visual record of selected steps, not necessarily sufficient to reconstruct
the entire evolving net or replay a solver witness.

## Bound memory, output and capture effort

Trace options include rule/name filters, step range, sampling interval, maximum
events, nodes/arcs per view, label bytes, neighborhood depth, inspected adjacency
entries, total emitted bytes, and optional capture-time limit. Establish finite
defaults when the feature is implemented. Limits apply before allocation and
during expansion, not only when writing the finished view.

Neighborhood selection starts from the rule's affected objects and adds nearby
context. Use existing adjacency where available. If a missing transpose would
require global work, either request and account for that trace-only work within
its limit, or emit a narrower view with unavailable context marked. Never scan
a whole high-degree neighborhood merely to discover it exceeds a display cap.
If the focus itself exceeds the cap, show a labeled partial/batch view and
counts already known cheaply. Mark omitted context explicitly; a clipped arc
must not make a place appear structurally isolated without explanation.

Stream completed events to a caller-owned sink/file and release each payload.
Working memory is bounded by one before/after pair and fixed stream buffers;
there is no unbounded event vector or background queue. After the event/byte
limit is reached, disable further capture and record one bounded truncation
notice. Reserve space for that notice. The reduction continues normally.
An output failure disables tracing with one diagnostic; it does not roll back
a valid reduction. Trace limits never become reduction expansion limits.

Tracing can add wall time, so under a shared deadline it may indirectly shorten
the reduction. Report capture time separately when tracing is enabled. Compare
traced/untraced outcomes using deterministic work limits when assessing
noninterference; do not claim identical deadline-limited runs.

## Local exports and the PDF

Start from the established concepts in `io/FlowPrinter.h`: DOT, highlighted
places/transitions, protected-place styling and limited neighborhoods. That
printer currently constructs full transposes for large-net neighborhood views
and uses a static file counter; the bounded trace adapter needs explicit output
streams/paths, per-run numbering and strict selection limits. Reuse suitable
rendering code without importing those costs into rule matching.

The future PDF renderer consumes the stream or individual local DOT artifacts
after reduction. It labels each page with rule, step, goal and relevant guards;
show before and after side by side, or on consecutive pages. Keep positions of
unchanged objects stable within a pair, using the union of the two captured
views for layout. Mark skipped applications and truncated neighborhoods so
pages do not imply an exhaustive or full-net history. Page navigation supplies
the requested animation without depending on embedded PDF animation support.

Graphviz/PDF tooling is an optional IO-side dependency, not a dependency of the
header-only reducer or libHSC's normal analysis path. Rendering gets its own
page/output/time limits and runs outside the reduction loop. A renderer failure
leaves the captured artifacts available for inspection.

## Acceptance checks when implemented

Verify that disabled tracing never evaluates a trace payload builder and makes
no trace allocations; measure its hot-path overhead against the same rules
without instrumentation. With fixed work limits, tracing must preserve rule
selection and reduced-net semantics. Check clear/append/compaction identities,
before/after markings and highlights on a trivial agglomeration, and bounded
capture on a high-degree focus and a long reduction sequence. Exercise sink
failure, truncation and renderer failure. Inspect the resulting pages for
readable names, weights, rule labels and explicit omitted context.
