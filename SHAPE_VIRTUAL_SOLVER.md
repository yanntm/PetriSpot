# Shape sweep: a virtual parallel solver

Proposal for review before implementation, following `CLAUDE.md`.

## Purpose and scope

Show the coverage and response time achievable if all 17 shape heuristics
started together, each with the resources of its measured run. For each
question, the first conclusive answer wins. TRUE and FALSE both count as
answers; UNKNOWN, a partial reachable set and an unanswered question do not.
This is an ideal parallel reference, with no contention or scheduling cost.
It does not predict a portfolio sharing a limited CPU and memory allocation.

First increment: add this reference to the sweep's cactus plots and expose
its totals and per-instance results. Subset controls and heuristic
contribution analysis follow separately.

## Existing page and data

The generator is `~/git/libHSC/experiments/order/sweep_pages.py`; its design
and folder guide are `experiments/order/SWEEP.md` and `README.md` there.
The page currently calls the sum of maximum answer counts per instance
"virtual best". This is a best-single-run coverage measure: it cannot
capture disjoint answers from two partial runs, or their relative speed.
Existing unique/marginal/dominance columns also use counts alone.

The recent input is
`/data/ythierry/MCC26logs/hsc/reduce20260912/results/reduce20260912.tsv`;
logs are in the adjacent `reduce20260912/` folder. The existing rebuild
entry point is `rebuild-order-pages.py` in the campaign directory. It also
rebuilds the original sweep and preserves separate page sets.

The current reduced pages report 1681 instances per examination, with every
heuristic recorded for 1651 CTLC, 1636 CTLF and 1182 SS instances. Collection
is documented as complete, but these roster gaps have not been audited.
They must not silently become failed solver attempts.

## Definition

Keep campaigns separate. A question key is (instance, examination, formula
name), or (instance, SS, metric name). For strategy h, let t(h,q) be its
recorded time to a conclusive answer to question q, infinity if unanswered.
Then:

    virtual_time(q) = min_h t(h,q)
    completion_time(instance, exam) = max_q virtual_time(q)

The second expression ranges over the full expected question set, not just
the questions somebody answered. An examination is complete only when all
its questions are answered. Thus several partial runs can jointly complete
it even when no single strategy does. SS retains all four expected metrics;
unsupported TRANSITIONS remains missing rather than reducing the target.

Union coverage uses question identities, never max(answered) or summed
counts. Repeated timestamps for the same question count once, at the earliest
valid time. Retain tied winners and links to their measured source runs.

## Evidence and timing

The TSV has named `answer_times`, but only aggregate correctness counts.
Read output verdicts and oracle values to associate timestamps with answers.
Flag wrong answers and conflicting verdicts explicitly; exclude disputed
questions from the successful curves rather than silently choosing a value.
Answers with no known oracle value remain separately labelled unverified.
Apply the same evidence policy to physical and virtual curves.

Verify the timestamp clock origin and its coverage of early reduction
answers and StateSpace output before implementation chooses a time field.
Document whether parsing, reduction and shape construction are included.
The existing `wall_s` measures process duration; `reach_s` is a different
measurement and must not substitute for answer latency.

If an emitted answer lacks a timestamp, retain it in coverage and report
the timing gap. A separate curve may use process exit as a conservative
upper bound, clearly labelled and applied equally to every strategy. Do
not invent an exact first-answer time. Resolve duplicate-run references
within a campaign, including chains; flag missing sources or cycles. Reused
times are estimates from the measured representative, not independent runs.

## Page changes

1. Put a prominent virtual-solver summary beside the cactus plots: questions
   answered, examinations complete, timed coverage and evidence exclusions.
2. Overlay a thick, distinctive "Virtual parallel (17 heuristics)" curve
   on answers obtained by time and examinations completed by time. Use the
   same population, timing definition and correctness policy for all curves.
   Keep seconds on x and cumulative count on y, with log/linear time choice;
   preserve zero-time observations and explain any log-axis display floor.
3. Default comparisons to instances with the explicit 17-strategy roster
   accounted for. Show included/excluded counts. An optional available-data
   view may include incomplete rosters, labelled as provisional observed
   coverage. Never interpret an absent row as a measured timeout.
4. Add virtual coverage and completion time to instance detail, with earliest
   answer sources. Keep the virtual solver separate from physical heuristic
   rankings and resource statistics.
5. Rename existing count-based "virtual best" and contribution labels to
   describe their best-single-run meaning until formula-based contribution
   analysis replaces them in a later increment.

## Implementation and validation

After review, implement in libHSC's experiment tooling, updating its README
and sweep design first. Factor answer extraction and aggregation into a small
typed Python module if needed to keep the page builder manageable. No solver
changes or new benchmark runs are needed. Regenerate the local pages through
the existing rebuild entry point; no publishing or cluster activity.

Check disjoint partial answers, repeated events, ties, duplicate references,
wrong/conflicting answers, missing timestamps and missing roster entries.
For a common population the virtual timed-answer curve must dominate every
member, and its coverage must equal the union of accepted question keys.
Compare a few actual instances against their logs and inspect the rendered
plots. Report any clock or data limitation before interpreting speedups.
