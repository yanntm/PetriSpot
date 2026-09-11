# Structural reduction validation

Isolated from walk/CTL development tests. The primary integration path is the
normal `petri64 --reduce` query driver, checked against original-net answers
and MCC consensus. Unknown is not an agreement or an error; conflicting known
verdicts are errors and are reported with both source and oracle.

`mcc.py` runs the real CLI on P/T MCC instances, all five supported XML
examinations plus deadlock, original and reduced. Each model has a shared
15-second hard wall limit across its subprocesses. A bounded walk is not a
complete decision procedure; compare only emitted verdicts to consensus and
to each other, and compare emitted bound lower bounds to known maxima.
It reads archives directly from the existing corpus, temporarily extracts
inside that corpus and removes its own extraction after each model. Results
go to a caller-selected path outside the repository. No cluster tools run.

Example (requires filesystem access to the corpus and output folder):

```
python3 Petri/test/reduction/mcc.py --binary build/petri64 \
  --inputs /home/ythierry/git/pnmcc-models-2026/website/INPUTS \
  --oracles /data/ythierry/MCC26deploy/MCC-drivers/oracle \
  --output /data/ythierry/MCC26logs/local/native-reduction/results.jsonl
```

Use `--model GLOB` or `--max-models N` for a bounded subset. The default covers
all P/T archives; colored unfolding is outside the native PNML importer.
Use `--exams RC,RF,UB --lp` to validate the LP path separately.
Use `--standalone` to exercise the real `reduce` use case: export a matched
model/formula pair, then solve it in a separate invocation. The same 15-second
model allowance includes transformation and both analysis invocations.
`summarize.py RESULTS.jsonl` reports oracle matches, unverified answers,
conflicts, errors, timeouts and rule edit counts.

`validate.cpp` retains the small finite-state check from module bring-up as a
supplementary diagnostic. It compares projected reachability, deadlocks and
raw counts on one bounded generated net. It is not the main validation suite
and is not being expanded. If needed, compile with
`c++ -std=c++23 -O2 -I Petri/src Petri/test/reduction/validate.cpp -o build/reduction-validate`
and run one example with `timeout 15s build/reduction-validate SEED`.

## One-shot GPPP CTLC 08 investigation

`gppp_probe.cpp` solves the AG child of MCC26 CTLCardinality index 08 with
its original property seed, prints the counterexample, and replays its enabled
transitions. `gppp_replay.py` independently checks that path against PNML using
Python integers and verifies the endpoint inequalities that refute the formula.
These are temporary diagnostics for this specific oracle disagreement; remove
both once oracle curation and the missing root-trace follow-up are complete.
The conclusions are in `PS_REDUCTIONS.md` and the open actions in `HANDOFF.md`.

Use the MCC26 model/formulas, not a previous edition:

```sh
c++ -std=c++23 -O1 -g -fsanitize=undefined -fno-sanitize-recover=all \
  -I Petri/src Petri/test/reduction/gppp_probe.cpp -lexpat -pthread -o build/gppp-probe
timeout 15s build/gppp-probe \
  /data/ythierry/MCC26deploy/MCC-drivers/INPUTS/GPPP-PT-C0010N1000000000 \
  > Petri/test/logs/gppp-probe.log 2>&1
python3 Petri/test/reduction/gppp_replay.py \
  /data/ythierry/MCC26deploy/MCC-drivers/INPUTS/GPPP-PT-C0010N1000000000 \
  Petri/test/logs/gppp-probe.log
```
