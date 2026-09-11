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

## Counting rule chains

`counting_chain.py --hsc HSC --petri PETRI [--places 2 --tokens 3 --components 1]`
checks one generated net of independent free cycles and an ordinary constant
place. Analytic binomial counts are compared with HSC before and after native
reduction, then after PNET export/reimport and a second reduction. The same
check covers PetriSpot's standalone export and its re-reduction. This catches
free-SCC fusion followed by constant-place removal losing its counting weight.
Each invocation has a shared 15-second limit and keeps diagnostics in
`Petri/test/logs/`; temporary model/export files are removed after the check.

## `check_statespace.sh`

One model through `petri64 reduce --goal STATESPACE`, its PNET counted by
`hsc-pn --states` (`HSC=` names the binary, default the libHSC build tree),
the four values against the deployed StateSpace oracle, 15 s each step; a
value the consumer leaves unanswered prints as missing, never as wrong.

```
bash Petri/test/reduction/check_statespace.sh bench/models/AutonomousCar-PT-01a
```

## `check_oracle.sh`

One model, one examination, `build/petri64` with the options given, against the
deployed oracle, 15 s: prints the reduction and initial-state lines, then one
line per FORMULA saying whether the oracle agrees; exit 1 on a disagreement.

```
bash Petri/test/reduction/check_oracle.sh bench/models/AirplaneLD-PT-0010 RC --reduce -q --totalTime=8 -t 8
```
