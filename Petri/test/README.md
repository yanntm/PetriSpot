# `test/` — test material

* `props/` — hand-written MCC-format property files for models that live
  outside the repository (see its README).
* `logs/` — run logs, git-ignored. `baseline_*.txt` are the invariant outputs
  of the example nets used as a regression reference.
* `bench_kers.sh` — KERS input/output benchmark of the invariant solver.
* `sexpr_roundtrip.sh` — MCC XML to s-expressions and back on the property
  files of the given model folders; the ASTs must agree.
* `probes/` — one-off analysis scripts kept for reference (place projection
  and the BridgeAndVehicles bound finding).
* `pnet_roundtrip.sh` — PNML versus PNET input on the given model folders:
  same invariant counts, same verdicts and step counts (step-bound walks).

Keep the development loop on the small models (Airplane, Angiogenesis, the
small Bridge); the challenge model goes through the MCC harness.

Models: extract MCC archives (`~/git/pnmcc-models-2026/website/INPUTS/*.tgz`)
into `bench/models/<model>/` at the repository root (git-ignored). Development
set used so far: AirplaneLD-PT-0010, Angiogenesis-PT-05,
ErlangenMainframeV1-PT-bP09C09 (the challenge), BridgeAndVehicles-PT-V*.

MCC harness: `~/git/MCC-drivers` (`BK_TOOL=petrispot ./run_test.pl
oracle/<model>-RF.out -t 300`, or `petrispotxred` to run the ITS-Tools reducer
first); its `petrispot/install.sh` takes `PETRISPOT_BIN=<path>/build/petri64`
to install a local build.
