# `test/` — test material

* `reduction/` — isolated structural-reduction validation: existing MCC models
  and formulas through original/reduced CLI analysis versus MCC oracles.

* `props/` — hand-written MCC-format property files for models that live
  outside the repository (see its README).
* `logs/` — run logs, git-ignored. `baseline_*.txt` are the invariant outputs
  of the example nets used as a regression reference.
* `bench_kers.sh` — KERS input/output benchmark of the invariant solver.
* `sexpr_roundtrip.sh` — MCC XML to s-expressions and back on the property
  files of the given model folders (reachability, bounds and CTL); the ASTs
  must agree.
* `ctl_oracle.sh` — the CTL checker on CTLCardinality and CTLFireability of
  the given model folders, several seeds, every verdict compared with the
  contest oracle (`~/git/MCC-drivers/oracle`); a WRONG line is a soundness
  bug. `props/Airplane-ctl.sexpr` is a hand-written CTL set for the Airplane
  example, run with `--trace` to read the evidence trees.
* `probes/` — one-off analysis scripts kept for reference (place projection
  and the BridgeAndVehicles bound finding; `perprop_payoff.py`,
  `reduction_resistance.py` and `dd_tail.py` measure the campaign claims of
  `PORTFOLIO.md` against a directory of MCC logs; `qla_props.py` writes the
  QuasiLivenessAll target set of a PNML as `(reach pI (fireable tI))` lines,
  the walker's yardstick on a large target set; `liveness_props.py` writes the
  Liveness examination as CTL, one `(ctl liveI (AG (EF (fireable tI))))` per
  transition or, with `--one`, the single conjunction).
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
