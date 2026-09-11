# `cli/` — command line and drivers

The `petri*` binaries are one `main` (`src/Petri.cpp`) over these headers.

* `Options.h` — the `Options` struct (every command line setting, with its
  default) and `addOptions(app, options)`, which declares them on a CLI11
  `App`. Options are grouped as they appear in `--help`: input, exports,
  invariants, reachability walk, CTL. Every option accepts both `--name=value` and
  `--name value`; flags keep the spelling the MCC driver and the ITS-Tools
  runner use.
* `InvariantDriver.h` — the invariant computations: from a KERS matrix
  (`--loadKERS`, program-to-program mode) or from a loaded net (P/T flows and
  semi-flows, printing, `--basisKERS` export).
* `WalkDriver.h` — the reachability side: property loading and printing,
  the round scheduler over open properties (`--totalTime`, `--escalate`), the portfolio run
  for one target and the `FORMULA` output, `--findDeadlock`. CTL properties
  of the file are handed to `CtlDriver.h` first.
  `--reduce` prepares the property/net pair through `reduction/Properties.h`
  before either engine compiles it; deadlock and LP queries use that adapter
  too. Original traces/hints bypass reduction until coordinate lifting exists.
* `CtlDriver.h` — the CTL properties: normalisation, the checker of `ctl/`
  in rounds of growing budgets (`--ctlSteps`, `--ctlRegion`, `--ctlRounds`)
  under a clock per property (`-t`, or a share of `--totalTime`),
  `--threads` properties checked at a time (each check single-threaded),
  `FORMULA` lines, evidence trees with `--trace`.

CLI11 (v2.7.2) is the documented all-in-one header, copied as is into
`src/cli11/CLI11.hpp` with its BSD-3 licence alongside; it is included only
from `Options.h` and `Petri.cpp`. Upgrading is copying a newer release file.

Adding an option: a field with its default in `Options`, one `add_option` or
`add_flag` line in `addOptions` in the right group, then use it in the driver
that needs it; document it in the root `README.md` table if it is
user-facing.
