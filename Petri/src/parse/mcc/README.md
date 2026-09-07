# `parse/mcc/` — MCC property XML

* `PropertyHandler.h` — expat SAX handler. Resolves place and transition ids
  against the net, builds linear atoms from `tokens-count` / `integer-*`
  comparisons, desugars `is-fireable` into pre-arc conditions, builds the
  path operators (`all-paths` / `exists-path` over `globally`, `finally`,
  `next`, `until` with `before` / `reach`) into a `CtlFormula`, and
  classifies each property: EF / AG over a state predicate and EF deadlock
  are the reachability kinds, `place-bound` a Bound, any other temporal
  formula a CTL property; unknown elements leave it Unsupported.
* `PropertyLoader.h` — `loadProperties(file, net)`.
