# `parse/mcc/` — MCC property XML

* `PropertyHandler.h` — expat SAX handler. Resolves place and transition ids
  against the net, builds linear atoms from `tokens-count` / `integer-*`
  comparisons, desugars `is-fireable` into pre-arc conditions, classifies each
  property as EF / AG / EF deadlock or Unsupported (nested temporal operators,
  bounds).
* `PropertyLoader.h` — `loadProperties(file, net)`.
