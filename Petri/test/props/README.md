# `test/props/` — hand-written property files

MCC-format property files for models that live outside the repository
(`bench/models/<model>/model.pnml`, extracted from the MCC archives).

* `BridgeAndVehicles-PT-*-Bridge.xml` — queries 00 and 01: all vehicles have
  crossed (`SORTI_A >= V && SORTI_B >= V`, and the side-A half); these are
  monotone and easy. Queries 02 and 03: the bridge is full
  (`SUR_PONT_A >= P`, `SUR_PONT_B >= P`, P = capacity, bounded by V): every
  vehicle of a batch must be authorised before any of them leaves, which an
  undirected walk essentially never does on the large instances.

* `Airplane.sexpr` — s-expression syntax sample for `Petri/examples/Airplane.pnml`:
  names, quoted names, indices, `fireable`, arithmetic, all three property
  kinds. Every property but the invariant on speeds has a witness.
