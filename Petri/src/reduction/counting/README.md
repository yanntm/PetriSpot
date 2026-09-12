# Transition counting after free-component fusion

[algorithm.md](algorithm.md) describes the in-memory reconstruction record used by
explicit diagram expansion and analytic counting. It supplements the existing
state-counting record in `../Counting.h`. `Enabling.h` preserves original
presets and follows free-group fusion, fixed totals and compaction. It is not
serialized.

The first witness is AutonomousCar-PT-01b: without reduction, 117338 states
and 521442 enabled transition occurrences. Free-SCC reduction retains the state
count but currently drops the transition-count certificate.

libHSC consumes it through `tools/enabling/count.hh` and an exact memoized DD
query. Export formats and general transformation tracing are outside its scope.
