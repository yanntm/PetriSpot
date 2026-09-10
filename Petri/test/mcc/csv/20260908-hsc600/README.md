# Campaign 20260908-hsc600 — libHSC alone on RD and OS, 600 s

`hsc-pn` as a tool of its own (`BK_TOOL=hsc`, `TAG=hsc600`), no ITS-Tools
beside it: the symbolic engine answering ReachabilityDeadlock and OneSafe
on 1954 instances at 600 s, one core. Logs in
`/data/ythierry/MCC26logs/hsc/20260908-hsc600/`.

|  | answered | ok | wrong | missed | none | bonus | at wall | failures | wall h |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ReachabilityDeadlock | 1036 | 1036 | **0** | 858 | 60 | 0 | 880 | 15 | 159.2 |
| OneSafe | 984 | 984 | **0** | 960 | 10 | 0 | 676 | 0 | 118.5 |

Read plainly: **sound and behind**. 2020 answers, not one of them wrong,
which is the credit; and **0 bonus** on either examination — there is no
instance where libHSC alone has the value. It misses 1818 that the field
holds, and 1343 of those are runs still going when the 600 s ran out, so
the budget is a real part of the gap and not the whole of it: 293 OneSafe
and the rest of the ReachabilityDeadlock misses come from runs that ended
early, which is the engine giving up rather than running out of clock.

Nothing here is a regression against the contest ITS-Tools: the report's
"missed values that only the contest ITS-Tools produced" is empty for both.

These are the numbers a decision to enter libHSC on these two examinations
has to beat, and at 600 s on one core it does not beat them yet. The
comparison that matters next is against ITS-Tools at the same budget on
the same instances, which is a set in the pages, not a number here.
