# What the Louvain decomposition is handed

## Before

Measured with `louvain-bench.sh` on a local product carrying `-louvainBench` (the
decomposition alone, no engine, so nothing answers a property first) and `GraphBuilder`
`DEBUG=2` (the files kept). `vars` is what the dependency matrix holds, `cstrs` the
comparisons of the properties whose support holds more than one variable, `maxcstr` the
largest of those supports, `edges` and `bytes` the graph written, `convert` and `louvain`
the seconds those two binaries need on it, `comms` the communities of the first level.

```
model                        vars  cstrs maxcstr      edges      bytes  convert  louvain  comms decomp_ms
Philosophers-PT-000005         25      3      10        210       2360        0        0      7      53
Philosophers-PT-000010         50    107      20      19070     219840        0        0      5      76
Philosophers-PT-000020        100     70      40      44620     573110        0        0      5      98
Philosophers-PT-000050        250     45     100     193000    2743369      .03        0      5     205
Philosophers-PT-000100        500    210     200    3814800   55260201      .56      .01      5    1256
SharedMemory-PT-000005         41     19      21       1951      22712        0        0      5      81
SharedMemory-PT-000010        131     38     100      92251    1195169      .01        0      5     152
SharedMemory-PT-000050       2651     23    2500   48802051  834636576     9.41      .22      5   16234
AirplaneLD-PT-0010             57     11      20       1417      16004        0        0      6      65
```

A comparison of support b contributes b*(b-1) edges, both directions, each weighted
10 * the number of variables; the net itself contributes `|control| x |write|` per
transition, weighted at most 1. So the properties, not the net, decide the size:

* 2651 variables produce 48.8 million edges and 835 MB, because one comparison covers
  2500 of them. `convert` needs 9.41 s of the 10 s it is given, on an idle machine.
* the same shape at every size: 500 variables and a 200 place comparison give 3.8 M edges.
* five communities whatever the model and whatever its size, which is what a clique of
  weight 10n against structural weights of 1/nbElts leaves of the modularity signal.

On a campaign model (VehicularWifi-COL-none, 8729 variables after reduction) the same
measurement gave 17.59 M edges of which 17.52 M carried the constraint weight: the net
contributed 64 418, the comparisons all the rest, and 4177 variables ended up with degree
above 8000 in a graph of 8729.

## After

A comparison over at most `MAX_CONSTRAINT` (20) variables contracts its variables into one
node of the graph, so that no partition can separate them; a wider one, and a group of
overlapping ones grown past half the net, is left to the partition instead. A transition
holding more than `MAX_CONTROL` (8) control places, or inducing more than `MAX_INDUCED`
(64) edges, is left out. Nothing weighted `10 * n` is emitted any more.

```
model                        vars  cstrs maxcstr      edges      bytes  convert  louvain  comms decomp_ms
Philosophers-PT-000005         25      3      10         62        600        0        0      4      46
Philosophers-PT-000010         50    107      20        160       1860        0        0     10      66
Philosophers-PT-000020        100     70      40        225       2230        0        0      2      73
Philosophers-PT-000050        250     45     100        800      10213        0        0     50     129
Philosophers-PT-000100        500    210     200       1600      21747        0        0    100     777
SharedMemory-PT-000005         41     19      21        162       2255        0        0      6      52
SharedMemory-PT-000010        131     38     100        631       9658        0        0      1      64
SharedMemory-PT-000050       2651     23    2500      17551     345036        0        0    100     930
AirplaneLD-PT-0010             57     11      20        193       1841        0        0      4      54
```

| model | edges | bytes | decomposition |
| --- | --- | --- | --- |
| Philosophers-PT-000100 | 3 814 800 -> 1 600 | 55 MB -> 21 kB | 1256 ms -> 777 ms |
| SharedMemory-PT-000050 | 48 802 051 -> 17 551 | 835 MB -> 345 kB | 16 234 ms -> 930 ms |
| AirplaneLD-PT-0010 | 1 417 -> 193 | 16 kB -> 1.8 kB | 65 ms -> 54 ms |

The communities follow the net now instead of the properties: Philosophers-PT-000100 gives
100 of them, one per philosopher, where every model of every size used to give 5. Two of
them give fewer than before (SharedMemory-PT-000010 gives 1, Philosophers-PT-000020 gives
2), which the cost does not tell us how to read: what the partition is worth to the
decision diagrams is measured in verdicts, and that measurement is still to do.
