# What the Louvain decomposition is handed, before the change

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
