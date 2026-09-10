#!/bin/bash
# probe_node.sh: what a compute node offers a job. The CPU model and the
# instruction sets a binary may assume, the cores and the memory the job's
# cgroup grants (and whether a large allocation survives), the Java found,
# then every deployed tool on AirplaneLD-PT-0010 through the harness itself,
# so a binary the node cannot run shows as it would in a campaign log.
# Submit from the harness root, one job per (class, core count) to compare:
#   cd ~/MCC26/MCC-drivers && mkdir -p probe && cd probe && \
#   oarsub -l "/nodes=1/core=6,walltime=0:20:0" -p "(host like 'small%')" \
#     "cd ~/MCC26/MCC-drivers && bash probe_node.sh"
# The log is OAR.<jobid>.stdout in the folder the job was submitted from.
# Every step is bounded (60 s at most); the whole probe takes a few minutes.
set -u
TREE=$(cd "$(dirname "$0")" && pwd)
cd "$TREE"
sec() { echo; echo "=== $1"; }

sec "job"
echo "host $(hostname)  date $(date -Is)  job ${OAR_JOB_ID:-none}"
echo "cores granted: $( [ -n "${OAR_NODE_FILE:-}" ] && wc -l < "$OAR_NODE_FILE" || echo '?')  nproc: $(nproc)  affinity: $(taskset -pc $$ 2>/dev/null | sed 's/.*: //')"

sec "cpu"
lscpu | grep -E "^(Model name|Architecture|CPU\(s\)|Thread\(s\) per core|Core\(s\) per socket|Socket\(s\)|CPU max MHz|L3 cache)" | sed 's/  */ /g'
flags=$(grep -m1 "^flags" /proc/cpuinfo)
for f in sse4_2 popcnt lzcnt movbe avx f16c fma bmi1 bmi2 avx2 avx512f avx512bw; do
  echo "$flags" | grep -qw "$f" && echo "  $f yes" || echo "  $f NO"
done
lvl=1; echo "$flags" | grep -qw sse4_2 && echo "$flags" | grep -qw popcnt && lvl=2
echo "$flags" | grep -qw avx2 && echo "$flags" | grep -qw bmi2 && echo "$flags" | grep -qw fma && lvl=3
echo "$flags" | grep -qw avx512f && lvl=4
echo "  x86-64 level: v$lvl"

sec "memory"
free -g | head -2
echo "ulimit -v: $(ulimit -v)  ulimit -m: $(ulimit -m)"
echo "cgroups of this shell:"; cat /proc/self/cgroup
for cg in $(cut -d: -f3 /proc/self/cgroup | sort -u); do
  for f in /sys/fs/cgroup$cg/memory.max /sys/fs/cgroup$cg/memory.high /sys/fs/cgroup/memory$cg/memory.limit_in_bytes /sys/fs/cgroup/memory$cg/memory.soft_limit_in_bytes /sys/fs/cgroup$cg/cpuset.cpus /sys/fs/cgroup/cpuset$cg/cpuset.cpus; do
    [ -f "$f" ] && echo "$f: $(cat $f)"
  done
done
if command -v python3 > /dev/null; then
  for gb in 4 12 20 28; do
    echo -n "allocate ${gb} GB: "
    timeout 60 python3 -c "
import sys
n=$gb*(1<<30); b=bytearray(n)
for i in range(0, n, 1<<12): b[i]=1
print('ok, touched')" 2>&1 | tail -1
    echo "  rc=${PIPESTATUS[0]}"
  done
else
  echo "no python3 on the node: allocation test skipped"
fi

sec "java"
command -v java && java -version 2>&1 | head -2 || echo "no java on PATH"

sec "its-tools-native alone"
timeout 60 itstools/itstools/its-tools-native --help 2>&1 | head -4; echo "rc=${PIPESTATUS[0]}"

for pair in "itstools OS" "itstools RC" "petrispot RC" "hsc SS"; do
  set -- $pair
  sec "harness: $1 on AirplaneLD-PT-0010-$2 (60 s)"
  BK_TOOL=$1 timeout 120 ./run_test.pl "oracle/AirplaneLD-PT-0010-$2.out" -t 60 2>&1 | grep -E "^FORMULA|^STATE_SPACE|CPU features|Illegal|Killed|error|Error|exit|BK_TIME|BK_STOP|BK_START|TIME LIMIT" | head -20
done
echo; echo "=== probe done $(date -Is)"
