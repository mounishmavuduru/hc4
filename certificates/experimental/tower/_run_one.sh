#!/bin/bash
# _run_one.sh <sing-file> <timeout-seconds>
# Runs one graded-tower Singular script detached, with a long budget, writing
# a result file next to it.  Detached (setsid) so the host's ~10-min job kill
# does not reap it; run at most TWO of these at once (memory).
d="/mnt/c/Users/mouni/hc4/certificates/experimental/tower"
cd "$d" || exit 3
sing="$1"; to="${2:-1500}"
base="${sing%.sing}"
res="${base}_result.txt"
setsid bash -c "cd '$d'; echo START \$(date +%H:%M:%S) budget=${to}s > '$res'; timeout $to Singular -q < '$sing' >> '$res' 2>&1; echo DONE rc=\$? \$(date +%H:%M:%S) >> '$res'" >/dev/null 2>&1 &
echo "launched $sing -> $res (pid $!)"
