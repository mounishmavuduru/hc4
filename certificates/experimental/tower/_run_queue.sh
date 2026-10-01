#!/bin/bash
# Serial queue runner for the compute-bound graded-tower Singular scripts.
# ONE Singular at a time (never OOMs the box), each capped at ~12.5 GB virtual
# memory (ulimit) so a runaway fails cleanly, plus a per-job time budget.
# Detached via setsid so the host's ~10-min job kill does not reap it.
# Writes <base>_result.txt per job and appends progress to _queue.log.
d="/mnt/c/Users/mouni/hc4/certificates/experimental/tower"
cd "$d" || exit 3

setsid bash -c "
  cd '$d'
  echo QUEUE START \$(date +%H:%M:%S) > _queue.log
  ulimit -v 12500000
  run () {
    sing=\$1; to=\$2; base=\${sing%.sing}; res=\${base}_result.txt
    echo \"--- \$sing budget=\${to}s START \$(date +%H:%M:%S) ---\" | tee -a _queue.log
    echo START \$(date +%H:%M:%S) budget=\${to}s > \$res
    timeout \$to Singular -q < \$sing >> \$res 2>&1
    rc=\$?
    echo DONE rc=\$rc \$(date +%H:%M:%S) >> \$res
    echo \"    \$sing DONE rc=\$rc [\$(grep -m1 -E '^DIM ' \$res || echo 'no DIM')]\" | tee -a _queue.log
  }
  run r2a_1_32003_lam0_x4incr.sing 600
  run r2a_1_32003_x4incr.sing 600
  run r1_np_1_32003_x4incr.sing 600
  run r2c_np1_1_32003.sing 1200
  run r2c_np0_1_32003.sing 1200
  echo QUEUE DONE \$(date +%H:%M:%S) | tee -a _queue.log
" >/dev/null 2>&1 &
echo "queue launched (pid $!) -> _queue.log"
