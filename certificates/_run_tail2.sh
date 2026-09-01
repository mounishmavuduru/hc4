#!/bin/bash
set -u
D=/mnt/c/Users/mouni/hc4/certificates/_lean/tail
for f in "$D"/radA_*.sing "$D"/punc_*.sing; do
  [ -e "$f" ] || continue
  echo "########## $(basename "$f") ##########"
  timeout 500 Singular -q < "$f" 2>&1 | grep -aE 'RADA|PUNC|RESULT|dim|GBsize|SURVIVOR'
  echo "[rc=$?]"
done
echo ALLDONE
