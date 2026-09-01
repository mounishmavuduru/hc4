#!/bin/bash
set -u
D=/mnt/c/Users/mouni/hc4/certificates/_lean/tail
for f in "$D"/tail_*.sing; do
  echo "########## $(basename "$f") ##########"
  timeout 400 Singular -q < "$f" 2>&1 | grep -aE 'COMPONENT|RESULT|dim|GBsize'
  echo "[rc=$?]"
done
echo ALLDONE
