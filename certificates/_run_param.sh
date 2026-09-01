#!/bin/bash
set -u
D=/mnt/c/Users/mouni/hc4/certificates/_lean/param
for f in "$D"/param_b*_32003.sing; do
  echo "### $(basename "$f")"
  timeout 480 Singular -q < "$f" 2>&1 | grep -aE 'RESULT|GBsize|std time|dim|reduced|branch'
  echo "[done]"
done
echo ALLDONE
