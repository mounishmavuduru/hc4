#!/bin/bash
set -u
D=/mnt/c/Users/mouni/hc4/certificates/_lean/fiber
for f in "$D"/fib_*.sing; do
  echo "### $(basename "$f")"
  timeout 150 Singular -q < "$f" 2>&1 | grep -aE 'RESULT|GBsize'
done
echo ALLDONE
