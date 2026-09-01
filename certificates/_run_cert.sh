#!/bin/bash
set -u
D=/mnt/c/Users/mouni/hc4/certificates/_lean/param
for k in 0 1; do
  echo "########## branch $k ##########"
  timeout 300 Singular -q < "$D/cert_b$k.sing" > "$D/cert_b$k.out" 2>&1
  echo "rc=$?  lines=$(wc -l < "$D/cert_b$k.out")"
  grep -aE 'NGEN|SELFCHECK' "$D/cert_b$k.out"
done
echo ALLDONE
