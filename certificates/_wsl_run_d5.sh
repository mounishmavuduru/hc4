#!/bin/bash
# Run the degree-5 branch decisions in local Singular (WSL).
set -u
cd /mnt/c/Users/mouni/hc4/certificates || exit 3

echo "=== print sanity ==="
printf 'ring r=0,(x),dp; ideal I=x2; size(I); quit;\n' | Singular -q 2>&1

for f in d5_web_branch1.singular.txt d5_web_branch0.singular.txt; do
  echo
  echo "=== $f ==="
  # strip CR robustly, append an explicit quit so batch mode exits cleanly
  python3 -c "s=open('$f').read().replace(chr(13),''); open('/tmp/clean.sing','w').write(s+'\nquit;\n')"
  timeout 540 Singular -q /tmp/clean.sing
  echo "[singular rc=$?]"
done
