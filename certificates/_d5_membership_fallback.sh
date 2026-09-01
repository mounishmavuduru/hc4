#!/bin/bash
# Fallback decision: instead of computing the full radical(J) (expensive over
# ~16 vars), test each y-coefficient t of det A for membership t in sqrt(J)
# by Rabinowitsch:  t in sqrt(J)  <=>  1 in J + <t*w - 1>  in R[w].
# Branch EMPTY  <=>  every such t is in sqrt(J)  <=>  every Rabinowitsch
# ideal reduces to <1>.  We build this transform in Singular from the same
# generated file by reusing its ring + ideal J + ideal T, then loop over T.
set -u
cd /mnt/c/Users/mouni/hc4/certificates || exit 3

for f in d5_web_branch1.singular.txt d5_web_branch0.singular.txt; do
  echo
  echo "=== membership test on $f ==="
  python3 - "$f" <<'PY'
import sys, re
src = open(sys.argv[1]).read().replace(chr(13),'')
# keep everything up to (and including) the ideal T = ... ; declaration,
# drop the original radical/reduce tail, and append a Rabinowitsch loop.
m = re.search(r'ideal T =.*?;', src, re.S)
head = src[:m.end()]
# extract the ring variable list to build the extended ring R[w]
rv = re.search(r'ring R = 0,\((.*?)\),dp;', head, re.S).group(1)
loop = r'''
int empty = 1;
for (int i = 1; i <= size(T); i++) {
  ring Rw = 0,(''' + rv + r''',rabw),dp;
  ideal J = imap(R, J);
  ideal T = imap(R, T);
  ideal Q = J, T[i]*rabw - 1;
  ideal g = std(Q);
  // 1 in Q  <=>  reduced Groebner basis is <1>
  if (size(g) != 1 || g[1] != 1) { empty = 0; "coefficient", i, "NOT in sqrt(J)"; }
  setring R;
  kill Rw;
}
"empty (1 = branch has NO counterexample):"; empty;
quit;
'''
open('/tmp/memb.sing','w').write(head + loop)
PY
  timeout 900 Singular -q /tmp/memb.sing
  echo "[singular rc=$?]"
done
