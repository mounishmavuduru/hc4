#!/bin/bash
# Fast mod-p decision for the degree-5 rank-3 pivot-free branches.
# Over F_p the Groebner engine is orders of magnitude faster than char 0.
# For each y-coefficient t of det A we test  t in sqrt(J)  by Rabinowitsch:
#   t in sqrt(J)  <=>  1 in <J, t*w - 1>  (over F_p).
# Branch EMPTY (mod p)  <=>  every t is in sqrt(J).
# NOTE: mod-p membership is strong EVIDENCE, not a char-0 proof (a prime could
# be unlucky).  We run several primes; agreement across primes + the earlier
# exact 400-point sampling is the evidence base.  Labelled EXPERIMENTAL.
set -u
cd /mnt/c/Users/mouni/hc4/certificates || exit 3

PRIMES="32003 32009 100003"

for f in d5_web_branch1.singular.txt d5_web_branch0.singular.txt; do
  for p in $PRIMES; do
    echo
    echo "=== $f   (mod $p) ==="
    python3 - "$f" "$p" <<'PY'
import sys, re
src = open(sys.argv[1]).read().replace(chr(13),'')
p   = sys.argv[2]
m = re.search(r'ideal T =.*?;', src, re.S)
head = src[:m.end()]
head = head.replace('ring R = 0,', 'ring R = %s,' % p)
rv = re.search(r'ring R = %s,\((.*?)\),dp;' % p, head, re.S).group(1)
loop = (
'\nint bad = 0;\n'
'for (int i = 1; i <= size(T); i++) {\n'
'  ring Rw = %s,(%s,rabw),dp;\n' % (p, rv) +
'  ideal J = imap(R, J);\n'
'  ideal T = imap(R, T);\n'
'  ideal Q = J, T[i]*rabw - 1;\n'
'  ideal g = std(Q);\n'
'  if (size(g) != 1 || g[1] != 1) { bad = bad + 1; "  coeff", i, "NOT in sqrt(J)"; }\n'
'  setring R;\n'
'  kill Rw;\n'
'}\n'
'"size(T) =", size(T), "   bad =", bad, "   EMPTY(mod p)=", (bad==0);\n'
'quit;\n')
open('/tmp/modp.sing','w').write(head + loop)
PY
    timeout 300 Singular -q /tmp/modp.sing
    echo "[rc=$?]"
  done
done
