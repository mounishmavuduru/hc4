import re, os, sys
os.chdir(os.path.dirname(os.path.abspath(__file__)))
f = sys.argv[1] if len(sys.argv) > 1 else "d5_web_branch1.singular.txt"
p = sys.argv[2] if len(sys.argv) > 2 else "32003"
s = open(f).read().replace(chr(13), "")
head = s[:re.search(r"ideal J =.*?;", s, re.S).end()]
head = head.replace("ring R = 0,", "ring R = %s," % p)
# drop the LIB line's radical machinery need; slimgb only.
out = head + (
    '\nshort=0;\n'
    'option(redSB);\n'
    'int t0 = timer;\n'
    'ideal G = slimgb(J);\n'
    'int t1 = timer;\n'
    '"slimgb time (1/100 s):", t1 - t0;\n'
    '"GB size:", size(G);\n'
    'if (size(G) == 1 && G[1] == 1) { "J = (1): residual INCONSISTENT -> branch EMPTY outright"; }\n'
    'else { "dim V(J):", dim(std(G)); }\n'
    'quit;\n')
open("/tmp/gbprobe.sing", "w").write(out)
print("wrote /tmp/gbprobe.sing for", f, "mod", p)
