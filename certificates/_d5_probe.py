import re, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
for f in ["d5_web_branch1.singular.txt", "d5_web_branch0.singular.txt"]:
    s = open(f).read().replace(chr(13), "")
    rv = re.search(r"ring R = 0,\((.*?)\),dp;", s, re.S).group(1)
    J = re.search(r"ideal J =(.*?);", s, re.S).group(1)
    T = re.search(r"ideal T =(.*?);", s, re.S).group(1)
    jg = [g for g in J.split(",\n") if g.strip()]
    tg = [g for g in T.split(",\n") if g.strip()]

    def maxdeg(block):
        md = 0
        for mono in re.split(r"[+\-]", block):
            if not mono.strip():
                continue
            d = 0
            # sum exponents; bare var = degree 1
            for m in re.finditer(r"[a-z][a-z0-9]*(\^([0-9]+))?", mono):
                d += int(m.group(2)) if m.group(2) else 1
            md = max(md, d)
        return md

    print(f)
    print("   vars:", len(rv.split(",")))
    print("   J: gens", len(jg), " maxdeg", maxdeg(J),
          " gen-char-sizes", sorted(len(g) for g in jg)[:3], "...", sorted(len(g) for g in jg)[-3:])
    print("   T: gens", len(tg), " maxdeg", maxdeg(T),
          " gen-char-sizes", sorted(len(g) for g in tg)[:3], "...", sorted(len(g) for g in tg)[-3:])
