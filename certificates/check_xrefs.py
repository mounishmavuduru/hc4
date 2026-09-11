# check_xrefs.py -- P2.4: tracked doc<->script cross-reference integrity.
#
# Greps every project document for cited certificate/script filenames
# (*.py, *.sing, *.sh under certificates/ and src/) and asserts each one exists
# on disk, so a doc edit cannot silently introduce a dangling reference.
# Fail-closed: exit 1 if any cited file is missing.
#
#   py -u check_xrefs.py
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DOCS = ['status.md', 'lemma_ledger.md', 'research_log.md', 'reviewer_packet.md',
        'prior_art.md', 'ROADMAP.md', 'main_result.tex', 'theorem_G_paper.tex',
        os.path.join('certificates', 'README.md')]
# a cited name: word chars/underscores/dashes/dots ending in one of the exts.
# LaTeX escapes underscores as \_ ; strip those before matching.
PAT = re.compile(r'\b([A-Za-z0-9_\-.]+?\.(?:py|sing|sh))\b')
# names that are templates/placeholders, not real files
PLACEHOLDERS = {'d5_branchK_famM.sing', 'hc4-backup.bundle'}
# search roots for existence
SEARCH = [HERE, os.path.join(ROOT, 'src'), ROOT]


def _all_scripts():
    names = set()
    for d in SEARCH:
        try:
            names.update(f for f in os.listdir(d) if f.endswith(('.py', '.sing', '.sh')))
        except OSError:
            pass
    return names


ALL_SCRIPTS = None


def exists(name):
    global ALL_SCRIPTS
    for d in SEARCH:
        if os.path.exists(os.path.join(d, name)):
            return True
    # Tolerate a name that is the TAIL of a line-wrapped longer filename
    # (e.g. "decision.py" from "d5_rank3_pivotfree_\ndecision.py" in the ledger):
    # accept if some real script ends with it and the match is at an underscore
    # boundary of that script.
    if ALL_SCRIPTS is None:
        ALL_SCRIPTS = _all_scripts()
    return any(s.endswith(name) and s != name and s[-len(name) - 1] == '_'
               for s in ALL_SCRIPTS)


def main():
    cited, missing = {}, {}
    for doc in DOCS:
        path = os.path.join(ROOT, doc)
        if not os.path.exists(path):
            continue
        txt = open(path, encoding='utf-8', errors='replace').read().replace('\\_', '_')
        for m in PAT.finditer(txt):
            name = m.group(1)
            if name in PLACEHOLDERS or name.startswith('.'):
                continue
            # skip obvious non-filenames (e.g. 'e.g.py' won't occur; but guard
            # against version-like tokens such as 'v1.2.py' being real anyway)
            cited.setdefault(name, set()).add(doc)
    for name, docs in sorted(cited.items()):
        if not exists(name):
            missing[name] = sorted(docs)
    print(f'{len(cited)} distinct cited script names across {len(DOCS)} docs')
    if missing:
        print('DANGLING REFERENCES:')
        for name, docs in missing.items():
            print(f'  {name}  <- cited in {", ".join(docs)}')
        print('FAILED')
        sys.exit(1)
    print('ALL CHECKS PASSED: every cited script exists on disk.')
    sys.exit(0)


if __name__ == '__main__':
    main()
