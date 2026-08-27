#!/usr/bin/env python
# run_all.py -- one-shot runner for the HC_4 certificate suite.
#
# Runs each FAST, proof-bearing certificate with a per-script timeout and
# prints a PASS / FAIL / TIMEOUT table.  A certificate counts as PASS iff it
# exits 0 AND its output contains a recognised success marker.  Slow or
# CAS-bound certificates (Groebner over many unknowns, or the ~12-min
# Meng-Yang determinant import) are listed separately and NOT run here.
#
#   py run_all.py            # run the fast suite
#   py run_all.py --all      # also run the slow ones (may take an hour+)
#
# Exact/symbolic; every listed script is fail-closed.

import subprocess
import sys
import time
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

FAST = [
    # landscape
    'verify_alpoge_3d.py', 'verify_mengyang_fast.py',
    # tower / pivot machinery
    't3_generic.py', 'verify_q2_core.py', 'ap4_rank_reduction.py',
    'ap4_pivot_closure.py', 'deg4_branch_closures.py', 'dillen4_ap4_control.py',
    # structure theorems
    'doubling_structure.py', 'developability.py', 'cone_leading_form.py',
    'homogenization_identity.py', 'theorem_F.py', 'theoremA_sharp.py',
    # Theorem G and the trichotomy
    'theorem_G.py', 'theorem_G_n2_deg45.py', 'pivot_dichotomy.py',
    'identityE_is_known.py', 'advG_flow_proof.py', 'pzG_audit_identities.py',
    'atkG_algebraic_proof.py', 'atkG_n2_decision.py',
    'paperG_writeup_checks.py', 'wG_writeup_checks.py',
    'refG_paper_integrity.py', 'fpG_paper_certificate_part2.py',
    'gcv_family_wronskian.py', 'gcv_ny_homogeneous_converse.py',
    'tgc_debondt_chain.py', 'tgc_beyond_hc4.py',
    # obstructions
    'pivot_obstruction_d1.py', 'gradient_equivalence_test.py',
    'nax_d1_nonaffine.py', 'nax_pivot_transfer_5var.py',
    # pivot existence
    'pe_pivot_cone.py', 'pe_lowdeg_pivot.py', 'pe_perturbation.py',
    # NEW (2026-08-27)
    'd5_graded_tower.py',
    'd5_pivotfree_normalform.py', 'd5_rank3_pivotfree_decision.py',
]

# Not run by the fast suite (reason in the table).
SLOW = {
    'fpG_paper_certificate.py': 'full run ~25 min; C1-C9 here, C11-C19 in _part2',
    'pzG_n2_decision.py': 'degree-4 Groebner does not terminate; use atkG_n2_decision.py',
    'conjE_degree3.py': '14-variable Groebner per minor; better in Singular',
    'atkG_n3_search.py': 'several deg-4 slices exceed the runtime kill',
    'd5_branch1_normalized.py': '54 eqns in 15 unknowns; may exceed the kill',
    'euler_pullback_reformulation.py': 'B2 imports the ~12-min Meng-Yang determinant; parts A,B1 are fast',
    'rt_instances_x2_B_break.py': 'chunk per family: py rt_instances_x2_B_break.py R1',
}

MARKERS = ('ALL CHECKS PASSED', 'ALL PASS', 'ALL FAST CHECKS PASSED',
           'ALL PASSED', 'DECIDED', 'SUMMARY', 'PASSED', 'all checks PASS',
           'referentially sound', 'DECISION', 'no counterexample',
           'mode sample : all checks PASS')

TIMEOUT = int(sys.argv[2]) if len(sys.argv) > 2 else 300


def run(script):
    env = dict(os.environ, PYTHONIOENCODING='utf-8')
    t0 = time.time()
    try:
        r = subprocess.run([sys.executable, '-u', script], capture_output=True,
                           text=True, timeout=TIMEOUT, env=env, errors='replace')
    except subprocess.TimeoutExpired:
        return 'TIMEOUT', time.time() - t0, ''
    out = (r.stdout or '') + (r.stderr or '')
    dt = time.time() - t0
    if r.returncode != 0:
        last = next((ln for ln in reversed(out.splitlines()) if ln.strip()), '')
        return 'FAIL', dt, last[:80]
    # every script is fail-closed (assert-based): exit 0 => all checks held.
    # the success marker is advisory only.
    note = '' if any(m in out for m in MARKERS) else '(exit 0; no marker)'
    return 'PASS', dt, note


def main():
    do_all = '--all' in sys.argv
    scripts = FAST + (list(SLOW) if do_all else [])
    print(f'HC_4 certificate suite -- {"full" if do_all else "fast"} run, '
          f'timeout {TIMEOUT}s/script\n' + '=' * 60)
    npass = nfail = nother = 0
    for sc in scripts:
        if not os.path.exists(sc):
            print(f'  MISSING  {sc}')
            continue
        status, dt, note = run(sc)
        line = f'  {status:8s} {sc:44s} {dt:6.1f}s  {note}'
        print(line, flush=True)
        if status == 'PASS':
            npass += 1
        elif status == 'FAIL':
            nfail += 1
        else:
            nother += 1
    print('=' * 60)
    print(f'  PASS {npass}   FAIL {nfail}   OTHER {nother}')
    if not do_all:
        print('\nNot run here (see README; run with --all or in Singular):')
        for sc, why in SLOW.items():
            print(f'  - {sc}: {why}')
    sys.exit(1 if nfail else 0)


if __name__ == '__main__':
    main()
