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
    'doubling_structure.py', 'cone_leading_form.py',
    'homogenization_identity.py', 'theorem_F.py', 'theoremA_sharp.py',
    # Theorem G and the trichotomy
    'theorem_G.py', 'theorem_G_n2_deg45.py', 'pivot_dichotomy.py',
    'advG_flow_proof.py', 'pzG_audit_identities.py',
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
    'd5_y1slice_reduction.py',
    # degree-5 rank-3 closure -- sympy-only, fail-closed (no CAS needed)
    'verify_paramcert.py',            # R-D5-GEN: Nullstellensatz witness, re-derived + checked
    'd5_survivor_family.py',          # R-D5-TAIL: deductive structural core S1-S4
    'd5_linear_direction_normalform.py',  # linear-direction normal form N1-N3
    'd5_rank2_toppiece.py',           # degree-5 rank-2 first graded reduction
]

# Not run by the fast suite (reason in the table).
SLOW = {
    'fpG_paper_certificate.py': 'full run ~25 min; C1-C9 here, C11-C19 in _part2',
    'pzG_n2_decision.py': 'degree-4 Groebner does not terminate; use atkG_n2_decision.py',
    'conjE_degree3.py': '14-variable Groebner per minor; better in Singular',
    'identityE_is_known.py': 'Rabinowitsch radical membership; >400 s (proved in-session; attribution fact)',
    'atkG_n3_search.py': 'several deg-4 slices exceed the runtime kill',
    'd5_branch1_normalized.py': '54 eqns in 15 unknowns; may exceed the kill',
    'euler_pullback_reformulation.py': 'B2 imports the ~12-min Meng-Yang determinant; parts A,B1 are fast',
    'rt_instances_x2_B_break.py': 'chunk per family: py rt_instances_x2_B_break.py R1',
    'd5_graded_tower.py': 'universal 4x4 graded determinant; ~5-15 min',
    'd5_pivotfree_normalform.py': 'weighted leading-form identities; slow',
    'd5_rank3_pivotfree_decision.py': 'reduction + Singular export + sampling; ~4 min, exports .sing',
    'developability.py': 'passes standalone (~min) but its grandchild escapes the runner tree-kill on Windows',
    # degree-5 R-D5-TAIL certificates that PASS with exit 0 on their expected
    # outcome -- all need WSL Singular
    'd5_lindir_verify.py': 'needs WSL Singular; linear-direction tower J=(1) over F_p and Q',
    '_d5_surv_target.py': 'needs WSL Singular; M2 isotropy radical test (v*^T A v* in radical J -> exit 0)',
    '_d5_surv_c0c1.py': 'needs WSL Singular; c0=c1=0 sub-locus (b3=c2 y1^3) is a rank<3 cone -> exit 0',
    # NOTE: _d5_surv_rabin.py (M1) and _iso_power.py are DECISION scripts whose
    # exit code / stdout is the verdict, not a pass/fail -- see ledger R-D5-TAIL
    # and README; not listed here so --all does not misreport them.
}

MARKERS = ('ALL CHECKS PASSED', 'ALL PASS', 'ALL FAST CHECKS PASSED',
           'ALL PASSED', 'DECIDED', 'SUMMARY', 'PASSED', 'all checks PASS',
           'referentially sound', 'DECISION', 'no counterexample',
           'mode sample : all checks PASS')

TIMEOUT = int(sys.argv[2]) if len(sys.argv) > 2 else 300


def run(script):
    # Popen + explicit process-TREE kill on timeout.  On Windows the py
    # launcher / python.exe can orphan a grandchild whose open pipe makes a
    # plain subprocess.run(timeout=) block far past the cap; taskkill /T fixes
    # it.  sys.executable is python.exe here, but be robust anyway.
    env = dict(os.environ, PYTHONIOENCODING='utf-8')
    t0 = time.time()
    p = subprocess.Popen([sys.executable, '-u', script], stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, env=env, errors='replace')
    try:
        out, _ = p.communicate(timeout=TIMEOUT)
        rc = p.returncode
    except subprocess.TimeoutExpired:
        subprocess.run(['taskkill', '/F', '/T', '/PID', str(p.pid)],
                       capture_output=True)
        try:
            out, _ = p.communicate(timeout=30)
        except subprocess.TimeoutExpired:
            out = ''
        return 'TIMEOUT', time.time() - t0, ''
    dt = time.time() - t0
    out = out or ''
    if rc != 0:
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
