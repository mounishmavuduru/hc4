# Cloud compute plan — degree-5 graded-tower closure (HC₄)

**Status: a plan, not executed.** Written 2026-10-01 from session 2b12535a after
Phase 28 confirmed all five open towers are compute-bound on the local 14 GB
machine. This document specifies exactly what to run on a high-RAM cloud box and
how to orchestrate it. Nothing here is run locally.

## Honest scope (read first)

This compute **closes or strengthens finite sub-cases and evidence** for the
degree-5 rank ≤ 2 / rank-3 pivot-free frontier. It does **not** prove HC₄: the
general obstruction is JC₂ (the planar Jacobian conjecture), a conceptual
problem no Gröbner computation settles. A successful run can turn a specific
tower from "evidence" into "no counterexample of that shape exists (mod p, then
char 0)", or — the headline — exhibit a pivot-free constant-Hessian quintic. It
cannot manufacture a theorem about HC₄ in general. All mod-p / single-instance
results remain EXPERIMENTAL until a char-0 lift is proved.

## Why cloud, and what spec

The binding resource is **RAM, not cores.** In Phase 28 a single `std` for
`r2c_np1` (34 024 determinant terms, 236 equations over F₃₂₀₀₃) exceeded a
12.5 GB `ulimit` and aborted with "no more memory"; the other four hit their
time budgets with partial progress. One `std` had reached 10.8 GB before being
killed.

- **Machine:** Linux, **≥ 64 GB RAM** (128 GB preferred for `r2c_np1` and the
  full, non-staged towers), 8+ cores (Singular's `std`/`slimgb` is largely
  single-threaded, so cores help only by running different towers in parallel).
- **Software:** Singular ≥ 4.3 native (`apt-get install singular`), Python 3 +
  sympy for the independent re-verification step.
- **Portability tweak:** `graded_tower.py` assumes a Windows+WSL host and prefixes
  Singular calls with `base.WSL`. On a native Linux box set that prefix to empty
  (invoke `Singular -q < file` directly, feed via **stdin**, never as a file
  argument). This is the only code change needed.

## The five towers (p = 32003, seed 1)

Each is launched with its driver command; raise `--timeout` to the per-job
budget below. The `.sing` inputs already exist under `tower/` and can be reused,
but regenerating from the driver is cleaner and records sizes.

| # | driver command | size | local wall | cloud budget (RAM / time) |
|---|---|---|---|---|
| 1 | `py -u graded_tower.py --mode r2c_np1 --seed 1 --timeout T` | 34024 terms, 236 eq | memory (>12.5 GB) | 128 GB / 12 h |
| 2 | `py -u graded_tower.py --mode r2c_np0 --seed 1 --timeout T` | 20112 terms, 176 eq | timeout 1200 s | 64 GB / 8 h |
| 3 | `py -u graded_tower.py --mode r2a --seed 1 --x4incr --timeout T` | 28126 terms, 294 eq | timeout (stages 4,3 done) | 64 GB / 12 h |
| 4 | `py -u graded_tower.py --mode r2a --seed 1 --lam0 --x4incr --timeout T` | 16222 terms, 195 eq | timeout (stage 3 done) | 64 GB / 8 h |
| 5 | `py -u graded_tower.py --mode r1_np --seed 1 --x4incr --timeout T` | 30163 terms, 293 eq | timeout (stage 4 done) | 64 GB / 12 h |

Priority order: **#1 (`r2c_np1`) first** — it is the most informative (N₄ = 0,
f₃″ a genuine cube by Hesse, a true no-pivot condition) and the only one whose
wall is memory rather than time, so a big-RAM box most directly unblocks it.
Then #2 (`r2c_np0`, the f₃″ = 0 companion). Then the three `--x4incr` towers,
which bank partial per-stage dims even if the full run times out.

Also try, per tower, the **non-staged full run** (drop `--x4incr`) on the
128 GB box — staging trades memory for time and may not be needed with headroom.

## Interpretation (verify the script's convention before trusting a label)

- Final **`DIM -1`** (empty ideal / whole ring with the Rabinowitsch unit) =
  **INCONSISTENT** = no constant-Hessian quintic of that shape exists for this
  prime. Proved-empty mod p (experimental).
- Final **`DIM ≥ 0`** = **CONSISTENT** = candidate survivor. Then: slice to a
  rational point, extract f, and **independently re-verify in sympy (exact
  arithmetic) that det Hess f is a nonzero constant**, and compute the pivot cone
  {v : D²_v f constant}. A **pivot-free** consistent point is the headline.
- **No `DIM`** (timeout / halt / out-of-memory) = still compute-bound; raise RAM
  or budget, or stage harder.

## Across primes and the char-0 lift

1. Run each tower at **p = 32003**, then **p = 40009** (`HC4_P=40009`, as in
   Phase 27's survivor work). Agreement across two primes is the evidence bar for
   an EMPTY verdict; a CONSISTENT point must reproduce at both primes.
2. For any CONSISTENT (ideally pivot-free) point: attempt a **characteristic-0
   lift** — re-solve the same tower over ℚ (rational coefficients) seeded at the
   mod-p point, or run the symbolic-survivor tower (`lift_tower_sym.py`,
   flagged external-machine in the HANDOFF). Only a char-0 construction with an
   independent exact det-Hess check is a real counterexample; mod-p alone is not.

## Orchestration (how I'd drive it)

- **One tower per worker, five workers**, independent and fail-closed. Each
  worker: regenerate its `.sing`, run Singular with its budget, write a
  `*_result.txt`, exit non-zero on no-DIM (so a missing verdict is never read as
  success).
- No shared mutable state between workers; each owns its mode. Memory isolation
  per worker (one heavy `std` per machine, or cgroup-capped if co-located).
- **Aggregate** into a single verdict table (mode × prime → DIM / class), plus
  the independent sympy re-verification artifact for any CONSISTENT point.
- Fold every verdict into `research_log.md` as the next phase, each line labeled
  EMPTY-mod-p / CONSISTENT / compute-bound, negatives stated plainly, no
  manufactured closure. Commit, push, refresh the backup bundle.

## Definition of done for this cloud run

- Every tower returns a real verdict (DIM or a documented still-compute-bound),
  at two primes, or
- a CONSISTENT pivot-free point is found, independently re-verified, reproduced
  at a second prime, and a char-0 lift is attempted.

Either outcome is a genuine, honestly-labeled contribution to the degree-5
case — and in neither case do we claim HC₄ itself is resolved.
