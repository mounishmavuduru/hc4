# experimental/ — probes, NOT certificates

Everything here is **experimental** in the project's labeling scheme: exact
arithmetic, but modular (F_32003) and/or instance-level (one witness, one
random survivor point), so nothing here is a theorem. The proof-bearing
certificates live one level up (`certificates/*.py`, run by `run_all.py`).
Each script writes its output next to itself (`*.txt`); the small JSON files in
`survpt/` are the explicit survivor points the tower runs used. Singular
scratch (`*.sing`, `probe8/`) is gitignored. Scripts import `_d5_close` from
the parent directory.

Degree-5 rank ≤ 2 (ROADMAP P5.A2; proved parts are in `../d5_rank2_deg9.py`):

- `d5r2_probe9.py` / `d5r2a_struct.py` — discovery probes for the [det]_9
  structure in the u-constant sub-case (superseded by the formal certificate).
- `d5r1_probe.py` — r = 1 first reduction, instance check.
- `d5r2a_probe8.py` — the [det]_8 = 0 system at the λ = 0 and λ ≠ 0 witnesses
  (mod p): consistent in both cases (dim 9 and dim 12).

Degree-5 rank-3 survivor lifting (ROADMAP P5.A1):

- `survivor_point_r3.py [seed]` (prime from env `HC4_P`, default 32003) — Rabinowitsch-localised rank-3 locus of V(J)
  (dimension 3 for each c tried), sliced to a rational point (`survpt/*.json`).
- `lift_rank9.py` — weight-9 lifting piece at the point: rank 13 of 22, null
  space 9 = 6 gauge + 3 (gauge directions verified).
- `lift_tower.py [timeout] [seed]` — the FULL tower (41 lower-weight unknowns,
  294 positive-weight coefficient equations + Rabinowitsch for the weight-0
  constant) at the point: **inconsistent** for all five points tried (seeds
  21–23 at p = 32003; seeds 24, 25 at p = 40009 via `HC4_P=40009`,
  outputs `*_40009.txt`).
- `lift_tower_sym.py [timeout] [seed]` — same with the survivor point left
  symbolic (fixed c): attempts a mod-p theorem for that c — std did NOT finish
  within the 570 s local budget (`lift_tower_sym_21.txt`); external-machine item.
- `lift_debug.py` — sanity checks on the point (det Hess F ≡ 0, J = 0, gauge).

Ordinary-degree towers for rank ≤ 2 leading structures (research_log Phase 27):

- `graded_tower.py --mode {r2a,r2c,r1,r2c_np,r2c_np1,r2c_np0,r1_np} [--seed N]
  [--p P] [--timeout T] [--lam0] [--incremental] [--x4incr]` — full tower
  (all nonconstant coefficients of det Hess f = 0, Rabinowitsch unit for the
  constant, optional second unit for a "no pivot" condition); on consistency
  extracts a point, re-checks det Hess f in sympy and computes the pivot cone.
  Outputs `tower_<mode>_<seed>[_incr|_lam0|...].txt`; Singular scratch in
  `tower/` (gitignored `*.sing`). Results: `r2c`, `r1`, `r2c_np` (seeds 1, 2)
  inconsistent (`r2c_np` caveat: its random-minor no-pivot condition may clash
  with the forced cube structure — not evidence); `r1` positive-degree part
  consistent (dim 21), fails only at the unit; `r2a`, `r2a --lam0`, `r1_np`
  timed out. `r2c_np1`, `r2c_np0` and the `--x4incr` runs were killed for low
  memory and have no output — rerun at most two at a time.
