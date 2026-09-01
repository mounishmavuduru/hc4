# HC_4 project — roadmap to a complete, hardened, maximized contribution

Provenance: this roadmap is the synthesis of an 8-dimension audit + code review of
the repository (workflow run wf_3db29f87-1f1, 2026-09-01: three independent
math/code reviews of the degree-5 work, a rigor/labeling audit, an open-problem
map, a novelty audit, a reproducibility audit, a publication-readiness audit, and
a completeness critic). Item locations below are `file:line` into the repo.

## 0. The honest ceiling — what "done / perfect / maximized" can mean

HC_4 is an **open** conjecture, and the bridge **HC_4 ⟹ JC_2** makes it the
*stronger* of the two survivors: it cannot be proved without simultaneously
proving JC_2 (famously open; any JC_2 counterexample has max component degree
≥ 108), and it can be refuted only via (a) a *pivot-free* 4-variable
constant-Hessian potential with non-injective gradient, or (b) a JC_2
counterexample lifted through the doubling. So **"the research is over" cannot
mean "HC_4 solved."**

"Done / perfect / optimized / maximized" is achieved when the project's
**partial-results contribution** is: correct, exactly (exact-arithmetic)
verified, **honestly labeled** (known / proved-here / machine-modular /
experimental / conjectured), **reproducible** from a clean checkout,
**novelty-checked**, **expert-reviewed**, and packaged as a submittable artifact —
with every *mechanically-closable* gap closed and the frontier pushed as far as
exact computation and proof allow. The audit's verdict on the current state: the
recent degree-5 work is mathematically **correct** (no false theorem; all traced
failure modes are fail-closed); the remaining work is labeling, reproducibility,
packaging, novelty, external review, plus a small number of bounded computations.

The work below is in execution order. P0–P2 are integrity/hygiene (cheap, high
expected value); P3–P4 are the artifact and its validation; P5 is honest frontier
extension. Nothing in P0–P4 needs new mathematics.

---

## P0 — Preserve & de-risk (do FIRST; near-free, highest expected loss)

The completeness critic's top item, flagged by no single dimension:

- **P0.1 Back up the repository.** `git remote -v` is empty; the entire
  R-D5-GEN / R-D5-TAIL closure exists only on the local branch
  `d5-rank3-pivotfree-closure` (master lags at 91fd6a6). This is the project's
  newest, most novel result with zero backup — the exact loss pattern
  `MEMORY.md` warns about twice. Create a **private** remote, push all
  branches/tags, verify commits `3f5d02d, 8f3013b, 877ebfc, c14db93` landed,
  then merge/fast-forward the reviewed branch into master.
  *(Outward-facing: creating a remote publishes the repo to a host — get the
  owner's go-ahead before pushing. Safe interim with no external exposure:
  `git bundle create hc4-backup.bundle --all` to a durable local path, and merge
  the branch to master locally.)*

---

## P1 — Correctness & honesty integrity (the audit's real defects)

- **P1.1 Fix the confirmed mod-p bug.** `d5_lindir_decide.py:63` `s_poly_modp`
  does `int(co) % p`, which *truncates* rational coefficients instead of taking
  the modular inverse. The E4=0 branch denominators (1/(9c9²), 1/(4c5), …) are
  not cleared by random integer c, so the F_p ideal in `d5_lindir_decide.py` and
  the F_p column of `d5_lindir_verify.py` is **corrupted** (verified: 865 wrong
  coefficients for sub-A/branch-0/p=32003; e.g. 103627404/15007 → 6905, correct
  2090). **Scope (important): isolated to the `d5_lindir_*` unification
  cross-check only** — the R-D5-TAIL certificates use the exact emitter `s_poly`
  on the integer-coefficient family `b3 = y1²·(…)` and are unaffected; the result
  still holds via the ℚ path + R-D5-GEN. Fix: `r = sp.Rational(co); if r.q % p ==
  0: skip prime; c = (int(r.p) * pow(int(r.q), p-2, p)) % p` (or clear
  denominators per generator). Rerun; then correct `lemma_ledger.md:768`
  "over F_p AND over Q" → "over Q (11 a-params)".

- **P1.2 Relabel R-D5-TAIL machine-modular / conditional.** The two decisive
  radical facts — M1 (det Hess₃a5 ∉ √J) and M2 (v*ᵀAv* ∈ √J) — are certified only
  by multi-prime, random-c, random-λ Rabinowitsch, i.e. exactly the modular
  sampling the project's rules forbid inferring a char-0 theorem from. The
  deductive **structural core S1–S4 is genuine char-0**; only M1/M2 are modular.
  Rewrite the R-D5-TAIL header to `[machine-modular; deductive structural core +
  multi-prime radical certificates]`; soften the bolded "is therefore EMPTY" to
  "empty conditional on the char-0 lift of M1/M2, currently certified only over
  F_p at sampled (c, λ)"; **delete the self-contradictory** phrase "the char-0
  power certificate … over F_p" (it is a single-prime, single-c F_p reduction).
  Mirror in `status.md` and `reviewer_packet.md`. Note the sharper true statement
  the reviewer found: emptiness needs only **M2** (a rank<3 cone also exits a
  rank-3 branch), so decouple the conclusion from M1.
  (`lemma_ledger.md:729-762`, `status.md:205-220`.)

- **P1.3 Close or disclose the two completeness gaps.**
  - **solve()-exhaustiveness (undisclosed fail-open):** `_d5_close.py:173`
    enumerates strata via `sp.solve(E4=0 ∧ S, C)`, which has no completeness
    guarantee; a missed component is silently never visited, yet the ledger
    credits only Noetherian termination. Either replace with an ideal-theoretic
    decomposition of the c-locus (Singular `minAssGTZ`/`primdecGTZ` over ℚ, or
    assert 0-dimensionality + solution-count = vdim per node), or add an explicit
    ledger caveat and downgrade to "EMPTY conditional on solve()-completeness".
    (Strongly mitigated by R-D5-GEN + the lindir unification, but currently
    unguarded.)
  - **c0=c1=0 sub-locus (v* = 0):** the excluding direction v* = (0,c1,−c0)
    vanishes on b3 = c2·y1³, so M2's argument is vacuous there and closes
    nothing. Add a dedicated decision for b3 = c2·y1³ (std(J)=(1), or a cone /
    isotropy test on the rank-1 kernel span{e2,e3}) as a distinct leaf; re-scope
    the tail claim if it survives. (`d5_survivor_family.py:66`,
    `_d5_surv_target.py:80-84`.)

- **P1.4 Fix the "pivot-free residual = JC_2" conflation.** `status.md:222` (and
  check `lemma_ledger.md:675`) says the *pivot-free* residual is JC_2-equivalent —
  the inverse of the truth. The **with-pivot** class is JC_2-equivalent (via the
  trichotomy); the pivot-free residual is precisely the part *not* reduced to
  JC_2 (the genuine heart of HC_4). JC_2 is the ultimate blocker simply because
  HC_4 ⟹ JC_2. One-line fix; no mathematics affected.

- **P1.5 Foreground R-D5-GEN as the standalone exact result.** Generic-b3
  emptiness is a real **char-0** Nullstellensatz certificate (`verify_paramcert.py`,
  Σ gᵢTᵢ = 1 over ℚ(c), independently re-checked in sympy; D0, D1 match Z0, Z1).
  Separate it visibly from the modular tail so the tail caveat does not bleed onto
  the exactly-proved part.

---

## P2 — Reproducibility (clean-clone integrity)

- **P2.1 Track the certificates.** `verify_paramcert.py` reads its
  Nullstellensatz witnesses from gitignored, untracked `certificates/_lean/param/`
  → it **errors on a fresh clone**, defeating the "CAS-free independent check".
  Commit the witnesses (and the R-D5-TAIL modular transcripts: surv_rabin /
  surv_iso / iso_pow) into a **tracked** `certificates/certs/`, carved out of the
  `_lean/` gitignore; repoint the loader. Also have `verify_paramcert.py`
  **re-derive** the J generators from the E-system and assert set-equality with
  the loaded ones (closes the last trust gap — currently it trusts the generator
  file).
- **P2.2 Extend the suite runner.** No degree-5 script is in `run_all.py`. Add
  the sympy-only proof-bearing ones (`d5_survivor_family.py`,
  `d5_linear_direction_normalform.py`, `verify_paramcert.py`) to FAST; add the
  WSL-Singular ones to SLOW with a "needs WSL Singular" reason. Regenerate
  `run_all.out.txt` at HEAD; fix README's stale "Suite result of record"
  (currently quotes a transcript that contradicts the committed runner).
- **P2.3 Mark superseded README notes.** Banner the "Singular NOT installed"
  notes (`README.md:361-363,404`) as superseded by the WSL-install section; add a
  "SUPERSEDED — false cone on this family" header to `_d5_surv_slice.py`.
- **P2.4 Keep a tracked doc↔script cross-reference check** (the audit's throwaway
  xref) in `run_all.py` so future edits cannot introduce a dangling reference.

---

## P3 — Publication packaging (editorial; no new mathematics)

- **P3.1 Make `main_result.tex` submittable.** It **does not compile**
  (`\Jac`:191, `\rea`:766 undefined) and — worse — its inline proof of the
  headline Theorem G still prints the **audit-rejected** analytic-continuation
  step (296-297), which its own Remark (334) disavows. Define the macros; excise
  the continuation paragraph; import the repaired proof from `theorem_G_paper.tex`
  (eigensplit/cone-chart or the atkG derivation proof) or cite it as the proof of
  record.
- **P3.2 Fix the packaging.** Two manuscripts redundantly re-prove Theorem G, the
  doubling, and the trichotomy, and **cross-cite each other**. Choose a clean
  one-directional split (recommended: `theorem_G_paper.tex` owns rigidity +
  trichotomy + localization; `main_result.tex` owns the pivot theorem case
  analysis, deg≤4, the obstructions, Lemma R, degree-5, and *cites* the former)
  or a single merge. Remove the duplication.
- **P3.3 One contribution block, partial-results framing leading.** A crisp
  numbered "Main results" statement; lead with "HC_4 open; HC_4 ⟹ JC_2, itself
  open (deg ≥ 108)". Tighten the run-on abstract.
- **P3.4 Write the degree-5 section, correctly scoped.** GEN = exact char-0;
  TAIL = machine-modular; **leading-form level only**; the lifting question
  explicitly **open**. Do not headline; do not call it "proved" in the hand-proof
  sense. Cite the certificate scripts.
- **P3.5 Reconcile cross-document inconsistencies.** Theorem-G proof count (3 vs
  4) unified; T5 (stale "in progress") vs Theorem B; Nagaoka–Yazawa attribution
  of Identity E propagated to `prior_art.md` / `lemma_ledger.md` /
  `reviewer_packet.md`; reviewer_packet's Theorem A summary sharpened (σ≠0 needs
  no degree bound); T3/T4 "planned" → "certified"; L-B provenance; stale header
  dates.
- **P3.6 Submission scaffolding.** Real author/affiliation; a referee summary
  keyed to the papers' theorem numbers (not the internal A/B/G names); complete
  `main_result.tex` bibliography (copy from companion); end-to-end `latexmk`
  build + integrity certificate; move `theorem_G_paper.prev_*.tex` out of the
  submission set.

---

## P4 — Novelty & external validation

- **P4.1 Finish the Theorem G novelty sweep.** Still open — only equiaffine
  sources (Reilly, Fox; both assume the bordered Hessian ≠ 0) were checked; three
  sweep agents died on quota. Because the homogeneous case collapses to the
  (known) Identity E, **all** of Theorem G's novelty is the *inhomogeneous* case,
  so the target is affine differential geometry (MSC 53A15) and degenerate/
  parabolic Monge–Ampère (35J96), **not** the projective vanishing-Hessian corpus.
  Concrete queries and sources in the audit; also re-run the arXiv "Hessian
  conjecture" full-text sweep for 2026-08…09. Keep "claimed new pending expert
  review" until it returns empty; then record the negative result in
  `prior_art.md`.
- **P4.2 Frame/novelty-check the unframed results.** Lemma R ("HC_4 ⟺ every
  constant-Hessian 4-var potential is a polynomial in its own four partials" +
  commuting framing I1–I4) — check the "polynomial-in-partials" equivalence (van
  den Essen's book Ch.1; de Bondt–van den Essen symmetric-Jacobian; Zhao LND
  papers; `lit/garbagnati-repetto/gr.tex:410`), then promote to a stated result.
  Frame the degree-5 method as an **application of comprehensive Gröbner bases**
  (cite Weispfenning 1992 / Suzuki–Sato) — the technique is classical; only the
  HC_4 application is new.
- **P4.3 Independent adversarial re-audit of the CORE.** This audit round only
  re-derived the *degree-5* material; Theorem A, Theorem B, the pivot trichotomy,
  and Q2 got no fresh from-scratch re-derivation. Passing certificates ≠ correct
  theorem statement/reduction (the degree-5 reviews found real logic slips exactly
  by re-deriving). Commission a disprove-intent re-derivation of the A/B/
  trichotomy/Q2 chains, checking case exhaustiveness and hypothesis usage.
- **P4.4 Human expert referee pass.** All verification so far is machine/agent
  side. Route `reviewer_packet.md` + both manuscripts to an affine-geometry /
  polynomial-automorphism specialist (the four questions in reviewer_packet §6).
  Prerequisite to any "submittable" claim.

---

## P5 — Frontier extension (honest classes; do not mis-scope)

**(A) Achievable / compute-bounded — legitimate new work:**
- **P5.A0 — the ≥108 deductive shortcut (try FIRST, before brute lifting).**
  HC_4 ⟹ JC_2 via an *explicit* doubling (Meng 2006 Prop 1.4; `doubling_structure.py`),
  and every JC_2 counterexample has degree ≥ 108 (Guccione et al., cited). A
  pivot-free degree-d constant-Hessian potential induces a JC_2 counterexample of
  a **computable** degree; if that induced degree is < 108 for d = 5 (and small
  d), **no such object can exist** — closing the degree-5 lifting problem by a
  *known theorem* rather than by infeasible elimination. This is degree
  bookkeeping against a cited bound, not an invented lemma — verify it before
  banking on computation.
- **P5.A1 Degree-5 LIFTING decision** (if A0 doesn't close it): impose det Hess₄f
  ≡ const on the full quintic f = f5+f4+f3+f2, solve for the lower-weight
  coefficients, decide solvability and pivot-freeness. Nominally finite
  elimination but likely beyond the local ceiling (export to a bigger CAS); a
  NEGATIVE answer closes degree-5 HC_4, a pivot-free YES is the first pivot-free
  4-var constant-Hessian potential. Do **not** assume it is mechanical.
- **P5.A2 Degree-5 rank ≤ 2 leading-form branches** — not yet addressed; by
  analogy with degree 4 (S1/S2/S3 force a pivot) they are expected to land in the
  with-pivot (JC_2) class, but the degree-5 graded formulas must be derived.
  Required before any "degree-5 leading-form classification complete" statement.
- **P5.A3 char-0 upgrade of M1/M2** — re-run the two full-variety radical tests
  over ℚ(c), or extend the exact lift certificate (as done for R-D5-GEN), on a
  machine without the ~10-min kill. Pure compute-ceiling limitation.

**(B) Needs a new idea — research directions, NOT scheduled compute:**
- Conjecture C2 (descent rigidity): Kronecker structure theory of corank-1
  symmetric pencils M(s,w) over C(w) to force doubling structure.
- Corrected non-affine D1: a route around the Perazzo obstruction (or accept the
  affine-only statement as final). Low priority — strengthens an obstruction on
  *known* HC_5 data, not HC_4.
- Lemma-R local-nilpotency attack: use the commuting framing (I3) + constant
  divergence (I1) to force LND in the pivot-free residual — the one class not
  reduced to JC_2. Speculative.

**(C) Equivalent to JC_2 — cannot be closed without solving JC_2; keep OUT of the
"achievable" bucket:**
- HC_4 restricted to potentials-with-a-pivot (the trichotomy makes this exact).
- General-degree AP4 = the Meng-doubling class. Use only to map which JC_2 degree
  ranges each HC_4 degree stratum touches (degree-5 needs only Moh-range planar
  Keller maps). **Do not list these as compute-bounded or mechanical.**

**(D) Genuinely open / hard — the true heart; never manufacture closure:**
- HC_4 itself.
- **Pivot-free existence:** does a 4-var constant-Hessian potential with *no*
  pivot exist? If NO, HC_4 ⟺ JC_2 outright. If YES, that is the only place an
  HC_4-specific counterexample lives. Known bounds (PE-2/3/4): degree ≥ 5, ≥ 2
  storeys of degree ≥ 3, dim span{N_α} ≥ 4; such potentials do exist one
  dimension up (Meng–Yang HC_5, empty pivot cone). This is where all
  degree-by-degree work is ultimately aimed.

---

## Definition of done (checklist)

- [ ] Repo backed up (remote + master merged).
- [ ] Zero confirmed code defects; every certificate fail-closed and, where
      claimed char-0, actually char-0 (or explicitly relabeled modular).
- [ ] Every claim's label matches what its certificate establishes; no bolded
      "proved" resting on modular sampling.
- [ ] Clean-clone reproducible: `run_all.py` green at HEAD incl. the sympy-only
      degree-5 proofs; all cited witnesses tracked.
- [ ] One coherent, compiling, partial-results artifact; degree-5 written up with
      honest scope.
- [ ] Theorem G novelty sweep complete (or scoped in-paper); attributions
      consistent; core results independently re-audited; expert referee pass done.
- [ ] Frontier documented in the four honest classes; the achievable computations
      (A0 first) attempted; open/JC_2-equivalent items clearly marked as such.

The end state is a **maximized partial-results artifact on HC_4**, not a claimed
resolution. Reaching it is editorial, computational (bounded), and validation
work — not a proof of HC_4, which remains open and reduces to JC_2.
