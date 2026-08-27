# Classical theorems for the HC_4 project: primary-source extracts

Verified 2026-08-08. All quotes are verbatim from the stated primary sources, which are
downloaded locally where indicated. TeX-ASCII transliteration of formulas is used; page
numbers refer to the original journal pagination.

Local primary sources used here:
- `lit/bcw/bcw1982.pdf` — H. Bass, E. H. Connell, D. Wright, *The Jacobian Conjecture:
  reduction of degree and formal expansion of the inverse*, Bull. Amer. Math. Soc. (N.S.)
  **7** (1982), no. 2, 287–330. (Open access at AMS.)
- `lit/bbr/bbr1962.pdf` — A. Białynicki-Birula, M. Rosenlicht, *Injective morphisms of
  real algebraic varieties*, Proc. Amer. Math. Soc. **13** (1962), 200–203.
- `lit/ax/ax1969-pjm.pdf` — J. Ax, *Injective endomorphisms of varieties and schemes*,
  Pacific J. Math. **31** (1969), no. 1, 1–7.
- `lit/cynk-rusek/cynk-rusek1991.pdf` — S. Cynk, K. Rusek, *Injective endomorphisms of
  algebraic and analytic sets*, Ann. Polon. Math. **56** (1991), no. 1, 29–35.
- `lit/moh/moh_00000144.jpg … moh_00000147.jpg` — page images (GDZ Göttingen) of
  T. T. Moh, *On the Jacobian conjecture and the configurations of roots*, J. Reine
  Angew. Math. **340** (1983), 140–212. DOI 10.1515/crll.1983.340.140, Zbl 0525.13011.

Throughout, BCW = Bass–Connell–Wright 1982.

---

## 1. Wang's theorem: quadratic Keller maps are invertible

**Original source.** S. S.-S. Wang, *A Jacobian criterion for separability*, J. Algebra
**65** (1980), no. 2, 453–494. DOI 10.1016/0021-8693(80)90233-1, MR 585736. (Elsevier;
full text not freely fetchable — statement verified below from BCW's published
reproduction, which cites this paper as [Wa].)

**Statement (BCW, p. 298, Theorem (2.4)), verbatim:**

> (2.4) THEOREM (S. WANG [Wa]). *Let k be a field of characteristic ≠ 2. Let
> F ∈ End(A^n_k) have invertible Jacobian J(F). Assume that deg(F) ≤ 2. Then F
> is invertible.*

Notes:
- deg(F) ≤ 2 means each component has total degree ≤ 2 — **not** necessarily
  homogeneous (BCW (1.3) Remark 1, p. 292: "let F ∈ End(A^n_k) be a quadratic map
  (deg(F) ≤ 2) with J(F) invertible. Then S. Wang [Wa] (see (2.4) below) has shown
  that F is invertible, and he conjectured that its inverse G has degree ≤ 2^{n-1}.").
- Any dimension n; any field of characteristic ≠ 2 (in particular all of char 0).
- Wang's conjecture deg(F^{-1}) ≤ 2^{n-1} is TRUE by BCW Corollary (1.4) (see §5 below).
- Per BCW footnote 5 (p. 301), A. V. Yagzhev ("Jagžev") and S. Oda each independently
  rediscovered Wang's theorem.

**The "midpoint" proof (due to S. Oda, reproduced verbatim in BCW, p. 298):**

> This result was rediscovered by S. Oda [O], with a much simpler proof than Wang's.
> A slight elaboration of Oda's argument is presented in [Wr2, Lemma 3.5]. The argument
> is so short that we reproduce it here. In view of Theorem (2.1), condition (c), it
> suffices to show that F is injective. Suppose, on the contrary, that F(a) = F(b) with
> a ≠ b in k^n. Replacing F(X) by G(X) = F(X + a) − F(a), we have 0 = G(0) = G(c) where
> c = b − a ≠ 0, and G is still quadratic; write G = G_(1) + G_(2) where the components
> of G_(d) are homogeneous of degree d (d = 1, 2). Then
>
>     0 = G(c) = G_(1)(c) + G_(2)(c)
>              = G_(1)(c) + 2 t_0 G_(2)(c)        (t_0 = 1/2)
>              = d/dt ( G_(1)(c) t + G_(2)(c) t^2 ) |_{t = t_0}
>              = d/dt ( G(tc) ) |_{t = t_0} = J(G)(t_0 c) · c.
>
> Since J(G), like J(F), is invertible, and c ≠ 0, this is a contradiction.

(The evaluation point t_0 = 1/2 is the midpoint of the segment [0, c] — hence "midpoint
proof". Injectivity plus Theorem (2.1)(c) ⇒ (a) then gives invertibility; over C one can
instead invoke any of the injective-⇒-automorphism results of §2.)

Oda's preprint: Susumu Oda, *The Jacobian problem and the simply-connectedness of A^n
over a field k of characteristic zero*, Osaka Univ., preprint, 1980 (BCW reference [O]).
Elaboration: D. Wright, *On the Jacobian Conjecture*, Illinois J. Math. **25** (1981),
423–440, Lemma 3.5 (BCW reference [Wr2]).

---

## 2. Injectivity implies automorphism

### 2.1 Białynicki-Birula–Rosenlicht 1962 (primary, full text in lit/bbr/)

Proc. Amer. Math. Soc. **13** (1962), 200–203. Received March 15, 1961.

**§3, p. 201, verbatim:**
> If f: R^n → R^n *is an injective morphism* (of real algebraic sets, e.g. a real
> polynomial map), *then f is also surjective.*

**§5, p. 203, verbatim (the algebraically closed case):**
> We remark finally the following easy proof that *if f: k^n → k^n is an injective
> polynomial map, where k is any algebraically closed field, then f is also surjective*:
> Here f must be birational or purely inseparable and a result of Chevalley [2, p. 195]
> shows that f is open, so f(k^n) = k^n − X, with X closed. Whether or not f^{-1} is
> defined at a particular point depends on the poles of rational functions, so either
> dim X = n − 1 or X is empty. If X were nonempty a nonconstant polynomial function on
> k^n with zero locus contained in X would give a similar function with no zero locus,
> which is impossible.

Note: BBR state surjectivity; over char 0 the argument shows f is birational, and
polynomiality of the inverse then follows (e.g. via ZMT, or from BCW Theorem (2.1)
(b) ⇒ (a)). Cynk–Rusek (below) attribute to [BBR] the statement "every injective
polynomial transformation of k^n is a polynomial automorphism".

### 2.2 Ax's theorem (primary, full text in lit/ax/)

J. Ax, *Injective endomorphisms of varieties and schemes*, Pacific J. Math. **31**
(1969), 1–7. (The variety case was first observed in J. Ax, *The elementary theory of
finite fields*, Ann. of Math. **88** (1968), 239–271, §14; the PJM paper says: "A
finiteness property (Corollary 1 to the theorem) of algebraic varieties observed in
[1, §14] is that every injective endomorphism of a variety is surjective.")

**PJM p. 1, verbatim:**
> THEOREM. *Let Y be a scheme of finite type over a scheme X. Let Y --φ--> Y be an
> X-morphism. If φ is injective then φ is surjective.*
>
> COROLLARY 1. *Let Y be an algebraic variety over an algebraically closed field k.
> Let Y --φ--> Y be a morphism. Assume that the induced mapping φ(k) of the k-valued
> points Y(k) of Y to Y(k) is injective. Then φ(k) is surjective.*

The proof is by model-theoretic transfer from finite fields ("Lefshetz Principle"),
p. 2–3. Ax also credits (p. 7): Borel's cohomological proof (char 0), Shimura's
reduction-mod-p argument, and "the first and only previous result of this kind was
obtained by A. Bialynicki-Birula and M. Rosenlicht who gave a simple proof in [2] of
the special case of Corollary 1 when Y is affine space A^n."

### 2.3 Cynk–Rusek 1991: injective ⇒ (biregular) automorphism (primary, in lit/cynk-rusek/)

Ann. Polon. Math. **56.1** (1991), 29–35.

**Abstract, verbatim:** "We prove that every injective endomorphism of an affine
algebraic variety over an algebraically closed field of characteristic zero is an
automorphism."

**Theorem 2.2, p. 32, verbatim:**
> THEOREM 2.2. *Let V be an affine algebraic variety and let F : V → V be a regular
> mapping. Then the following statements are equivalent:*
> (i) *F is injective.*
> (ii) *F is a bijection of V.*
> (iii) *F is an automorphism of V.*

(Standing hypothesis of their §2: all varieties are over a field k algebraically closed
of characteristic zero. (i) ⇒ (ii) is "the Ax Theorem"; (ii) ⇒ (iii) is their
contribution — the inverse is regular, i.e. for V = k^n the set-theoretic inverse is a
polynomial map.)

**Grothendieck's scheme version, as quoted by Cynk–Rusek (p. 29), verbatim:**
> Grothendieck proved in [5, Prop. 17.9.6] a counterpart of the Zariski Main Theorem for
> the category of S-preschemes, where S is a fixed prescheme. His theorem says that
> every injective S-endomorphism of an S-prescheme of finite presentation is an
> automorphism.

([5] = A. Grothendieck, EGA IV (Étude locale des schémas et des morphismes de schémas),
Inst. Hautes Études Sci. Publ. Math. **32** (1967). Quoted here via Cynk–Rusek; the EGA
text itself was not fetched.)

### 2.4 Rudin 1995 (NOT verified from primary source)

W. Rudin, *Injective polynomial maps are automorphisms*, Amer. Math. Monthly **102**
(1995), no. 6, 540–543 — an elementary (no algebraic geometry) proof that an injective
polynomial map C^n → C^n is bijective with polynomial inverse. Paywalled (JSTOR); exact
statement not independently verified here. Cite Cynk–Rusek or BCW (2.1) instead if a
verified source is required.

### 2.5 BCW Theorem (2.1): the equivalence package over a field (primary, in lit/bcw/)

**BCW p. 287 (Introduction), verbatim:** "How do we recognize when such an F is
invertible? The question is unambiguous since, once F is bijective, its set theoretic
inverse is automatically polynomial (see Theorem 2.1)."

**BCW p. 288, verbatim:** "In contrast, if F: C^n → C^n is polynomial, det J(F) = 1,
and F is injective, then F is invertible (see Theorem 2.1)."

**BCW pp. 294–295, Theorem (2.1), verbatim:**
> (2.1) THEOREM. *Let k be a field and let F: A^n_k → A^n_k be a morphism with J(F)
> invertible. Consider the following conditions.*
> (a) *F is invertible, i.e. k[F] = k[X].*
> (b) *F is birational, i.e. k(F) = k(X).*
> (c) *F: k^n → k^n is injective.*
> (d) *The integral closure \bar{k[F]} of k[F] in k[X] is unramified over k[F].*
> (e) *k[X] is a finitely generated k[F]-module.*
> (e′) *F is proper.*
> (f) *k[X] is a projective k[F]-module.*
> (g) *k(X) is galois over k(F).*
> *If char(k) = 0 the conditions are all equivalent. In general we have*
> (g) ⇐ (a) ⇔ (b) ⇐ (c),  (a) ⇒ (e),  (d) ⇐ (e) ⇔ (e′) ⇐ (f)
> *except that k should be infinite for (c) ⇒ (b).*

(The implications (d) ⇒ (a) and (g) ⇒ (a) "invoke the simple connectivity of C^n";
(g) ⇒ (a) was first proved by Campbell [Math. Ann. 205 (1973), 243–248].)

---

## 3. Bass–Connell–Wright / Yagzhev reduction to degree 3 (with the dimension cost)

Primary source: BCW 1982, Chapter II (lit/bcw/bcw1982.pdf).

**Summary statement (BCW p. 288, Introduction), verbatim:**
> In Chapter II it is shown that the Jacobian Conjecture will follow once it is shown
> for all F = (F_1,...,F_n) of the form F_i = X_i − H_i where each H_i is a cubic
> homogeneous polynomial, and the matrix J(H) = (∂H_i/∂X_j) is nilpotent. (This
> contrasts temptingly with Wang's proof of the conjecture for quadratic F.)

**The Reduction Theorem (BCW p. 304, Theorem (2.1) of Ch. II), verbatim:**
> (2.1) THEOREM. *Let k be a commutative ring, and let F ∈ MA^1_n(k) have invertible
> Jacobian J(F). There exist an integer m ≥ 0, elements G, H ∈ EA^0_{n+m}(k), and
> \tilde{F}(T) ∈ MA_{n+m}(k[T]), where T is an indeterminate, with the following
> properties.*
> (a) *For T = 1 we have \tilde{F}(1) = G ∘ F^{[m]} ∘ H. Thus if \tilde{F} is invertible
> then so also is F.*
> (b) *The k[T]-algebra endomorphism φ_{\tilde F} of k[T]^{[n+m]} defined by \tilde F
> can be viewed as a k-algebra endomorphism of k[X_1,...,X_{n+m},T] = k^{[r]}, where
> r = n + m + 1. As such it defines an element L ∈ MA_r(k), which is invertible iff
> \tilde F is. We have L = X_r + N where N is cubic homogeneous (N = N_(3)), and linear
> in each variable except quadratic in T, and J(N) is nilpotent.*

**Corollary (BCW p. 304, (2.2)), verbatim:**
> (2.2) COROLLARY. *Suppose, for all n and all F ∈ MA_n(k) of the form F = X + N with N
> cubic homogeneous and J(N) nilpotent, that F is invertible. Then for all n and all
> F ∈ MA_n(k) with J(F) invertible, F is invertible.*

**The dimension cost is explicit:**
- F^{[m]} = (F_1,...,F_n, X_{n+1},...,X_{n+m}) is the *stabilization* of F to n + m
  variables (BCW p. 304): "This is implemented here as follows: Let F = (F_1,...,F_n) ∈
  End(A^n_k) have invertible Jacobian J(F). For any m ≥ 0 consider F^{[m]} = (F_1,...,
  F_n, X_{n+1},...,X_{n+m}) ∈ End(A^{n+m}_k)... F is invertible iff F^{[m]} is
  invertible." (p. 300).
- Degree-lowering step (Proposition (3.1), p. 305): each elimination of one monomial of
  top degree d ≥ 4 costs **2 extra variables** (the elementary maps G, H live in
  EA^1_{n+2}(k)); iterating brings deg ≤ 3.
- Unipotent/homogenization step (Proposition (5.2), p. 308), verbatim:
  > (5.2) PROPOSITION. *Let F ∈ MA^1_n(k) be of degree ≤ d + 1. There exist elements
  > G, H ∈ EA^0_{dn}(k) such that the element F′ = G ∘ F^{[(d−1)n]} ∘ H has the
  > following form. Writing X_{dn} = (X_n^{(1)},...,X_n^{(d)}) as in Proposition (5.1),
  > we have F′ = X_{dn} + N where N = (F_(2) − X_n^{(d)}, F_(d+1), F_(d) − X_n^{(2)},
  > ..., F_(3) − X_n^{(d−n)}). If J(F) is invertible then J(N) is nilpotent.*
  So a map of degree ≤ d+1 in n variables becomes X + N (N as displayed, J(N)
  nilpotent) in **dn variables**. BCW's history footnote 4 (p. 292): John Tyrrell "also
  knew, but never published, the 'reduction to degree 3' of Proposition (3.1) in
  Chapter II below."

**Yagzhev's independent proof (BCW footnote 5, p. 301), verbatim:**
> Since writing this paper L. Avramov has kindly brought to our attention the very
> interesting paper of Jagžev [Ja] (1980). Jagžev proves a theorem substantially
> equivalent to our Reduction Theorem (2.1) of Chapter II, though with some differences
> of detail, and by somewhat different methods. He rediscovers the implication
> (c) ⇒ (a) of Theorem (2.1) above, and uses it, as in [Wr2], to prove Wang's Theorem
> ((2.4) above), which, he, like Oda [O], rediscovered.

Reference [Ja]: A. V. Jagžev (Yagzhev), *On Keller's problem*, Siberian Math. J. **21**
(1980), 141–150 (Russian original: Sibirsk. Mat. Zh. 21 (1980), no. 5, 141–150).

**Degree-2 in char p caveat (BCW p. 298), verbatim:** "This theorem tempts one to
conjecture, as does Oda, that the Jacobian Conjecture might be true in characteristic
p > 0 for morphisms F of degree d < p. However this cannot be so since we show below
that the Jacobian Conjecture (in any characteristic) follows once it is known for
morphisms of degree ≤ 3."

---

## 4. Moh: JC_2 holds up to degree 100

Primary source: T. T. Moh, *On the Jacobian conjecture and the configurations of
roots*, J. Reine Angew. Math. (Crelle) **340** (1983), 140–212. Verified from GDZ page
images (lit/moh/moh_00000144.jpg = p. 140, ..., moh_00000147.jpg = p. 143). Article =
GDZ PPN243919689_0340, structure LOG_0011 (images 144–216); full-article PDF:
https://gdz.sub.uni-goettingen.de/download/pdf/PPN243919689_0340/LOG_0011.pdf

**p. 140 (Introduction), verbatim:**
> The well-known Jacobian conjecture states "let k be an algebraically closed field of
> characteristic zero and f(x, y), g(x, y) two polynomials over k. If
> f_x g_y − f_y g_x = 1 then k[x, y] = k[f, g]." ... Throughout this article the field
> k is assumed to be algebraically closed and of characteristic zero.

**p. 143, verbatim (the main theorem and the degree-100 computation):**
> Then we compile the results in Propositions 5.4 and 6.1 to be the main theorem of
> this article.
>
> In Appendix I we collect all needed propositions of the polynomial solutions of the
> differential equations involved. There is a simple computer program essentially based
> on the discussion at the end to determine the possible counterexamples to the
> Jacobian conjecture for all polynomials of degrees less than or equal to 100. In fact
> there are only very few possibilities. However the number of coefficients of those
> possible counterexamples range from 3370 to 7328. In Appendix II we give some
> theoretic arguments to reduce the number of coefficients to at most 12. Thus it
> becomes a computable problem to show that the Jacobian conjecture is valid for all
> polynomials of degrees less than or equal to 100. The direct computation is presented
> in Appendix II.
>
> The computer program is a routine one and will be furnished to anyone upon request.
> It uses 147.148 sec of a CDC 6500 and costs only $2.41.

**Cross-check (BCW p. 294), verbatim:** "Moreover Moh [Mo], using characteristic pairs
and a computer search, has proved the conjecture when d_1, d_2 are ≤ 100." (d_i =
deg(F_i), n = 2. BCW p. 288: "Notable among these is Moh's proof of the Jacobian
Conjecture for n = 2 when deg(F) ≤ 100 [Mo].")

So the citable statement is: **JC_2 holds for all Keller maps (f, g) over an
algebraically closed field of characteristic zero with max(deg f, deg g) ≤ 100**, by
theory (Props. 5.4, 6.1) plus a finite computer verification (Appendix II).

---

## 5. Keller-map basics: formal inverse, étale, degree of the inverse

All from BCW (lit/bcw/bcw1982.pdf) unless stated.

### 5.1 Keller maps; the Jacobian Property

**BCW p. 289, verbatim:** Chain rule "(1) J(G(F)) = J(G)(F) · J(F)"; hence
"(2) F ∈ GA_n(k) ⇒ J(F) ∈ GL_n(k^{[n]}). In this case det J(F) is a unit of k^{[n]}.
When k is reduced (without nonzero nilpotent elements) the units of k^{[n]} are just
the (constant) units k^× of k. The Jacobian Problem asks about the converse of (2),
i.e. about the validity of the following 'Jacobian Property':
JP_n(k): if F ∈ End(A^n_k) and J(F) is invertible then F is invertible."

Keller's original paper: O. H. Keller, *Ganze Cremona-Transformationen*, Monats. Math.
Physik **47** (1939), 299–306 (BCW reference [K]; BCW p. 288: "The Jacobian Conjecture
seems first to have been formulated by O. H. Keller in 1939.").

Char p failure (BCW (1.1) Remark 5, p. 290): "If k is a field of characteristic p > 0
then JP_n(k) fails for all n ≥ 1. Just take F_1 = X_1 + X_1^p and F_i = X_i for i ≥ 2."

Lefschetz principle (BCW (1.1) Remark 4, p. 290): "JP_n(C) ⇒ JP_n(k) for every integral
domain k of characteristic 0."

### 5.2 The formal inverse

**BCW (1.1) Remark 2, pp. 289–290, verbatim:**
> Suppose that F(0) = 0 and that J(F)(0) is invertible. Then (Implicit Function
> Theorem) F is formally invertible at the origin. In other words there is a formal
> inverse G = (G_1,...,G_n) defined by the conditions, G_i ∈ k^{[[n]]} =
> k[[X_1,...,X_n]], G_i(F) = X_i (i = 1,...,n). ... *In order that F ∈ GA_n(k) it is
> necessary and sufficient that the power series G_i be polynomials.*
>
> In any case the map φ_F: f ↦ f(F) defines an automorphism of k^{[[n]]}, and hence an
> *injective* endomorphism of k^{[n]}, i.e. F_1,...,F_n are "algebraically independent"
> over k. Invertibility of F is equivalent to the condition k[F_1,...,F_n] =
> k[X_1,...,X_n].

**Abhyankar's inversion formula (BCW Ch. III, p. 312), verbatim:**
> (2.2) COROLLARY. *Let G = (G_1,...,G_n) be the inverse of F: G_i(F) = X_i for
> i = 1,...,n. Then* G_i = Σ_{p ∈ N^n} D^{[p]} ( X_i · j(F) · H^p ),
where H = X − F, j(F) = det J(F), D^{[p]} = D_1^{p_1}···D_n^{p_n}/(p_1!···p_n!)
(Theorem (2.1), attributed to Abhyankar [Ab2, Purdue lecture notes 1974]; the formula
goes back through Gurjar/Goursat, and in the diagonal case to Good 1960 and Jacobi
1830).

**Tree expansion (BCW Ch. III, Theorem (4.1), p. 321):** for F = X − H with H_i
homogeneous of degree δ ≥ 2 and j(F) = 1: G_i = Σ_{d≥0} G_i^{(d)}, G_i^{(0)} = X_i,
G_i^{(1)} = H_i, and for d ≥ 2, "(12) G_i^{(d)} = Σ_T (1/α(T)) Σ_f P_{T,f}" where T
runs over isomorphism classes of rooted trees with d vertices, α(T) = Card Aut(T), f
over i-rooted labelings, P_{T,f} = Π_{v∈V(T)} (D_{f_{v_+}} H_{f_v}). "G_i^{(d)} is a
homogeneous polynomial of degree d(δ − 1) + 1, and G_i is a polynomial iff G_i^{(d)} = 0
for all d ≫ 0."

### 5.3 A Keller map is étale

**BCW p. 295–296 (proof of Theorem (2.1)), verbatim:**
> In the exact sequence Ω_{k[F]/k} ⊗_{k[F]} k[X] --J--> Ω_{k[X]/k} → Ω_{k[X]/k[F]} → 0
> the matrix of J with respect to the bases (dF_i ⊗ 1) and (dX_i), respectively, is
> J(F). Thus
> (1) J(F) invertible ⇔ Ω_{k[X]/k[F]} = 0 ⇔ k[X]/k[F] is unramified.
> If x ∈ k^n and y = F(x) the inclusion of local rings k[F]_y → k[X]_x is unramified
> and they have the same residue class field k, hence the same completions [A+K, VI,
> (3.7)]. Consequently F is étale at x [A+K, VI, (4.5)]. Passing to the algebraic
> closure of k we see that F, being étale at all closed points, is étale [A+K, VI,
> (4.6)]. Thus
> **(2) J(F) invertible ⇔ k[X]/k[F] is étale (= flat and unramified).**
> Flatness has the following useful consequence. (3) k[X] ∩ k(F) = k[F].

Also (p. 296): "(e) ⇔ (e′): A theorem of Chevalley asserts the equivalence of
properness and finiteness for any affine morphism... Note that when k = C properness in
the topological and algebraic-geometric senses coincide."

### 5.4 Étale + (universally) injective ⇒ open immersion

**Stacks Project, Tag 025G (Theorem 41.14.1), verbatim** (fetched 2026-08-08 from
https://stacks.math.columbia.edu/tag/025G):
> Theorem 41.14.1. Let f : X → Y be a morphism of schemes. The following are
> equivalent: (1) f is an open immersion, (2) f is universally injective and étale, and
> (3) f is a flat monomorphism, locally of finite presentation.

(For a Keller map F: C^n → C^n, which is étale by §5.3, injectivity on C-points implies
universal injectivity — F is injective with separable residue field extensions — so an
injective Keller map is an open immersion; combined with Ax surjectivity (§2.2) this
gives an isomorphism. Related: an étale monomorphism is an open immersion, EGA IV
17.9.1.)

### 5.5 Degree bound for the inverse (Gabber)

**BCW p. 292, verbatim:**
> (1.4) COROLLARY. *If k is a field and F ∈ GA_n(k) then* deg(F^{-1}) ≤ deg(F)^{n-1}.
>
> The following result was communicated to us by Ofer Gabber... He attributed it to an
> unrecalled colloquium lecturer at Harvard. [Footnote 4: Since writing this paper John
> Tyrrell (Kings College, Univ. of London) has indicated that this result was "well
> known" to the classical geometers.]
>
> (1.5) THEOREM. *Let f: P^n --→ P^n be a birational map with inverse f^{-1}. Then*
> deg(f^{-1}) ≤ deg(f)^{n-1}.

Consequences (BCW p. 292–294): Wang's conjecture deg(G) ≤ 2^{n-1} for quadratic F is
true ((1.3) Remark 1); the group G_(d) of degree-≤d automorphisms is a *closed*
subvariety of E_(d) (Corollary (1.6)), via: "According to Corollary (1.4) we have
F ∈ G_(d) ⇔ H_(r) = 0 for all r ≥ d^{n-1}" (p. 294).

---

## 6. One-line summary table for citation

| # | Statement | Cite | Verified from |
|---|-----------|------|---------------|
| 1 | Keller map of degree ≤ 2 (any n, char ≠ 2, not nec. homogeneous) is invertible | Wang, J. Algebra 65 (1980) 453–494; statement + Oda's midpoint proof in BCW Thm (2.4) | BCW pdf (primary reproduction) |
| 2 | Injective polynomial self-map of k^n (k alg. closed) is surjective | Białynicki-Birula–Rosenlicht, Proc. AMS 13 (1962) 200–203, §3 and §5 | BBR pdf |
| 3 | Injective endomorphism of a finite-type scheme is surjective (Ax–Grothendieck) | Ax, Pacific J. Math. 31 (1969) 1–7 (also Ann. of Math. 88 (1968) §14); Grothendieck EGA IV 17.9.6 | Ax pdf; EGA via Cynk–Rusek quote |
| 4 | Injective endomorphism of affine variety over alg. closed char-0 field is an automorphism (regular inverse) | Cynk–Rusek, Ann. Polon. Math. 56 (1991) 29–35, Thm 2.2 | C–R pdf |
| 5 | JC reduces to F = X + N, N cubic homogeneous, J(N) nilpotent — at the cost of increasing dimension (stabilization; deg ≤ d+1 in n vars ⇒ dn vars) | BCW Bull. AMS 7 (1982), Ch. II Thm (2.1), Cor (2.2), Props (3.1), (5.2); Yagzhev, Siberian Math. J. 21 (1980) 141–150 | BCW pdf |
| 6 | JC_2 true for deg ≤ 100 | Moh, Crelle 340 (1983) 140–212 | GDZ page images |
| 7 | J(F) invertible ⇔ F étale; injective ⇔ invertible over char 0; bijective ⇒ polynomial inverse; formal inverse exists and F invertible iff it is polynomial | BCW Thm (2.1), Remark (1.1)2, eq. (2) p. 296 | BCW pdf |
| 8 | Universally injective + étale ⇔ open immersion | Stacks Project Tag 025G (Thm 41.14.1) | stacks.math.columbia.edu |
| 9 | deg(F^{-1}) ≤ deg(F)^{n-1} | BCW Cor (1.4)/Thm (1.5) (Gabber) | BCW pdf |
