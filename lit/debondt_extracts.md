# Verbatim extracts: de Bondt, "Polynomials with constant Hessian determinants in dimension three"

**Source (primary):** arXiv:1203.6605 e-print tarball (v3, 20 Jan 2015), downloaded 2026-08-08 to
`C:\Users\mouni\hc4\lit\debondt\` (files `consthess3.tex`, `consthess3.bbl`). Published as
J. Pure Appl. Algebra 219 (2015), DOI 10.1016/j.jpaa.2014.12.020. Versions on arXiv: v1 (29 Mar 2012),
v2 (28 Mar 2013), v3 (20 Jan 2015). All quotes below are verbatim TeX from `consthess3.tex`
(line numbers refer to that file).

Throughout, `\hess` = Hessian matrix, `\jac` = Jacobian matrix, `\grad` = gradient; K is an
arbitrary field of characteristic zero.

---

## 1. Main theorem (HC_2 and HC_3), verbatim (lines 263-271)

```tex
\begin{theorem}[Main result] \label{dillenext}
Let $K$ be a field of characteristic zero and $f \in K[x] = K[x_1,x_2,\ldots,x_n]$
be a polynomial of degree $d$. If $n \le 3$, then $\grad f$ satisfies the Jacobian conjecture.

If $n \le 3 \le d$ and $\det \hess f \in K$, then there exists a $T \in \GL_n(K)$ such that all
entries below the anti-diagonal of the Hessian of $f(Tx)$ are zero. In particular, the quadratic
part of $f$ vanishes somewhere at $K^n \setminus \{0\}^n$ when $2 \le n \le 3$, namely at the last
column of $T$.
\end{theorem}
```

Notes on the exact form of the invertibility claim:
- "\grad f satisfies the Jacobian conjecture" means: if the Keller condition
  det Hess f in K^* holds, then grad f is invertible. Part iii) of the proof (lines 420-424)
  reduces to `F := \grad f(Tx)` with anti-triangular `\jac F = \hess(f(Tx))` and invokes
  Lemma \ref{antitri}.
- The inverse is genuinely POLYNOMIAL: Lemma \ref{antitri} and its proof (lines 429-443, verbatim):

```tex
\begin{lemma} \label{antitri}
Suppose that $A$ is a commutative $\Q$-algebra and $F \in A[x]^n = A[x_1,x_2,\ldots,x_n]^n$. If all entries
below the anti-diagonal of $\jac F$ are zero and $\det \jac F \in A^{*}$, then $F$ is invertible.
\end{lemma}

\begin{proof}
Since the entries below the anti-diagonal
of $\jac F$ are zero, we see that $F_{n+1-i} \in A[x_1, x_2, \ldots, x_i]$
for each $i$. Hence it follows from $\det \jac F \in A^{*}$, that the entries
on the anti-diagonal of $\jac F$ are nonzero constants. Say that these constants are
$c_1, c_2, \ldots, \allowbreak c_n$ from left to right. Then $F_n - c_1 x_1 \in A$ and
$F_{n+1-i} - c_i x_i \in A[x_1, x_2, \ldots, x_{i-1}]$ for all $i \ge 2$. By induction on $i$,
we obtain that $x_i \in A[F_{n+1-i},F_{n+2-i},\ldots,F_n]$ for each $i$,
so $F$ is invertible. \qedhere
\end{proof}
```

So x_i lies in the polynomial ring A[F_{n+1-i},...,F_n]: the inverse map is polynomial.

Abstract statement (lines 62-68, verbatim):

```tex
In this paper, we show that the Jacobian conjecture holds for gradient maps
in dimension $n \le 3$ over a field $K$ of characteristic zero. We do this by
extending the following result for $n \le 2$ by F. Dillen to $n \le 3$:
if $f$ is a polynomial of degree larger than two in $n \le 3$ variables such that the
Hessian determinant of $f$ is constant, then after a suitable linear transformation
(replacing $f$ by $f(Tx)$ for some $T \in \GL_n(K)$), the Hessian
matrix of $f$ becomes zero below the anti-diagonal. The result does not hold for larger $n$.
```

Important precision: the anti-triangular normal form is proved for det Hess f in K
(constant, INCLUDING zero), with n <= 3 <= d; the JC-for-gradient-maps consequence is the
case det Hess f in K^*.

## 2. Anti-triangular normal form: scope and FAILURE for n >= 4

- Proved for: n <= 3 <= d = deg f, det Hess f constant (in K). The d >= 3 hypothesis is
  necessary over non-algebraically-closed fields (lines 273-278: `f = \frac12 x_1^2 + ... +
  \frac12 x_n^2` over R admits no such T).
- Stated to FAIL for n >= 4. Abstract: "The result does not hold for larger $n$."
  Introduction (lines 216-221, verbatim):

```tex
As mentioned above, our main result (in theorem \ref{dillenext}) is a similar statement for polynomials
with constant Hessian determinant in dimension three:
if the degree of the polynomial is larger than two, then after a suitable linear transformation,
every entry below the anti-diagonal of its Hessian matrix becomes zero.
This result does not hold in dimensions larger than three and neither for quadratic polynomials over
$\R$ in dimensions two and three.
```

### Explicit 4-variable counterexample (lines 242-256, verbatim)

```tex
Now that we know how linear transformations influence the Hessian, we are able to see that Dillen's result
cannot be extended to dimension four. Take for instance
\begin{equation} \label{dillen4}
f = (x_1 + x_2^2) x_3 + (x_2 + (x_1 + x_2^2)^2) x_4
\end{equation}
Then the cubic part of $f$ is equal to $x_2^2 x_3 + x_1^2 x_4$, and the rows of its Hessian, whose
entries are linear forms, are independent over $\C$. This is maintained after a linear transformation,
so if we could obtain by way of a transformation that
the Hessian of $f$ becomes zero below the anti-diagonal, the lower left corner of $\hess f$ would
get a nontrivial linear part, because the last row of the Hessian of the cubic part of
$f$ cannot become zero. This however contradicts that the Hessian determinant is a
nonzero constant. The polynomial $f$ in (\ref{dillen4}) was made by applying the reduction
of the Jacobian conjecture to gradient maps, as in \cite[Prop.\@ 1.4]{MR2221506}, on the planar
invertible map $F = (x_1 + x_2^2, x_2 + (x_1 + x_2^2)^2)$.
```

([Men] = MR2221506 = Guowu Meng, "Legendre transform, Hessian conjecture and tree formula",
Appl. Math. Lett. 19(6):503-510, 2006.)

After Lemma antitri (lines 445-456) the paper also shows this same f defeats Theorem
\ref{weight} in n = 4 (verbatim):

```tex
We showed earlier that $f$ in \eqref{dillen4} is a counterexample to theorem \ref{dillenext} with
$n = 4$. Using some techniques in the above proof, one can show that $f$ is a counterexample
to theorem \ref{weight} with $n = 4$ as well. For that purpose, take any $T \in \GL_4(K)$
and any weight function $w$ such that $0 < w(x_1) \le w(x_2) \le w(x_3) \le w(x_4)$.

Since the rows of the Hessian of the cubic part of $f$ are independent, there exists
a $j$ such that $(\hess f(Tx))_{nj} \notin K$. Furthermore, $\det \hess f(Tx) \ne 0$
tells us that there exists an $i \ge 3$ and a $k \ge 2$ such that $(\hess f(Tx))_{ik} \ne 0$.
From $0 < w(x_1) \le w(x_2) \le w(x_3) \le w(x_4)$, we deduce that
$n w(f(Tx)) - 2 w(x_1x_2x_3x_4) \ge n w(f(Tx)) - 2 w(x_jx_kx_ix_n) \ge
w\big((\hess f(Tx))^2_{nj}(\hess f(Tx))^2_{ik}\big) > 0$, which contradicts \eqref{wdegeq}.
```

## 3. The key weight lemma (Theorem \ref{weight}) and its failure for n >= 4

Verbatim statement (lines 296-301):

```tex
\begin{theorem} \label{weight}
Let $K$ be a field of characteristic zero and assume that $f \in K[x] = K[x_1,x_2,\ldots,x_n]$
satisfies $\det \hess f \ne 0$. If $n \le 3$, then there exist a $T \in \GL_n(K)$ and
a weight function $0 < w(x_1) \le w(x_2) \le \cdots \le w(x_n)$, such that the Hessian determinant
of the $w$-leading part of $f(Tx)$ is nonzero.
\end{theorem}
```

Abstract already flags the failure (lines 70-76): "This result does not hold for larger $n$
either (even if we replace `positive' by `nontrivial' above)."

Counterexample for every n >= 4 (Example \ref{weightcounter}, lines 306-329, verbatim):

```tex
\begin{example} \label{weightcounter}
Let $n \ge 4$ and
$$
f := x_1 x_2 + t x_1 x_2^2 + (x_2 + x_1 x_3)^3 + x_1^4 (1 + x_4) +
     (x_5^7 + \cdots + x_n^{n+2})
$$
Then $\det \hess f =  t g$, where
$$
g = -\tfrac1{450} (n+1)!\, (n+2)!\, x_1^9(x_2 + x_1 x_3) x_5^5 \cdots x_n^n
\in \Z[x] \setminus \{0\}
$$
In particular,
$$
\det \hess (f|_{t=0}) = 0 \ne g = \det \hess (f|_{t=1})
$$
Consequently, for each
$T \in \GL_n(\C)$ and each $(w(x_1),w(x_2),\ldots,w(x_n)) \in \R^n \setminus \{0\}^n$, the
$w$-leading part of $f(Tx)|_{t=0}$ has Hessian determinant zero. We will show in section \ref{second}
that the same holds for $f(Tx)|_{t=1}$, although its Hessian determinant is nonzero.
Hence the condition $n \le 3$ in theorem \ref{weight} is necessary.

This example was inspired by formula (9) in \cite[Th.\@ 3.5]{MR2095579}, and $f|_{t=0}$ is of the form
of this formula in dimension $n = 4$.
```

### Newton-polytope obstruction remark for n >= 4 (introduction, lines 130-135, verbatim)

```tex
Example \ref{weightcounter} makes clear that theorem \ref{weight} is no longer true in dimensions
larger than three. In fact, in dimension four and up, it is possible that $f$ has the following
property, in such a way that it cannot be undone by applying linear transformations: the
property of $f$ that $\det \hess \bar{f} \ne 0$ is entirely encapsulated by the Newton polytope of $f$.
More precisely, the part of $f$ consisting of monomials whose supports lie on the boundary of
the Newton polytope of $f$ is a polynomial whose Hessian determinant is zero.
```

## 4. Input from the singular-Hessian classification (Gordan-Noether line)

Theorem \ref{zerohess} (lines 335-350, verbatim), imported from [dBvdE2] = de Bondt &
van den Essen, "Singular Hessians", J. Algebra 282(1):195-204, 2004 (MR2095579), whose
techniques go back to Gordan-Noether 1876 (MR1509898):

```tex
\begin{theorem} \label{zerohess}
Let $K$ be a field of characteristic zero. Suppose that $h \in K[x] = K[x_1,x_2,\ldots,x_n]$
has no terms of degree less than two and that $\det \hess h = 0$. If $n \le 3$, then there exists a
$T \in \GL_n(K)$ such that all entries below the anti-diagonal of the Hessian of $h(Tx)$ are zero.
\begin{enumerate}

\item[i)] If $n = 2$, then $h \in K[l_1]$ for some linear form $l_1 \in K[x]$.

\item[ii)] If $n = 3$, then either $h \in K[l_1,l_2]$ for some linear forms
$l_1, l_2 \in K[x]$, or $h = g_1(l_1) x_1 + g_2(l_1) x_2 + g_3(l_1)x_3$ for some
linear form $l_1 \in K[x]$ and polynomials $g_1, g_2, g_3 \in K[x_1]$.
Furthermore, the leading homogeneous part of $h$ is of the form $l_1^{\deg h - 1} l_4$ for
some linear form $l_4 \in K[x]$ in the latter case.

\end{enumerate}
\end{theorem}
```

Gordan-Noether context in the introduction (lines 155-163, verbatim):

```tex
But the converse may not be true. However, in \cite{MR1509898},
the authors show that $h$ can indeed be written as a polynomial in $n-1$ (or less) linear forms
over $\C$ in case $n \le 4$ and $h$ has Hessian determinant zero,
and give counterexamples for all $n \ge 5$ and all $d \ge 3$.
In \cite{MR2095579}, A. van den Essen and the author classify all (not necessarily homogeneous)
polynomials $h \in K[x_1,x_2,\ldots,x_n]$ with $n \le 3$, such that the Hessian determinant of $h$ is zero,
where $K$ is a field of characteristic zero, using techniques of \cite{MR1509898}.
We shall use these results to prove our main lemma (theorem \ref{weight}) and the
case $\det \hess f = 0$ of our main theorem (theorem \ref{dillenext}).
```

## 5. Method of proof of the main theorem (outline)

Proof of Theorem \ref{dillenext} (lines 380-427):
1. Case det Hess f = 0: anti-triangular T exists directly by Theorem \ref{zerohess}
   (singular-Hessian classification, Gordan-Noether techniques).
2. Case det Hess f in K^*: apply Theorem \ref{weight} to get T in GL_n(K) and a positive
   weight w with det Hess of the w-leading part of f(Tx) nonzero, hence a nonzero CONSTANT;
   homogeneity of the determinant in weights forces (line 395-397, eq. \eqref{wdegeq})
   `n w(f(Tx)) - 2w(x_1 x_2 \cdots x_n) = 0`. Weight bookkeeping (lines 399-410) then kills
   entry (n,n) (and (3,2) when n = 3) of Hess f(Tx), i.e. everything below the anti-diagonal.
3. Invertibility (with polynomial inverse) from anti-triangularity: Lemma \ref{antitri}.
4. Theorem \ref{weight} itself (Section 3) is proved by induction-like case analysis on the
   leading homogeneous part fbar: if det Hess fbar != 0, done with trivial weights; otherwise
   Theorem \ref{zerohess} constrains fbar (in K[l_1,l_2] or the l_1^{d-1} l_4 shape), and one
   perturbs weights one variable at a time (Lemmas \ref{weight2}, \ref{weight1}, \ref{xkdiv} -
   "increase w(x_3) as much as possible while a shared term survives" arguments over the
   Newton polytope). No degree bounds are used; the machinery is weight functions/Newton
   polytopes plus the n <= 3 singular-Hessian classification.

## 6. What obstructs dimension 4 (explicit statements in the paper)

1. The normal-form theorem itself fails at n = 4: the explicit counterexample \eqref{dillen4},
   `f = (x_1 + x_2^2) x_3 + (x_2 + (x_1 + x_2^2)^2) x_4`, built by pulling the planar Keller map
   `F = (x_1 + x_2^2, x_2 + (x_1 + x_2^2)^2)` through Meng's dimension-doubling reduction
   [Men, Prop. 1.4]. Obstruction mechanism: the cubic part x_2^2 x_3 + x_1^2 x_4 has a Hessian
   with linearly independent rows of linear forms, a property invariant under GL_4, incompatible
   with anti-triangularity + constant nonzero Hessian determinant. (Lines 242-256.)
2. The weight lemma fails at n >= 4 in a strong sense: Example \ref{weightcounter} - the ONLY
   weight vector making the leading part's Hessian determinant nonzero is w = 0 (its proof,
   lines 728-786, shows w(x_i) = 0 for all i), so no nontrivial (let alone positive) weight
   works for any T in GL_n(C).
3. Newton-polytope encapsulation (lines 130-135, quoted above): in dimension >= 4 the
   nonvanishing of det Hess can live strictly inside the Newton polytope - the boundary part
   has Hessian determinant zero - and no linear transformation undoes this.
4. Stakes remark (lines 186-189, verbatim): "We shall show that the Jacobian conjecture holds
   for gradient maps in dimension three. In dimension two, this problem has already been solved
   in 1991 by F. Dillen, see below. Notice that an affirmative answer to the same problem in
   dimension four would imply the planar Jacobian conjecture."
   [Verbatim TeX: `Notice that an affirmative answer to the same problem in dimension four
   would imply the planar Jacobian conjecture.`] - i.e. HC_4 => JC_2, consistent with the
   project's framing.
5. Open problem the paper leaves for n > 3 (lines 864-866, verbatim):

```tex
\begin{problem}
Does theorem \ref{definite} also hold for all fields $K$ when $n > 3$.
\end{problem}
```

where Theorem \ref{definite} (lines 812-816, verbatim) is:

```tex
\begin{theorem} \label{definite}
Let $K$ be a field of characteristic zero and $f \in K[x] = K[x_1,x_2,\allowbreak \ldots,x_n]$
such that $\det \hess f \in K^{*}$ and $f$ is anisotropic over $K$ at $\lambda$ for
some $\lambda \in K^n$. If $n \le 3$ or $K = \R$, then $\deg f = 2$.
\end{theorem}
```

(For K = R this rests on Pogorelov [Pog, MR0319126]; corollary, lines 850-854: real gradient
Keller maps with identity linear part are translations.)

## 7. Bibliography key resolution (from consthess3.bbl)

- MR2095579 [dBvdE2]: de Bondt & van den Essen, "Singular Hessians", J. Algebra 282(1):195-204, 2004.
- MR1509898 [GN]: Gordan & Noether, "Ueber die algebraischen Formen, deren Hesse'sche Determinante identisch verschwindet", Math. Ann. 10(4):547-568, 1876.
- MR2221506 [Men]: Meng, "Legendre transform, Hessian conjecture and tree formula", Appl. Math. Lett. 19(6):503-510, 2006.
- MR2138860 [dBvdE4]: de Bondt & van den Essen, "A reduction of the Jacobian conjecture to the symmetric case", Proc. AMS 133(8):2201-2205, 2005.
- MR3179983 [dB]: de Bondt, "Symmetric Jacobians", Cent. Eur. J. Math. 12(6):787-800, 2014.
- MR1107649 [Dil]: Dillen, "Polynomials with constant Hessian determinant", J. Pure Appl. Algebra 71(1):13-18, 1991.
- MR2110519 [dBvdE3]: de Bondt & van den Essen, "Nilpotent symmetric Jacobian matrices and the Jacobian conjecture II", J. Pure Appl. Algebra 196(2-3):135-148, 2005.
- MR0319126 [Pog]: Pogorelov, "On the improper convex affine hyperspheres", Geom. Dedicata 1(1):33-46, 1972.
