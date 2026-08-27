# Extracts from Gao, "Counterexamples to the Jacobian conjecture in dimensions greater than two" (arXiv:2608.00222)

Source file: `C:\Users\mouni\hc4\lit\gao\Jacobian_CE.tex` (2267 lines, dated July 31, 2026).
All LaTeX below is verbatim from the source (line references given). Transcription note: line numbers refer to the .tex file; nothing has been altered inside code fences.

Author/date block (lines 24–31): title "Counterexamples to the Jacobian conjecture in dimensions greater than two", author Shuhong Gao (Clemson), date July 31, 2026. AI disclosure footnote: "The main idea and framework are due to the author, and Claude Fable 5 assisted in the proofs and in the writing up of the paper."

Global conventions (lines 59–71): F(x) is a COLUMN vector; the Jacobian matrix is `J(F) = (∂_1 F, ..., ∂_n F)` — columns are partial derivatives, variables in displayed order. "Geometric degree" = number of points in the generic fiber (lines 105–108).

The six maps of the paper: F (Alpöge's, n=3, deg 3 generic fiber), G (n=3, geometric degree 4), F_4 (n=4, degree 5), F_5 (n=4, degree 10), F_6 (n=5, degree 6), F_7 (n=5, degree 12).

---

## 1. The two DIMENSION-FOUR counterexamples F_4 and F_5

### 1.1 F_4 (subsec:specI, lines 838–1045)

Direction field and syzygy matrix (lines 841–858):

```latex
The column $\Delta=(1,w_1,w_1^2)^{\mathsf T}$ is unimodular, with the natural
syzygy matrix of \S\ref{subsec:forms} (rows
$-\delta_{i+1}e_1^{\mathsf T}+e_{i+1}^{\mathsf T}$):
\[
N=\begin{pmatrix}-w_1&1&0\\[1pt] -w_1^2&0&1\end{pmatrix},
\qquad N\Delta=0,\qquad
\det(\Delta,V,W)=\det\bigl(N\,(V,W)\bigr)\ \ \forall V,W .
\]
```

Ruled-surface ansatz and branch relations (lines 859–900): for `X=A(w_1)+w_2 B(w_1)`, expanding `B=b_1\Delta+b_2\Delta'+b_3\Delta''`, `A'=a_1\Delta+a_2\Delta'+a_3\Delta''`, with `c_2=a_2+w_2(b_1+b_2')`, `c_3=a_3+w_2(b_2+b_3')`:

```latex
P=\begin{pmatrix}c_2 & b_2\\ 2(w_1c_2+c_3) & 2(w_1b_2+b_3)\end{pmatrix},
\qquad
\widetilde P=N\,J(\Delta)=\begin{pmatrix}1&0\\ 2w_1&0\end{pmatrix},
\]
and therefore
\[
\det\bigl(P+\gamma\widetilde P\bigr)\;=\;2\,(c_2b_3-c_3b_2)\;+\;2\,\gamma\,b_3 .
```

The M-coefficient vanishes automatically (L-branch); the normalization L≡2 and tangency reduce to (line 884):

```latex
b_3=1,\qquad b_1=b_2^2-b_2',\qquad a_2=a_3b_2 .
```

Chosen data (lines 902–917): `b_2=w_1`, `a_3=1` (so `a_2=w_1`, `b_1=w_1^2-1`), `a_1=\tfrac13-\tfrac83w_1`:

```latex
A(w_1)=\Bigl(\tfrac13w_1-\tfrac43w_1^2,\;\;
\tfrac23w_1^2-\tfrac89w_1^3,\;\;
2w_1+\tfrac79w_1^3-\tfrac23w_1^4\Bigr)^{\mathsf T},\qquad
B(w_1)=\bigl(w_1^2-1,\;w_1^3,\;w_1^4+w_1^2+2\bigr)^{\mathsf T},
\]
so the surface $X=A+w_2B=(p,q,r)^{\mathsf T}$ and the tangent field $\Delta$ are, explicitly,
\begin{align*}
p(w_1,w_2)&=\tfrac13w_1-\tfrac43w_1^2+w_2\,(w_1^2-1),\\
q(w_1,w_2)&=\tfrac23w_1^2-\tfrac89w_1^3+w_2\,w_1^3,\\
r(w_1,w_2)&=2w_1+\tfrac79w_1^3-\tfrac23w_1^4+w_2\,(w_1^4+w_1^2+2),\\
\Delta(w_1,w_2)&=(1,\;w_1,\;w_1^2)^{\mathsf T}
\qquad(\text{independent of }w_2).
\end{align*}
```

Definition of F_4 (lines 921–931):

```latex
Define on $\C^4\ni(x,y,z,t)$:
\[
\gamma=1+xy+x^2t,\quad u_1=1+xz,\quad
u_2=-\tfrac79+\tfrac{22}9xy+\tfrac83xz,\quad
w_1=\gamma u_1,\quad w_2=\gamma u_2,\quad C=\gamma x,
\]
\[
S_i=A_i(w_1)+w_2B_i(w_1)+\gamma\,\Delta_i(w_1)\quad(i=1,2,3),
\]
\[
\boxed{\;F_4\;=\;\Bigl(C,\;\;\frac{S_1}{C},\;\;\frac{S_2}{C^2},\;\;\frac{S_3}{C}\Bigr)^{\mathsf T}.\;}
\]
```

**Theorem thm:F4 (lines 951–955), verbatim:**

```latex
\begin{theorem}\label{thm:F4}
$F_4$ is a polynomial map $\C^4\to\C^4$ with components of degrees $4,11,12,21$;
its Jacobian determinant is identically $-\tfrac{44}{9}$; and its generic fiber
consists of exactly $5$ points.
\end{theorem}
```

Chain factorization of the Jacobian (lines 1030–1039, verbatim excerpt): "The proof of the Jacobian statement is the verified chain factorization `x^4\cdot(-\tfrac{22}9)\cdot\gamma^2\cdot2\gamma^2\cdot(\gamma x)^{-4} =-\tfrac{44}9`, in which the factor `2\gamma^2` is Lemma lem:JL for the surface above (L=2, M=0; a polynomial identity), the factor `\gamma^2` comes from the scaling `(w_1,w_2)=(\gamma u_1,\gamma u_2)`, `x^4` from the monomial substitution `(v_1,v_2,v_3)=(xy,xz,x^2t)` with x-degrees (1,1,2), and `-\tfrac{22}9` is the determinant of the affine stage. The five side conditions (divisibility of S_1,S_3 by C and of S_2 by C^2) determine the constants `(\gamma_0,u_1^0,u_2^0,\alpha_0,\alpha_1)=(1,1,-\tfrac79,\tfrac13,-\tfrac83)`."

Tangency quintic controlling the fibers of F_4 (from the proof sketch, lines 976–995): with `(Y_1,Y_2,Y_3)=(v_1v_2, v_1^2v_3, v_1v_4)`, `\gamma=Y_1-A_1-w_2B_1`:

```latex
R(w_1;Y)\;=\;\det\bigl(A(w_1)-Y,\;B(w_1),\;\Delta(w_1)\bigr)
\;=\;\tfrac29w_1^5+\tfrac29w_1^4+\bigl(\tfrac89+Y_1\bigr)w_1^3
-\bigl(\tfrac43+2Y_2\bigr)w_1^2+(2Y_1+Y_3)\,w_1-2Y_2 ,
```

a quintic with constant leading coefficient 2/9. (Line 1041–1045: the raw tangency equation of the swept ruled family has degree 7, always divisible by the spurious factor `w_1^2-1`, leaving the squarefree quintic; witness target `v=(\tfrac23,-\tfrac15,\tfrac12,3)`; numerical target (0.7,-0.4,0.3,0.5) has five preimages, three real and one conjugate pair.)

**Full expanded polynomials of F_4 (lines 998–1027), verbatim:**

```latex
F_{4,1} ={}& x^3t+x^2y+x\\[2pt]
F_{4,2} ={}&\tfrac{22}{9}x^6yz^2t^2+\tfrac{8}{3}x^6z^3t^2+\tfrac{44}{9}x^5y^2z^2t+\tfrac{16}{3}x^5yz^3t+\tfrac{44}{9}x^5yzt^2+\tfrac{41}{9}x^5z^2t^2+\tfrac{22}{9}x^4y^3z^2\\
&+\tfrac{8}{3}x^4y^2z^3+\tfrac{88}{9}x^4y^2zt+14x^4yz^2t+\tfrac{16}{3}x^4z^3t+\tfrac{22}{9}x^4yt^2+\tfrac{10}{9}x^4zt^2+\tfrac{44}{9}x^3y^3z+\tfrac{85}{9}x^3y^2z^2\\
&+\tfrac{16}{3}x^3yz^3+\tfrac{44}{9}x^3y^2t+12x^3yzt+\tfrac{70}{9}x^3z^2t-\tfrac{7}{9}x^3t^2+\tfrac{22}{9}x^2y^3+\tfrac{98}{9}x^2y^2z+\tfrac{92}{9}x^2yz^2+\tfrac{8}{3}x^2z^3\\
&+\tfrac{10}{3}x^2yt-\tfrac{4}{9}x^2zt+\tfrac{37}{9}xy^2+\tfrac{40}{9}xyz+\tfrac{29}{9}xz^2-\tfrac{26}{9}xt-\tfrac{26}{9}y-\tfrac{35}{9}z\\[2pt]
F_{4,3} ={}& \tfrac{22}{9}x^6yz^3t^2+\tfrac{8}{3}x^6z^4t^2+\tfrac{44}{9}x^5y^2z^3t+\tfrac{16}{3}x^5yz^4t+\tfrac{22}{3}x^5yz^2t^2+\tfrac{65}{9}x^5z^3t^2+\tfrac{22}{9}x^4y^3z^3\\
&+\tfrac{8}{3}x^4y^2z^4+\tfrac{44}{3}x^4y^2z^2t+\tfrac{58}{3}x^4yz^3t+\tfrac{16}{3}x^4z^4t+\tfrac{22}{3}x^4yzt^2+\tfrac{17}{3}x^4z^2t^2+\tfrac{22}{3}x^3y^3z^2\\
&+\tfrac{109}{9}x^3y^2z^3+\tfrac{16}{3}x^3yz^4+\tfrac{44}{3}x^3y^2zt+26x^3yz^2t+\tfrac{122}{9}x^3z^3t+\tfrac{22}{9}x^3yt^2+\tfrac{1}{3}x^3zt^2+\tfrac{22}{3}x^2y^3z\\
&+\tfrac{61}{3}x^2y^2z^2+16x^2yz^3+\tfrac{8}{3}x^2z^4+\tfrac{44}{9}x^2y^2t+\tfrac{46}{3}x^2yzt+\tfrac{26}{3}x^2z^2t-\tfrac{7}{9}x^2t^2+\tfrac{22}{9}xy^3+15xy^2z\\
&+16xyz^2+\tfrac{19}{3}xz^3+\tfrac{10}{3}xyt-2xzt+\tfrac{37}{9}y^2+\tfrac{16}{3}yz+\tfrac{11}{3}z^2-\tfrac{22}{9}t\\[2pt]
F_{4,4} ={}& \tfrac{22}{9}x^{12}yz^4t^4+\tfrac{8}{3}x^{12}z^5t^4+\tfrac{88}{9}x^{11}y^2z^4t^3+\tfrac{32}{3}x^{11}yz^5t^3+\tfrac{88}{9}x^{11}yz^3t^4+\tfrac{89}{9}x^{11}z^4t^4\\
&+\tfrac{44}{3}x^{10}y^3z^4t^2+16x^{10}y^2z^5t^2+\tfrac{352}{9}x^{10}y^2z^3t^3+\tfrac{148}{3}x^{10}yz^4t^3+\tfrac{32}{3}x^{10}z^5t^3+\tfrac{88}{9}x^9y^4z^4t\\
&+\tfrac{32}{3}x^9y^3z^5t+\tfrac{44}{3}x^{10}yz^2t^4+\tfrac{116}{9}x^{10}z^3t^4+\tfrac{176}{3}x^9y^3z^3t^2+\tfrac{266}{3}x^9y^2z^4t^2+32x^9yz^5t^2\\
&+\tfrac{22}{9}x^8y^5z^4+\tfrac{8}{3}x^8y^4z^5+\tfrac{176}{3}x^9y^2z^2t^3+\tfrac{272}{3}x^9yz^3t^3+\tfrac{350}{9}x^9z^4t^3+\tfrac{352}{9}x^8y^4z^3t+\tfrac{620}{9}x^8y^3z^4t\\
&+32x^8y^2z^5t+\tfrac{88}{9}x^9yzt^4+6x^9z^2t^4+88x^8y^3z^2t^2+\tfrac{584}{3}x^8y^2z^3t^2+\tfrac{394}{3}x^8yz^4t^2+16x^8z^5t^2\\
&+\tfrac{88}{9}x^7y^5z^3+\tfrac{59}{3}x^7y^4z^4+\tfrac{32}{3}x^7y^3z^5+\tfrac{352}{9}x^8y^2zt^3+\tfrac{248}{3}x^8yz^2t^3+\tfrac{440}{9}x^8z^3t^3+\tfrac{176}{3}x^7y^4z^2t\\
&+\tfrac{1520}{9}x^7y^3z^3t+146x^7y^2z^4t+32x^7yz^5t+\tfrac{22}{9}x^8yt^4-\tfrac{4}{9}x^8zt^4+\tfrac{176}{3}x^7y^3zt^2+212x^7y^2z^2t^2\\
&+\tfrac{616}{3}x^7yz^3t^2+\tfrac{172}{3}x^7z^4t^2+\tfrac{44}{3}x^6y^5z^2+52x^6y^4z^3+\tfrac{482}{9}x^6y^3z^4+16x^6y^2z^5+\tfrac{88}{9}x^7y^2t^3\\
&+\tfrac{112}{3}x^7yzt^3+20x^7z^2t^3+\tfrac{352}{9}x^6y^4zt+200x^6y^3z^2t+264x^6y^2z^3t+\tfrac{1120}{9}x^6yz^4t+\tfrac{32}{3}x^6z^5t\\
&-\tfrac{7}{9}x^7t^4+\tfrac{44}{3}x^6y^3t^2+\tfrac{344}{3}x^6y^2zt^2+\tfrac{1354}{9}x^6yz^2t^2+\tfrac{655}{9}x^6z^3t^2+\tfrac{88}{9}x^5y^5z+\tfrac{194}{3}x^5y^4z^2\\
&+\tfrac{968}{9}x^5y^3z^3+\tfrac{604}{9}x^5y^2z^4+\tfrac{32}{3}x^5yz^5+\tfrac{20}{3}x^6yt^3-\tfrac{40}{9}x^6zt^3+\tfrac{88}{9}x^5y^4t+\tfrac{1040}{9}x^5y^3zt\\
&+\tfrac{2168}{9}x^5y^2z^2t+\tfrac{554}{3}x^5yz^3t+\tfrac{338}{9}x^5z^4t+\tfrac{74}{3}x^5y^2t^2+\tfrac{452}{9}x^5yzt^2+\tfrac{287}{9}x^5z^2t^2+\tfrac{22}{9}x^4y^5\\
&+\tfrac{116}{3}x^4y^4z+\tfrac{994}{9}x^4y^3z^2+\tfrac{1007}{9}x^4y^2z^3+40x^4yz^4+\tfrac{8}{3}x^4z^5-\tfrac{34}{9}x^5t^3+\tfrac{236}{9}x^4y^3t+\tfrac{1024}{9}x^4y^2zt\\
&+\tfrac{382}{3}x^4yz^2t+\tfrac{454}{9}x^4z^3t+\tfrac{52}{9}x^4yt^2-\tfrac{47}{9}x^4zt^2+9x^3y^4+\tfrac{532}{9}x^3y^3z+\tfrac{859}{9}x^3y^2z^2+\tfrac{542}{9}x^3yz^3\\
&+\tfrac{83}{9}x^3z^4+\tfrac{206}{9}x^3y^2t+\tfrac{346}{9}x^3yzt+\tfrac{250}{9}x^3z^2t-\tfrac{17}{3}x^3t^2+\tfrac{40}{3}x^2y^3+\tfrac{131}{3}x^2y^2z+\tfrac{404}{9}x^2yz^2\\
&+\tfrac{41}{3}x^2z^3+\tfrac{10}{3}x^2yt+\tfrac{10}{9}x^2zt+9xy^2+\tfrac{142}{9}xyz+\tfrac{89}{9}xz^2-\tfrac{28}{9}xt+\tfrac{20}{3}y+\tfrac{29}{3}z
```

### 1.2 The rigidity / non-equivalence statement after thm:F4 (lines 957–966), verbatim:

```latex
Consequently $F_4$ is a non-injective Keller map
in dimension four which is not equivalent, under composition with polynomial
automorphisms on either side, to $\Phi\times\mathrm{id}$ for any
three-dimensional Keller map $\Phi$ of geometric degree $\neq5$: geometric degree
is invariant under such compositions, and $\deg(\Phi\times\mathrm{id})
=\deg\Phi$. In particular $F_4$ is not a trivial extension of Alp\"oge's example
or of the map $G$. Whether $F_4$ is equivalent to $\Phi\times\mathrm{id}$ for a
degree-\emph{five} member $\Phi$ of Gallagher's family \cite{Gallagher2026} is
left open; the fiber stratification of a product has product form, so the
stratification deferred to Appendix~\ref{app:F4} is expected to decide this.
```

(Analogous statement for F_5, lines 1218–1221: "Consequently $F_5$ is a four-dimensional non-injective Keller map arising from the $M$-branch, of geometric degree ten, not equivalent to $\Phi\times\mathrm{id}$ for any three-dimensional Keller map $\Phi$ of geometric degree $\neq10$.")

### 1.3 F_5 (subsec:specII, lines 1047–1284)

Direction field `\Delta=(1,w_1,w_2)^{\mathsf T}` (unimodular), syzygy matrix (lines 1050–1055):

```latex
N=\begin{pmatrix}-w_1&1&0\\ -w_2&0&1\end{pmatrix},
\qquad N\Delta=0 ,
```

Here `\widetilde P=N\,J(\Delta)=I_2` and (lines 1072–1079):

```latex
P=N\,J(X)=\begin{pmatrix}a_2 & b_2\\ a_3 & b_3\end{pmatrix},
\qquad
\det\bigl(P+\gamma\widetilde P\bigr)
\;=\;(a_2b_3-a_3b_2)\;+\;\gamma\,(a_2+b_3)\;+\;\gamma^2 .
\]
Thus $M\equiv1$ identically and $L=a_2+b_3$, the trace of $P$: this $\Delta$
lives on the $M$-branch, with normalization $L\equiv0$ and
$\det J(F_0)=\gamma^3$.
```

Branch solved via rank-one trace-free parametrization `(a_2,b_2,a_3,b_3)=(-\rho\lambda,\lambda,-\rho^2\lambda,\rho\lambda)` and the linear wave-type equation (eq:wave, lines 1110–1112):

```latex
\lambda_{11}+2(\rho\lambda)_{12}+(\rho^2\lambda)_{22}\;=\;0 ,
```

For ρ=w_1 the substitution `\eta=w_2-\tfrac12w_1^2` turns eq:wave into the heat equation; polynomial solutions are spanned by Rosenbloom–Widder heat polynomials `v_k(x,t)=\sum_{j\le k/2}\frac{k!}{j!\,(k-2j)!}x^{k-2j}t^j`, giving one new polynomial solution `v_k(w_1,\tfrac12 w_1^2-w_2)` per degree (lines 1131–1154). Chosen: `\rho=w_1`, `\lambda=\nu_1+\nu_2w_1+\nu_3(w_2-w_1^2)`, and the side conditions at base point `(\gamma_0,u_1^0,u_2^0)=(1,1,1)` give `(\nu_1,\nu_2,\nu_3)=(\tfrac{51}{29},-\tfrac{60}{29},-\tfrac{240}{29})`, with surface (lines 1173–1177):

```latex
X_1&=\tfrac{160}{29}w_1^3-\tfrac{60}{29}w_1^2-\tfrac{240}{29}w_1w_2+\tfrac{51}{29}w_1+\tfrac{60}{29}w_2,\\
X_2&=\tfrac{60}{29}w_1^4-\tfrac{20}{29}w_1^3-\tfrac{120}{29}w_2^2+\tfrac{51}{29}w_2,\\
X_3&=-\tfrac{48}{29}w_1^5+\tfrac{15}{29}w_1^4+\tfrac{240}{29}w_1^3w_2-\tfrac{17}{29}w_1^3-\tfrac{60}{29}w_1^2w_2-\tfrac{240}{29}w_1w_2^2+\tfrac{51}{29}w_1w_2+\tfrac{30}{29}w_2^2 .
```

Unipotent nonlinear stage (lines 1179–1189):

```latex
with $v=(v_1,v_2,v_3)=(xy,x^2z,x^3t)$, set
\[
\gamma=1+\tfrac{20}{29}v_1,\qquad
u_1=1-\tfrac{49}{29}v_1+8v_2,\qquad
u_2=1-\tfrac{69}{29}v_1+9v_2+v_3+\tfrac{4197}{1682}\,v_1^2 ,
\]
where the $v_1$-column lies in $\ker(\nabla E_2)\cap\ker(\nabla E_3)$, the
$v_2$-column in $\ker(\nabla E_3)$, and the coefficient $\tau=\tfrac{4197}{1682}$
kills the $v_1^2$-obstruction while leaving the stage Jacobian constant
($=\tfrac{160}{29}$), because the corresponding cofactor vanishes by design.
```

Definition (lines 1189–1193):

```latex
With $w_1=\gamma u_1$, $w_2=\gamma^2u_2$, $C=\gamma x$ and
$S_i=X_i(w_1,w_2)+\gamma\Delta_i(w_1,w_2)$, define
\[
\boxed{\;F_5\;=\;\Bigl(C,\;\;\frac{S_1}{C},\;\;\frac{S_2}{C^2},\;\;\frac{S_3}{C^3}\Bigr)^{\mathsf T}.\;}
\]
```

**Theorem thm:F5 (lines 1212–1216), verbatim:**

```latex
\begin{theorem}\label{thm:F5}
$F_5$ is a polynomial map $\C^4\to\C^4$ with components of degrees $3,12,14,16$;
its Jacobian determinant is identically $\tfrac{160}{29}$; and its generic fiber
consists of exactly $10$ points.
\end{theorem}
```

Fiber-degree-10 resultant (proof sketch, lines 1227–1244): with `(Y_1,Y_2,Y_3)=(v_1v_2, v_1^2v_3, v_1^3v_4)`, `\gamma=Y_1-X_1(w_1,w_2)`, resultant w.r.t. w_2 of `X_2+(Y_1-X_1)w_1=Y_2`, `X_3+(Y_1-X_1)w_2=Y_3`, after rescaling by `-29^4/30`:

```latex
R(w_1;Y)={}&20480w_1^{10}-51200w_1^{9}+57600w_1^{8}+\bigl(222720Y_1-126720\bigr)w_1^{7}\\
&+\bigl(-204160Y_1-742400Y_2+234480\bigr)w_1^{6}
+\bigl(-119712Y_1+890880Y_2+890880Y_3-155448\bigr)w_1^{5}\\
&+\bigl(201840Y_1^2+59160Y_1-417600Y_2-1113600Y_3+65025\bigr)w_1^{4}\\
&+\bigl(-67280Y_1^2-538240Y_1Y_2+41412Y_1+306240Y_2+556800Y_3-44217\bigr)w_1^{3}\\
&+\bigl(-146334Y_1^2+201840Y_1Y_2+807360Y_1Y_3+75429Y_1-88740Y_2-459360Y_3\bigr)w_1^{2}\\
&+\bigl(97556Y_1^3-42891Y_1^2-403680Y_1Y_3+177480Y_3\bigr)w_1\\
&-97556Y_1^2Y_2+42891Y_1Y_2+171564Y_1Y_3-25230Y_2^2+201840Y_2Y_3-403680Y_3^2-75429Y_3 ,
```

with constant leading coefficient. Witness target `v=(\tfrac12,\tfrac23,-\tfrac13,1)`: R squarefree, coprime to the γ=0 resultant and to the first subresultant, so exactly 10 fiber points.

**Full expanded polynomials of F_5 (lines 1262–1284), verbatim:**

```latex
F_{5,1} ={}& \tfrac{20}{29}x^2y+x\\[2pt]
F_{5,2} ={}& \tfrac{32768000}{24389}x^7y^2z^3-\tfrac{602112000}{707281}x^6y^3z^2-\tfrac{768000}{24389}x^6y^2zt+\tfrac{3276800}{841}x^6yz^3+\tfrac{2076288000}{20511149}x^5y^4z\\
&+\tfrac{4704000}{707281}x^5y^3t-\tfrac{54835200}{24389}x^5y^2z^2+\tfrac{80752000}{20511149}x^4y^5-\tfrac{76800}{841}x^5yzt+\tfrac{81920}{29}x^5z^3+\tfrac{152428800}{707281}x^4y^3z\\
&+\tfrac{374400}{24389}x^4y^2t-\tfrac{1044480}{841}x^4yz^2+\tfrac{169140800}{20511149}x^3y^4-\tfrac{1920}{29}x^4zt+\tfrac{515520}{24389}x^3y^2z+\tfrac{3360}{841}x^3yt+\tfrac{9600}{29}x^3z^2\\
&+\tfrac{907520}{707281}x^2y^3-\tfrac{108960}{841}x^2yz-\tfrac{180}{29}x^2t+\tfrac{14050}{24389}xy^2-\tfrac{252}{29}xz+y\\[2pt]
F_{5,3} ={}& \tfrac{98304000}{24389}x^8y^2z^4-\tfrac{2408448000}{707281}x^7y^3z^3+\tfrac{9830400}{841}x^7yz^4+\tfrac{22127616000}{20511149}x^6y^4z^2-\tfrac{191692800}{24389}x^6y^2z^3\\
&-\tfrac{90354432000}{594823321}x^5y^5z-\tfrac{48000}{24389}x^6y^2t^2+\tfrac{245760}{29}x^6z^4-\tfrac{201456000}{20511149}x^5y^4t+\tfrac{45158400}{24389}x^5y^3z^2\\
&-\tfrac{73022484000}{17249876309}x^4y^6-\tfrac{864000}{24389}x^5y^2zt-\tfrac{1310720}{841}x^5yz^3-\tfrac{5316643200}{20511149}x^4y^4z-\tfrac{4800}{841}x^5yt^2-\tfrac{13521600}{707281}x^4y^3t\\
&-\tfrac{25906560}{24389}x^4y^2z^2-\tfrac{4696088400}{594823321}x^3y^5-\tfrac{86400}{841}x^4yzt+\tfrac{112640}{29}x^4z^3+\tfrac{69664320}{707281}x^3y^3z-\tfrac{120}{29}x^4t^2\\
&+\tfrac{62760}{24389}x^3y^2t-\tfrac{1613760}{841}x^3yz^2+\tfrac{40879390}{20511149}x^2y^4-\tfrac{2160}{29}x^3zt+\tfrac{3661320}{24389}x^2y^2z+\tfrac{240}{29}x^2yt+\tfrac{9480}{29}x^2z^2\\
&-\tfrac{516820}{707281}xy^3-\tfrac{105360}{841}xyz-\tfrac{189}{29}xt-\tfrac{1437}{1682}y^2-z\\[2pt]
F_{5,4} ={}& -\tfrac{629145600}{24389}x^9y^2z^5+\tfrac{19267584000}{707281}x^8y^3z^4+\tfrac{49152000}{24389}x^8y^2z^3t-\tfrac{62914560}{841}x^8yz^5-\tfrac{132882432000}{20511149}x^7y^4z^3\\
&-\tfrac{903168000}{707281}x^7y^3z^2t+\tfrac{1975910400}{24389}x^7y^2z^4-\tfrac{15504384000}{20511149}x^6y^5z^2-\tfrac{768000}{24389}x^7y^2zt^2+\tfrac{4915200}{841}x^7yz^3t\\
&-\tfrac{1572864}{29}x^7z^5+\tfrac{2308608000}{20511149}x^6y^4zt-\tfrac{15174451200}{707281}x^6y^3z^3+\tfrac{3799290048000}{17249876309}x^5y^6z+\tfrac{4704000}{707281}x^6y^3t^2\\
&-\tfrac{85708800}{24389}x^6y^2z^2t+\tfrac{54312960}{841}x^6yz^4+\tfrac{8448384000}{594823321}x^5y^5t-\tfrac{11697561600}{20511149}x^5y^4z^2+\tfrac{2437443220800}{500246412961}x^4y^7\\
&-\tfrac{76800}{841}x^6yzt^2+\tfrac{122880}{29}x^6z^3t+\tfrac{195724800}{707281}x^5y^3zt-\tfrac{496404480}{24389}x^5y^2z^3+\tfrac{184263724800}{594823321}x^4y^5z\\
&+\tfrac{374400}{24389}x^5y^2t^2-\tfrac{1873920}{841}x^5yz^2t+\tfrac{184320}{29}x^5z^4+\tfrac{484262400}{20511149}x^4y^4t+\tfrac{1251989760}{707281}x^4y^3z^2\\
&+\tfrac{135717082080}{17249876309}x^3y^6-\tfrac{1920}{29}x^5zt^2+\tfrac{2238720}{24389}x^4y^2zt-\tfrac{844800}{841}x^4yz^3-\tfrac{4465704480}{20511149}x^3y^4z+\tfrac{2760}{841}x^4yt^2\\
&+\tfrac{7680}{29}x^4z^2t-\tfrac{6755640}{707281}x^3y^3t-\tfrac{9897600}{24389}x^3y^2z^2-\tfrac{2476692498}{594823321}x^2y^5-\tfrac{145200}{841}x^3yzt+\tfrac{123776}{29}x^3z^3\\
&-\tfrac{61829880}{707281}x^2y^3z-\tfrac{210}{29}x^3t^2-\tfrac{126390}{24389}x^2y^2t-\tfrac{1938744}{841}x^2yz^2+\tfrac{128225745}{41022298}xy^4-\tfrac{2412}{29}x^2zt\\
&+\tfrac{5740350}{24389}xy^2z+\tfrac{11001}{841}xyt+\tfrac{9318}{29}xz^2+\tfrac{1629343}{1414562}y^3-\tfrac{99879}{841}yz-\tfrac{160}{29}t
```

---

## 2. The general tangent-sweep construction

### 2.1 The plane sweep and monomial twist (sec:sweep, subsec:twist, lines 200–374)

Plane sweep: `\mathcal K = (p(w), q(w))^T` with `\deg q \le \deg p + 1`. `T_0(\gamma,w)=(p(w)+\gamma p'(w), q(w)+\gamma q'(w))^T`, `\det J(T_0)=\gamma(p'q''-p''q')`. Normalizing to direction field `\delta=(2,w)^T` (tangent exactly when `q'=\tfrac{w}{2}p'`, i.e. `q(w)=\int_0^w \tfrac{s}{2}p'(s)ds` — eq:normalization) gives the inflection-free normal form `S(\gamma,w)=(p(w)+2\gamma, q(w)+\gamma w)^T` with `\det J(S)=2\gamma`. Tangency polynomial (eq:tangency): `W_{X,Y}(w)=q(w)+\tfrac{w}{2}(X-p(w))-Y`; key derivative identity (eq:Wprime): `W'_{X,Y}(w)=\tfrac12(X-p(w))=\gamma`, `W''_{X,Y}(w)=-\tfrac12 p'(w)`. Lemma lem:fibroots: roots of W_{X,Y} biject with fiber of S; deg = d+1 with constant leading coefficient `-\mathrm{lc}(p)/(2(d+1))`; generically (d+1)-to-one (= the class of K, Plücker `m=\deg(\deg-1)-2\delta-3\kappa`). Prop prop:fibsizes: sweep fiber sizes d+1 / d (smooth point) / d−1 (cusp or node); after the twist: d+1 (generic), d−1 (smooth curve points), d−2 (cusps), d−3 (nodes), plus possible extra twist strata.

**The monomial twist (subsec:twist, lines 336–374), verbatim core:**

```latex
The sweep is not a Keller map ($\det J(S)=2\gamma$), but the factor $\gamma$ can be
cancelled by conjugation with monomial maps. Pad the sweep to three variables and
compose:
\begin{multline*}
\begin{pmatrix}x\\ y\\ z\end{pmatrix}
\;\xrightarrow{\ \text{monomial}\ }\;
\begin{pmatrix}x\\ xy\\ x^2z\end{pmatrix}
\;\xrightarrow{\ \text{affine}\ }\;
\begin{pmatrix}x\\ \gamma\\ u\end{pmatrix}
\;\xrightarrow{\ w=\gamma u\ }\;
\begin{pmatrix}x\\ \gamma\\ w\end{pmatrix}\\
\xrightarrow{\ \text{padded sweep}\ }\;
\begin{pmatrix}\gamma x\\ P\\ Q\end{pmatrix}
\;\xrightarrow{\ \text{twist}\ }\;
\begin{pmatrix}C\\ P/C\\ Q/C^2\end{pmatrix},
\end{multline*}
where $\gamma=\gamma_0+a\,xy+b\,x^2z$, $u=1+xy$, $P=p(w)+2\gamma$,
$Q=q(w)+\gamma w$ and $C=\gamma x$. The Jacobians multiply to
\[
C^{-3}\cdot 2\gamma^2\cdot\gamma\cdot x^3\cdot(-b)\;=\;-2b,
\]
a nonzero constant. ... The twist requires $C\mid P$ and $C^2\mid Q$; these
divisibilities amount to finitely many \emph{linear} conditions on the
coefficients of $p$ (and on $\gamma_0,a$), the ``side conditions'' of
\cite{Speyer2026}; \S\ref{subsec:stage} systematizes these conditions in every
dimension. Because the
ramification locus $\gamma=0$ of the sweep is exactly where the twisted coordinate
$C=\gamma x$ degenerates, the composite is everywhere unramified; the
$d+1$ sheets of the sweep survive, but sheets that would have merged over the
curve now escape to infinity. \emph{The construction converts ramification into
non-properness.}
```

### 2.2 Sweeping a direction field on a hypersurface (subsec:sweephyp, lines 540–600)

**Tangency criterion (eq:tangcrit), verbatim:**

```latex
Let $w=(w_1,\dots,w_{n-2})$, let $X\in\C[w]^{\,n-1}$ be a column vector
of $n-1$ polynomials, viewed as a parametrized hypersurface
$X\colon\C^{n-2}\to\C^{n-1}$ ... A \emph{polynomial tangent field} is a column vector $\Delta\in\C[w]^{\,n-1}$
with $\Delta(w)$ in the column span of $J(X)(w)$ for all $w$, equivalently the
\emph{tangency criterion}
\begin{equation}\label{eq:tangcrit}
\det\bigl(\Delta,\,J(X)\bigr)\;\equiv\;0 ,
\end{equation}
... Under the standard identification
$\Lambda^{n-1}\C^{\,n-1}\cong\C$, the left-hand side is the exterior product
$\Delta\wedge\partial_1X\wedge\dots\wedge\partial_{n-2}X$ ... The padded sweep is
\[
F_0(x,\gamma,w)\;=\;\begin{pmatrix}\gamma x\\ X(w)+\gamma\Delta(w)\end{pmatrix}\colon
\ \C^n\to\C^n .
\]
```

**Lemma lem:JL (lines 564–575), verbatim:**

```latex
\begin{lemma}\label{lem:JL}
Let $\Delta$ be a polynomial tangent field of $X$, i.e.\
$\det\bigl(\Delta,J(X)\bigr)\equiv0$. Then, with the variables ordered
$(x,\gamma,w_1,\dots,w_{n-2})$,
\[
\det J(F_0)\;=\;\gamma\,\det\bigl(\Delta,\;J(X)+\gamma J(\Delta)\bigr)
\;=\;\sum_{k=1}^{n-2}L_k\,\gamma^{\,k+1},
\]
where $L_k$ is the sum of the determinants $\det(\Delta,C_1,\dots,C_{n-2})$ in
which exactly $k$ of the columns $C_i$ are taken from $J(\Delta)$ and the rest
from $J(X)$.
\end{lemma}
```

Explicit branch coefficients: for n=4 (`L:=L_1`, `M:=L_2`): `\det J(F_0)=\gamma^2 L + \gamma^3 M`, `L=\det(\Delta,\partial_1X,\partial_2\Delta)+\det(\Delta,\partial_1\Delta,\partial_2X)`, `M=\det(\Delta,\partial_1\Delta,\partial_2\Delta)`. For n=5: `\det J(F_0)=\gamma^2L_1+\gamma^3L_2+\gamma^4L_3` with L_1, L_2, L_3 the sums with 1, 2, 3 columns from J(Δ) respectively (lines 586–596).

### 2.3 Prescribed direction fields (subsec:forms, lines 602–713)

Roles reversed: PRESCRIBE Δ, solve for X. Assume Δ unimodular; complete to `U\in\mathrm{SL}_{n-1}(\C[w])`; let N = last n−2 rows of U^{-1} (so NΔ=0, rows of N a free basis of syzygies of Δ). Key pairing (eq:pairing): `\det(\Delta,V_1,\dots,V_{n-2})=\det(N(V_1,\dots,V_{n-2}))`. Set `P:=N J(X)`, `\widetilde P:=N J(\Delta)`; then (eq:wedge): `\det J(F_0)=\gamma\det(P+\gamma\widetilde P)` with `\det P=0` (tangency); L_k = coefficient of γ^k in det(P+γP̃).

**The branches (lines 656–674), verbatim core:**

```latex
The normalization conditions now split into $n-2$ \emph{branches}, one for each
index $1\le k\le n-2$:
\[
\text{the \emph{$L_k$-branch}:}\qquad
L_j\equiv0\ (j\neq k),\qquad L_k\equiv\text{const}\neq0,
\qquad\text{giving}\qquad \det J(F_0)=L_k\,\gamma^{\,k+1} .
\]
... since every term of $L_k$ contains $k$ columns
of $\widetilde P$, one has $L_j\equiv0$ automatically for all
$j>\operatorname{rank}\widetilde P$, so the rank of $\widetilde P$ bounds the
highest available branch. At the two extremes, a direction field depending on a
single parameter forces $\operatorname{rank}\widetilde P\le1$ (only the bottom
branch survives), while the tautological field with $\widetilde P$ the identity
opens the top branch $k=n-2$.
```

Potentials (when δ_1=1, syzygy rows `-\delta_{i+1}e_1^T+e_{i+1}^T`): `G_i:=X_i-\delta_i X_1` (i=2..n−1), `\Gamma=(G_2,\dots,G_{n-1})^T`, `\widehat\Delta=(\delta_2,\dots,\delta_{n-1})^T`; then `P=J(\Gamma)+X_1 J(\widehat\Delta)`, `\widetilde P=J(\widehat\Delta)`, and the branch analysis lives in the pencil `P+\gamma\widetilde P=J(\Gamma)+\mu J(\widehat\Delta)`, `\mu=X_1+\gamma` (lines 682–713).

### 2.4 The twist in dimension n: the stage equations (subsec:stage, lines 715–834)

Output shape: `F=(C, S_1/C^{e_1}, \dots, S_{n-1}/C^{e_{n-1}})^T`, `C=\gamma x`. Three layers:

(i) **Discrete data** (eq:discrete, lines 733–739), verbatim:

```latex
Positive integers: weights $d_j$ for the parameters
($\operatorname{wt}w_j=d_j$), twist exponents $e_i$ for the target components,
and $x$-degrees $m_j$ for the source monomials. They must satisfy
\begin{equation}\label{eq:discrete}
\operatorname{wt}(X_i)\ge e_i\ \ (\text{every monomial}),\qquad
1+\operatorname{wt}(\Delta_i)\ge e_i,\qquad
\sum_j m_j=\sum_i e_i,\qquad
\sum_j d_j+(k+1)=\sum_i e_i .
\end{equation}
The first two conditions make
$E_i:=S_i\bigl(\gamma,\;\gamma^{d_1}u_1,\dots,\gamma^{d_{n-2}}u_{n-2}\bigr)
/\gamma^{\,e_i}$ a \emph{polynomial} in $(\gamma,u)$; the last two make the
$x$- and $\gamma$-powers cancel in the determinant chain below.
```

(ii) **The stage**: `v_j=x^{m_j}y_j`; `(\gamma,u_1,\dots,u_{n-2})^T = \theta^0 + Lv + \sum_\mu q_\mu \mu(v)` — base point, invertible linear part, correction vectors on monomials of degree ≥2. Lemma lem:stage: if every correction vector lies in `V=\mathrm{span}\{\mathrm{col}_j : j\notin J\}` (J = indices of variables occurring in correction monomials), the stage is a polynomial automorphism with constant Jacobian det L.

(iii) **Side conditions** (eq:side): `[\mu](E_i \circ \text{stage}) = 0` for every monomial μ(v) of x-degree < e_i. Triangular hierarchy: r=0 gives `E_i(\theta^0)=0` (hypersurface through origin); r=1 gives `\nabla E_i(\theta^0)\cdot\mathrm{col}_j=0`; r≥2 gives kernel conditions plus inhomogeneous linear equations on correction vectors. "Every unknown ... enters (eq:side) for the first time linearly, so the hierarchy is solved by finite-dimensional linear algebra, order by order."

**Chain factorization (eq:chain, lines 808–816), verbatim:**

```latex
\det J(F)\;=\;
\underbrace{x^{\sum m_j}}_{\text{monomials}}\cdot
\underbrace{\det L}_{\text{stage}}\cdot
\underbrace{\gamma^{\sum d_j}}_{\text{scalings}}\cdot
\underbrace{c\,\gamma^{\,k+1}}_{\det J(F_0)}\cdot
\underbrace{(\gamma x)^{-\sum e_i}}_{\text{twist}}
\;=\;c\,\det L ,
```

**Discrete-data table for the six maps (lines 818–831), verbatim:**

```latex
\begin{tabular}{lcccccc}
\hline
map & $n$ & $\det J(F_0)$ & $(d_j)$ & $(e_i)$ & $(m_j)$ & $\det J(F)=\pm c\det L$\\
\hline
$F$   & $3$ & $2\gamma^2$ & $(1)$     & $(1,2)$     & $(1,2)$     & $-2$\\
$G$   & $3$ & $2\gamma^2$ & $(1)$     & $(1,2)$     & $(1,2)$     & $2$\\
$F_4$ & $4$ & $2\gamma^2$ & $(1,1)$   & $(1,2,1)$   & $(1,1,2)$   & $-\tfrac{44}{9}$\\
$F_5$ & $4$ & $\gamma^3$  & $(1,2)$   & $(1,2,3)$   & $(1,2,3)$   & $\tfrac{160}{29}$\\
$F_6$ & $5$ & $\gamma^2$  & $(1,3,4)$ & $(1,2,3,4)$ & $(1,2,3,4)$ & $-290$\\
$F_7$ & $5$ & $\gamma^2$  & $(1,3,4)$ & $(1,2,3,4)$ & $(1,2,3,4)$ & $119377$\\
\hline
\end{tabular}
```

---

## 3. Theorem collapse5 and ALL open problems

### 3.1 Theorem thm:collapse5 (lines 1642–1669), verbatim:

```latex
\begin{theorem}[collapse on the characteristic-invariance locus]\label{thm:collapse5}
Let $c\neq0$, let $\varphi(w_1,w_2)$ be arbitrary, and suppose the datum $G_2$ is
constant along the characteristic field, $\mathcal WG_2=0$, i.e.\
$G_2=G_2(w_1,\theta)$. Set
\[
K(w_1,\theta)\;=\;c\,\frac{\partial G_2}{\partial w_1}\Big|_{\theta}
\;+\;\frac{\partial\,(G_2^{\,2})}{\partial\theta}.
\]
Then the master equation \eqref{eq:master5} factors as a perfect square,
\[
4c\,D_0-D_1^{\,2}\;=\;-\bigl(\mathcal WG_4-K\bigr)^2 ,
\]
so it linearizes to the transport equation $\mathcal WG_4=K$, whose general
polynomial solution is $G_4=(w_2/c)\,K(w_1,\theta)+h(w_1,\theta)$ with $h$
arbitrary. For every such solution the components $X_1,X_2,X_3$ of the swept
threefold depend on $(w_1,\theta)$ alone, $X_4=h(w_1,\theta)$, and in the
parameters $(\gamma,w_1,\theta,w_2)$ the sweep is the triangular suspension
\[
S\;=\;\bigl(T(\gamma,w_1,\theta),\ \gamma w_2+h(w_1,\theta)\bigr),
\qquad
T\;=\;\widehat X(w_1,\theta)+\gamma\,(1,w_1,w_1^2)^{\mathsf T},
\qquad
\det J(T)=c\,\gamma:
\]
a four-dimensional bottom-branch sweep with direction field $(1,w_1,w_1^2)^{\mathsf T}$,
extended by one shear-type coordinate. In particular no datum with
$\mathcal WG_2=0$ produces a genuinely five-dimensional example.
\end{theorem}
```

Context (middle branch, lines 1621–1640): `D_2\equiv c` integrates to `G_3=2w_1G_2-c\,w_3+\varphi(w_1,w_2)`; `X_1=-D_1/(2c)`; the master equation (eq:master5) is `4c\,D_0=D_1^{\,2}`; the characteristic field is `\mathcal W=c\,\partial_2+(\partial_2\varphi)\,\partial_3` with invariant `\theta=w_3-\varphi(w_1,w_2)/c`. The equation is "a Hamilton--Jacobi-type equation degenerate along $\mathcal W$."

Follow-up Remark (lines 1693–1700), verbatim:

```latex
\begin{remark}
The collapse is not confined to the locus $\mathcal WG_2=0$. The
complementary explicit families we have solved --- e.g.\ $\varphi=0$,
$G_2$ independent of $w_3$ and $G_4$ linear in $w_3$ with constant coefficient ---
also collapse, in a different mode ($w_3$ enters all components linearly through
constant coefficients, and a linear change of target coordinates splits off a
graph coordinate).
\end{remark}
```

### 3.2 Open Problem prob:rigid5 (lines 1702–1709), verbatim — the ONLY formal `\begin{problem}` environment in the paper:

```latex
\begin{problem}[rigidity of the middle branch]\label{prob:rigid5}
Does the master equation \eqref{eq:master5} admit a polynomial solution that is
not equivalent to a suspension of a four-dimensional sweep? We conjecture that
it does not. If so, this is a genuine rigidity phenomenon: for the mixed field
$\Delta=(1,w_1,w_1^2,w_2)^{\mathsf T}$ the middle escape $\det J(F_0)=c\gamma^3$ is
unrealizable by irreducibly five-dimensional data, even though it is
unobstructed at the level of branch counting.
\end{problem}
```

### 3.3 All other open questions stated in the paper (verbatim):

(a) After thm:F4 (lines 963–966): "Whether $F_4$ is equivalent to $\Phi\times\mathrm{id}$ for a degree-\emph{five} member $\Phi$ of Gallagher's family \cite{Gallagher2026} is left open; the fiber stratification of a product has product form, so the stratification deferred to Appendix~\ref{app:F4} is expected to decide this."

(b) Concluding remarks (lines 1960–1973), verbatim:

```latex
Several questions suggest themselves. Which fiber-size sets are realizable by
Keller maps in dimension $n$? (For $F$ no fiber has size one less than the
generic size, yet $G$ does attain size $4-1=3$ --- over the degenerate hyperplane
rather than over a contact stratum --- so the gap phenomenon is not universal and
seems to be contact-geometric in nature.) Do
non-graph surfaces in $\C^4$ admit hyperbolic sweep normalizations, which would
give four-dimensional examples with two independent escape divisors? Is the
middle branch of the mixed direction field $(1,w_1,w_1^2,w_2)^{\mathsf T}$ genuinely
rigid --- that is, does the master equation \eqref{eq:master5} admit any
polynomial solution that is not a suspension of a four-dimensional sweep
(Theorem~\ref{thm:collapse5} and Problem~\ref{prob:rigid5})? More generally, for
which pairs (direction field, branch index) in dimension $n$ is the branch
realizable by irreducibly $n$-dimensional data? And what is
the minimal geometric degree of a Keller counterexample in dimension $n\ge4$?
```

(c) Unexplored items flagged in text: the M-branch heat-polynomial family ("every polynomial solution $(\rho,\lambda)$ of the linear equation \eqref{eq:wave} is a candidate surface, and the fiber geometry of the resulting degree-ten (and higher) maps remains to be understood", lines 1955–1958); complete stratifications of F_4–F_7 deferred to a subsequent version (lines 1948–1955, 2165–2167); propagation remark (lines 1816–1828): for n≥6 one may take H to be a CONSTANT family equal to a Keller counterexample in n−3 variables ("Alp\"oge's $F$ for $n=6$, the maps $F_4$ or $F_5$ for $n=7$, $F_6$ or $F_7$ for $n=8$, and so on --- the construction feeds on its own output, three dimensions at a time. ... For such non-invertible $H$ the fiber system is no longer triangular, and its geometry is unexplored.").

---

## 4. The five-dimensional maps F_6 and F_7

### 4.1 F_6 (bottom branch of the mixed field Δ=(1,w_1,w_1^2,w_2)^T; subsubsec:bottom5)

Branch solution (lines 1394–1423): `D_2\equiv0` gives `G_3=2w_1G_2+\varphi(w_1,w_2)`; taking `\varphi=\tau w_2`, `D_1\equiv c` reads `\partial_3(\tau G_4+G_2^2-c\,w_3)=0`, so the complete solution is `G_2` arbitrary trivariate, `\psi\in\C[w_1,w_2]` arbitrary, `G_3=2w_1G_2+\tau w_2`, `G_4=(c\,w_3-G_2^2)/\tau+\psi`, then `X_1=-D_0/c` — the branch is FLEXIBLE.

Chosen data (lines 1427–1442): `\tau=c=1`, `G_2=w_2+w_1w_3+w_1^2-w_1^3`, `\varphi=w_2`, `\psi=2w_1^4-2w_1^5`; explicitly:

```latex
X_1&=-w_3+2w_2-2w_1+2w_1w_3+5w_1^2-2w_1^3+8w_1^4-10w_1^5,\\
X_2&=w_2+2w_1w_2-w_1^2+2w_1^2w_3+4w_1^3-2w_1^4+8w_1^5-10w_1^6,\\
X_3&=w_2+2w_1w_2+w_1^2w_3+2w_1^2w_2+2w_1^3w_3+3w_1^4-2w_1^5+8w_1^6-10w_1^7,\\
X_4&=w_3-w_2w_3+w_2^2-2w_1w_2+3w_1^2w_2-w_1^2w_3^2-2w_1^3w_3+w_1^4+2w_1^4w_3+8w_1^4w_2-10w_1^5w_2-w_1^6 ,
```

of degrees 5,6,7,6, with `\det J(S)=\gamma`. Weights (1,3,4) on (w_1,w_2,w_3); scalings `w_1=\gamma u_1, w_2=\gamma^3u_2, w_3=\gamma^4u_3`; gradients at base point `\nabla E_2=(-16,-17,3,2)`, `\nabla E_3=(-17,-18,5,3)`, `\nabla E_4=(-2,-2,0,1)`. Stage on `v=(xy, x^2z_1, x^3z_2, x^4z_3)` (lines 1479–1485):

```latex
\gamma&=1-29v_1+999v_1^2+355v_1v_2-41553v_1^3+v_4,\\
u_1&=1+27v_1-5v_2,\\
u_2&=v_1-12v_2+\tfrac{2128}{5}v_1^2+v_3,\\
u_3&=-4v_1-10v_2 .
```

(elementary polynomial automorphism, constant Jacobian −290). With `S_1=X_1+\gamma`, `S_2=X_2+\gamma w_1`, `S_3=X_3+\gamma w_1^2`, `S_4=X_4+\gamma w_2`, `C=\gamma x`:

```latex
\boxed{\;F_6\;=\;\Bigl(C,\;\;\frac{S_1}{C},\;\;\frac{S_2}{C^2},\;\;
\frac{S_3}{C^3},\;\;\frac{S_4}{C^4}\Bigr)^{\mathsf T}\;:\;\C^5\to\C^5 .\;}
```

**Theorem thm:F6 (lines 1534–1538), verbatim:**

```latex
\begin{theorem}\label{thm:F6}
$F_6$ is a polynomial map $\C^5\to\C^5$ with components of degrees
$7,38,40,42,44$; its Jacobian determinant is identically $-290$; and its
generic fiber consists of exactly $6$ points.
\end{theorem}
```

Additional properties (lines 1540–1547): F_6 fixes the axis `{y=z_1=z_2=z_3=0}` pointwise; for every C_0≠0 the fiber over (C_0,0,0,0,0) has exactly 4 points; F_6 is "not equivalent to $\Phi\times\mathrm{id}$ for any four-dimensional Keller map $\Phi$ of geometric degree $\neq6$". Chain: `x^{10}\cdot(-290)\cdot\gamma^{8}\cdot\gamma^{2}\cdot(\gamma x)^{-10}=-290`. Fiber polynomial (lines 1584–1586):

```latex
R_Y(w_1)\;=\;2w_1^6-2w_1^5-w_1^3+(Y_1+1)\,w_1^2
+\bigl(Y_4+Y_2^{\,2}-Y_1Y_3-2Y_2+Y_1\bigr)\,w_1+\bigl(Y_3-Y_2\bigr)\;=\;0 ,
```

degree exactly 6, constant leading coefficient 2. Over Y=0: `R_0(w_1)=w_1^2(w_1-1)(2w_1^3-1)`. Only the first component is printed (line 1610): `F_{6,1}=x-29x^2y+999x^3y^2-41553x^4y^3+355x^4yz_1+x^5z_3`; the remaining components have 342, 421, 507, 904 terms, "reproduced, together with all verification scripts, in the ancillary files" (NOT in this .tex).

### 4.2 F_7 (curve-type field, subsec:curvefam)

**Theorem thm:curvered (reduction for the curve-type field, lines 1723–1749), verbatim:**

```latex
\begin{theorem}[reduction for the curve-type field]\label{thm:curvered}
Let $n\ge3$ and $\Delta=(1,w_1,\dots,w_1^{\,n-2})^{\mathsf T}$. With the
potentials $\Gamma=(G_2,\dots,G_{n-1})$ of \S\ref{subsec:forms}, set
\[
H_{i+1}\;:=\;G_{i+1}-i\,w_1^{\,i-1}G_2\qquad(2\le i\le n-2).
\]
Then $\widetilde P=\widehat\Delta{}'\,e_1^{\mathsf T}$ is of rank one, with
$\widehat\Delta{}'=(1,2w_1,\dots,(n-2)w_1^{\,n-3})^{\mathsf T}$, the pencil
is linear in $\mu$,
\[
\det\bigl(J(\Gamma)+\mu\widetilde P\bigr)\;=\;D_0+\mu D_1,
\qquad
D_0=\det J(\Gamma),
\qquad
D_1\;=\;\det J_{(w_2,\dots,w_{n-2})}\bigl(H_3,\dots,H_{n-1}\bigr),
\]
and the (unique) branch normalization $D_1\equiv c\neq0$ says precisely:
\begin{itemize}
\item[(i)] the family
$H(w_1;\,\cdot\,)=(H_3,\dots,H_{n-1})^{\mathsf T}\colon\C^{n-3}\to\C^{n-3}$ has constant
Jacobian determinant $c$ in the transverse parameters $(w_2,\dots,w_{n-2})$,
for every value of the spectator parameter $w_1$;
\item[(ii)] $G_2\in\C[w]$ is completely free.
\end{itemize}
Given such data, $X_1=-D_0/c$ and $X_{i+1}=G_{i+1}+w_1^{\,i}X_1$ yield
$\det J(F_0)=c\,\gamma^2$.
\end{theorem}
```

**Proposition prop:curvefib (univariate fiber theory, lines 1794–1814), verbatim:**

```latex
\begin{proposition}[univariate fiber theory]\label{prop:curvefib}
Assume in addition that each $H(w_1;\,\cdot\,)$ is a polynomial automorphism
of $\C^{n-3}$ whose inverse depends polynomially on $w_1$ (e.g.\ any
triangular family). Then for every target $Y\in\C^{n-1}$ the fiber system
$S(\gamma,w)=Y$ triangularizes: $\gamma=Y_1-X_1$, the potential equations
become $G_{i+1}(w)=Y_{i+1}-w_1^{\,i}Y_1$, hence
\[
H\bigl(w_1;\,w_2,\dots,w_{n-2}\bigr)\;=\;h(w_1),
\qquad
h_{i+1}\;=\;Y_{i+1}-i\,w_1^{\,i-1}Y_2+(i-1)\,w_1^{\,i}Y_1 ,
\]
so $(w_2,\dots,w_{n-2})=H^{-1}(w_1;h(w_1))$ are polynomial functions of $w_1$,
and the last remaining equation is the univariate \emph{fiber polynomial}
\[
R_Y(w_1)\;:=\;G_2\bigl(w_1,\,H^{-1}(w_1;h(w_1))\bigr)-\bigl(Y_2-w_1Y_1\bigr)
\;=\;0 .
\]
The fibers of the twisted map over targets with nonzero first coordinate are
in bijection with the roots of $R_Y$ at which $\gamma\neq0$. For $n=3$, $R_Y$
is the tangency polynomial of \S\ref{sec:sweep}.
\end{proposition}
```

Data for F_7 (lines 1830–1845): n=5, c=1, `H_3=w_2+w_1^3+w_1^2w_3`, `H_4=w_3+w_1^4`, `G_2=w_2^{\,2}+w_2+w_1w_3+w_1^3`; `G_3=2w_1G_2+H_3`, `G_4=3w_1^2G_2+H_4`, `X_1=-\det J(G_2,G_3,G_4)`, components X_1..X_4 of degrees 7,8,9,10 (not printed in .tex), `\det J(S)=\gamma`, base point `S(1,1,0,-1)=0`, `D_0(1,0,-1)=1`. Weights (1,3,4), gradients `\nabla E_2=(-15,0,-3,4)`, `\nabla E_3=(-20,3,-1,6)`, `\nabla E_4=(-19,8,-1,7)`. Stage (lines 1857–1864):

```latex
\gamma&=1-61v_1+9012v_1^2+\tfrac{238}{19}v_1v_2-\tfrac{35823670}{19}v_1^3+8v_3+v_4,\\
u_1&=1+59v_1+v_2-5178v_1^2+19v_3,\\
u_2&=-7v_1-27v_2,\\
u_3&=-1-234v_1-5v_2 ,
```

(elementary polynomial automorphism, constant Jacobian 119377). With `C=\gamma x`, `S_i=X_i(w_1,w_2,w_3)+\gamma w_1^{i-1}`:

```latex
\boxed{\;F_7\;=\;\Bigl(C,\;\;\frac{S_1}{C},\;\;\frac{S_2}{C^2},\;\;
\frac{S_3}{C^3},\;\;\frac{S_4}{C^4}\Bigr)^{\mathsf T}\;:\;\C^5\to\C^5 .\;}
```

**Theorem thm:F7 (lines 1878–1882), verbatim:**

```latex
\begin{theorem}\label{thm:F7}
$F_7$ is a polynomial map $\C^5\to\C^5$ with components of degrees
$7,86,89,92,95$; its Jacobian determinant is identically $119377$; and its
generic fiber consists of exactly $12$ points.
\end{theorem}
```

Additional (lines 1884–1891): chain `x^{10}\cdot119377\cdot\gamma^{8}\cdot\gamma^{2}\cdot(\gamma x)^{-10}`; F_7 fixes the axis `{y=z_1=z_2=z_3=0}` pointwise; fiber over (C_0,0,0,0,0), C_0≠0, has exactly 7 points; not equivalent to Φ×id for any 4-dim Keller map Φ of geometric degree ≠12. Fiber polynomial (lines 1898–1904):

```latex
p_2(w_1)=Y_4-3Y_2w_1^2+2Y_1w_1^3-w_1^4,\qquad
p_1(w_1)=Y_3-2Y_2w_1+Y_1w_1^2-w_1^3-w_1^2\,p_2(w_1),
\]
\[
R_Y(w_1)\;=\;p_1(w_1)^2+p_1(w_1)+w_1\,p_2(w_1)+w_1^3+Y_1w_1-Y_2 ,
```

degree exactly 12, constant leading coefficient 1. Over Y=0: `R_0(w_1)=w_1^{\,5}(w_1-1)(w_1^6+w_1^5+w_1^4-w_1^3-w_1^2-w_1+1)`. First component (lines 1923–1926): `F_{7,1}=x-61x^2y+9012x^3y^2-\tfrac{35823670}{19}x^4y^3+\tfrac{238}{19}x^4yz_1+8x^4z_2+x^5z_3`; the remaining components have between 14,000 and 26,000 terms, in ancillary files (NOT in this .tex).

---

## 5. Symmetric Jacobians / gradient maps / Hessians / dimension two

**Hessians, symmetric Jacobians, gradient maps: NOT MENTIONED ANYWHERE in the paper.** Verified by grep over the full .tex: zero occurrences of "Hessian"/"hessian"/"symmetric". The words "potential(s)" and "gradient(s)" occur only in the internal sense of the construction — the scalar potentials `G_i = X_i - \delta_i X_1` (integration devices for the tangency PDE, subsec:forms lines 680–713) and gradients `\nabla E_i` of the side-condition functions — never in the sense of gradient maps F = ∇f or symmetric Jacobian matrices. The paper contains no statement bearing on the Hessian conjecture HC_n, and none of its counterexamples is claimed to have symmetric Jacobian. (Nearest structural connection, for what it is worth, is unstated: eq:wave and the heat-polynomial solutions arise from equality of mixed partials, i.e., integrability, but the paper draws no symmetric-Jacobian consequence.)

**The two-dimensional case (lines 178–182), verbatim — the only remark on dimension two:**

```latex
In dimension two, Moh
\cite{Moh1983} verified the conjecture for maps of degree at most $100$; the
two-dimensional case remains open and is untouched by the counterexamples
discussed here, which exist only in dimension $\ge3$ (by Wang's theorem, degree
$2$ examples are impossible, and the known constructions produce degree $\ge3$).
```

Related framing (lines 184–198): Białynicki-Birula–Rosenlicht (injective polynomial self-map of C^n is bijective with polynomial inverse) forces a Keller counterexample to fail injectivity; étale-ness forces non-injectivity to occur only through failure of properness (Jelonek non-properness set is empty or a hypersurface); Pinchuk's real 2-dim example (nonvanishing, NONconstant Jacobian) cited as an early warning — noted as real category, nonconstant Jacobian, hence not a 2-dim complex Keller counterexample.

---

## 6. Fiber stratification theorems for the 3-dim maps (context)

- Theorem thm:alpoge (lines 395–399): Alpöge's F has component degrees 7,6,4, `\det J(F)\equiv-2`, generic fiber exactly 3 points. Triple collision: `F(0,0,-\tfrac14)=F(1,-\tfrac32,\tfrac{13}{2})=F(-1,\tfrac32,\tfrac{13}{2})=(-\tfrac14,0,0)`.
- Theorem thm:strat3 (lines 434–452): full stratification of F via `c_3(v)=27v_1^2v_3^2-18v_1v_2v_3+16v_1+v_2^3v_3-v_2^2`; fiber sizes 3/1/0, never 2; image = C^3 minus the curve `\mathcal C=\{(\tfrac4{27}t^{-2},\tfrac43t^{-1},t):t\ne0\}`; Jelonek set = irreducible surface {c_3=0}.
- Theorem thm:deg4 (lines 478–482): G has component degrees 4,11,12, `\det J(G)\equiv2`, generic fiber exactly 4 points. Data: `p(w)=w^3-6w^2+6w`, `q(w)=\tfrac38w^4-2w^3+\tfrac32w^2`, `\gamma=2-4xy-x^2z`, `u=1+xy`, `w=\gamma u`.
- Theorem thm:fibG (lines 2110–2134): G attains every fiber size in {4,3,2,1,0}; reconstruction formulas `x=v_1/\gamma`, `y=(w-\gamma)/v_1`, `z=\gamma(6\gamma-4w-\gamma^2)/v_1^2`; over v_1=0 the fiber has 3 points if `v_2^2\ne24v_3`, else 1.
- Appendix app:F (lines 1997–2040): the six-element reduced lex Gröbner basis h_1..h_6 of the graph ideal of F, certified (S-polynomials, substitution checks, mod-p reproduction at p=2147483629).
- Appendix app:G: reduced lex basis of J_G has 174 elements (139 MB); elimination quartic `R_G(x;v)=E(X,Y)x^4-24v_1^2(3X^2+8X+16Y-24)x^2-64v_1^3(X+4)x-16v_1^4`, with E the implicit equation of the swept quartic; z-side quartic Z_G with leading coefficient v_1^2.
- Appendix app:F4 (lines 2149–2167): Gröbner bases for F_4–F_7 out of reach; fiber theory via the univariate polynomials above; complete stratifications left to a subsequent version.

## 7. Verification-relevant notes for symbolic re-checkers

- Column-vector convention: Jacobian COLUMNS are ∂_i F; variable order matters for signs (e.g., `(x,\gamma,w_1,\dots)` in lem:JL; `(\gamma,w)` for plane sweep).
- F_4 uses source variables (x,y,z,t) with monomials `(v_1,v_2,v_3)=(xy,xz,x^2t)` of x-degrees (1,1,2) — note e-vector (1,2,1) is NOT sorted: `F_4=(C,S_1/C,S_2/C^2,S_3/C)`.
- F_5 uses `(v_1,v_2,v_3)=(xy,x^2z,x^3t)`, scalings `(w_1,w_2)=(\gamma u_1,\gamma^2 u_2)`, exponents (1,2,3), and a unipotent NONLINEAR stage (the `\tfrac{4197}{1682}v_1^2` correction).
- F_6/F_7 use source (x,y,z_1,z_2,z_3), `v=(xy,x^2z_1,x^3z_2,x^4z_3)`, scalings `(\gamma u_1,\gamma^3u_2,\gamma^4u_3)`, exponents (1,2,3,4).
- All identities claimed verified in exact rational arithmetic; scripts "available from the author" / ancillary files (not in this .tex). Full F_6 components (342–904 terms) and F_7 components (14k–26k terms) are NOT in this file — reconstruct from the recipe data above.
- Paper's own consistency cross-checks that a re-verifier can replicate: chain factorizations (eq:chain and per-map versions), tangency identities `\det(\Delta,J(X))\equiv0`, branch normalizations (L=2,M=0 for F_4; L=0,M=1 for F_5; D_1=c,D_2=0 for F_6; D_1=c for F_7), fiber polynomials and their Y=0 factorizations.
