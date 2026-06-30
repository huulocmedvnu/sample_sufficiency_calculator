# Anisotropic Sample-Sufficiency for Single-Cell Perturbation-Direction Estimation

**A rigorous derivation of the cell quota `n*` under a general (non-isotropic) noise covariance.**

The first release used the isotropic simplification `Σ = σ²I`. That is inadmissible for single-cell
data, whose per-cell covariance is strongly anisotropic and low-rank-dominated (cell cycle, lineage,
ambient/technical axes). This document replaces it with the exact quadratic-form result and states the
estimator, its full sampling distribution, finite-sample (tail) guarantees, the regime of validity, and
the plug-in procedure.

---

## 1. Model and assumptions

Work in a fixed `d`-dimensional embedding (e.g. the shared PCA space; the argument is coordinate-free).

- **(A1) Sampling.** Treated cells $x^t_1,\dots,x^t_{n_t}$ are i.i.d. with mean $\mu_t$ and covariance
  $\Sigma_t$; control (DMSO) cells $x^c_1,\dots,x^c_{n_c}$ i.i.d. with mean $\mu_c$, covariance
  $\Sigma_c$; the two pools are independent. Second moments are finite. (Heteroscedasticity allowed:
  $\Sigma_t \neq \Sigma_c$.)
- **(A2) Estimands.** The perturbation vector $v := \mu_t-\mu_c$, magnitude $m:=\lVert v\rVert>0$,
  and **direction** $u := v/m \in \mathbb S^{d-1}$. The scientific quantity of interest in a
  similarity/MoA screen is $u$ (cosine geometry); $m$ is a nuisance (potency).
- **(A3) Estimator.** $\hat v := \hat\mu_t-\hat\mu_c$ with $\hat\mu_t,\hat\mu_c$ the sample means;
  $\hat u := \hat v/\lVert\hat v\rVert$.

By the multivariate CLT and independence,

$$\hat v \;\sim\; \mathcal N\!\big(v,\; S\big),\qquad
  S \;:=\; \frac{\Sigma_t}{n_t}+\frac{\Sigma_c}{n_c}. \tag{1}$$

`S` is the sampling covariance of the *measured perturbation vector*; it is `O(1/n)`. Write the error
$e:=\hat v-v\sim\mathcal N(0,S)$.

Let $P := I-uu^\top$ be the orthogonal projector onto the tangent space $T_u\mathbb S^{d-1}$
(rank $d-1$). Note $P=P^\top=P^2$, $Pu=0$.

---

## 2. The angular error, to first order (Delta method)

**Definition.** The angular error is $\theta := \arccos(\hat u^\top u)\in[0,\pi]$.

**Lemma 1 (Jacobian of normalization).** The map $g(x)=x/\lVert x\rVert$ has Jacobian
$Dg(x)=\frac{1}{\lVert x\rVert}\big(I-\tfrac{xx^\top}{\lVert x\rVert^2}\big)$. At $x=v$ this is
$Dg(v)=\tfrac1m P$. The radial direction $u$ lies in $\ker Dg(v)$: changes in $\hat v$ **along** $u$
do not move the direction; only the **perpendicular** part does.

*Proof.* Differentiate $g_i=x_i\lVert x\rVert^{-1}$:
$\partial g_i/\partial x_j=\delta_{ij}\lVert x\rVert^{-1}-x_ix_j\lVert x\rVert^{-3}
=\lVert x\rVert^{-1}(\delta_{ij}-\hat x_i\hat x_j)$. At $x=v$, $\hat x=u$. $\square$

**Theorem 1 (anisotropic mean-squared angular error).** Under (A1)–(A3), in the high-SNR regime
(Section 5),

$$\boxed{\;\mathbb E[\theta^2]\;=\;\frac{\operatorname{tr}\!\big(P S P\big)}{m^2}\;+\;O(\rho^{-4})
   \;=\;\frac{\operatorname{tr}(S)-u^\top S\,u}{m^2}\;+\;O(\rho^{-4})\;}\tag{2}$$

with $S=\Sigma_t/n_t+\Sigma_c/n_c$ and $\rho^2:=m^2/(u^\top S u)$ the along-signal SNR.

*Proof.* To first order $\hat u-u = Dg(v)\,e+O(\lVert e\rVert^2)=\tfrac1m Pe+O(\lVert e\rVert^2)$.
Because $Pe\perp u$, for small angles $\theta=\lVert\hat u-u\rVert+O(\theta^3)=\tfrac1m\lVert Pe\rVert
+O(\lVert e\rVert^2)$. Hence
$\mathbb E[\theta^2]=\tfrac1{m^2}\mathbb E\lVert Pe\rVert^2+\dots
=\tfrac1{m^2}\mathbb E[e^\top P e]+\dots=\tfrac1{m^2}\operatorname{tr}(PSP)+\dots$,
using $\mathbb E[e^\top Pe]=\operatorname{tr}(P\,\mathbb E[ee^\top])=\operatorname{tr}(PS)$ and
$\operatorname{tr}(PSP)=\operatorname{tr}(PS)=\operatorname{tr}(S)-u^\top S u$ (idempotence of $P$).
The remainder is computed in Section 5. $\square$

**Reading.** $\operatorname{tr}(S)$ is the *total* sampling noise; $u^\top S u$ is the noise **along the
signal**, which is harmless to the direction and is *subtracted out*. Only the perpendicular noise
$\operatorname{tr}(PSP)$ rotates $\hat u$. The isotropic formula keeps all $d-1$ perpendicular
dimensions at the *same* variance; the anisotropic truth weights each axis by its own variance and by
its alignment with $u$.

---

## 3. The cell quota `n*`

Fix a target RMS angular tolerance $\theta_\star$ and set $\mathbb E[\theta^2]=\theta_\star^2$.

**Heteroscedastic, general arms.** With $S=\Sigma_t/n_t+\Sigma_c/n_c$, solving (2) gives the feasible
$(n_t,n_c)$ frontier $\operatorname{tr}\!\big(P(\Sigma_t/n_t+\Sigma_c/n_c)P\big)=m^2\theta_\star^2$.

**Equal arms** ($n_t=n_c=n$, $\Sigma_t\approx\Sigma_c\approx\Sigma$): $S=2\Sigma/n$ and

$$\boxed{\;n^\star_{\text{aniso}}
   =\frac{2\,\operatorname{tr}(P\Sigma P)}{m^2\,\theta_\star^2}
   =\frac{2\big(\operatorname{tr}\Sigma-u^\top\Sigma u\big)}{m^2\,\theta_\star^2}\;}\tag{3}$$

**Large shared DMSO pool** ($n_c\to\infty$): the control arm is noiseless, factor $2\to1$:
$n^\star=\operatorname{tr}(P\Sigma_t P)/(m^2\theta_\star^2)$.

**Isotropic reduction (consistency check).** If $\Sigma=\sigma^2I$ then
$\operatorname{tr}(P\Sigma P)=\sigma^2\operatorname{tr}P=\sigma^2(d-1)$ and (3) collapses to
$n^\star=2\sigma^2(d-1)/(m^2\theta_\star^2)$ — exactly the v1 formula. The anisotropic result thus
**contains** the isotropic one as the special case $\Sigma\propto I$.

**Effective dimension.** Writing $\bar\sigma^2_\perp:=\operatorname{tr}(P\Sigma P)/(d-1)$ (mean noise
variance in directions $\perp u$), (3) reads $n^\star=2(d-1)\bar\sigma^2_\perp/(m^2\theta_\star^2)$:
the *form* is unchanged, but $\sigma^2$ is replaced by the **signal-orthogonal mean variance**.

---

## 4. Full sampling distribution, variance, and a finite-sample (tail) guarantee

Controlling the *mean* angular error is not enough for a reviewer: we want
$\Pr(\theta>\theta_\star)\le\delta$. Decompose $PSP=\sum_{i=1}^{d-1}\nu_i\,w_iw_i^\top$
(spectral, $\nu_i\ge0$). Then $\lVert Pe\rVert^2\stackrel{d}{=}\sum_i\nu_i z_i^2$, $z_i\stackrel{iid}\sim
\mathcal N(0,1)$ — a **generalized chi-square**, so

$$\theta^2\;\approx\;\frac1{m^2}\sum_{i=1}^{d-1}\nu_i z_i^2,\qquad
  \mathbb E[\theta^2]=\frac{\operatorname{tr}(PSP)}{m^2},\qquad
  \operatorname{Var}(\theta^2)=\frac{2\operatorname{tr}\!\big((PSP)^2\big)}{m^4}.\tag{4}$$

**Effective degrees of freedom** (participation ratio):

$$d_{\text{eff}}:=\frac{\big(\operatorname{tr}\,PSP\big)^2}{\operatorname{tr}\big((PSP)^2\big)}
   \;\le\;d-1,$$

which is $\ll d-1$ when a few noise axes dominate — the angle "sees" far fewer effective dimensions
than $d$, and the quota is correspondingly smaller than a naive $d$-scaling suggests.

**Proposition 1 (confidence quota, Hanson–Wright).** There is an absolute $c>0$ such that for $t>0$,
$\Pr\!\big(\lVert Pe\rVert^2-\operatorname{tr}(PSP)\ge t\big)\le
\exp\!\big(-c\min\{t^2/\lVert PSP\rVert_F^2,\;t/\lVert PSP\rVert_{\mathrm{op}}\}\big).$
Hence, to guarantee $\Pr(\theta>\theta_\star)\le\delta$ it suffices that

$$m^2\theta_\star^2\;\ge\;\operatorname{tr}(PSP)
  +\sqrt{\tfrac{2}{c}\,\lVert PSP\rVert_F^2\,\log\tfrac1\delta}
  +\tfrac{2}{c}\,\lVert PSP\rVert_{\mathrm{op}}\log\tfrac1\delta.$$

Because $S=2\Sigma/n$ makes every term $\propto 1/n$, this inverts to a closed-form confidence quota

$$n^\star_\delta
  =\frac{2}{m^2\theta_\star^2}\Big[\operatorname{tr}(P\Sigma P)
  +\sqrt{\tfrac{2}{c}\,\lVert P\Sigma P\rVert_F^2\log\tfrac1\delta}
  +\tfrac{2}{c}\,\lVert P\Sigma P\rVert_{\mathrm{op}}\log\tfrac1\delta\Big]
  \;=\;n^\star_{\text{aniso}}\cdot\Big(1+O\big(\sqrt{\log(1/\delta)/d_{\text{eff}}}\big)\Big).\tag{5}$$

So the tail-controlled quota is the mean quota inflated by a factor governed by $d_{\text{eff}}$ and
$\log(1/\delta)$ — it is *not* a different scaling law, just a multiplicative safety margin.

---

## 5. Regime of validity and the second-order bias

The Delta expansion is exact to leading order; the correction comes from the fluctuating lever arm.
From $\tan\theta=\lVert Pe\rVert/(m+u^\top e)$,

$$\mathbb E[\theta^2]=\frac{1}{m^2}\,\mathbb E\!\left[\frac{\lVert Pe\rVert^2}{(1+u^\top e/m)^2}\right]
  =\frac{\operatorname{tr}(PSP)}{m^2}\Big(1+\frac{3}{m^2}\big(u^\top S u
  +2\,\tfrac{u^\top S P S u}{\operatorname{tr}(PSP)}\big)+\dots\Big),$$

using Isserlis' theorem (the $O(e^3)$ cross term vanishes by symmetry). The relative correction is
$O\big(u^\top S u/m^2\big)=O(\rho^{-2})$. **Validity:** the formula holds when the along-signal SNR
$\rho^2=m^2/(u^\top S u)\gg1$, i.e. when $n\gg u^\top\Sigma u/m^2$ (per arm). At the quota itself this is
automatically satisfied whenever $\theta_\star\ll1$, since $\theta_\star^2=\operatorname{tr}(PSP)/m^2
\ge u^\top SPS u/\dots$ — concretely, $\rho^{-2}\approx\theta_\star^2\cdot(u^\top\Sigma u)/\operatorname{tr}(P\Sigma P)$,
which is small for tight tolerances. Outside this regime (very weak signals, $\rho\lesssim1$) the
angle is near-uniform on the sphere and no finite quota "clears" the noise — the correct report is
"undetectable at this depth," not a number.

---

## 6. Detection mode (the "clear over DMSO" reading) is a *different* anisotropic functional

If "enough" means *detecting* a shift rather than *orienting* it, the relevant statistic is the
two-sample Hotelling $T^2$, with non-centrality

$$\lambda=\frac{n_tn_c}{n_t+n_c}\,\Delta^\top\Sigma^{-1}\Delta,\qquad \Delta=\mu_t-\mu_c,$$

so the controlling scalar is the **Mahalanobis** magnitude $m_M^2:=v^\top\Sigma^{-1}v$, and for equal
arms $n^\star_{\text{detect}}=2\lambda_0(d,\alpha,\text{power})/m_M^2$. Note the contrast: detection
uses $\Sigma^{-1}$ (it *rewards* directions where noise is small), whereas direction-tolerance uses
$\operatorname{tr}(P\Sigma P)$ (it is *hurt* by perpendicular noise). A complete tool exposes both;
they share the $n^\star\propto(\text{noise})/(\text{signal}^2)$ skeleton but with different anisotropic
noise functionals.

---

## 7. Plug-in estimation from a pilot

In practice $\Sigma,u$ are unknown. Use the pooled **within-condition** sample covariance
$\hat\Sigma$ (centre each condition by its own mean, then pool — this removes the perturbation and
estimates residual heterogeneity, not total variance) and $\hat u=\hat v/\lVert\hat v\rVert$ from the
pilot. Then $\hat n^\star=2\operatorname{tr}(\hat P\hat\Sigma\hat P)/(\hat m^2\theta_\star^2)$.

- **Consistency.** $\hat\Sigma\to\Sigma$, $\hat u\to u$; the plug-in is consistent with leading error
  $O(\sqrt{d/n_{\text{pilot}}})$ from $\hat\Sigma$ and $O(\theta_{\text{pilot}})$ from $\hat u$
  (the trace functional is 1-Lipschitz in $\hat u$ to first order, so direction error enters at second
  order — a convenient robustness).
- **High dimension.** When $d\not\ll n_{\text{pilot}}$, regularize $\hat\Sigma$ (Ledoit–Wolf shrinkage)
  before taking the trace, or — most cheaply and exactly in a PCA basis — use the per-component
  variances as the diagonal of $\hat\Sigma$ (PCA decorrelates, so $\hat\Sigma\approx\operatorname{diag}
  (\ell_1,\dots,\ell_d)$ with $\ell_k$ the explained variances), giving
  $\operatorname{tr}(P\Sigma P)=\sum_k\ell_k-\sum_k u_k^2\ell_k$.

---

## 8. When isotropy lies (quantified) — the reviewer's objection, made precise

Let $\bar\sigma^2:=\operatorname{tr}\Sigma/d$. The isotropic quota mis-estimates the truth by

$$\frac{n^\star_{\text{aniso}}}{n^\star_{\text{iso}}}
  =\frac{\operatorname{tr}(P\Sigma P)}{(d-1)\bar\sigma^2}
  =\frac{\operatorname{tr}\Sigma-u^\top\Sigma u}{(d-1)\bar\sigma^2}.$$

This ratio is **$<1$ when the signal aligns with high-variance axes** (the big variance is along $u$,
harmless → the drug is *easier* than isotropy predicts) and **$>1$ when the signal lies in a quiet
direction while large nuisance axes sit perpendicular** (those nuisance axes dominate
$\operatorname{tr}(P\Sigma P)$ → *harder*). Single-cell $\Sigma$ is dominated by a handful of large
PCs (cell cycle, lineage), so this factor routinely departs from 1 in both directions — which is
exactly why an isotropic quota is not publishable, and why the diagonal $\sum_k(1-u_k^2)\ell_k$ form is
the minimal honest upgrade.

---

## 9. Summary of the operational object

$$
n^\star_{\text{aniso}}=\frac{2\operatorname{tr}(P\Sigma P)}{m^2\theta_\star^2}
=\frac{2\sum_k(1-u_k^2)\ell_k}{m^2\theta_\star^2}
\quad\text{(PCA basis)},\qquad
n^\star_\delta=n^\star_{\text{aniso}}\Big(1+O\big(\sqrt{\log(1/\delta)/d_{\text{eff}}}\big)\Big).
$$

Inputs from one pilot: $\{\ell_k\}$ (per-component residual variances), the perturbation direction
$u$ and magnitude $m$, a tolerance $\theta_\star$ (or detection power), and optional confidence
$\delta$. The isotropic `σ²` of v1 is the degenerate case $\ell_k\equiv\sigma^2$. Empirical magnitude
of the correction on the Tahoe-100M atlas is reported in `src/calibrate.py` (anisotropic vs isotropic
quota per drug).
