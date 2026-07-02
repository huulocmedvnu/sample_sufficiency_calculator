---
title: "How Many Cells? A Plain-Math Derivation of the Cell Quota"
subtitle: "A companion to THEORY.md, using only basic matrix algebra and one integral"
geometry: margin=1in
header-includes: |
  \usepackage{amsmath}
  \setlength{\emergencystretch}{3em}
---

<!-- Plain-language companion to docs/THEORY.md. Build (from repo root):
       pandoc docs/THEORY_PRIMER.md -o theory_primer.pdf --pdf-engine=tectonic
     This file is written by hand (not generated); it carries no [@cite] keys and no bibliography. -->

## What this note is

`docs/THEORY.md` derives the sample-size rule at full rigour (Delta method, generalized
chi-square tails, a second-order expansion, a Lean proof). This note derives the **same headline
formula** using nothing beyond: adding and multiplying matrices, the *transpose* and *trace*, the idea
of an eigenvalue, and one basic integral. Everything else is bookkeeping with the expectation symbol
$\mathbb{E}[\cdot]$, which you can read as "the average over many repeats of the experiment."

The rule we will reach is

$$
\boxed{\;n^\star \;=\; \frac{2\,\operatorname{tr}(P\Sigma P)}{m^2\,\theta_\star^2}\;}
$$

the number of cells per arm needed to measure a perturbation's *direction* to within an angle
$\theta_\star$. Every symbol is defined below.

---

## 1. The setup

A drug screen compares **treated** cells to **control** (vehicle) cells. Each cell is a point in
$d$-dimensional space (here $d = 50$ coordinates; where they come from is spelled out at the end of this
section). Write a treated cell as $x^t \in \mathbb{R}^d$ and a control cell as $x^c$.

Stacking the cells row by row, the two groups are really two **matrices**, $X_t \in \mathbb{R}^{n_t\times d}$
and $X_c \in \mathbb{R}^{n_c\times d}$. The two groups almost never have the same number of cells
($n_t \neq n_c$), so you **cannot subtract the matrices**. Instead you first collapse each group to a
single point — its **centroid** (the average cell, a vector in $\mathbb{R}^d$):

$$
\hat\mu_t = \frac{1}{n_t}\sum_{i=1}^{n_t} x^t_i,
\qquad
\hat\mu_c = \frac{1}{n_c}\sum_{j=1}^{n_c} x^c_j ,
$$

and only then subtract the two $d$-vectors. The **perturbation vector** is that shift from control to
treated,

$$
\hat v = \hat\mu_t - \hat\mu_c ,
$$

an arrow in $\mathbb{R}^d$. Its true (infinite-data) value is $v$, with

- **length** $m = \lVert v\rVert$ — *how strongly* the drug acts (its potency), and
- **direction** $u = v / m$ — a unit vector saying *which way* the transcriptome moved (which
  biological program was engaged).

**The key modelling fact:** two drugs are called "similar mechanism" when their arrows point the same
way, i.e. when their **directions** $u$ agree. The length $m$ is a side issue. So the quantity we must
pin down is the *direction*, and the thing that can go wrong is that our measured direction $\hat u =
\hat v / \lVert\hat v\rVert$ is **tilted** away from the true $u$ by some angle $\theta$. We want enough
cells that this tilt stays below a tolerance $\theta_\star$.

**Where the coordinates come from — and why the order matters.** The $d = 50$ coordinates are a shared
"map" (a PCA space) fitted once on *individual* cells pooled across the whole atlas. The steps then run
in a fixed order: the per-cell scatter $\Sigma$ (Section 2) is measured from single cells *on that map*,
and only *afterwards* is each group averaged into its centroid. This order is essential — if you averaged
the raw gene counts *before* building the map, each group would collapse to one point with no scatter,
$\Sigma$ would be zero, and every formula below would fall apart. (From here on we take equal arms,
$n_t = n_c = n$, to keep the algebra light; unequal arms change nothing essential.)

---

## 2. Averaging beats down noise (the one integral)

Individual cells are noisy: a single treated cell scatters around $\mu_t$ with some per-cell
covariance matrix $\Sigma$ (a $d\times d$ symmetric matrix; its diagonal holds the variance of each
coordinate). Averaging $n$ of them shrinks that scatter.

**One-dimensional reminder.** If a single measurement $X$ has variance $\sigma^2$, the average of $n$
independent copies has variance $\sigma^2/n$. This is the only place probability enters, and it rests
on one integral: for a standard bell curve $\varphi(z) = \tfrac{1}{\sqrt{2\pi}}e^{-z^2/2}$,

$$
\mathbb{E}[z^2] = \int_{-\infty}^{\infty} z^2\,\varphi(z)\,dz = 1 .
$$

*(Integrate by parts with $u=z$, $dv = z\,e^{-z^2/2}dz$: the boundary term vanishes and what remains is
the total area $\int \varphi = 1$.)* This "variance of a standard coordinate is 1" is all we need;
scaling and summing gives everything else.

**In $d$ dimensions.** The centroid $\hat\mu_t$ therefore has covariance $\Sigma / n$, and likewise
$\hat\mu_c$. Because the two arms are independent, their difference adds the two covariances:

$$
\hat v = v + e,
\qquad
e \text{ has mean } 0 \text{ and covariance }
S = \frac{\Sigma}{n} + \frac{\Sigma}{n} = \frac{2\Sigma}{n}.
$$

Here $e$ is the measurement error — a small random vector that nudges the true arrow $v$ to the
observed arrow $\hat v$. Note $S$ shrinks like $1/n$: more cells, less wiggle. **This factor of 2 (two
noisy arms) is where the 2 in the final formula comes from.**

---

## 3. Only *sideways* noise tilts a direction

Now the geometry — this is the heart of the whole method, and it needs no calculus, just a picture.

The true arrow $v$ has length $m$ and points along $u$. The noise $e$ pushes the tip around. Split the
push into two pieces:

- the part **along** $u$: it lengthens or shortens the arrow but does **not** rotate it;
- the part **across** $u$ (perpendicular): it swings the tip sideways and **does** rotate it.

Only the second piece changes the direction. If the sideways push has length $\ell_\perp$ and the arrow
has length $m$, then by the right-angle triangle (opposite over adjacent),

$$
\tan\theta = \frac{\ell_\perp}{m} \;\approx\; \frac{\ell_\perp}{m}
\quad\Longrightarrow\quad
\theta \approx \frac{\ell_\perp}{m}
\qquad(\text{small angles: } \tan\theta \approx \theta).
$$

So the tilt angle is just *sideways displacement divided by arrow length*. Squaring,

$$
\theta^2 \approx \frac{\ell_\perp^2}{m^2}.
$$

This small-angle step is the only approximation in the derivation; it is excellent whenever the arrow
is much longer than the noise (a strong signal), which is exactly the regime a well-designed screen
operates in.

---

## 4. The projection matrix picks out the sideways part

How do we extract "the part of $e$ perpendicular to $u$"? With one matrix. Define

$$
P = I - u u^\top .
$$

Applied to any vector $e$, this **removes the component along $u$**:

$$
P e = e - u\,(u^\top e).
$$

The term $u^\top e$ is a single number — how much of $e$ points along $u$ — and $u(u^\top e)$ is that
"along" piece; subtracting it leaves exactly the sideways piece. Three facts about $P$, each a
one-line check using $u^\top u = 1$:

| Property | Check | Meaning |
|:--|:--|:--|
| $P^\top = P$ | $(uu^\top)^\top = uu^\top$ | symmetric |
| $Pu = 0$ | $u - u(u^\top u) = u - u = 0$ | it deletes the along-direction |
| $P^2 = P$ | $(I-uu^\top)^2 = I - 2uu^\top + u(u^\top u)u^\top = I - uu^\top$ | applying it twice = once |

So the sideways displacement is $\ell_\perp = \lVert P e\rVert$, and its square is

$$
\ell_\perp^2 = \lVert Pe\rVert^2 = (Pe)^\top(Pe) = e^\top P^\top P\, e = e^\top P\, e ,
$$

using $P^\top P = P^2 = P$. Combining with Section 3,

$$
\theta^2 \approx \frac{e^\top P e}{m^2}.
$$

---

## 5. Averaging the squared angle

The angle $\theta$ is random because $e$ is. We want its **typical** squared size, $\mathbb{E}[\theta^2]$.
We need one identity, and it follows from plain linearity of averaging.

**Identity.** For a symmetric matrix $A$ and a mean-zero random vector $e$ with covariance $S$,

$$
\mathbb{E}\!\left[e^\top A\, e\right] = \operatorname{tr}(A S).
$$

*Proof.* Write the quadratic form as a double sum, $e^\top A e = \sum_{i,j} A_{ij}\, e_i e_j$. Average
term by term; by definition of covariance, $\mathbb{E}[e_i e_j] = S_{ij}$. Hence

$$
\mathbb{E}[e^\top A e] = \sum_{i,j} A_{ij} S_{ij}
= \sum_i \Big(\sum_j A_{ij} S_{ji}\Big)
= \sum_i (A S)_{ii} = \operatorname{tr}(A S),
$$

using $S_{ij}=S_{ji}$ (covariance is symmetric). $\square$

Apply it with $A = P$ and $S = 2\Sigma/n$:

$$
\mathbb{E}[\theta^2]
= \frac{1}{m^2}\,\mathbb{E}[e^\top P e]
= \frac{1}{m^2}\operatorname{tr}\!\Big(P\cdot\frac{2\Sigma}{n}\Big)
= \frac{2}{n\,m^2}\operatorname{tr}(P\Sigma).
$$

Finally, a cosmetic rewrite makes the formula symmetric. Because $P^2 = P$ and the trace is unchanged
by cycling factors ($\operatorname{tr}(XY)=\operatorname{tr}(YX)$),

$$
\operatorname{tr}(P\Sigma) = \operatorname{tr}(P^2\Sigma) = \operatorname{tr}(P\Sigma P).
$$

So the **expected squared tilt** is

$$
\boxed{\;\mathbb{E}[\theta^2] = \frac{2\,\operatorname{tr}(P\Sigma P)}{n\,m^2}\;}
$$

More cells $n$: smaller tilt. Stronger signal $m$: much smaller tilt (it enters squared).

---

## 6. Solve for the number of cells

We want the typical tilt to sit at the tolerance: set $\mathbb{E}[\theta^2] = \theta_\star^2$ and solve
for $n$:

$$
\frac{2\,\operatorname{tr}(P\Sigma P)}{n\,m^2} = \theta_\star^2
\qquad\Longrightarrow\qquad
\boxed{\;n^\star = \frac{2\,\operatorname{tr}(P\Sigma P)}{m^2\,\theta_\star^2}\;}
$$

That is the cell quota. Read it as: *cells needed $=$ (twice the sideways noise) $\div$ (signal
strength squared $\times$ tolerance squared).*

---

## 7. What $\operatorname{tr}(P\Sigma P)$ means, and the simplest special case

The number $\operatorname{tr}(P\Sigma P)$ is the **total noise pointing sideways to the signal**.
There is a tidy way to see this. Using the identity from Section 5 in reverse, or just multiplying out
$P = I - uu^\top$,

$$
\operatorname{tr}(P\Sigma P) = \operatorname{tr}(\Sigma) - u^\top\Sigma\, u .
$$

In words: take the total per-cell noise $\operatorname{tr}(\Sigma)$ (the sum of the variances of all
$d$ coordinates), then **subtract the noise that lies along the signal**, $u^\top\Sigma u$ — because
that part only stretches the arrow and never tilts it. Only what's left can hurt the direction.

**Equal noise in every direction (isotropic case).** Suppose the noise is the same size in all
directions, $\Sigma = \sigma^2 I$. Then $\operatorname{tr}(\Sigma) = d\sigma^2$ and $u^\top\Sigma u =
\sigma^2$, so

$$
\operatorname{tr}(P\Sigma P) = \sigma^2(d-1),
\qquad
n^\star = \frac{2(d-1)\,\sigma^2}{m^2\,\theta_\star^2}.
$$

The $d-1$ (not $d$) is the honest bookkeeping: of the $d$ directions, the one *along* the signal is
free, leaving $d-1$ directions of harmful noise. This is the same $d-1$ that shows up as $P$ having one
zero eigenvalue (the $u$ direction) and $d-1$ eigenvalues equal to one.

---

## 8. Putting in real numbers

The Tahoe-100M atlas gives a per-coordinate noise of $\sigma^2 = 2.406$ in the $d=50$ PCA space. At the
standard tolerance $\theta_\star = 0.1$ radian (about $5.7^\circ$),

$$
n^\star = \frac{2\,(50-1)\,(2.406)}{m^2\,(0.1)^2}
= \frac{235.8}{0.01\,m^2}
\approx \frac{23{,}577}{m^2}\ \text{cells per arm}.
$$

So the whole design question collapses to the drug's signal strength $m$:

| signal strength $m$ | cells needed per arm $n^\star \approx 23{,}577/m^2$ |
|:--|:--|
| $m = 6$ (strong cytotoxic) | $\approx 650$ |
| $m = 3$ (moderate) | $\approx 2{,}600$ |
| $m = 1.3$ (typical/weak) | $\approx 14{,}000$ |
| $m = 1$ (very weak) | $\approx 23{,}600$ |

A strong perturbation is resolved in a few hundred cells; a weak one may need tens of thousands. Since
the median condition in Tahoe-100M has only $\sim$1,300 cells, most weak perturbations are
**under-sampled** — the single most important practical consequence of the formula.

---

## 9. What this rule does and does not promise

- **It controls the *average* tilt, not a worst case.** $\theta_\star$ is a root-mean-square (typical)
  angle. `THEORY.md` §4 adds a *confidence* version $n^\star_\delta$ that also bounds the chance of a
  large tilt; that step needs the distribution of $e^\top P e$ (a weighted sum of squared bell curves)
  and a tail inequality, which is why it lives in the full derivation.
- **It assumes the tilt is small** (the $\tan\theta\approx\theta$ step of Section 3). For very weak
  signals the arrow is barely longer than the noise, the angle is large, and *no* finite number of
  cells "clears" it — the honest report there is "undetectable at this depth," not a cell count.
- **It governs the direction of the *average* cell only** — not the fine structure of individual cells
  (rare subpopulations, local UMAP neighbourhoods). Those are set by different quantities.
- **The noise $\Sigma$ must be measured for your own platform.** It is a plug-in number, not a
  universal constant; the $\sigma^2 = 2.406$ above is specific to this atlas and pipeline.

Everything above is exact to leading order; `docs/THEORY.md` supplies the higher-order corrections, the
probability guarantee, and a machine-checked proof of the algebra. But the one-line takeaway is the
boxed formula of Section 6, and it follows from a triangle, a projection, and a single integral.
