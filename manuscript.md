---
title: "How fast does an epidemic die out on a network? Sharp decay envelopes for mean-field SIS epidemics and a proof of the Persoons–Sensi–Prasse–Van Mieghem conjectures"
author: "Anonymous"
date: "1 November 2026 · Version 0.1.0-candidate · [doi:10.5281/zenodo.23856194](https://doi.org/10.5281/zenodo.23856194)"
abstract: |
  **Proved.** For the N-intertwined mean-field (NIMFA) model of susceptible–infectious–susceptible (SIS) epidemics on networks, Persoons, Sensi, Prasse and Van Mieghem (Phys. Rev. E 109, 034308, 2024) conjectured two universal decay bounds. At or below the epidemic threshold the average infection level obeys $y(t)\le1/(1+t)$. Above it, the excess over the endemic level, starting from above, obeys the same bound. We prove both, the second with the sharper prefactor $1-\bar p$. The proof applies Jensen's inequality to the Volterra-type function $h(s)=-s-\log(1-s)$. The same argument covers time-varying networks, node-dependent recovery and, with Perron weights, directed networks. Via a known correlation inequality, the bounds also hold for the expected prevalence of the exact stochastic SIS process with a common recovery rate.

  **Sharpness.** Regular networks attain the envelope. At the threshold $t\,y(t)\to\kappa(G)\le1$, an explicit constant of the leading eigenvector equal to 1 only for regular graphs. On standard synthetic networks of 1000 nodes, $\kappa$ ranges from about 0.9 to about 0.04, decreasing as the eigenvector localises. Individual nodes can exceed the population envelope by a factor of order $n^{1/3}$.

  **Application.** For the timescale-separation question that motivated the conjectures, we give unconditional transition-time bounds; the growth part is exact on regular graphs but very conservative on strongly heterogeneous ones. High-precision numerics show that the paper's growth estimate (Lemma 5) is not a valid general bound away from $R_0\to1$.

  AI-assisted, unrefereed candidate. No external peer review, independent reproduction or formal verification is claimed.
---

# Plain-English summary

*The question.* An infection spreads through a network of contacts (people, farms, computers) but is not quite strong enough to persist. How quickly does it die out? Does the shape of the network matter, or only one summary number, the reproduction number $R_0$?

*The setting.* We study the most widely used deterministic network model of an infection that can recur, SIS (susceptible, infectious, susceptible again). Each node $i$ has an infection level $x_i(t)$ between 0 and 1, and $y(t)$ is the average over nodes. Time is measured in units of the average recovery time. With weekly recovery, 52 time units is a year.

*What was conjectured.* In 2024 Persoons, Sensi, Prasse and Van Mieghem ran extensive simulations and conjectured a universal speed limit. At or below the tipping point ($R_0\le1$), whatever the network and whatever the starting state, $y(t)\le1/(1+t)$. A similar limit holds for settling to the endemic level above the tipping point when the epidemic starts above that level. They did not prove either statement.

*What this paper shows.*

1. **Both conjectures are true** (Theorems 3 and 4). The proof rests on one "entropy-like" quantity that can only decrease at a controlled rate.
2. **The argument goes further.** It covers networks that change over time, nodes that recover at different speeds and, after a natural reweighting, one-way (directed) contacts. On directed networks the plain average can break the bound. Because the deterministic model over-estimates the exact random (stochastic) model on average, the speed limit also applies to the *expected* prevalence of the stochastic SIS model when all nodes recover at the same rate (Proposition 12).
3. **The limit is sharp, and shape matters at the tipping point.** Regular networks, where every node has the same number of contacts, attain $1/(1+t)$ exactly. Exactly at $R_0=1$, the average eventually behaves like $\kappa(G)/t$. Here $\kappa(G)\le1$ is computed from the network's leading eigenvector (Theorem 7) and equals 1 only for regular networks. On standard synthetic networks of 1000 nodes, $\kappa$ is about 0.8–0.9 for random and small-world networks. It is about 0.04–0.13 for scale-free (Barabási–Albert, Holme–Kim) and spatial (random geometric) networks, where infection concentrates in a small core. The time for the average to fall below 0.1% ranges from about 40 to about 920 recovery times at the same $R_0=1$.
4. **Small averages can hide persistent cores.** A small $\kappa$ means the *average* falls quickly while a well-connected core stays infected. A hub can carry a multiple of $1/t$ that grows like $n^{1/3}$ with network size $n$ (Proposition 9). This is not good news for elimination.
5. **$R_0=1$ is a knife-edge.** Just below the tipping point, decay turns exponential after roughly $1/(1-R_0)$ recovery times. The slow $1/t$ decay is specific to the threshold itself.
6. **The motivating question.** The 2024 paper used these bounds to decide when a changing contact network can be treated as frozen ("timescale separation"). We make the decay part unconditional, show that an assumption made there about mixed starting states is never needed, and give rigorous transition-time bounds for every $R_0$ (Section 9). Our rigorous growth bound is exact for regular networks but very conservative for strongly heterogeneous ones. Separately, high-precision numerical integration shows that one of the paper's growth estimates is not a valid general bound away from $R_0\to1$.

*What this paper does not show.* The sharpness results ($\kappa$, attainment) concern the deterministic mean-field model, not a particular real epidemic. The mean-field model is known to be least accurate near the threshold (Van Mieghem, Omic and Kooij 2009; Prasse and Van Mieghem 2020). The numerical examples use synthetic networks only.

# 1. Introduction

## 1.1 The model and the conjectures

Let $A$ be the adjacency matrix of a simple undirected graph on $n$ nodes, $\tau>0$ the effective infection rate and $B=\tau A$. With the recovery rate scaled to one, the NIMFA SIS equations (Lajmanovich and Yorke 1976; Van Mieghem, Omic and Kooij 2009) are
$$
\dot x_i=-x_i+(1-x_i)(Bx)_i,\qquad x_i(0)\in[0,1],\qquad y(t)=\frac1n\sum_{i=1}^n x_i(t). \tag{1.1}
$$
The basic reproduction number is $R_0=\rho(B)$, the spectral radius. For $R_0\le1$ every solution tends to $0$. For $R_0>1$ on a connected graph there is a unique endemic equilibrium $p\in(0,1)^n$, which attracts every nonzero initial state (Lajmanovich and Yorke 1976; Fall et al. 2007; Khanafer et al. 2016).

Persoons et al. (2024, Section III A) state the following. Here $u$ is the all-ones vector, $\tau_c=1/\rho(A)$, $K_N$ is the complete graph and $y(t;G,\tau,V(0))$ is the prevalence.

> **Conjecture 1** (their eq. 9). If $\tau\le\tau_c(G)$, then $y(t;G,\tau,V(0))\le y(t;G,\tau_c(G),u)\le y(t;K_N,\tau_c(K_N),u)=1/(1+t)$.
>
> **Conjecture 2** (their eq. 10). If $\tau>\tau_c(G)$ and $V(0)\ge V_\infty$ componentwise, then $|y(t;G,\tau,V(0))-y_\infty|\le|y(t;G,\tau,u)-y_\infty|\le 1/(1+t)$.

Their Appendix E outlines an argument for Conjecture 1 but leaves two inequalities, (E8) and (E9), unproved. The conjectures feed into bounds on an *upper-transition time* $T(r)$: the time after which a static network's prevalence is within $r$ of its steady state. These bounds decide whether a slowly changing network may be treated as quenched (frozen). Section 9 returns to that question.

## 1.2 What is new here, and what is not

| Result | Where | Status |
|-----------------------------------------------------|------------|------------------------------------------------|
| Conjecture 1, with temporal, heterogeneous-recovery and weighted directed versions | Thm 1, Cor 2, Thm 3 | new proof; the quantitative envelope is new as far as we found |
| Conjecture 2, with prefactor $1-\bar p$ and an exponential form | Thm 4 | new proof |
| Envelope attained by regular graphs; subcritical rate $1-\rho$ exact for all connected graphs | Prop 5, Prop 6 | elementary |
| Critical constant: $t\,x_i(t)$ converges for every node and initial state, with explicit limits (directed and heterogeneous too); $\kappa(G)\le1$ with equality iff regular; $r\,T(r)\to\kappa(G)$ | Thm 7, Cor 8 | the $1/t$ mechanism is the standard centre-manifold one; the global statement, constants and characterisation are new as far as we found |
| Node-wise: universal $\sqrt{n/(t(1+t))}$ bound; worst asymptotic node constant of order exactly $n^{1/3}$ | Prop 9 | new as far as we found |
| Directed unweighted bound and arithmetic-mean parameter are false | Prop 10 | counterexamples |
| SIR/SIRS infectious fractions, and the expected prevalence of the exact stochastic SIS process with common recovery rate, obey the envelopes | Prop 11, Prop 12 | direct corollaries (the stochastic step is Simon and Kiss 2018) |
| Mixed initial states never converge slower than both extremes; unconditional transition-time bounds for all $R_0$ (growth part exact on regular graphs, very conservative on strongly heterogeneous ones); numerical failures of the growth formula of Persoons et al. away from $R_0\to1$ | Prop 13, Thm 14, Sec 10 | new; the failures are numerical, not interval-certified |

**Known tools.** The function $h(s)=-s-\log(1-s)$ is the classical Volterra (Goh, Korobeinikov) Lyapunov function $\sigma-1-\log\sigma$, written in the susceptible fraction $\sigma=1-s$ at the disease-free state. That family underlies the global-stability theory of multigroup and network epidemic models (Shuai and van den Driessche 2013). The averages $H_w$ below are separable Lyapunov functions for a monotone system (Dirr et al. 2015). The diagonal weights of Lemma 2 are the classical diagonal stability of Metzler matrices (Rantzer 2015), used for SIS on directed graphs by Khanafer et al. (2016, Lemma A.1). Comparison for cooperative systems is classical (Smith 2008).

**What is new.** We did not find, in a bounded search recorded with this release, the Jensen step that turns the separable function into a *quantitative* Riccati rate valid at the critical point, where linear copositive Lyapunov functions degenerate. Below threshold, Mei et al. (2017, Theorem 4, which they attribute to Lajmanovich and Yorke) give exponential decay of a left-Perron-weighted average at rate $1-\rho$. That bound degenerates at $R_0=1$, the case the conjectures are about. We found no proof of either conjecture. That includes the authors' own October 2026 follow-up on temporal clique graphs (Persoons and Van Mieghem 2026), whose abstract does not mention the conjectures; we did not read its full text.

## 1.3 Roadmap

Section 2 fixes notation and three tools. Section 3 proves the envelope theorem and Section 4 the two conjectures. Section 5 asks how sharp the envelope is and introduces $\kappa(G)$. Section 6 compares nodes with averages and Section 7 shows which stronger statements fail. Section 8 extends the bounds to SIR/SIRS and to the stochastic model. Section 9 applies everything to the timescale-separation question, and Section 10 gives numerical illustrations. Section 11 lists limitations and open questions.

# 2. Setting and three tools

## 2.1 The general model

Let $B(t)\in\mathbb R^{n\times n}_{\ge0}$ and $\delta_i(t)>0$ be piecewise continuous on $[0,\infty)$, with finitely many pieces and finite one-sided limits on every bounded interval. In particular they are bounded on bounded intervals, and $\inf_i\delta_i>0$ there. Consider
$$
\dot x_i=-\delta_i(t)\,x_i+(1-x_i)\,(B(t)x)_i=:F_i(t,x). \tag{2.1}
$$
On each piece $F$ is continuous in $t$ and locally Lipschitz in $x$, so solutions exist and are unique. Concatenating over the pieces gives a continuous, piecewise $C^1$ solution, as in the switching model of Persoons et al. (their eq. 5).

The cube $[0,1]^n$ is forward invariant: at $x_i=0$, $\dot x_i\ge0$, and at $x_i=1$, $\dot x_i=-\delta_i<0$. The half-open cube $[0,1)^n$ is also forward invariant. On bounded intervals $x$ stays a positive distance from $\{x_i=1\}$, because $\frac{d}{dt}(1-x_i)\ge-(1-x_i)\,\|B(t)\|_\infty$ gives $1-x_i(t)\ge(1-x_i(0))\exp(-\int_0^t\|B\|_\infty)$. The system is *cooperative*: $\partial F_i/\partial x_j=(1-x_i)B_{ij}(t)\ge0$ for $j\neq i$ on the cube.

We write $y_w=\sum_iw_ix_i/\sum_iw_i$ for a weighted average ($y=y_w$ with $w=u$), $\rho(\cdot)$ for spectral radius, $\lambda_{\max}$ for the largest eigenvalue of a symmetric matrix, $\operatorname{diag}(w)$ for a diagonal matrix and $\circ$ for the entrywise product.

## 2.2 Comparison

**Lemma 1 (comparison).** Let $F^X$ and $F^Y$ be right-hand sides of the form (2.1) with data $(B^X,\delta^X)$ and $(B^Y,\delta^Y)$, where $B^X(t)\le B^Y(t)$ entrywise and $\delta^X(t)\ge\delta^Y(t)$. Let $X:[0,T]\to[0,1]^n$ be continuous and piecewise $C^1$ with $\dot X\le F^X(t,X)$ componentwise, using one-sided derivatives at break points. Let $Y:[0,T]\to[0,1]^n$ solve $\dot Y=F^Y(t,Y)$. If $X(0)\le Y(0)$, then $X(t)\le Y(t)$ for all $t\in[0,T]$.

This is the classical Kamke comparison theorem for cooperative systems (Smith 2008, Chapter 3). Appendix A gives a short proof in exactly this form. Two consequences are used throughout. First, solutions are monotone in the initial state, in $B$ (entrywise) and in $-\delta$. This is the first inequality in each conjecture (Persoons et al.'s Theorem 1, taken from earlier work). Second, a *subsolution* $X$ (with $\dot X\le F(X)$) stays below the true solution.

## 2.3 Diagonal weights

The bound we need is a quadratic-form inequality $x^{\mathsf T}WBx\le\theta\,x^{\mathsf T}W\Delta x$ with a positive diagonal $W$. For symmetric $B$ it is the Rayleigh quotient. For directed $B$ it needs weights.

**Lemma 2 (weights).** Let $S$ be symmetric with nonnegative off-diagonal entries, and let $Sv=0$ for some $v>0$. Then for all $x$,
$$x^{\mathsf T}Sx=-\tfrac12\sum_{i\neq j}S_{ij}v_iv_j\Big(\frac{x_i}{v_i}-\frac{x_j}{v_j}\Big)^2\le0 .$$
Consequently, let $B\ge0$ be irreducible, $\Delta=\operatorname{diag}(\delta_i)$ with $\delta_i>0$, and $\theta=\rho(\Delta^{-1}B)$. Let $v>0$ and $u>0$ satisfy $Bv=\theta\Delta v$ and $u^{\mathsf T}B=\theta u^{\mathsf T}\Delta$. Then $W=\operatorname{diag}(u_i/v_i)$ satisfies
$$x^{\mathsf T}WBx\le\theta\,x^{\mathsf T}W\Delta x\qquad(x\in\mathbb R^n). \tag{2.2}$$

*Proof.* With $a=x/v$, expanding gives $\sum_{i,j}S_{ij}v_iv_j(a_i-a_j)^2=2\sum_ia_i^2v_i(Sv)_i-2x^{\mathsf T}Sx$, and $Sv=0$. The diagonal terms vanish. For the consequence, $v$ and $\tilde u=\Delta u$ are right and left Perron vectors of the irreducible nonnegative matrix $\Delta^{-1}B$, with eigenvalue $\theta$. Put $M=B-\theta\Delta$ and $S=WM+M^{\mathsf T}W$. Then $S$ is symmetric with nonnegative off-diagonal entries, and $Sv=WMv+M^{\mathsf T}u=0$ because $Wv=u$. Hence $2x^{\mathsf T}WMx=x^{\mathsf T}Sx\le0$. $\square$

## 2.4 The Volterra-type functional

Let
$$h(s)=-s-\log(1-s),\qquad h'(s)=\frac{s}{1-s},\qquad 0\le s<1, \tag{2.3}$$
which is $\sigma-1-\log\sigma$ in the susceptible fraction $\sigma=1-s$. For $r\ge0$ and $d>0$ define $q_{r,d}(s)=d\,\dfrac{s^2}{1-s}-r\,s^2$.

**Lemma 3.** (a) $h$ is strictly increasing and strictly convex, maps $[0,1)$ onto $[0,\infty)$, and satisfies $s^2/2\le h(s)\le s^2/(2(1-s))$.
(b) $\Phi_{r,d}=q_{r,d}\circ h^{-1}$ is convex and $C^1$ on $[0,\infty)$, with $\Phi_{r,d}(0)=0$ and right derivative $\Phi_{r,d}'(0)=2(d-r)$. It is Lipschitz on every bounded interval, uniformly for $(r,d)$ in bounded sets.
(c) Along (2.1), $\frac{d}{dt}h(x_i)=-\delta_i\dfrac{x_i^2}{1-x_i}+x_i(Bx)_i$.

*Proof.* (a) On $(0,1)$, $h'>0$ and $h''=(1-s)^{-2}>0$. Also $\frac{d}{ds}(h-s^2/2)=s^2/(1-s)\ge0$ and $\frac{d}{ds}\big(s^2/(2(1-s))-h\big)=s^2/(2(1-s)^2)\ge0$, and all three functions vanish at $0$.

(b) For $u=h(s)$ with $s\in(0,1)$, $\Phi'(u)=q'(s)/h'(s)=d\frac{2-s}{1-s}-2r(1-s)$. Its $s$-derivative is $d/(1-s)^2+2r>0$, so $\Phi'$ is increasing in $u$ and $\Phi$ is convex. The right limit at $s=0$ is $2(d-r)$. The Lipschitz bounds hold because $\Phi'$ is bounded for $s$ in compact subsets of $[0,1)$.

(c) is the chain rule. All identities are checked symbolically (`exact_checks.py`, E1–E3, E13). $\square$

The key point is that $\Phi_{r,d}$ is convex *even when* $r>d$, where $q_{r,d}$ itself takes negative values. This is what lets one envelope cover the regimes below, at and above threshold.

# 3. The envelope theorem

In words: a single number, the average $H$ of $h(x_i)$, can only fall at a controlled rate. That rate makes it obey a scalar logistic (Riccati) equation, and Jensen's inequality turns $H$ back into the ordinary average.

**Theorem 1 (spectral Riccati envelope).** In the setting (2.1), let $W=\operatorname{diag}(w)$ with constant $w_i>0$. Let $r(t)\ge0$ and $d(t)>0$ be piecewise continuous, with finitely many pieces on bounded intervals, bounded on bounded intervals and with $\inf d>0$ there. Suppose
$$
x^{\mathsf T}WB(t)x-\sum_iw_i\delta_i(t)\frac{x_i^2}{1-x_i}\;\le\;-\sum_iw_i\,q_{r(t),d(t)}(x_i)\qquad\text{for all }x\in[0,1)^n,\ t\ge0. \tag{D}
$$
Let $H_w(x)=\sum_iw_ih(x_i)/\sum_iw_i$ and let $z$ solve
$$\dot z=(r(t)-d(t))\,z-r(t)\,z^2,\qquad z(0)=a. \tag{3.1}$$
(i) If $x(0)\in[0,1)^n$ and $a=h^{-1}(H_w(x(0)))$, then for all $t\ge0$,
$$y_w(t)\;\le\;h^{-1}\big(H_w(x(t))\big)\;\le\;z(t). \tag{3.2}$$
(ii) For any $x(0)\in[0,1]^n$, with $a=1$, $y_w(t)\le z(t)$ for $t\ge0$ and $H_w(x(t))\le h(z(t))$ for $t>0$.

*Proof.* (i) By Lemma 3(c) and (D), wherever the solution is differentiable,
$$\dot H_w=\frac{1}{\sum w}\Big[x^{\mathsf T}WBx-\sum_iw_i\delta_i\frac{x_i^2}{1-x_i}\Big]\le-\frac{\sum_iw_i\,\Phi_{r,d}(h(x_i))}{\sum_iw_i}\le-\Phi_{r,d}(H_w),$$
where the last step is Jensen's inequality for the convex $\Phi_{r,d}$ (Lemma 3b). The function $K=h(z)$ satisfies $\dot K=h'(z)\dot z=-q_{r,d}(z)=-\Phi_{r,d}(K)$ (checked as E3), with $K(0)=H_w(x(0))$. Also $z<1$ for all $t$, since $\dot z=-d<0$ at $z=1$. On $[0,T]$ both $H_w$ and $K$ stay in a bounded interval (Section 2.1), on which $\Phi$ is Lipschitz in its argument uniformly in $t$. Scalar comparison (Gronwall applied to $(H_w-K)^+$) gives $H_w\le K=h(z)$. Since $h$ is increasing, $h^{-1}(H_w)\le z$. Jensen's inequality for the convex increasing $h$ gives $h(y_w)\le H_w$, hence $y_w\le h^{-1}(H_w)$.

(ii) Let $x^{(\varepsilon)}$ be the solution from $(1-\varepsilon)u$, and $x^{u}$ the solution from $u$. Then $H_w(x^{(\varepsilon)}(0))=h(1-\varepsilon)$, and (i) applies with $a=1-\varepsilon$: $y_w(x^{(\varepsilon)})\le z_{1-\varepsilon}$ and $H_w(x^{(\varepsilon)})\le h(z_{1-\varepsilon})$. As $\varepsilon\downarrow0$, continuous dependence gives $x^{(\varepsilon)}\to x^{u}$ and $z_{1-\varepsilon}\to z_1$. For $t>0$, $x^{u}(t)\in[0,1)^n$: each coordinate leaves $1$ immediately, since $\dot x_i=-\delta_i<0$ there, and $[0,1)^n$ is forward invariant. $H_w$ is continuous on $[0,1)^n$. So $y_w(x^u)\le z_1$ and $H_w(x^u(t))\le h(z_1(t))$. By Lemma 1, any solution from $x(0)\le u$ satisfies $x\le x^{u}$, and $y_w$, $H_w$ are coordinatewise nondecreasing. $\square$

**Closed form.** For $a>0$, $w=1/z$ solves the linear equation $\dot w=(d-r)w+r$. So
$$
z(t)=\Big[a^{-1}e^{\int_0^t(d-r)}+\int_0^tr(s)\,e^{\int_s^t(d-r)}ds\Big]^{-1},\qquad
z(t)=\begin{cases}\dfrac{a\,e^{(r-d)t}}{1+\frac{ra}{r-d}(e^{(r-d)t}-1)},&r\neq d,\\[1.2ex]\dfrac{a}{1+rat},&r=d,\end{cases} \tag{3.3}
$$
the second for constant $r,d$ (checked as E4). Because $\dot z=-dz+rz(1-z)$, scalar comparison shows the envelope is increasing in $r$ and decreasing in $d$.

**Corollary 2 (when (D) holds).** Condition (D) holds in each of the following cases.

(a) $\delta_i(t)\ge d(t)$ for all $i$, and $x^{\mathsf T}WB(t)x\le r(t)\,x^{\mathsf T}Wx$ for all $x\ge0$.

(b) $x^{\mathsf T}WB(t)x\le\theta(t)\sum_iw_i\delta_i(t)x_i^2$ for all $x\ge0$, with $\theta(t)\in[0,1]$, $d(t)=\min_i\delta_i(t)$ and $r(t)=\theta(t)d(t)$.

In particular the envelope applies:

1. **Undirected networks** ($B(t)$ symmetric, $\delta_i\ge d$): $W=I$, $r(t)=\rho(B(t))$.
2. **Reversible networks** ($WB(t)$ symmetric for a fixed $W$): $r(t)=\rho(B(t))$, since $W^{1/2}BW^{-1/2}$ is symmetric and similar to $B$.
3. **Directed networks** (static irreducible $B$, $\delta_i=1$): $W=\operatorname{diag}(u_i/v_i)$ with $u,v$ the left and right Perron vectors of $B$, and $r=\rho(B)$, by Lemma 2.
4. **Heterogeneous recovery** (static $B$, $\Delta=\operatorname{diag}(\delta_i)$, $\theta=\rho(\Delta^{-1}B)\le1$). Take $W=I$ if $B$ is symmetric (then $\theta=\lambda_{\max}(\Delta^{-1/2}B\Delta^{-1/2})$), or $W$ from Lemma 2 if $B$ is only irreducible. Then $d=\min_i\delta_i$ and $r=\theta d$.
5. **Temporal networks**: any of the above, provided one fixed $W$ serves for all $t$. This always holds for undirected snapshots, with $W=I$.

*Proof.* (a) The left side of (D) is at most $r\sum w_ix_i^2-d\sum w_ix_i^2/(1-x_i)=-\sum w_iq_{r,d}(x_i)$. (b) The left side is at most $-\sum_iw_i\delta_i\big(\frac{x_i^2}{1-x_i}-\theta x_i^2\big)$. The bracket is nonnegative because $\theta\le1\le(1-x_i)^{-1}$, so replacing $\delta_i$ by $d=\min\delta_i$ can only increase the expression. The result is $-\sum w_iq_{\theta d,d}(x_i)$. The listed cases follow from the Rayleigh quotient and Lemma 2. $\square$

# 4. The two conjectures

## 4.1 Conjecture 1

**Theorem 3 (Conjecture 1).** Let $B\ge0$ be symmetric (the graph may be disconnected), $\delta=1$ and $\rho=\rho(B)\le1$. Then for every $x(0)\in[0,1]^n$ and $t\ge0$,
$$
y(t)\;\le\;y(t;u)\;\le\;\frac{1-\rho}{e^{(1-\rho)t}-\rho}\;\le\;\frac{1}{1+t}, \tag{4.1}
$$
where the middle expression is read as $1/(1+t)$ when $\rho=1$. In the notation of Persoons et al., for $\tau\le\tau_c(G)$,
$$y(t;G,\tau,V(0))\le y(t;G,\tau_c(G),u)\le y(t;K_N,\tau_c(K_N),u)=\frac1{1+t}.$$
The same bound $1/(1+t)$ holds in three further settings:

- **temporal networks** with symmetric snapshots satisfying $\rho(B(t))\le1$, with no condition on how fast or how often the network changes;
- **directed irreducible networks** with $\rho(B)\le1$, for the weighted average $y_w$ with $w_i=u_i/v_i$;
- **heterogeneous recovery**, for the plain average in the form $y(t)\le1/(1+\delta_{\min}t)$, when $B$ is symmetric and $\rho(\Delta^{-1}B)\le1$.

*Proof.* The first inequality in (4.1) is Lemma 1 ($x(0)\le u$), and the same lemma gives monotonicity in $\tau$. By Corollary 2(1), Theorem 1(ii) with $W=I$, $r=\rho$, $d=1$ and $a=1$ gives the middle bound (formula (3.3), E4d). The envelope increases with $r$, giving $1/(1+t)$ at $r=1$ (E4e). On $K_N$ at $\tau_c$, every row of $B$ sums to $1$, so the solution from $u$ is $x_i=1/(1+t)$ (Proposition 5). The three further settings are Corollary 2(5), (3) and (4), with $r=\theta d\le d$. $\square$

## 4.2 Conjecture 2

**Theorem 4 (Conjecture 2, sharpened).** Let $B\ge0$ be symmetric and irreducible with $\rho(B)>1$ and $\delta=1$. Let $p\in(0,1)^n$ be the endemic equilibrium, so $(Bp)_i=p_i/(1-p_i)$. Write $s_i=1-p_i$, $\bar s=1-\bar p$ and $\kappa_p=\min_ip_i/s_i>0$. If $p\le x(0)\le u$, then
$$
0\;\le\;y(t)-\bar p\;\le\;\frac{(1-\bar p)\,\kappa_p}{(1+\kappa_p)e^{\kappa_pt}-1}\;\le\;\frac{1-\bar p}{1+t}\;\le\;\frac1{1+t}. \tag{4.2}
$$
In particular $|y(t;G,\tau,V(0))-y_\infty|\le|y(t;G,\tau,u)-y_\infty|\le1/(1+t)$ for $V(0)\ge V_\infty$. On a regular graph started from $u$, the first upper bound is an equality.

For a symmetric $B$ on a disconnected graph, the algebraic bound $(1-\bar p)/(1+t)$ holds with $p$ the maximal equilibrium. On components with spectral radius at most 1 the equilibrium is $0$ and Theorem 3 applies. Componentwise, with components $G_k$ of sizes $n_k$: $y(t)-\bar p\le\frac1n\sum_kn_k\,\bar s_k\,E_k(t)$, where $E_k(t)=\kappa_k/((1+\kappa_k)e^{\kappa_kt}-1)$ on supercritical components ($\kappa_k$ the component's $\kappa_p$, $\bar s_k$ its mean of $1-p_i$) and $E_k(t)=1/(1+t)$, $\bar s_k=1$ on the others.

*Proof.* Lemma 1 and the equilibrium property give $p\le x(t)\le u$. Set $S=\operatorname{diag}(s)$ and $\zeta_i=(x_i-p_i)/s_i\in[0,1]$. Expanding (1.1) with $Bp=p/s$ gives (checked as E5)
$$\dot\zeta_i=-\frac{\zeta_i}{s_i}+(1-\zeta_i)(C\zeta)_i,\qquad C=BS. \tag{4.3}$$
This is system (2.1) with recovery rates $1/s_i$ and infection matrix $C$.

Now $SBp=p$ with $p>0$, so $\rho(SB)=1$, since a positive eigenvector of an irreducible nonnegative matrix belongs to the spectral radius. The symmetric matrix $S^{1/2}BS^{1/2}$ is similar to $SB$, so $S^{1/2}BS^{1/2}\preceq I$. Hence, with $W=S$,
$$x^{\mathsf T}WCx=(S^{1/2}x)^{\mathsf T}\big(S^{1/2}BS^{1/2}\big)(S^{1/2}x)\le x^{\mathsf T}Sx=x^{\mathsf T}Wx .$$
Corollary 2(a) applies with $r=1$ and $d=d_*=\min_i1/s_i=1+\kappa_p$. Theorem 1(ii) with $a=1$ gives $\sum_is_i\zeta_i/\sum_is_i\le(d_*-1)/(d_*e^{(d_*-1)t}-1)$ (formula (3.3), E16). Here $\sum_is_i\zeta_i=\sum_i(x_i-p_i)=n(y-\bar p)$ and $\sum_is_i=n\bar s$.

The second inequality in (4.2) is $\kappa(1+t)\le(1+\kappa)e^{\kappa t}-1$, true at $t=0$ with a larger derivative on the right. On a regular graph $p$ is constant, $\kappa_p=\rho-1$, and (4.3) from $\zeta(0)=u$ is exactly the envelope equation (E15). For disconnected graphs, apply the argument to each component and add. $\square$

*Remarks.* The prefactor $1-\bar p$ makes (4.2) strictly stronger than the conjecture. The bound concerns convergence *from above*. From below no bound free of the initial state is possible, since the disease-free state is an equilibrium (Persoons et al., Lemma 1). Section 9 treats growth under the paper's lower bound $V(0)\ge ru$.

# 5. How sharp is the envelope?

In words: regular networks attain the envelope. Below threshold its exponential rate is right for every connected network. At the threshold the long-run constant is $\kappa(G)$, which equals 1 only for regular networks.

**Proposition 5 (attainment).** If $B(t)u=r(t)u$ (constant row sums), $\delta_i(t)=d(t)$ and $x(0)=au$, then $x(t)=z(t)u$ with $z$ from (3.1), so (3.2) holds with equality. Hence, among bounds depending only on $(r(\cdot),d(\cdot),a)$, the envelope is the smallest possible. In particular the constant $1$ in Conjecture 1 cannot be lowered: every regular graph at its threshold, started fully infected, has $y(t)=1/(1+t)$.

*Proof.* Substitute $x=zu$ into (2.1). $\square$

**Proposition 6 (exact subcritical rate).** Let $B$ be symmetric and irreducible with $\rho=\rho(B)<1$, $\delta=1$ and $x(0)\neq0$. Then $L=\lim_{t\to\infty}e^{(1-\rho)t}y(t)$ exists and is positive. So the exponent $1-\rho$ in (4.1) is optimal for every connected graph, not only regular ones.

*Proof.* Let $v>0$ be the unit Perron vector, $c=v^{\mathsf T}x$ and $v_{\min}=\min_iv_i$. Since $x\in[0,1]^n$, $\|x\|_\infty\le\|x\|_2$, so $\|x\circ Bx\|_2\le\|x\|_\infty\|Bx\|_2\le\rho\|x\|_2^2\le\rho\|x\|_1^2\le\rho\,c^2/v_{\min}^2$, using $c\ge v_{\min}\|x\|_1$.

Now $\dot c=(\rho-1)c-v^{\mathsf T}(x\circ Bx)$, so $g=e^{(1-\rho)t}c$ satisfies $-\rho v_{\min}^{-2}e^{-(1-\rho)t}g^2\le\dot g\le0$. Thus $g$ is nonincreasing, and $1/g$ increases by at most $\rho v_{\min}^{-2}\int_0^\infty e^{-(1-\rho)s}ds<\infty$, so $g\to L_c>0$.

Write $x=cv+w$ with $w\perp v$. Then $\dot w=(B-I)w-P(x\circ Bx)$, where $P=I-vv^{\mathsf T}$ commutes with $B$ and $\|e^{(B-I)s}P\|=e^{-(1-\lambda_2)s}$. Here $1-\lambda_2=(1-\rho)+\gamma$, with spectral gap $\gamma=\rho-\lambda_2>0$. Variation of constants with $\|x\|_2^2=O(e^{-2(1-\rho)t})$ (from (4.1)) gives $e^{(1-\rho)t}\|w(t)\|\le e^{-\gamma t}\|w(0)\|+C\int_0^te^{-\gamma(t-s)}e^{-(1-\rho)s}ds\to0$. Hence $e^{(1-\rho)t}y=e^{(1-\rho)t}(c\sum_iv_i+\sum_iw_i)/n\to L_c\sum_iv_i/n>0$. $\square$

**Theorem 7 (the critical constant).** Let $B\ge0$ be irreducible and let $\Delta=\operatorname{diag}(\delta_i)$ be constant with $\delta_i>0$ and $\rho(\Delta^{-1}B)=1$, the threshold. Let $v,u>0$ satisfy $Bv=\Delta v$, $u^{\mathsf T}B=u^{\mathsf T}\Delta$ and $u^{\mathsf T}v=1$. Then for every $x(0)\in[0,1]^n\setminus\{0\}$ and every node $i$,
$$
\lim_{t\to\infty}t\,x_i(t)=\frac{v_i}{\sigma},\qquad \sigma=\sum_ju_j\delta_jv_j^2,\qquad |t\,x_i(t)-v_i/\sigma|=O\Big(\frac{\log t}{t}\Big). \tag{5.1}
$$
For symmetric $B$ the left and right vectors are proportional, and we may choose $u=v$. In the homogeneous symmetric case ($\delta=1$, $B=B^{\mathsf T}$), the population average satisfies
$$
\lim_{t\to\infty}t\,y(t)=\kappa(B):=\frac{\big(\frac1n\sum_jv_j\big)\big(\frac1n\sum_jv_j^2\big)}{\frac1n\sum_jv_j^3}\in(0,1], \tag{5.2}
$$
with $\kappa(B)=1$ if and only if every row sum of $B$ equals 1. For $B=\tau_cA$ this means the graph is regular.

The limit does not depend on the initial state, although the $O(\log t/t)$ constant does, through $u^{\mathsf T}x(0)$. Note the contrast: Theorem 3 says no network does worse than $1/(1+t)$, while (5.2) says every non-regular network does *strictly better* in the long run, by the factor $\kappa(G)$.

*Relation to known theory.* At the threshold the Jacobian at $0$ has a simple zero eigenvalue with eigenvector $v$. The $1/t$ law with coefficient $\sigma$ is exactly what a one-dimensional centre-manifold reduction $\dot c=-\sigma c^2+O(c^3)$ predicts (Carr 1981). The same cubic Perron sum appears in the near-threshold steady state of Van Mieghem (2012) and Prasse and Van Mieghem (2020, Theorem 1), as $R_0\downarrow1$. What Theorem 7 adds is a global, elementary statement: every nonzero initial state, every node, directed and heterogeneous networks, an explicit rate, and the characterisation $\kappa\le1$ with equality only for regular graphs. A local reduction alone does not give these.

*Proof.* Put $M=B-\Delta$. This matrix is irreducible and Metzler, with $Mv=0$ and $u^{\mathsf T}M=0$. By Perron–Frobenius, $0$ is a simple eigenvalue and all other eigenvalues have negative real part. The projector $Q=I-vu^{\mathsf T}$ commutes with $M$, and $\|e^{Ms}Q\|\le Ke^{-\gamma s}$ for some $K\ge1$, $\gamma>0$. Write $c=u^{\mathsf T}x$, $w=Qx=x-cv$, $u_{\min}=\min_iu_i$, and use $\|\cdot\|=\|\cdot\|_2$. Constants $C,C_1,\dots$ depend on $B$, $\Delta$ and $x(0)$ only.

*Step 1: $c\asymp1/t$.* Since $u^{\mathsf T}M=0$, $\dot c=-u^{\mathsf T}(x\circ Bx)$. As $0\le u^{\mathsf T}(x\circ Bx)\le\|u\|_\infty\|x\|_\infty\|B\|_1\|x\|_1\le\|u\|_\infty\|B\|_1\|x\|_1^2$ and $\|x\|_1\le c/u_{\min}$, we get $\dot c\ge-\beta c^2$ with $\beta=\|u\|_\infty\|B\|_1/u_{\min}^2$. Since $c(0)>0$, $c(t)\ge c(0)/(1+\beta c(0)t)\ge c_1/(1+t)$ with $c_1=\min\{c(0),1/\beta\}$. For the upper bound, Lemma 2 applied to $M$ gives the weights $W=\operatorname{diag}(u/v)$. By Corollary 2(4) and Theorem 1(ii), $y_w\le1/(1+\delta_{\min}t)$. Hence $c=\sum_iv_iw_ix_i\le\max_iv_i\,(\sum_iw_i)\,y_w\le C_1/(1+t)$ and $\|x\|_1\le(\sum_iw_i/\min_iw_i)\,y_w\le C_1'/(1+t)$.

*Step 2: $\|w\|=O(t^{-2})$.* Here $\dot w=Mw-Q(x\circ Bx)$ and $\|x\circ Bx\|\le\|B\|\,\|x\|_1^2\le C_2/(1+t)^2$. Splitting the variation-of-constants integral at $t/2$ gives
$$\|w(t)\|\le Ke^{-\gamma t}\|w(0)\|+KC_2\|Q\|\Big(\tfrac t2e^{-\gamma t/2}+\tfrac{1}{\gamma(1+t/2)^2}\Big)\le\frac{C_w}{(1+t)^2}.$$

*Step 3: the reduced equation.* Since $Bv=\Delta v$,
$$\dot c=-u^{\mathsf T}\big((cv+w)\circ B(cv+w)\big)=-\sigma c^2+E,\qquad|E|\le C_3\,(c\|w\|+\|w\|^2).$$
For $g=1/c$, $\dot g=\sigma-E/c^2$, with $|E|/c^2\le C_3(\|w\|/c+(\|w\|/c)^2)\le C_4/(1+t)$ by Steps 1–2. Integrating, $|g(t)-g(0)-\sigma t|\le C_4\log(1+t)$, so $t\,c(t)=1/\sigma+O(\log t/t)$.

*Step 4: conclusion.* $t\,x_i=t\,c\,v_i+t\,w_i\to v_i/\sigma$, with the stated rate. (Here $w_i$ in $t\,w_i$ is the component of $w=Qx$; in Step 1 $w_i$ denoted the weights.) In the homogeneous symmetric case take $u=v$ with $\|v\|_2=1$. Then $\sigma=\sum v_j^3$, and averaging gives (5.2), which is invariant under rescaling $v$. Finally (E14),
$$n\sum v_j^3-\sum v_j\sum v_j^2=\tfrac12\sum_{i,j}(v_i-v_j)(v_i^2-v_j^2)\ge0,$$
which is Chebyshev's sum inequality. Each term is nonnegative, and all vanish iff $v$ is constant, that is, iff every row sum of $B$ equals 1. $\square$

**Corollary 8 (transition time at the threshold).** Under (5.2), let $T(r)$ be the upper-transition time of Persoons et al. (their eq. 8): the time after which $y(t)\le r$ for all $V(0)\in[r,1]^n$. Then $T(r)\le1/r-1$ and $r\,T(r)\to\kappa(B)$ as $r\downarrow0$.

*Proof.* By Lemma 1 the worst case is $V(0)=u$. Since $x(h;u)\le u$, Lemma 1 and autonomy show $y(t;u)$ is nonincreasing, so $T(r)$ solves $y(T(r);u)=r$ and $T(r)\to\infty$. Then $rT(r)=T\,y(T;u)\to\kappa$, and $T(r)\le1/r-1$ by Theorem 3. $\square$

**Why $\kappa$ can be small: localisation.** $\kappa(G)$ is small when the Perron vector is *localised* on a few nodes (Goltsev et al. 2012; Pastor-Satorras and Castellano 2016). Then $\frac1n\sum v_j^3$ is large relative to $\frac1n\sum v_j\cdot\frac1n\sum v_j^2$. At $R_0=1$ only that core is near-critical, the rest of the network is effectively subcritical, and the *average* collapses quickly while the core persists. Localisation grows with size for Barabási–Albert and random geometric graphs, so their $\kappa$ decreases with $n$ (Table 1b). Small $\kappa$ therefore signals a persistent core, not easy elimination.

# 6. Individual nodes versus the population average

**Proposition 9.** Let $B$ be symmetric with $\rho(B)\le1$ and $\delta=1$.
(a) Every node satisfies $x_i(t)\le\min\{1,\sqrt{n/(t(1+t))}\}$ for all $t>0$ and all initial states.
(b) Star $K_{1,m}$ at threshold, $x(0)\neq0$: $t\,x_{\rm hub}(t)\to2\sqrt m/(\sqrt m+1)$, which exceeds 1 for $m>1$, while $\kappa=2\sqrt m/(m+1)$.
(c) Let $1\le d\le N-1$ with $Nd$ even. Join one hub to every vertex of a $d$-regular graph on $N$ vertices, set $B=A/\rho(A)$, and take $x(0)\neq0$. Then $t\,x_{\rm hub}(t)\to R_N$ with
$$R_N\;\ge\;\frac{d^2/N}{1+(d+N/d)^3/N^2}. \tag{6.1}$$
Choosing $N^{2/3}/2\le d\le N^{2/3}$ (with $d$ even if $N$ is odd; such $d$ exists for $N\ge 64$) gives $R_N\ge N^{1/3}/112$.
(d) If $B$ is symmetric and irreducible with $\rho(B)=1$ (and $\delta=1$), then $\max_i\lim_{t\to\infty}t\,x_i(t)\le n^{1/3}$.

Together, (c) and (d) show that the worst asymptotic node constant over $n$-node networks is of order exactly $n^{1/3}$. So no constant independent of the network gives $x_i(t)\le C/(1+t)$ for all nodes.

*Proof.* (a) Theorem 1(ii) with $W=I$ and $r=\rho\le1$ gives $H(x(t))\le h(z(t))\le h(1/(1+t))$. Since $h\ge0$, $h(x_i)\le nH$. Lemma 3(a) gives $x_i^2/2\le h(x_i)\le n\,h(1/(1+t))\le n/(2t(1+t))$.

(b) The Perron vector is $(\sqrt m,1,\dots,1)$. Insert it in (5.1) and (5.2) (E9).

(c) A $d$-regular graph on $N$ vertices exists when $Nd$ is even, for example a circulant. By symmetry the Perron vector is $(1,\varepsilon,\dots,\varepsilon)$, with $\lambda=N\varepsilon$ and $\lambda\varepsilon=1+d\varepsilon$. So $\lambda^2-d\lambda-N=0$ and $\lambda-d=N/\lambda\le N/d$ (E10). By (5.1), $R_N=(1+N\varepsilon^2)/(1+N\varepsilon^3)$, where $N\varepsilon^2=\lambda^2/N\ge d^2/N$ and $N\varepsilon^3=\lambda^3/N^2\le(d+N/d)^3/N^2$. For the stated $d$, $d^2/N\ge N^{1/3}/4$ and $(d+N/d)^3/N^2\le(3N^{2/3})^3/N^2=27$.

(d) Scale $v$ so that $\max_iv_i=1$. By (5.1) the node constants are $v_i\sum v_j^2/\sum v_j^3\le\sum v_j^2/\sum v_j^3$. Hölder's inequality gives $\sum v_j^2\le n^{1/3}(\sum v_j^3)^{2/3}$, so the constant is at most $n^{1/3}(\sum v_j^3)^{-1/3}\le n^{1/3}$, because $\sum v_j^3\ge1$ (E18). $\square$

For $N=2000$ and $d=158$, the hub's long-run constant is $6.93$ (Section 10), so the hub carries about seven times the population envelope. The order of the worst node *uniformly in time and initial state* lies between $n^{1/3}$ (from (c)) and $n^{1/2}$ (from (a)); the exact order is open.

# 7. Stronger statements that are false

**Proposition 10.** (a) *Directed networks need weights.* Take $B=\begin{pmatrix}0&4\\1/16&0\end{pmatrix}$, so $\rho(B)=1/2$. From $x(0)=u$, $y(0)=1$, $\dot y(0)=-1$ and $\ddot y(0)=97/32>2$, so $y(t)>1/(1+t)$ for small $t>0$. With the Perron weights $w=(1/8,8)$ of Corollary 2(3), $y_w(t)\le(1/2)/(e^{t/2}-1/2)$ holds.
(b) *The parameter $h^{-1}(H(x(0)))$ cannot be replaced by the initial arithmetic mean.* For the star $K_{1,2}$ with $B=A/\sqrt2$ ($\rho=1$) and $x(0)=(1/2,0,0)$ (centre first), $\dot y(0)=(\sqrt2-1)/6>0$. So $y(t)\le y(0)/(1+y(0)t)$ fails for small $t$.

*Proof.* Direct computation (E11, E12). The weighted inequality follows from Lemma 2 and Theorem 1. $\square$

Numerically, on 25 random strongly connected directed graphs at $\rho=1$, the unweighted average exceeds $1/(1+t)$ by up to $0.057$, while the weighted average never does (`numerical_checks.py`, N2, C1). The identity weight $W=I$ fails the quadratic-form inequality for such graphs (control C1c).

# 8. SIR, SIRS and the stochastic SIS model

**Proposition 11 (SIR/SIRS).** Consider the network SIRS model
$$\dot s_i=-s_i(Bx)_i+\xi_i\varrho_i,\qquad \dot x_i=-\delta_ix_i+s_i(Bx)_i,\qquad \dot\varrho_i=\delta_ix_i-\xi_i\varrho_i,$$
with $s_i+x_i+\varrho_i=1$ and waning rates $\xi_i\ge0$ ($\xi=0$ is SIR). The infectious fractions satisfy $x(t)\le X(t)$, where $X$ solves the SIS system (2.1) with the same $B,\delta$ and $X(0)=x(0)$. Hence every bound of Theorems 1 and 3 and Corollary 2 applies to the infectious fractions.

*Proof.* The simplex constraints are invariant, so $s_i\le1-x_i$. Since $(Bx)_i\ge0$, $\dot x_i\le-\delta_ix_i+(1-x_i)(Bx)_i$, and Lemma 1 applies. $\square$

**Proposition 12 (exact stochastic SIS).** Consider the Markovian SIS process on a weighted, possibly directed network, with infection rates $B_{ij}$ ($i\neq j$) and a common curing rate $\delta$. Suppose the nodal states at time 0 are independent, for example a deterministic initial state. Let $\mathbb E[X(t)]$ be the vector of infection probabilities. Then $\mathbb E[X(t)]\le x(t)$, where $x$ solves NIMFA (2.1) with $x(0)=\mathbb E[X(0)]$. Consequently, for symmetric $B$ with $\rho(B)\le\delta$, the expected prevalence satisfies $\mathbb E[Y(t)]\le1/(1+\delta t)$, and the weighted directed version of Theorem 3 transfers in the same way. (We do not claim the heterogeneous-recovery case here: the cited correlation result assumes a common recovery rate.)

*Proof.* Independent initial states give $\mathbb E[S_i]\mathbb E[I_j]-\mathbb E[S_iI_j]=0$ at $t=0$. By Simon and Kiss (2018, Theorem 3), $\mathbb E[S_i]\mathbb E[I_j]-\mathbb E[S_iI_j]\ge0$ for all $t$, that is, $S$–$I$ pairs have nonpositive covariance (equivalently, infection states are non-negatively correlated; see also Cator and Van Mieghem 2014). Their Corollary 1 then gives $\mathbb E[X_i(t)]\le x_i(t)$. Apply Theorem 3 to $x$ and use linearity of expectation. $\square$

These bounds are far from tight for the stochastic process: on a finite network it dies out exponentially fast, so its expected prevalence eventually decays exponentially even at $R_0=1$. The sharpness results of Section 5 are specific to the deterministic model. The value of Proposition 12 is that the universal limit $1/(1+\delta t)$ is rigorous for the exact process as well.

# 9. Consequences for timescale separation

Persoons et al. define the upper-transition time $T(r)$ as the first time after which $|y(t)-y_\infty|\le r$ for all $V(0)\in[r,1]^n$. Under timescale separation, a network that changes every $\Delta t\ge T(r)$ can be treated as a sequence of static networks. Their bounds use the two conjectures, an assumption on mixed initial states, and a growth estimate derived as $R_0\downarrow1$. We revisit each.

**Proposition 13 (mixed initial states).** For any nonnegative $B$, positive $\delta$ and $V(0)\in[r,1]^n$,
$$y(t;ru)-y_\infty\le y(t;V(0))-y_\infty\le y(t;u)-y_\infty,$$
so $|y(t;V(0))-y_\infty|\le\max\{|y(t;ru)-y_\infty|,|y(t;u)-y_\infty|\}$.

*Proof.* Lemma 1 with $ru\le V(0)\le u$. $\square$

Their Assumption 1 reads, literally, that a mixed initial condition "converges slower to the steady state prevalence $y_\infty$ than the same process with the initial condition either $V(0)=u$ or $V(0)=ru$". Their surrounding text explains the intended property: mixed cases "converge faster than at least one of the extremal cases", so the extremal cases bound the error. Proposition 13 proves this intended property, which is all their bounds use, for every network. The literal wording, slower than both extremes, is never true.

**Theorem 14 (unconditional transition-time bounds).**

(a) *$R_0\le1$.* Let $B$ be symmetric with $\rho=\rho(B)\le1$. Then $y(t)\le r$ for all $V(0)\in[0,1]^n$ and all $t\ge T_\le(r)$, where
$$T_\le(r)=\frac{\log\big(1+(1-\rho)(1/r-1)\big)}{1-\rho}\quad(\rho<1),\qquad T_\le(r)=\frac1r-1\quad(\rho=1). \tag{9.1}$$
$T_\le(r)$ is at most $\min\{(1-r)/r,\ \log(1/r)/(1-\rho)\}$. It makes their eq. (11) unconditional, improves their Lemma 2, and needs no crossover point between the two regimes.

(b) *$R_0>1$.* Let $B$ be symmetric and irreducible with $R_0=\rho(B)>1$ and endemic state $p$. Let $v$ be the Perron vector scaled to $\max_iv_i=1$, and set $c_*=1-1/R_0$, $\kappa_p=\min_ip_i/(1-p_i)$ and $m_1=\tfrac{c_*}2\min_iv_i/p_i$ (always $\le1/2$). Define
$$
\begin{aligned}
T_D(r)&=\frac1{\kappa_p}\log^+\frac{1+\kappa_p\bar s/r}{1+\kappa_p},\qquad
T_G(r)=\frac{1}{R_0-1}\log^+\frac{c_*-r}{r}+t_2,\\
t_2&=\begin{cases}0,&\bar p(1-m_1)\le r,\\[0.5ex]\dfrac1{\kappa_p}\log\dfrac{(1-m_1)(\bar p-r)}{m_1r},&\text{otherwise,}\end{cases}
\end{aligned} \tag{9.2}
$$
with $\log^+a=\max\{0,\log a\}$ and the first term of $T_G$ read as $0$ when $r\ge c_*$. Then $|y(t)-y_\infty|\le r$ for all $V(0)\in[r,1]^n$ and all $t\ge\max\{T_D(r),T_G(r)\}$. The growth part $T_G$ remains valid for directed irreducible $B$, with $v$ the right Perron vector.

(c) *Exactness and asymptotics.* On a regular graph, $T_G(r)=\frac{2}{R_0-1}\log\frac{p-r}{r}$ for $0<r<p/2$ and $T_G(r)=0$ for $r\ge p/2$, which is exactly the growth time from $ru$. For $r<p/2$ this coincides with the first form of Persoons et al.'s Lemma 5, $\frac{2}{R_0-1}\log\frac{\max_ip_i-r}{r}$. On regular graphs, therefore, their growth formula is exact and now proved for every $R_0>1$. In general, $T_G(r)\le\big(\frac1{R_0-1}+\frac1{\kappa_p}\big)\log\frac1r+O(1)$ as $r\to0$, against the lower bound $\log\frac{y_\infty-r}{r}/(R_0-1)$ of their Lemma 6.

*Proof.* (a) By Theorem 3, $y(t)\le z(t)$ with $z$ decreasing, and $z(T_\le)=r$ (E4f). The comparison with $(1-r)/r$ and $\log(1/r)/(1-\rho)$ uses $\log(1+a)\le a$ and $1+(1-\rho)(1/r-1)\le1/r$.

(b) *Decay side.* By Proposition 13 and Theorem 4, $y(t;V(0))-y_\infty\le y(t;u)-\bar p\le\bar s\kappa_p/((1+\kappa_p)e^{\kappa_pt}-1)$. This bound is decreasing in $t$ and at most $r$ for $t\ge T_D$.

*Growth side, phase 1.* (For directed irreducible $B$, $x(t;ru)\to p$ by Lajmanovich and Yorke 1976.) Let $c$ solve $\dot c=c\,(R_0-1-R_0c)$ with $c(0)=r$. If $r<c_*$, then $c$ increases to $c_*$, and $cv\in[0,1]^n$. With $Bv=R_0v$ and $v_i\le1$,
$$F_i(cv)=c\,v_i\,(R_0-1-R_0cv_i)\ge\dot c\,v_i$$
(E7). So $cv$ is a subsolution, and since $c(0)v\le ru$, Lemma 1 gives $x(t;ru)\ge c(t)v$. As $t\to\infty$, $x(t;ru)\to p$ and $c\to c_*$, so $c_*v\le p$ and $m_1\le1/2$. Let $t_1=\frac1{R_0-1}\log^+\frac{c_*-r}{r}$. If $r<c_*/2$, then $c(t_1)=c_*/2$ (E17). If $r\ge c_*/2$, then $t_1=0$ and $x(0)\ge ru\ge\frac{c_*}{2}v$. Either way, $x(t_1;ru)\ge\frac{c_*}2v\ge m_1p$.

*Phase 2.* Let $m$ solve $\dot m=\kappa_pm(1-m)$ with $m(0)=m_1$. Using $Bp=p/s$, $F_i(mp)=m(1-m)p_i^2/s_i\ge\dot m\,p_i$ (E6). So $m(\cdot)p$ is a subsolution, and $x(t_1+\tau;ru)\ge m(\tau)p$. Hence $y_\infty-y(t;ru)\le\bar p\,(1-m(t-t_1))$, which is decreasing and at most $r$ for $t-t_1\ge t_2$ (E17). Proposition 13 combines the two sides. Neither phase used symmetry.

(c) On a regular graph, $p=c_*u$, $v=u$, $\kappa_p=R_0-1$ and $m_1=1/2$. If $r<p/2$, both phases are exact logistic solutions, and each term of (9.2) equals $\frac1{R_0-1}\log\frac{p-r}{r}$. The solution from $ru$ is the scalar logistic $x/(p-x)=\frac{r}{p-r}e^{(R_0-1)t}$, which reaches $p-r$ at exactly the total time. If $r\ge p/2$, then $p-r\le r$, so the growth time and $T_G$ are both $0$. The asymptotic statement follows from (9.2). `numerical_checks.py` (N12) checks (b) end to end on random graphs for $r$ on both sides of $c_*/2$, and (c) on a regular graph for $r$ below $p/2$, between $p/2$ and $p$, at $p$ and above $p$. $\square$

**What this says about the motivating question.** Four points follow.

1. *At and below the threshold, the upper bound depends only on $R_0$, but the true behaviour depends on shape.* The envelope depends only on $R_0$. This is consistent with the observation of Persoons et al. that $T(r)$ is "almost entirely determined by $R_0$", and turns it into a theorem as an *upper* bound. At $R_0=1$, however, Corollary 8 gives $T(r)\approx\kappa(G)/r$. On the synthetic networks of Section 10, $T(10^{-3})$ ranges from 39 (random geometric) to 922 (small world). For their $N=50$ Erdős–Rényi ensembles we find a median $\kappa$ of $0.96$ (10% quantile $0.81$; Table 1b), which is consistent with the dependence not being visible there. This is a population-average effect: asymptotically the worst node can exceed the envelope by a factor of order $n^{1/3}$, but no more (Proposition 9).
2. *Temporal networks need no timescale separation for die-out bounds.* By Theorem 3, if every snapshot has $R_0\le1$, then $y(t)\le1/(1+t)$ across all switches, whatever the inter-update times. More generally, by Corollary 2(5), $y(t)\le z(t)$ with $r(t)=R_0(t)$ in (3.3). This envelope is computable from the sequence of snapshot reproduction numbers alone, and holds in the annealed, intermediate and quenched regimes alike.
3. *Their growth formula is exact for regular graphs but is not a general bound away from $R_0\to1$.* Persoons et al. derive their eq. (20) in the limit $R_0\downarrow1$ and suggest, from simulations on $N=50$ Erdős–Rényi graphs, that it holds for all $R_0>1$. In high-precision integration with two independent solvers (Table 3), it fails within its stated hypotheses ($r\le\min_ip_i$, $V(0)=ru$). Examples: a 2-D random geometric graph at $R_0=2$ (gap $82\,r$ at the formula's time), a 16-node lollipop graph at $R_0=3$ ($67\,r$), and the same 1000-node Erdős–Rényi graph at $R_0=3$ and $5$ ($1.9\,r$ and $3.6\,r$). The cause is leaving the near-threshold regime in which (20) was derived, which heterogeneity amplifies. These are numerical findings, not interval-certified proofs.
4. *Their combined bound (23) for $R_0>1$ survives our tests.* It takes the larger of $(1-r)/r$ and (20), and $(1-r)/r$ dominates in every case we ran. Its decay part is now proved (Theorem 4 gives the sharper $\bar s/r-1$). Its growth part is not established by our results. Our rigorous replacement (9.2) is exact on regular graphs, close on near-regular ones, and orders of magnitude too large when some node has an endemic level near zero ($\kappa_p\approx0$; Table 2). A sharp general growth bound, perhaps the natural conjecture that growth from $ru$ takes at most $(1-r)/r$, remains open.

# 10. Numerical illustration on standard synthetic networks

We integrated (1.1) on eight standard synthetic network models, using fixed seeds and the largest connected component (`illustrations.py`, `eq20_check.py`). All rigorous bounds evaluated held. Integration used LSODA with relative tolerance $10^{-10}$ and absolute tolerance $10^{-15}$. Table 3 also uses Radau with the analytic Jacobian. These are illustrations, not proofs.

**Table 1.** At the threshold $R_0=1$, started fully infected (seed 11). $\kappa(G)$ from (5.2); $t\,y(t)$ at $t=10^4$; the largest node constant $\max_iv_i\sum v^2/\sum v^3$; and the transition time $T(10^{-3})$ (Corollary 8 gives $T\approx\kappa/r$, the envelope $T\le999$).

| Graph | $n$ | $\kappa(G)$ | $t\,y(t)$ at $t=10^4$ | max node constant | $T(10^{-3})$ |
|----------------------------------------------|-------:|---------:|-----------:|-----------:|---------:|
| Watts–Strogatz small world ($k=10$, $p=0.1$) | 1000 | 0.920 | 0.920 | 1.62 | 922 |
| Stochastic block model (4 communities) | 1000 | 0.822 | 0.822 | 2.05 | 822 |
| Erdős–Rényi (mean degree 10) | 1000 | 0.809 | 0.809 | 1.79 | 810 |
| Hub joined to 158-regular graph (denser) | 2001 | 0.592 | 0.593 | 6.93 | 603 |
| Star, 199 leaves | 200 | 0.141 | 0.141 | 1.87 | 143 |
| Holme–Kim (scale-free, clustered) | 1000 | 0.121 | 0.121 | 2.68 | 127 |
| Barabási–Albert ($m=2$) | 1000 | 0.050 | 0.050 | 2.16 | 57 |
| Random geometric (2-D, radius 0.06) | 1000 | 0.035 | 0.035 | 1.28 | 39 |

**Table 1b.** Seed and size dependence of $\kappa$, with localisation. Median and range of $\kappa$ over 10 further seeds at $n=1000$; median $n\cdot\mathrm{IPR}$, where $\mathrm{IPR}=\sum v_j^4$ for unit $v$ and larger values mean stronger localisation; and $\kappa$ for one seed at $n=250,500,1000,2000$.

| Family | $\kappa$ median (range), 10 seeds | median $n\cdot\mathrm{IPR}$ | $\kappa$ at $n=250$ | 500 | 1000 | 2000 |
|-----------------------------------|----------------:|---------:|-------:|-------:|-------:|-------:|
| Erdős–Rényi (mean degree 10) | 0.799 (0.775–0.807) | 1.5 | 0.821 | 0.812 | 0.812 | 0.804 |
| Watts–Strogatz small world | 0.908 (0.879–0.933) | 1.2 | 0.947 | 0.934 | 0.916 | 0.903 |
| Random geometric (2-D) | 0.036 (0.029–0.050) | 37.2 | 0.187 | 0.151 | 0.065 | 0.025 |
| Holme–Kim (scale-free, clustered) | 0.126 (0.112–0.162) | 41.2 | 0.252 | 0.183 | 0.128 | 0.084 |
| Barabási–Albert ($m=2$) | 0.079 (0.060–0.106) | 79.0 | 0.195 | 0.136 | 0.087 | 0.053 |

For comparison, Persoons et al.'s ensembles ($N=50$ Erdős–Rényi graphs, link density uniform on $(0,1)$; 200 connected samples) have median $\kappa=0.96$, 10% quantile 0.81 and minimum 0.57.

**Table 2.** Above threshold, $r=10^{-3}$, seed 11. Each cell gives the observed time, then our rigorous bound from (9.2) after the slash: decay is from $u$, growth is from $ru$. The eq. (20) columns give the growth formula of Persoons et al. (their eq. 20, first form). $^\dagger$: observed growth exceeds the formula within its hypotheses. $^\ddagger$: it exceeds the formula, but $r>\min_ip_i$, so Lemma 5's hypotheses fail.

| Graph | $\min p_i$ at $R_0=1.2$ | $R_0=1.2$ decay | $R_0=1.2$ growth | eq. (20) | $R_0=2$ decay | $R_0=2$ growth | eq. (20) |
|--------------------|---------:|-------------:|-------------------:|---------:|-------------:|-------------:|---------:|
| Watts–Strogatz | 1.0e-1 | 25.1 / 39 | 51.5 / 81 | 54.4 | 5.7 / 8.1 | 12.6 / 18 | 12.8 |
| Stochastic block model | 3.3e-2 | 24.9 / 98 | 50.7 / 212 | 56.3 | 6.1 / 23 | 12.9 / 47 | 13.0 |
| Erdős–Rényi | 2.7e-2 | 24.8 / 115 | 50.5 / 251 | 55.7 | 6.2 / 27 | 12.9 / 56 | 12.9 |
| Hub + 158-regular | 1.3e-1 | 28.4 / 31 | 53.6 / 77 | 64.8 | 6.2 / 6.2 | 13.2 / 17 | 13.6 |
| Star | 2.4e-2 | 18.4 / 130 | 35.0 / 201 | 56.6 | 7.8 / 43 | 10.9 / 65 | 13.2 |
| Holme–Kim | 2.2e-3 | 23.4 / 516 | 41.1 / 2571 | 59.2 | 11.1 / 87 | 16.2 / 237 | 13.4 $^\dagger$ |
| Barabási–Albert | 3.1e-4 | 22.5 / 862 | 36.0 / 16112 | 58.2 | 13.2 / 207 | 17.7 / 829 | 13.3 $^\dagger$ |
| Random geometric | 1.8e-9 | 65.5 / 990 | 67.0 / $7.8\times10^{9}$ | 53.3 $^\ddagger$ | 20.3 / 256 | 37.9 / 4531 | 12.7 $^\dagger$ |

**Table 3.** Persoons et al.'s eq. (20) at $r=10^{-3}$, $V(0)=ru$, all within Lemma 5's hypotheses. Gap is $\bar p-y(L_5)$ at the formula's time $L_5$, in units of $r$, from LSODA and Radau, which agree to all digits shown. Gap $>r$ means the formula under-estimates the growth time.

| Graph | $n$ | $R_0$ | $\min p_i$ | $L_5$ | gap/$r$, LSODA | gap/$r$, Radau | fails |
|------------------------------------------|-----:|-----:|--------:|-------:|---------:|---------:|------:|
| Random geometric (2-D), n=1000, seed 11 | 1000 | 2 | 0.007 | 12.67 | 82.27 | 82.27 | yes |
| Holme–Kim, n=1000, seed 11 | 1000 | 2 | 0.039 | 13.37 | 4.05 | 4.05 | yes |
| Barabási–Albert m=2, n=1000, seed 11 | 1000 | 2 | 0.012 | 13.35 | 5.57 | 5.57 | yes |
| Erdős–Rényi mean degree 10, n=1000, seed 11 | 1000 | 3 | 0.255 | 6.66 | 1.91 | 1.91 | yes |
| Erdős–Rényi mean degree 10, n=1000, seed 11 | 1000 | 5 | 0.415 | 3.39 | 3.57 | 3.57 | yes |
| Lollipop L(6,10) (16 nodes) | 16 | 3 | 0.055 | 6.53 | 67.28 | 67.28 | yes |
| Erdős–Rényi mean degree 10, n=1000, seed 11 | 1000 | 1.2 | 0.027 | 55.67 | 0.37 | 0.37 | no |
| Watts-Strogatz, n=1000, seed 11 | 1000 | 2 | 0.395 | 12.78 | 0.83 | 0.83 | no |

![At the threshold every network obeys $y(t)\le1/(1+t)$ (dashed line, proved in Theorem 3 and attained by regular graphs). The long-run value of $(1+t)\,y(t)$ is the structural constant $\kappa(G)$ of Theorem 7, which ranges from 0.92 to 0.035 across these synthetic networks.](fig_threshold_decay.pdf){width=92%}

Table 2 also shows the limits of our growth bound. It is close for near-regular graphs and exact for regular ones. It is useless (around $10^{10}$) for the random geometric graph at $R_0=1.2$, where some peripheral nodes have endemic levels near zero and $\kappa_p\approx0$.

# 11. Limitations and open questions

- **Model scope.** The sharpness results concern the deterministic NIMFA equations on finite networks. For the exact stochastic SIS process, Proposition 12 transfers the upper bounds only, and these are conservative.
- **Directed temporal networks** require one fixed weight matrix $W$ valid for all snapshots. Without it we have no envelope, and the unweighted bound fails (Proposition 10).
- **Growth from below.** Theorem 14(b) is rigorous but conservative on strongly heterogeneous graphs. A bound governed by the decay rate at $p$, rather than by $\min_ip_i/(1-p_i)$, is the natural target, as is the conjecture that growth from $ru$ takes at most $(1-r)/r$.
- **Node-wise constants uniformly in time.** Asymptotically the worst node has constant of order $n^{1/3}$ (Proposition 9). Uniformly in time and initial state we only know the range $n^{1/3}$ to $n^{1/2}$.
- **Size and localisation.** We have not proved how $\kappa$ scales with $n$ for random graph families. Numerically it decreases for Barabási–Albert and random geometric graphs (Table 1b).
- **Priority.** The prior-art search was bounded and is recorded in the release. We claim no historical priority beyond "not found". We read the abstract, but not the full text, of Persoons and Van Mieghem (2026).

# Verification and reproducibility

- `exact_checks.py` (SymPy; 37 checks, 4 negative controls) re-derives every algebraic identity used in a proof.
- `numerical_checks.py` (NumPy/SciPy/NetworkX; 23 checks, 9 negative controls) tests each theorem on random instances. This includes Theorem 14's formulas end to end, for $r$ on both sides of $c_*/2$ and with mixed initial states, and the directed and heterogeneous limits of Theorem 7. The negative controls are false statements, and each must be *observed to fail*. They include the directed unweighted bound, the arithmetic-mean parameter, node-wise $1/(1+t)$, a constant 1% below $\kappa$, an exponent $10^{-3}$ above $1-\rho$, and identity weights on a directed graph.
- `illustrations.py` produces Tables 1, 1b and 2 and Figure 1, and exits nonzero if any rigorous bound it evaluates is violated.
- `eq20_check.py` produces Table 3 with two independent solvers.

No script uses `assert`. Each records every check and exits with status 1 on any failure, so all behave identically under `python -O`. A fault-injection switch (`NIMFA_INJECT_FAULT=<check>`) confirms that a flipped check makes the run fail. Numerical checks supplement, and do not replace, the written proofs. They test the stated results on selected instances; they do not test them exhaustively.

# AI-use statement

This candidate was produced with AI assistance. A scouting note proved the two conjectures. The extension work, the scripts and this manuscript were then developed and checked within an AI-assisted workflow. That workflow included AI-run multi-perspective reviews and an AI-assisted adversarial review. No human peer review, independent reproduction or formal (machine-checked) verification is claimed. The work is released so that specialists can check it.

# Appendix A. Proof of Lemma 1

Fix $\varepsilon>0$. Let $Y_\varepsilon$ solve $\dot Y_\varepsilon=F^Y(t,Y_\varepsilon)+\varepsilon u$, $Y_\varepsilon(0)=Y(0)+\varepsilon u$, with $F^Y$ extended to $\mathbb R^n$ by the same (locally Lipschitz) formula. For small $\varepsilon$, $Y_\varepsilon$ exists on $[0,T]$, and $Y_\varepsilon\to Y$ uniformly as $\varepsilon\downarrow0$ (continuous dependence).

Let $t_0=\sup\{t\le T:\ X\le Y_\varepsilon\text{ on }[0,t]\}$. Then $t_0>0$ because $X(0)<Y_\varepsilon(0)$, and $X(t_0)\le Y_\varepsilon(t_0)$ by continuity. Suppose $t_0<T$, and let $I$ be the set of coordinates with $X_i(t_0)=Y_{\varepsilon,i}(t_0)$. For $j\notin I$ the strict inequality persists on some $[t_0,t_0+\eta]$. For $i\in I$, $0\le X_i(t_0)\le1$. With $B^X\le B^Y$, $X(t_0)\le Y_\varepsilon(t_0)$ and $\delta^X\ge\delta^Y$,
$$F^X_i(t_0,X(t_0))=-\delta^X_iX_i+(1-X_i)(B^XX)_i\le-\delta^Y_iX_i+(1-X_i)(B^YY_\varepsilon)_i=F^Y_i(t_0,Y_\varepsilon(t_0)),$$
using $X_i=Y_{\varepsilon,i}$, $X_i\ge0$ and $1-X_i\ge0$. Hence the right derivative of $X_i-Y_{\varepsilon,i}$ at $t_0$ is at most $-\varepsilon<0$; at a break point, the right-hand limits of the data are used. So $X_i<Y_{\varepsilon,i}$ on some $(t_0,t_0+\eta]$.

Thus $X\le Y_\varepsilon$ on $[0,t_0+\eta]$, contradicting the definition of $t_0$. So $X\le Y_\varepsilon$ on $[0,T]$, and letting $\varepsilon\downarrow0$ gives $X\le Y$. $\square$

# References

Carr, J. (1981). *Applications of Centre Manifold Theory.* Applied Mathematical Sciences 35, Springer. <https://doi.org/10.1007/978-1-4612-5929-9>

Cator, E., and Van Mieghem, P. (2014). Nodal infection in Markovian susceptible-infected-susceptible and susceptible-infected-removed epidemics on networks are non-negatively correlated. *Physical Review E* 89(5), 052802. <https://doi.org/10.1103/PhysRevE.89.052802>

Dirr, G., Ito, H., Rantzer, A., and Rüffer, B. S. (2015). Separable Lyapunov functions for monotone systems: Constructions and limitations. *Discrete and Continuous Dynamical Systems – B* 20(8), 2497–2526. <https://doi.org/10.3934/dcdsb.2015.20.2497>

Fall, A., Iggidr, A., Sallet, G., and Tewa, J. J. (2007). Epidemiological models and Lyapunov functions. *Mathematical Modelling of Natural Phenomena* 2(1), 62–83. <https://doi.org/10.1051/mmnp:2008011>

Goltsev, A. V., Dorogovtsev, S. N., Oliveira, J. G., and Mendes, J. F. F. (2012). Localization and spreading of diseases in complex networks. *Physical Review Letters* 109(12), 128702. <https://doi.org/10.1103/PhysRevLett.109.128702>

Khanafer, A., Başar, T., and Gharesifard, B. (2016). Stability of epidemic models over directed graphs: A positive systems approach. *Automatica* 74, 126–134. <https://doi.org/10.1016/j.automatica.2016.07.037>

Lajmanovich, A., and Yorke, J. A. (1976). A deterministic model for gonorrhea in a nonhomogeneous population. *Mathematical Biosciences* 28(3–4), 221–236. <https://doi.org/10.1016/0025-5564(76)90125-5>

Mei, W., Mohagheghi, S., Zampieri, S., and Bullo, F. (2017). On the dynamics of deterministic epidemic propagation over networks. *Annual Reviews in Control* 44, 116–128. <https://doi.org/10.1016/j.arcontrol.2017.09.002>

Pastor-Satorras, R., and Castellano, C. (2016). Distinct types of eigenvector localization in networks. *Scientific Reports* 6, 18847. <https://doi.org/10.1038/srep18847>

Persoons, R., Sensi, M., Prasse, B., and Van Mieghem, P. (2024). Transition from time-variant to static networks: Timescale separation in N-intertwined mean-field approximation of susceptible-infectious-susceptible epidemics. *Physical Review E* 109(3), 034308. <https://doi.org/10.1103/PhysRevE.109.034308> (preprint arXiv:2305.12446)

Persoons, R., and Van Mieghem, P. (2026). Susceptible-infectious-susceptible epidemics on temporal contact graphs with clique structure. *Physical Review E* 114(4), 044304. <https://doi.org/10.1103/6jqp-bz2r>

Prasse, B., and Van Mieghem, P. (2020). Time-dependent solution of the NIMFA equations around the epidemic threshold. *Journal of Mathematical Biology* 81(6–7), 1299–1355. <https://doi.org/10.1007/s00285-020-01542-6>

Rantzer, A. (2015). Scalable control of positive systems. *European Journal of Control* 24, 72–80. <https://doi.org/10.1016/j.ejcon.2015.04.004>

Shuai, Z., and van den Driessche, P. (2013). Global stability of infectious disease models using Lyapunov functions. *SIAM Journal on Applied Mathematics* 73(4), 1513–1532. <https://doi.org/10.1137/120876642>

Simon, P. L., and Kiss, I. Z. (2018). On bounding exact models of epidemic spread on networks. *Discrete and Continuous Dynamical Systems – B* 23(5), 2005–2020. <https://doi.org/10.3934/dcdsb.2018192>

Smith, H. (2008). *Monotone Dynamical Systems: An Introduction to the Theory of Competitive and Cooperative Systems.* Mathematical Surveys and Monographs 41, American Mathematical Society (first published 1995). <https://doi.org/10.1090/surv/041>

Van Mieghem, P. (2012). Epidemic phase transition of the SIS type in networks. *EPL (Europhysics Letters)* 97(4), 48004. <https://doi.org/10.1209/0295-5075/97/48004>

Van Mieghem, P., Omic, J., and Kooij, R. (2009). Virus spread in networks. *IEEE/ACM Transactions on Networking* 17(1), 1–14. <https://doi.org/10.1109/TNET.2008.925623>
