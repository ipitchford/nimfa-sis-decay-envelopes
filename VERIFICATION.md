# Verification

Numerical and symbolic checks supplement, and do not replace, the written proofs in `nimfa_envelopes.pdf`. They
test the stated results on selected instances; they do not test them exhaustively. No formal (machine-checked)
verification is claimed.

## Scripts

| Script | What it checks | Runtime | Output |
|---|---|---|---|
| `exact_checks.py` | Exact SymPy re-derivation of every algebraic identity used in a proof: derivative of $h$, convexity identity of Lemma 3, envelope ODE, closed forms (3.3), subcritical and critical limits, residual transform (4.3), both subsolution identities of Theorem 14, the quadratic-form identity of Lemma 2, star and hub-plus-regular eigenvector formulas, the directed counterexample ($\ddot y(0)=97/32$) and its Perron weights, the arithmetic-mean counterexample, the bounds on $h$, Chebyshev's identity, the regular-graph endemic exactness, logistic hitting times, and the Hölder step of Proposition 9(d) on 200 rational vectors. Negative controls: false variants evaluated at an exact rational point must be nonzero. | seconds | `exact_checks.json` |
| `numerical_checks.py` | Theorems 1, 3, 4 on random weighted graphs below, at and above threshold, including random initial states with the entropy parameter. Directed Perron weights (quadratic form and envelope); heterogeneous recovery, symmetric and directed; noncommuting temporal snapshots with time-varying recovery; disconnected graphs. Theorem 7 limits (undirected, directed, heterogeneous) and $r\,T(r)\to\kappa$. Proposition 9(a); Theorem 14 end to end for $r$ on both sides of $c_*/2$ and mixed initial states, with exact regular-graph growth times at $r<p/2$, $p/2\le r<p$, $r=p$, $r>p$; SIR and SIRS; regular-graph equality; exact subcritical rate. Negative controls (each must be observed to fail): unweighted directed bound, two-node counterexample, identity weights on a directed graph, constant 1 at the threshold for a Barabási–Albert graph, node-wise $1/(1+t)$, constant $0.99/(1+t)$ on a regular graph, exponent $10^{-3}$ above $1-\rho$, a constant 1% below $\kappa$, the arithmetic-mean parameter. | seconds | `numerical_checks.json` |
| `illustrations.py` | Tables 1, 1b, 2 and Figure 1. Exits nonzero if any rigorous bound it evaluates (envelope, node bound, Theorem 14 decay and growth bounds) is violated. | about 5 minutes | `illustrations.json`, `fig_threshold_decay.*` |
| `eq20_check.py` | Table 3: integrates growth from $ru$ with LSODA and with Radau (analytic Jacobian), at relative tolerance $10^{-11}$. Records whether Persoons et al.'s eq. (20) under-estimates the growth time within its hypotheses, and fails if the two solvers disagree. | about 2 minutes | `eq20_check.json` |

## Assurance properties

- No `assert` statements. Every check is recorded, and the scripts exit with status 1 on any failure, so they behave
  identically under `python -O`. All shipped outputs were produced with `python3 -I -O -B`.
- Fault injection: `NIMFA_INJECT_FAULT=<check name>` flips that check; the run then fails. This was confirmed for
  `exact_checks.py` (E3) and `numerical_checks.py` (N5a).
- A clean replay of the checkers from a fresh copy reproduced the recorded JSON exactly. An independent reviewer re-ran
  both fast checkers on the revised version, also with identical outputs.

## Limits

- Integration tolerances (LSODA, relative $10^{-10}$–$10^{-11}$) bound the meaning of one-sided numerical checks. A
  slack of $10^{-8}$ is allowed.
- The eq. (20) failures in Table 3 are high-precision numerical findings with two independent solvers. They are not
  interval-certified.
- Synthetic networks only. No real contact network is analysed.
