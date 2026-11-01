# How fast does an epidemic die out on a network? Sharp decay envelopes for mean-field SIS epidemics and a proof of the Persoons–Sensi–Prasse–Van Mieghem conjectures

Evidence Press research release, 1 November 2026. Version 0.1.0-candidate. DOI: [10.5281/zenodo.23856194](https://doi.org/10.5281/zenodo.23856194) (concept DOI, all versions: [10.5281/zenodo.23856193](https://doi.org/10.5281/zenodo.23856193)). Release page: <https://evidencepress.org/releases/nimfa-sis-decay-envelopes/>. Author: Anonymous. UK English. Text, data and diagrams CC0-1.0; code MIT. Unrefereed candidate; internal AI-assisted review only.

## In plain words

When an infection spreads through a network of contacts but is not quite strong enough to persist, how fast does it die out, and does the shape of the network matter? In the standard deterministic network model of a recurring infection (SIS, in its N-intertwined mean-field form, NIMFA), Persoons, Sensi, Prasse and Van Mieghem (2024) conjectured from simulations a universal speed limit. At or below the epidemic threshold, the average infection level satisfies $y(t)\le1/(1+t)$ (time in units of the mean recovery time), whatever the network and starting state. A matching limit holds above the threshold for settling to the endemic level from above.

This paper proves both conjectures. The same short argument extends them to time-varying networks, node-dependent recovery and (with Perron weights) one-way contacts. Via a known correlation inequality, it also covers the expected prevalence of the exact stochastic model. The limit is attained by regular networks. Exactly at the threshold, network shape matters through one number, $\kappa(G)\le1$, computed from the leading eigenvector: the average decays like $\kappa(G)/t$. On standard synthetic networks $\kappa$ ranges from about 0.9 (random, small-world) to about 0.04 (spatial, scale-free), where infection concentrates in a persistent core. Individual hubs can exceed the population limit by a factor of order $n^{1/3}$.

## Main results

- **Theorem 1 (envelope).** Under a quadratic-form condition, a weighted average of $h(x_i)=-x_i-\log(1-x_i)$ obeys a scalar Riccati comparison. Jensen's inequality transfers it to the weighted average of $x_i$.
- **Theorems 3–4.** Conjectures 1 and 2 (the latter with prefactor $1-\bar p$ and an exponential form), plus temporal, heterogeneous-recovery and weighted-directed versions.
- **Propositions 5–6, Theorem 7, Corollary 8.** Attainment on regular graphs, exact subcritical rate, and the critical constant $t\,x_i(t)\to v_i/\sum_ju_j\delta_jv_j^2$ (directed and heterogeneous). The constant satisfies $\kappa(G)\le1$ with equality iff the graph is regular, and $r\,T(r)\to\kappa(G)$.
- **Propositions 9–12.** A node-wise bound $\sqrt{n/(t(1+t))}$ and asymptotic worst-node order $n^{1/3}$; counterexamples to stronger statements; SIR/SIRS; the exact stochastic SIS process (common recovery rate).
- **Proposition 13, Theorem 14.** Timescale separation: mixed initial states never converge slower than both extremes, and unconditional transition-time bounds hold for all $R_0$. The growth part is exact on regular graphs and conservative on strongly heterogeneous ones. Two-solver numerics show that the growth formula of Persoons et al. (their eq. 20) is not a general bound away from $R_0\to1$.

Not done: a sharp growth bound on heterogeneous graphs; the uniform-in-time node-wise order; scaling laws for $\kappa$; interval-certified versions of the eq. (20) failures.

## Files

| File | Contents |
|---|---|
| `nimfa_envelopes.pdf` | The paper, with complete proofs |
| `manuscript.md` | Source of the paper; `build_pdf.sh` and `paper_header.tex` rebuild the PDF (`nimfa_envelopes.tex` is the generated LaTeX) |
| `novelty_report.md` | Prior-art search record and its limits |
| `REVIEW_SUMMARY.md` | Internal AI-assisted review: objections and how each was handled |
| `VERIFICATION.md` | What each script checks, how to run it, and what it does not establish |
| `exact_checks.py` | SymPy re-derivation of every algebraic identity used in a proof (37 checks, 4 negative controls) |
| `numerical_checks.py` | Numerical tests of every theorem on random instances (23 checks, 9 negative controls) |
| `illustrations.py` | Tables 1, 1b, 2 and Figure 1 on standard synthetic networks |
| `eq20_check.py` | Table 3: two-solver check of the growth formula of Persoons et al. |
| `*.json`, `fig_threshold_decay.*` | Recorded outputs and the figure |
| `tools/replay_compare.py` | Clean replay in a temporary copy, with output comparison and fault injection (`--quick`: the two fast checkers) |
| `refs.bib` | Bibliography, copied from Crossref records |
| `STATUS.md`, `ASSURANCE.md`, `PROVENANCE.md`, `CLAIMS.json`, `RESEARCH_GATES.json`, `AI_INDEX.md` | Status, assurance limits, provenance, machine-readable claims, prior-art gates, agent navigation |
| `REPLAY_RECEIPT.md`, `ENVIRONMENT.txt`, `SHA256SUMS`, `.github/workflows/replay.yml` | Replay record, environment, integrity sums, CI |
| `CITATION.cff`, `.zenodo.json`, `LICENSE`, `LICENSE-MIT.txt`, `LICENSES.md`, `PUBLIC_DOMAIN.md` | Citation and licences |

## Reproduce

Python 3.11 with the pinned packages in `requirements.txt`:

```bash
uv run --python 3.11 --with-requirements requirements.txt bash build.sh
```

`build.sh` runs every checker under `python -O` (about 8 minutes in total, mostly `illustrations.py`) and then `build_pdf.sh`, which needs Pandoc and Tectonic. No script uses `assert`. Each exits with status 1 on any failed check, and `NIMFA_INJECT_FAULT=<check name>` flips one check to confirm the harness fails.

Clean replay with output comparison (no PDF build): `python -I -B tools/replay_compare.py .` (add `--quick` for the two fast checkers only).
