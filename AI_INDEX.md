# Agent navigation

## Identity

- Release slug: `nimfa-sis-decay-envelopes`. Version `0.1.0-candidate`, tag `v0.1.0-candidate`. DOI 10.5281/zenodo.23856194; concept DOI 10.5281/zenodo.23856193.
- Repository: https://github.com/ipitchford/nimfa-sis-decay-envelopes. Release page: https://evidencepress.org/releases/nimfa-sis-decay-envelopes/
- Scholarly creator: Anonymous. Status: unrefereed candidate ([STATUS.md](STATUS.md)).
- Manuscript: [PDF](nimfa_envelopes.pdf) (19 pages), Markdown source [manuscript.md](manuscript.md), generated LaTeX [nimfa_envelopes.tex](nimfa_envelopes.tex), builder [build_pdf.sh](build_pdf.sh) with [paper_header.tex](paper_header.tex); bibliography [refs.bib](refs.bib) (Crossref records).
- Machine-readable claims: [CLAIMS.json](CLAIMS.json). Overview: [README.md](README.md).
- Assurance limits: [ASSURANCE.md](ASSURANCE.md). Review record (rating-free): [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). Provenance: [PROVENANCE.md](PROVENANCE.md).
- Prior art: [RESEARCH_GATES.json](RESEARCH_GATES.json), [novelty_report.md](novelty_report.md).
- Citation: [CITATION.cff](CITATION.cff), [.zenodo.json](.zenodo.json). Licences: [LICENSES.md](LICENSES.md), [LICENSE](LICENSE), [LICENSE-MIT.txt](LICENSE-MIT.txt), [PUBLIC_DOMAIN.md](PUBLIC_DOMAIN.md).

## Scientific map

Object: the N-intertwined mean-field (NIMFA) SIS equations dx_i/dt = -delta_i x_i + (1 - x_i)(Bx)_i on [0,1]^n, with average infection level y(t) and R0 = rho(B) for unit recovery. Question: the two decay conjectures of Persoons, Sensi, Prasse and Van Mieghem (Phys. Rev. E 109, 034308, 2024) — y(t) <= 1/(1+t) at or below threshold, and the same bound for the excess over the endemic level above threshold from above — and how sharp they are.

- Plain-English summary and Section 1: the model, the conjectures (quoted), known versus new (table in Section 1.2).
- Section 2: Lemma 1 (comparison), Lemma 2 (diagonal weights), Lemma 3 (the Volterra-type function h(s) = -s - log(1 - s)).
- Section 3: Theorem 1 (spectral Riccati envelope), Corollary 2 (when condition (D) holds: symmetric, temporal, heterogeneous recovery, directed with Perron weights).
- Section 4: Theorem 3 (Conjecture 1), Theorem 4 (Conjecture 2 with prefactor 1 - p_bar).
- Section 5: Proposition 5 (attainment on regular graphs), Proposition 6 (exact subcritical rate), Theorem 7 (critical constant kappa(G) <= 1, equality iff regular), Corollary 8 (r T(r) -> kappa).
- Section 6: Proposition 9 (node-wise bounds, order n^(1/3)). Section 7: Proposition 10 (counterexamples). Section 8: Propositions 11-12 (SIR/SIRS, exact stochastic SIS upper bound).
- Section 9: Proposition 13, Theorem 14 (timescale separation, unconditional transition-time bounds). Section 10: Tables 1, 1b, 2, 3 and Figure 1 (synthetic networks; eq. (20) of Persoons et al. checked with two solvers). Section 11: limitations and open questions. Appendix A: proof of Lemma 1.

## Evidence and replay

- Script guide and what each check does not establish: [VERIFICATION.md](VERIFICATION.md). One-command build: [build.sh](build.sh). Comparer: [tools/replay_compare.py](tools/replay_compare.py) (temporary copy; every checker with -I -O -B; regenerated outputs compared within a relative tolerance of 1e-7; fault injection; `--quick` runs only the two fast checkers). Environment: [requirements.txt](requirements.txt), [ENVIRONMENT.txt](ENVIRONMENT.txt). CI: [.github/workflows/replay.yml](.github/workflows/replay.yml). Record: [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). Integrity: [SHA256SUMS](SHA256SUMS).
- Scripts and outputs: [exact_checks.py](exact_checks.py) -> [exact_checks.json](exact_checks.json) (SymPy identities, 37 checks + 4 negative controls); [numerical_checks.py](numerical_checks.py) -> [numerical_checks.json](numerical_checks.json) (23 checks + 9 negative controls); [illustrations.py](illustrations.py) -> [illustrations.json](illustrations.json), [fig_threshold_decay.pdf](fig_threshold_decay.pdf), [fig_threshold_decay.png](fig_threshold_decay.png) (Tables 1, 1b, 2, Figure 1); [eq20_check.py](eq20_check.py) -> [eq20_check.json](eq20_check.json) (Table 3).

## Trust boundary

All review was internal and AI-assisted. There was no peer review, no independent reproduction and no formal verification. Symbolic checks verify identities; numerical checks test the theorems on selected instances; the theorems rest on the written proofs. The eq. (20) failures (Table 3) are numerical with two solvers, not interval-certified. The results concern the deterministic NIMFA model (stochastic transfer is an upper bound only, common recovery rate) and synthetic networks. The conjectures are due to Persoons et al.; the Volterra-type function, separable Lyapunov functions, diagonal weights and the centre-manifold mechanism are classical and credited. Priority holds only within the documented bounded searches.
