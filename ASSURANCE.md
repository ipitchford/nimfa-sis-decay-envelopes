# Assurance boundary

## What was done

1. **Written proofs** in the manuscript ([nimfa_envelopes.pdf](nimfa_envelopes.pdf); source [manuscript.md](manuscript.md)).
2. **Computational checks** ([VERIFICATION.md](VERIFICATION.md)): `exact_checks.py` re-derives every algebraic identity used in a proof with SymPy (37 checks and 4 negative controls evaluated at exact rational points); `numerical_checks.py` tests every theorem on random instances below, at and above threshold, including temporal, heterogeneous-recovery, directed and disconnected cases (23 checks and 9 negative controls that must be observed to fail); `illustrations.py` (Tables 1, 1b, 2 and Figure 1; exits nonzero if any rigorous bound it evaluates is violated); `eq20_check.py` (Table 3; two solvers, fails if they disagree). No script uses `assert`; each exits with status 1 on any failed check, so the checks stay active under `python -O`. `NIMFA_INJECT_FAULT=<check name>` flips one check to confirm the harness fails.
3. **Internal AI-assisted review.** Five reviewer roles (editor, methodology/proof, domain, perspective, devil's advocate) as Claude Code subagents; an adversarial review and, after revision, an audit of the new proofs by a separate model (OpenAI Codex, through a local companion tool); a separate proof audit; a clean replay of the checkers by the editor (outputs identical); one revision round; a re-review. A rating-free summary is in [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).
4. **Replay.** From a fresh extraction of the frozen archive every checker was rerun with `-I -O -B` in a temporary copy by [tools/replay_compare.py](tools/replay_compare.py), the regenerated JSON outputs were compared with the retained ones (numbers within a relative tolerance of 1e-7; package-version fields ignored), and fault injection was confirmed to make the exact and numerical checkers fail; see [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). CI reruns the exact and numerical checkers on Linux as a blocking job and the two slower scripts as a non-blocking job.
5. **Citations.** Every reference was copied from its Crossref record ([refs.bib](refs.bib)); DOIs were checked against the registries at release.

## What was not done

- No human specialist or journal peer review, no independent reproduction by an unaffiliated party, no proof-assistant formalisation.
- Symbolic and numerical checks corroborate identities and test the theorems on selected instances; the theorems rest on the written proofs. The eq. (20) failures (Table 3) are high-precision numerical findings with two solvers, not interval-certified.
- **Model scope.** The results concern the deterministic NIMFA model (and, for Proposition 12, an upper bound on the expected prevalence of the stochastic SIS model with a common recovery rate). The mean-field model is least accurate near the threshold. The examples use synthetic networks only.
- **Priority is bounded** by the documented searches in [novelty_report.md](novelty_report.md). The conjectures are due to Persoons, Sensi, Prasse and Van Mieghem (2024). The function $h$ is the classical Volterra-type Lyapunov function; separable Lyapunov functions, diagonal weights (Khanafer et al. 2016; Rantzer 2015), comparison for cooperative systems and the centre-manifold mechanism behind the $1/t$ decay are classical and credited.
- Agreement between AI reviewers is not independent evidence.
