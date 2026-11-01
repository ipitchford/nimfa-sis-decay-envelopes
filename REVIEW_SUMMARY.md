# Internal review summary (AI-assisted; not external peer review)

Dates: 1 November 2026. Routes:

- Five reviewer roles (editor-in-chief, methodology/proof, domain, cross-disciplinary perspective, devil's
  advocate), run as independent Claude Code subagents following the academic-paper-reviewer protocol.
- An adversarial review by a separate model (OpenAI Codex "Astra", via a local companion tool) and, after revision,
  its audit of the new proofs.
- A separate proof audit of the new and rewritten proofs.
- A clean replay of the checkers from a fresh copy.
- A re-review of the revision against every round-one objection.

This is internal checking, not independent reproduction or human peer review.

**Outcome.** No reviewer or audit found a critical issue. The proofs of the two conjectures (Theorems 1, 3, 4) were
re-derived independently by several reviewers and survived targeted attempts to break them. These included regimes
with infection rate above recovery, temporal switching, directed weights and heterogeneous recovery. Major issues
concerned one domain gap in a transition-time theorem, the scope of the numerical refutation, novelty framing and
readability. All were handled in one revision round plus a short verification fix.

| Objection | Response |
|---|---|
| Transition-time theorem: phase-1 time is negative when $r>c_*/2$; a logarithm could have a negative argument; the regular-graph exactness formula failed for $r\ge p/2$ | Rewritten with clipped phase times, piecewise second phase and $m_1\le1/2$ proved; exact regular-graph formula stated piecewise; tested end to end across $r$, including $r\ge p$ |
| The claimed failure of Persoons et al.'s eq. (20) included one case outside that lemma's hypotheses ($r>\min_ip_i$); failure diagnosed as heterogeneity only | That case relabelled; a $\min p_i$ column added; new two-solver Table 3 with cases inside the hypotheses, including Erdős–Rényi at $R_0=3,5$ and a 16-node lollipop; cause stated as leaving the $R_0\to1$ regime; labelled numerical, not certified |
| "Assumption 1 always holds" misquotes the source, whose literal wording says "slower" | Quoted verbatim; the paper now says the intended property (error bounded by the worse extreme) always holds, and the literal wording never does |
| The rigorous growth bound is orders of magnitude loose on heterogeneous graphs, while the source's combined bound survived all tests | Stated in the abstract, introduction and Section 9 |
| The $1/t$ law at the threshold is the standard centre-manifold mechanism | Stated with a citation; novelty restricted to the global statement, explicit constants and the characterisation $\kappa\le1$ iff regular |
| $h$ is the classical Volterra/Goh function; separable Lyapunov functions and diagonal stability are classical | Credited (Shuai–van den Driessche; Dirr et al.; Rantzer; Khanafer et al.) |
| "Network shape matters" is really eigenvector localisation, depends on size and seed, and small $\kappa$ is not reassuring | Localisation explained and cited; Table 1b with medians over 10 seeds, $n\cdot\mathrm{IPR}$ and size dependence; plain-English text says small $\kappa$ signals a persistent core |
| The claim that the source's simulations had $\kappa$ near 1 was not checked | Computed: their $N=50$ Erdős–Rényi ensembles have median $\kappa=0.96$ |
| The stochastic model was disclaimed although NIMFA upper-bounds it | New Proposition 12 via Simon and Kiss (2018), restricted to a common recovery rate after the audit noted that the cited result assumes one |
| Theorem 7 could be extended to directed and heterogeneous recovery | Done, with a numerical check |
| The node-wise gap between $n^{1/3}$ and $\sqrt n$ | Hölder argument shows the asymptotic worst-node order is exactly $n^{1/3}$; only the uniform-in-time order remains open |
| Disconnected exponential form ($\kappa_p=0$), $x(0)\neq0$ in node limits, local boundedness in Theorem 1, notation clashes, redundant minima | Fixed |
| Weak negative controls; verification text overstated coverage; stale numbering in scripts | Near-boundary controls added; symbolic controls evaluated at exact rational points; text narrowed; numbering fixed |
| Title and abstract too long; scale and knife-edge missing for lay readers | Shorter title; three-part abstract; plain-English scale ("52 time units is a year with weekly recovery") and knife-edge sentences |

Recorded as open: a sharp growth bound on heterogeneous graphs (perhaps at most $(1-r)/r$); the uniform-in-time
node-wise order; scaling of $\kappa$ with network size; interval-certified eq. (20) counterexamples; real contact
networks.
