# Prior-art and novelty record

Dates: 31 October – 1 November 2026. This is a bounded search, recorded so that readers can check it. It is not
a proof of priority. Permitted language: "we did not find"; not "first".

## Target

Conjectures 1 and 2 of Persoons, Sensi, Prasse and Van Mieghem, *Physical Review E* 109, 034308 (2024),
arXiv:2305.12446v3, their eqs. (9)–(10). The full text was read, including Section III (Lemma 5, eq. (20),
Assumption 1, eq. (23)) and Appendix E, which leaves inequalities (E8) and (E9) unproved.

## Evidence Press corpus

The live feed (210 releases, 1 November 2026) was scanned for NIMFA, SIS, epidemic, contagion, spreading and
mean-field. Related but distinct releases: algebraic extinction in acyclic compartment models (two releases) and
reconstruction and continuation in network SIS (Keliger–Horváth question). None treats NIMFA decay envelopes or
these conjectures.

## Citing works and follow-ups

- Semantic Scholar listed eight works citing the 2024 paper on 1 November 2026: cyber-risk modelling; community
  quarantine; fractional Markov epidemics; simulation of coevolving networks; finding patient zero (two records);
  individual decision-making; discrete-time delayed SIS. From titles and abstracts, none proves or addresses the
  conjectures.
- Persoons and Van Mieghem, *Phys. Rev. E* 114, 044304 (October 2026), on temporal clique graphs: the abstract
  (Crossref) concerns exact and approximate solutions on temporal clique graphs and does not mention the
  conjectures. Full text not read.

## Related work read in full

- Prasse and Van Mieghem (2020), *J. Math. Biol.* 81, 1299: near-threshold approximation as $R_0\downarrow1$ for growth
  from below. The cubic Perron sum also appears in their steady-state formula. Their results are limit
  approximations, not all-time bounds, and do not treat decay at $R_0=1$ exactly.
- Khanafer, Başar and Gharesifard (2016), *Automatica* 74: global stability on directed graphs. Lemma A.1 gives the
  diagonal Lyapunov weights used in our Lemma 2. Their result is qualitative, with no decay rate.
- Mei, Mohagheghi, Zampieri and Bullo (2017), *Annu. Rev. Control* 44: Theorem 4, exponential decay of a
  left-Perron-weighted average below threshold, which degenerates at $R_0=1$.
- Van Mieghem (2016), *PRE* 93, 052312, and arXiv:1612.01386: tanh/Riccati bounds for the Markovian SIS prevalence,
  conditional on bounding a remainder by constants. Not an unconditional NIMFA envelope.
- Simon and Kiss (2018), *DCDS-B* 23, 2005 (arXiv:1704.01726): NIMFA upper-bounds the exact stochastic SIS node
  probabilities (common recovery rate), used in our Proposition 12.
- Chang and Van Mieghem (2026), *Commun. Phys.*: generation-time distributions; unrelated to decay envelopes.

## Searches

Web and Semantic Scholar queries included: "Persoons Sensi Prasse Van Mieghem conjecture proof NIMFA 1/(1+t)";
"networked SIS critical spectral radius one convergence rate 1/t"; "NIMFA convergence"; "SIS mean-field decay rate";
formula-level queries on $-x-\log(1-x)$ Lyapunov functions in SIS; plus arXiv title searches for the works above.
An earlier scouting search (31 October) is summarised in the same terms. One Semantic Scholar snippet search hit a
rate limit.

## Adjudication

| Contribution | Finding |
|---|---|
| Proofs of Conjectures 1 and 2 | No prior proof found |
| Quantitative Riccati envelope via Jensen on the Volterra-type function | Not found. The function itself is classical (Volterra/Goh/Korobeinikov; Shuai and van den Driessche 2013), and so are separable Lyapunov functions and diagonal weights |
| Critical constant $\kappa(G)$, global statement, characterisation $\kappa\le1$ iff regular | Not found as stated. The $1/t$ mechanism is the standard centre-manifold reduction |
| Node-wise order $n^{1/3}$; unconditional transition-time bounds; numerical failures of eq. (20) | Not found |
