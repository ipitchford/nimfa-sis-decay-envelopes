"""Two-solver numerical check of where Persoons et al.'s growth formula (their eq. 20, first form)
fails as an upper bound on the growth time (Section 9 of the paper).

Their Lemma 5 (derived as R0 -> 1 from above) states, for V(0) with r <= v_i(0) <= p_i and |y(0)-y_inf| > r,
    T(r) <= L5 := 2/(R0-1) * log((p_max - r)/r).
For each case we check the lemma's hypotheses (r <= min_i p_i), integrate the growth from V(0) = r u with
two independent solvers (LSODA, and Radau with the analytic Jacobian), and report the gap
pbar - y(L5) at the formula's time. A case is recorded as a counterexample to the formula only if both
solvers give a gap exceeding r by more than a factor 1.2 (the solvers agree to far better than this).
We also record the time at which the combined bound (23) of Persoons et al., max((1-r)/r, L5), would apply.

This is a high-precision numerical finding, not an interval-certified proof. Exit status 1 only if the two
solvers disagree materially (relative gap difference > 1e-6) or a listed expected counterexample is not found.
No `assert`; safe under python -O.

Usage: python3 -I -O eq20_check.py [--out eq20_check.json]
"""
import json
import sys

import networkx as nx
import numpy as np
from scipy.integrate import solve_ivp

R_ACC = 1e-3


def endemic(B):
    p = np.ones(len(B))
    for _ in range(10**6):
        bp = B @ p
        nxt = bp / (1 + bp)
        if np.max(np.abs(nxt - p)) < 1e-15:
            return nxt
        p = nxt
    raise RuntimeError("no convergence")


def lollipop(m, k):
    return nx.lollipop_graph(m, k)


def cases():
    return [
        ("Random geometric (2-D), n=1000, seed 11", nx.random_geometric_graph(1000, 0.06, seed=11), 2.0, True),
        ("Holme-Kim, n=1000, seed 11", nx.powerlaw_cluster_graph(1000, 3, 0.3, seed=11), 2.0, True),
        ("Barabasi-Albert m=2, n=1000, seed 11", nx.barabasi_albert_graph(1000, 2, seed=11), 2.0, True),
        ("Erdos-Renyi mean degree 10, n=1000, seed 11", nx.erdos_renyi_graph(1000, 0.01, seed=11), 3.0, True),
        ("Erdos-Renyi mean degree 10, n=1000, seed 11", nx.erdos_renyi_graph(1000, 0.01, seed=11), 5.0, True),
        ("Lollipop L(6,10) (16 nodes)", lollipop(6, 10), 3.0, True),
        ("Erdos-Renyi mean degree 10, n=1000, seed 11", nx.erdos_renyi_graph(1000, 0.01, seed=11), 1.2, False),
        ("Watts-Strogatz, n=1000, seed 11", nx.watts_strogatz_graph(1000, 10, 0.1, seed=11), 2.0, False),
    ]


def gap_at(B, p, t_end, method):
    n = len(p)

    def f(_t, x):
        return -x + (1 - x) * (B @ x)

    def jac(_t, x):
        return np.diag(-1 - B @ x) + (1 - x)[:, None] * B

    kw = {"jac": jac} if method == "Radau" else {}
    sol = solve_ivp(f, (0, t_end), R_ACC * np.ones(n), method=method, rtol=1e-11, atol=1e-14,
                    dense_output=False, t_eval=[t_end], **kw)
    if not sol.success:
        raise RuntimeError(sol.message)
    return float(p.mean() - sol.y[:, -1].mean())


def main():
    out = "eq20_check.json"
    if "--out" in sys.argv:
        out = sys.argv[sys.argv.index("--out") + 1]
    rows, problems = [], []
    for name, G, R0, expect_fail in cases():
        G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
        A = nx.to_numpy_array(G)
        lam = np.linalg.eigvalsh(A)[-1]
        B = R0 * A / lam
        p = endemic(B)
        L5 = 2 / (R0 - 1) * np.log((p.max() - R_ACC) / R_ACC)
        hyp = bool(R_ACC <= p.min() and abs(R_ACC - p.mean()) > R_ACC)
        g1 = gap_at(B, p, L5, "LSODA")
        g2 = gap_at(B, p, L5, "Radau")
        agree = abs(g1 - g2) <= 1e-6 * max(1.0, abs(g1))
        fails = bool(hyp and g1 > 1.2 * R_ACC and g2 > 1.2 * R_ACC)
        if not agree:
            problems.append(f"{name} R0={R0}: solvers disagree ({g1} vs {g2})")
        if expect_fail and not fails:
            problems.append(f"{name} R0={R0}: expected counterexample not found")
        rows.append({"graph": name, "n": len(p), "R0": R0, "p_min": float(p.min()), "p_max": float(p.max()),
                     "L5": float(L5), "lemma5_hypotheses_hold": hyp, "gap_LSODA": g1, "gap_Radau": g2,
                     "gap_over_r": g1 / R_ACC, "formula_fails": fails,
                     "persoons_eq23_time": float(max((1 - R_ACC) / R_ACC, L5))})
        print(f"{name[:44]:44s} R0={R0:3.1f} pmin={p.min():.2e} hyp={hyp!s:5s} L5={L5:6.2f} "
              f"gap/r LSODA={g1 / R_ACC:8.2f} Radau={g2 / R_ACC:8.2f} fails={fails}")
    summary = {"script": "eq20_check.py", "r": R_ACC, "numpy": np.__version__, "rows": rows,
               "problems": problems, "status": "PASS" if not problems else "FAIL"}
    with open(out, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")
    print(json.dumps({"status": summary["status"], "problems": problems}))
    if problems:
        sys.exit(1)


if __name__ == "__main__":
    main()
