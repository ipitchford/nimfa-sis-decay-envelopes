"""Numerical illustration on standard synthetic network models (Section 10 of the paper).

Synthetic graphs only (fixed seeds). For each graph at the threshold R0 = 1 the script reports
the structural constant kappa(G), the observed t*y(t), the largest node constant, and how
tight 1/(1+t) is. Above threshold it compares observed transition times with the rigorous
bounds of Theorem 14 and with the (non-rigorous, R0 -> 1) bound of Persoons et al.
It writes illustrations.json and fig_threshold_decay.pdf/.png.

Illustrations are not proofs. The script exits with status 1 if any rigorous bound it
evaluates is violated (no `assert`; safe under python -O).

Usage: python3 -I -O illustrations.py [--outdir .]
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

FAIL = []


def integrate(B, x0, ts):
    f = lambda _t, x: -x + (1 - x) * (B @ x)  # noqa: E731
    sol = solve_ivp(f, (ts[0], ts[-1]), x0, t_eval=ts, method="LSODA", rtol=1e-10, atol=1e-15)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y


def perron(A):
    w, V = np.linalg.eigh(A)
    return w[-1], np.abs(V[:, -1])


def endemic(B):
    p = np.ones(len(B))
    for _ in range(10**6):
        bp = B @ p
        nxt = bp / (1 + bp)
        if np.max(np.abs(nxt - p)) < 1e-15:
            return nxt
        p = nxt
    raise RuntimeError("no convergence")


def graphs():
    n = 1000
    G = {
        "Erdos-Renyi (mean degree 10)": nx.erdos_renyi_graph(n, 0.01, seed=11),
        "Watts-Strogatz small world": nx.watts_strogatz_graph(n, 10, 0.1, seed=11),
        "Stochastic block model (4 communities)": nx.stochastic_block_model(
            [250] * 4, [[0.036 if i == j else 0.002 for j in range(4)] for i in range(4)], seed=11),
        "Random geometric (2-D)": nx.random_geometric_graph(n, 0.06, seed=11),
        "Holme-Kim (scale-free, clustered)": nx.powerlaw_cluster_graph(n, 3, 0.3, seed=11),
        "Barabasi-Albert (m=2)": nx.barabasi_albert_graph(n, 2, seed=11),
        "Star (1 hub, 199 leaves)": nx.star_graph(199),
    }
    hub = nx.random_regular_graph(158, 2000, seed=11)
    hub.add_edges_from((2000, j) for j in range(2000))
    G["Hub joined to 158-regular graph (2001 nodes)"] = hub
    out = {}
    for k, g in G.items():
        cc = max(nx.connected_components(g), key=len)
        out[k] = nx.convert_node_labels_to_integers(g.subgraph(cc).copy())
    return out


def transition_bounds(B, p, rho, r):
    """Rigorous bounds of Theorem 14(b): T_D (decay side) and T_G (growth side, two phases)."""
    s = 1 - p
    kap = float(np.min(p / s))
    sbar, pbar = float(s.mean()), float(p.mean())
    # decay side: exponential form of (4.2) is always the smaller one
    TD = 0.0 if r >= sbar else max(0.0, np.log((1 + sbar * kap / r) / (1 + kap)) / kap)
    lam, v = perron(B)
    v = v / v.max()
    g = rho - 1
    cstar = (rho - 1) / rho
    t1 = max(0.0, np.log((cstar - r) / r) / g) if r < cstar else 0.0
    m1 = (cstar / 2) * float(np.min(v / p))  # always <= 1/2
    if pbar * (1 - m1) <= r:
        t2 = 0.0
    else:
        t2 = np.log((1 - m1) * (pbar - r) / (m1 * r)) / kap
    TG = t1 + t2
    persoons = max((1 - r) / r, 2 / (rho - 1) * np.log((p.max() - r) / r))
    return {"kappa_endemic": kap, "p_min": float(p.min()), "T_decay_bound": TD, "T_growth_bound": TG,
            "T_bound": max(TD, TG), "persoons_eq23": persoons}


def observed_T(B, p, r, horizon):
    """Last time the decay (from u) or growth (from r u) prevalence error exceeds r."""
    ts = np.linspace(0, horizon, int(horizon / 0.005) + 1)
    Yu = integrate(B, np.ones(len(p)), ts).mean(0) - p.mean()
    Yr = p.mean() - integrate(B, r * np.ones(len(p)), ts).mean(0)

    def last(err):
        bad = np.nonzero(np.abs(err) > r)[0]
        if not len(bad):
            return 0.0
        return float(ts[bad[-1] + 1]) if bad[-1] + 1 < len(ts) else np.inf

    return last(Yu), last(Yr)


def main():
    outdir = "."
    if "--outdir" in sys.argv:
        outdir = sys.argv[sys.argv.index("--outdir") + 1]
    rows = []
    tl = np.r_[0.0, np.geomspace(0.01, 1e4, 300)]
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    for name, G in graphs().items():
        A = nx.to_numpy_array(G)
        lam, v = perron(A)
        B = A / lam
        n = len(A)
        X = integrate(B, np.ones(n), tl)
        y = X.mean(0)
        kap = v.sum() * (v**2).sum() / (n * (v**3).sum())
        node = v * (v**2).sum() / (v**3).sum()
        ratio = y * (1 + tl)
        sq = np.sqrt(n / (tl[1:] * (1 + tl[1:])))
        if np.max(ratio) > 1 + 1e-8:
            FAIL.append(f"{name}: y(t) exceeded 1/(1+t)")
        if np.max(X[:, 1:].max(0) - sq) > 1e-8:
            FAIL.append(f"{name}: node bound violated")
        row = {"graph": name, "n": n, "edges": G.number_of_edges(), "kappa": float(kap),
               "t_y_at_1e4": float(tl[-1] * y[-1]), "max_ratio_y_times_1_plus_t": float(np.max(ratio)),
               "max_node_constant": float(node.max()), "hub_t_x_at_1e4": float(tl[-1] * X[:, -1].max()),
               "T_at_threshold_r_1e-3": float(np.interp(-1e-3, -y, tl))}
        # above threshold transition times, r = 1e-3, R0 = 1.2 and 2
        for R0 in (1.2, 2.0):
            BR = R0 * B
            p = endemic(BR)
            bnd = transition_bounds(BR, p, R0, 1e-3)
            lemma5 = 2 / (R0 - 1) * np.log((p.max() - 1e-3) / 1e-3)
            TD_obs, TG_obs = observed_T(BR, p, 1e-3, horizon=4 * lemma5 + 300)
            obs = max(TD_obs, TG_obs)
            if not (TD_obs <= bnd["T_decay_bound"] * (1 + 1e-9) and TG_obs <= bnd["T_growth_bound"] * (1 + 1e-9)):
                FAIL.append(f"{name}: R0={R0} observed T {obs} exceeds a rigorous bound")
            row[f"R0_{R0}"] = {"observed_T": obs, "observed_T_decay": TD_obs, "observed_T_growth": TG_obs,
                               "persoons_lemma5_growth": float(lemma5),
                               "persoons_lemma5_exceeded": bool(TG_obs > lemma5 + 0.01),
                               "lemma5_hypotheses_hold": bool(1e-3 <= bnd["p_min"]),
                               **{k: float(val) for k, val in bnd.items()}}
        rows.append(row)
        ax.semilogx(tl[1:], ratio[1:], label=f"{name}  ($\\kappa$={kap:.2f})")
    ax.axhline(1.0, color="k", lw=1.2, ls="--", label="proved envelope (Theorem 3): $(1+t)\\,y=1$")
    ax.set_xlabel("time $t$ (units of the mean recovery time)")
    ax.set_ylabel("$(1+t)\\,y(t)$ at the threshold $R_0=1$")
    ax.set_ylim(0, 1.08)
    ax.set_title("At the threshold every network obeys $y\\leq 1/(1+t)$;\nthe long-run constant $\\kappa(G)$ depends on its shape", fontsize=10)
    ax.legend(fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig_threshold_decay.pdf"))
    fig.savefig(os.path.join(outdir, "fig_threshold_decay.png"), dpi=200)
    # kappa across seeds and sizes, with localisation (N * IPR), and Persoons et al.'s N=50 ER ensembles
    def kap_ipr(G):
        G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
        A = nx.to_numpy_array(G)
        _, v = perron(A)
        v = v / np.linalg.norm(v)
        return v.sum() * (v**2).sum() / (len(v) * (v**3).sum()), len(v) * float((v**4).sum())

    fam = {
        "Erdos-Renyi (mean degree 10)": lambda n, sd: nx.erdos_renyi_graph(n, 10 / (n - 1), seed=sd),
        "Watts-Strogatz small world": lambda n, sd: nx.watts_strogatz_graph(n, 10, 0.1, seed=sd),
        "Random geometric (2-D)": lambda n, sd: nx.random_geometric_graph(n, 0.06 * np.sqrt(1000 / n), seed=sd),
        "Holme-Kim (scale-free, clustered)": lambda n, sd: nx.powerlaw_cluster_graph(n, 3, 0.3, seed=sd),
        "Barabasi-Albert (m=2)": lambda n, sd: nx.barabasi_albert_graph(n, 2, seed=sd),
    }
    seeds = {}
    for name, mk in fam.items():
        ks = [kap_ipr(mk(1000, sd)) for sd in range(101, 111)]
        sizes = {n: kap_ipr(mk(n, 7))[0] for n in (250, 500, 1000, 2000)}
        seeds[name] = {"kappa_median": float(np.median([k for k, _ in ks])),
                       "kappa_min": float(min(k for k, _ in ks)), "kappa_max": float(max(k for k, _ in ks)),
                       "NIPR_median": float(np.median([q for _, q in ks])), "kappa_by_size": sizes}
    rng = np.random.default_rng(50)
    ens = []
    while len(ens) < 200:
        G = nx.erdos_renyi_graph(50, float(rng.uniform(0, 1)), seed=int(rng.integers(1 << 30)))
        if G.number_of_edges() and nx.is_connected(G):
            ens.append(kap_ipr(G)[0])
    persoons_ens = {"n_graphs": len(ens), "kappa_median": float(np.median(ens)),
                    "kappa_q10": float(np.quantile(ens, 0.1)), "kappa_min": float(min(ens))}
    summary = {"script": "illustrations.py", "numpy": np.__version__, "networkx": nx.__version__,
               "rows": rows, "kappa_across_seeds_and_sizes": seeds,
               "persoons_N50_ER_ensemble": persoons_ens, "failures": FAIL, "status": "PASS" if not FAIL else "FAIL"}
    with open(os.path.join(outdir, "illustrations.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")
    for row in rows:
        print(f"{row['graph'][:44]:44s} n={row['n']:5d} kappa={row['kappa']:.3f} ty={row['t_y_at_1e4']:.3f} "
              f"node={row['max_node_constant']:.2f} T(1e-3)={row['T_at_threshold_r_1e-3']:.0f} | "
              + " | ".join(f"R0={R0}: TD {row[f'R0_{R0}']['observed_T_decay']:.1f}/{row[f'R0_{R0}']['T_decay_bound']:.0f} "
                           f"TG {row[f'R0_{R0}']['observed_T_growth']:.1f}/{row[f'R0_{R0}']['T_growth_bound']:.0f} "
                           f"L5 {row[f'R0_{R0}']['persoons_lemma5_growth']:.1f}" for R0 in (1.2, 2.0)))
    print(json.dumps({"status": summary["status"], "failures": FAIL}))
    if FAIL:
        sys.exit(1)


if __name__ == "__main__":
    main()
