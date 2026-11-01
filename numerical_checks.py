"""Numerical checks of every theorem in the NIMFA envelope paper (NumPy/SciPy/NetworkX).

Numerical integration supplements, and does not replace, the written proofs.
No `assert` statements are used; results are recorded by `check()` / `control()`
and the script exits with status 1 if any check fails or any negative control is
NOT detected, so it behaves identically under `python -O`.

Negative controls are claims the paper says are FALSE (directed unweighted bound,
arithmetic-mean entropy parameter, node-wise 1/(1+t), constants below the sharp
ones); each must be observed to fail.

Fault injection: NIMFA_INJECT_FAULT=<check name> flips that check (harness test).

Usage: python3 -I -O numerical_checks.py [--out numerical_checks.json]
"""
import json
import os
import sys

import networkx as nx
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

INJECT = os.environ.get("NIMFA_INJECT_FAULT", "")
RESULTS = []
TOL = 1e-8  # absolute slack for integration error in one-sided inequalities


def _num(v):
    if isinstance(v, str):
        return v
    if np.isscalar(v):
        return float(v)
    return [float(x) for x in v]


def check(name, ok, **detail):
    ok = bool(ok) and all(np.isfinite(v) for v in detail.values() if np.isscalar(v) and not isinstance(v, str))
    if name == INJECT:
        ok = not ok
    RESULTS.append({"name": name, "kind": "positive", "pass": ok, **{k: _num(v) for k, v in detail.items()}})


def control(name, rejected, **detail):
    rejected = bool(rejected)
    if name == INJECT:
        rejected = not rejected
    RESULTS.append({"name": name, "kind": "negative_control", "pass": rejected, **{k: _num(v) for k, v in detail.items()}})


def wmax(*vals):
    """max() that propagates NaN instead of silently dropping it."""
    vals = [float(v) for v in vals]
    return float("nan") if any(np.isnan(v) for v in vals) else max(vals)


def h(x):
    x = np.asarray(x, float)
    return -x - np.log1p(-x)


def hinv(u):
    if u <= 0:
        return 0.0
    return brentq(lambda x: float(h(x)) - u, 0.0, 1 - 1e-15, xtol=1e-16)


def envelope(t, r, d, a=1.0):
    """Closed form (3.3) of z' = (r-d) z - r z^2, z(0)=a, constant r,d."""
    t = np.asarray(t, float)
    g = r - d
    if abs(g) < 1e-14:
        return a / (1 + r * a * t)
    return a * np.exp(g * t) / (1 + r * a * np.expm1(g * t) / g)


def integrate(rhs, x0, ts, rtol=1e-10, atol=1e-14):
    sol = solve_ivp(rhs, (ts[0], ts[-1]), x0, t_eval=ts, method="LSODA", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y


def sis(B, delta=None):
    if delta is None:
        return lambda _t, x: -x + (1 - x) * (B @ x)
    return lambda _t, x: -delta * x + (1 - x) * (B @ x)


def perron_sym(B):
    w, V = np.linalg.eigh(B)
    return w[-1], np.abs(V[:, -1])


def perron_pair(M):
    """Right/left Perron vectors of an irreducible nonnegative (or Metzler) matrix."""
    ev, Rv = np.linalg.eig(M)
    i = int(np.argmax(ev.real))
    ev2, Lv = np.linalg.eig(M.T)
    j = int(np.argmax(ev2.real))
    return ev[i].real, np.abs(Rv[:, i].real), np.abs(Lv[:, j].real)


def endemic(B):
    p = np.ones(len(B))
    for _ in range(10**6):
        bp = B @ p
        nxt = bp / (1 + bp)
        if np.max(np.abs(nxt - p)) < 1e-15:
            return nxt
        p = nxt
    raise RuntimeError("equilibrium iteration did not converge")


def kappa(v):
    return v.sum() * (v**2).sum() / (len(v) * (v**3).sum())


def rand_sym(rng, n, density=0.4, weighted=True):
    while True:
        A = (rng.random((n, n)) < density) * (rng.uniform(0.2, 1.5, (n, n)) if weighted else 1.0)
        A = np.triu(A, 1)
        A = A + A.T
        if nx.is_connected(nx.from_numpy_array(A)):
            return A


def rand_dir(rng, n, density=0.3):
    while True:
        A = (rng.random((n, n)) < density) * rng.uniform(0.1, 3.0, (n, n))
        np.fill_diagonal(A, 0)
        if nx.is_strongly_connected(nx.DiGraph(A)):
            return A


def main():
    rng = np.random.default_rng(20261101)
    ts = np.linspace(0, 40, 801)

    # N1. Theorem 1 / Theorem 3, static symmetric, below/at/above threshold, entropy parameter and a = 1.
    worst1 = worst1b = -np.inf
    for n in (3, 8, 20):
        for rho in (0.3, 0.9, 1.0, 1.2, 3.0):
            A = rand_sym(rng, n)
            B = rho * A / perron_sym(A)[0]
            X = integrate(sis(B), np.ones(n), ts)
            worst1 = wmax(worst1, np.max(X.mean(0) - envelope(ts, rho, 1.0)))
            x0 = rng.uniform(0, 0.97, n)
            a = hinv(float(h(x0).mean()))
            X2 = integrate(sis(B), x0, ts)
            Hs = h(np.clip(X2, 0, 1 - 1e-16)).mean(0)
            worst1b = wmax(worst1b, np.max(X2.mean(0) - envelope(ts, rho, 1.0, a)),
                          np.max(np.array([hinv(u) for u in Hs]) - envelope(ts, rho, 1.0, a)))
    check("N1a_static_envelope_all_infected", worst1 < TOL, max_excess=worst1)
    check("N1b_static_envelope_entropy_parameter", worst1b < TOL, max_excess=worst1b)

    # N2. Directed irreducible: weights u/v (Perron left/right) give the weighted envelope.
    worst2 = -np.inf
    worst_nsd = -np.inf
    unweighted_excess = -np.inf
    for _ in range(25):
        n = int(rng.integers(4, 14))
        A = rand_dir(rng, n)
        lam, v, u = perron_pair(A)
        for rho in (0.6, 1.0, 1.5):
            B = rho * A / lam
            w = u / v
            Msym = (w[:, None] * B + (w[:, None] * B).T) / 2 - rho * np.diag(w)
            worst_nsd = wmax(worst_nsd, np.linalg.eigvalsh(Msym).max() / np.abs(Msym).max())
            X = integrate(sis(B), np.ones(n), ts)
            worst2 = wmax(worst2, np.max(w @ X / w.sum() - envelope(ts, rho, 1.0)))
            if rho == 1.0:
                unweighted_excess = wmax(unweighted_excess, np.max(X.mean(0) - 1 / (1 + ts)))
    check("N2a_directed_weight_nsd", worst_nsd < 1e-10, max_rel_eig=worst_nsd)
    Ad = rand_dir(rng, 10)
    lam_d = perron_pair(Ad)[0]
    Bd1 = Ad / lam_d
    control("C1c_identity_weight_not_enough_directed",
            np.linalg.eigvalsh((Bd1 + Bd1.T) / 2 - np.eye(10)).max() > 1e-6)
    check("N2b_directed_weighted_envelope", worst2 < TOL, max_excess=worst2)
    control("C1_directed_unweighted_fails", unweighted_excess > 1e-4, max_excess=unweighted_excess)
    Bc = np.array([[0, 4], [1 / 16, 0]], float)
    tc = np.array([0, 1e-3, 1e-2, 5e-2])
    Xc = integrate(sis(Bc), np.ones(2), tc)
    control("C1b_two_node_counterexample", np.all(Xc.mean(0)[1:] - 1 / (1 + tc[1:]) > 0),
            excess_at_0_05=Xc.mean(0)[-1] - 1 / (1 + tc[-1]))
    wc = np.array([1 / 8, 8.0])
    Xc2 = integrate(sis(Bc), np.ones(2), ts)
    check("N2c_two_node_weighted_ok", np.max(wc @ Xc2 / wc.sum() - envelope(ts, 0.5, 1.0)) < TOL)

    # N3. Heterogeneous recovery: symmetric B with rho(D^-1 B)=1 gives y <= 1/(1+dmin t);
    # directed heterogeneous with weights from the null vectors of B - D.
    worst3 = worst3d = -np.inf
    for _ in range(12):
        n = int(rng.integers(4, 12))
        dl = rng.uniform(0.3, 3.0, n)
        A = rand_sym(rng, n)
        lam = np.max(np.linalg.eigvals(A / dl[:, None]).real)
        B = A / lam
        X = integrate(sis(B, dl), np.ones(n), ts)
        worst3 = wmax(worst3, np.max(X.mean(0) - 1 / (1 + dl.min() * ts)))
        Ad = rand_dir(rng, n)
        lam_d, v, u = perron_pair(Ad / dl[:, None])
        Bd = Ad / lam_d
        _, vv, uu = perron_pair(Bd - np.diag(dl))
        w = uu / vv
        Xd = integrate(sis(Bd, dl), np.ones(n), ts)
        worst3d = wmax(worst3d, np.max(w @ Xd / w.sum() - 1 / (1 + dl.min() * ts)))
    check("N3a_heterogeneous_recovery", worst3 < TOL, max_excess=worst3)
    check("N3b_heterogeneous_directed", worst3d < TOL, max_excess=worst3d)

    # N4. Temporal, noncommuting symmetric snapshots, time-varying recovery d(t).
    n = 7
    worst4 = -np.inf
    state = np.r_[np.ones(n), 1.0]
    rates = [0.4, 1.0, 1.7, 0.8, 2.5, 0.5, 1.0]  # snapshot spectral radii r_k
    recov = [1.0, 0.7, 1.3, 1.1, 0.9, 1.5, 0.6]  # snapshot recovery rates d_k (each segment has length 2)
    comm = 0.0
    prev = None
    for k, (rr, dk) in enumerate(zip(rates, recov)):
        A = rand_sym(rng, n)
        B = rr * A / perron_sym(A)[0]
        if prev is not None:
            comm = max(comm, np.linalg.norm(B @ prev - prev @ B))
        prev = B
        seg = np.linspace(2.0 * k, 2.0 * (k + 1), 81)

        def rhs(_t, y, B=B, rr=rr, dk=dk):
            x, z = y[:-1], y[-1]
            return np.r_[-dk * x + (1 - x) * (B @ x), (rr - dk) * z - rr * z * z]

        Y = integrate(rhs, state, seg)
        worst4 = wmax(worst4, np.max(Y[:-1].mean(0) - Y[-1]))
        state = Y[:, -1]
    check("N4_temporal_envelope", worst4 < TOL and comm > 1e-3, max_excess=worst4, commutator=comm)

    # N5. Theorem 4 (Conjecture 2): decay above the endemic state, algebraic and exponential forms.
    worst5a = worst5e = -np.inf
    for n in (4, 10, 25):
        for R0 in (1.05, 1.5, 3.0):
            A = rand_sym(rng, n)
            B = R0 * A / perron_sym(A)[0]
            p = endemic(B)
            s = 1 - p
            kap = np.min(p / s)
            for x0 in (np.ones(n), p + (1 - p) * rng.random(n)):
                X = integrate(sis(B), x0, ts)
                err = X.mean(0) - p.mean()
                lo = np.min(X - p[:, None])
                worst5a = wmax(worst5a, np.max(err - s.mean() / (1 + ts)), -lo)
                worst5e = wmax(worst5e, np.max(err - s.mean() * kap / ((1 + kap) * np.exp(kap * ts) - 1)))
    check("N5a_endemic_algebraic", worst5a < TOL, max_excess=worst5a)
    check("N5b_endemic_exponential", worst5e < TOL, max_excess=worst5e)

    # N6. Disconnected graph: one supercritical and one subcritical component.
    A1 = rand_sym(rng, 6)
    A2 = rand_sym(rng, 5)
    B = np.zeros((11, 11))
    B[:6, :6] = 2.0 * A1 / perron_sym(A1)[0]
    B[6:, 6:] = 0.7 * A2 / perron_sym(A2)[0]
    p = np.r_[endemic(B[:6, :6]), np.zeros(5)]
    X = integrate(sis(B), np.ones(11), ts)
    check("N6_disconnected", np.max(X.mean(0) - p.mean() - (1 - p.mean()) / (1 + ts)) < TOL)

    # N7. Critical asymptotics: t y(t) -> kappa(G), t x_i(t) -> v_i sum v^2/sum v^3, r T(r) -> kappa.
    tl = np.array([0.0, 10.0, 100.0, 1000.0, 20000.0])
    graphs = {
        "ER": nx.erdos_renyi_graph(200, 0.05, seed=7),
        "BA": nx.barabasi_albert_graph(200, 2, seed=7),
        "star": nx.star_graph(49),
        "cycle": nx.cycle_graph(30),
    }
    asym = {}
    for name, G in graphs.items():
        G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
        A = nx.to_numpy_array(G)
        lam, v = perron_sym(A)
        B = A / lam
        X = integrate(sis(B), np.ones(len(v)), tl, rtol=1e-11, atol=1e-16)
        k = kappa(v)
        ty = tl[-1] * X[:, -1].mean()
        node_const = v * (v**2).sum() / (v**3).sum()
        node_err = np.max(np.abs(tl[-1] * X[:, -1] - node_const) / node_const)
        asym[name] = {"kappa": k, "t_y_at_2e4": ty, "rel_err": abs(ty - k) / k, "node_rel_err": node_err,
                      "max_node_const": node_const.max()}
    check("N7a_kappa_limit", all(a["rel_err"] < 5e-3 for a in asym.values()), detail_json=json.dumps(asym))
    check("N7b_node_limits", all(a["node_rel_err"] < 5e-3 for a in asym.values()))
    check("N7c_kappa_regular_is_one", abs(asym["cycle"]["kappa"] - 1) < 1e-12)
    control("C2_conjectured_constant_not_sharp_for_BA", abs(asym["BA"]["t_y_at_2e4"] - 1.0) > 0.5,
            t_y=asym["BA"]["t_y_at_2e4"])
    control("C3_nodewise_one_over_t_fails", asym["star"]["max_node_const"] > 1.5 and
            asym["BA"]["max_node_const"] > 1.5)
    # r T(r) -> kappa, using monotone decay from u on the BA graph.
    G = graphs["BA"]
    A = nx.to_numpy_array(G)
    lam, v = perron_sym(A)
    B = A / lam
    tt = np.geomspace(1, 2e4, 400)
    tt = np.r_[0.0, tt]
    Y = integrate(sis(B), np.ones(len(v)), tt, rtol=1e-11, atol=1e-16).mean(0)
    mono = np.all(np.diff(Y) <= 1e-15)
    rT = []
    for rtol_acc in (1e-2, 1e-3, 1e-4):
        i = int(np.argmax(Y <= rtol_acc))
        # linear interpolation for the crossing time
        T = tt[i - 1] + (tt[i] - tt[i - 1]) * (Y[i - 1] - rtol_acc) / (Y[i - 1] - Y[i])
        rT.append(rtol_acc * T)
    check("N7d_rT_to_kappa", mono and abs(rT[-1] - kappa(v)) / kappa(v) < 2e-2, rT=rT, kappa=kappa(v))

    # N8. Node-wise universal bound x_i <= sqrt(n/(t(1+t))) at/below threshold.
    worst8 = -np.inf
    for name, G in graphs.items():
        G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
        A = nx.to_numpy_array(G)
        B = A / perron_sym(A)[0]
        tq = np.linspace(0.05, 50, 400)
        X = integrate(sis(B), np.ones(len(A)), np.r_[0.0, tq])[:, 1:]
        worst8 = wmax(worst8, np.max(X.max(0) - np.sqrt(len(A) / (tq * (1 + tq)))))
    check("N8_nodewise_bound", worst8 < TOL, max_excess=worst8)

    # N9. Growth from below (two-phase subsolution), sandwich for mixed initial states,
    # and the combined transition-time bound, symmetric and directed.
    worst9 = -np.inf
    tg = np.linspace(0, 600, 3001)
    for directed in (False, True):
        for R0 in (1.1, 1.5, 3.0):
            n = 12
            A = rand_dir(rng, n) if directed else rand_sym(rng, n)
            lam, v, _u = perron_pair(A)
            B = R0 * A / lam
            p = endemic(B)
            s = 1 - p
            kap = np.min(p / s)
            v = v / v.max()
            acc = 1e-3
            cstar = (R0 - 1) / R0
            g = R0 - 1
            c = cstar / (1 + (cstar / acc - 1) * np.exp(-g * tg))  # overflow-free logistic
            Xr = integrate(sis(B), acc * np.ones(n), tg)
            phase1 = np.max(c[None, :] * v[:, None] - Xr)
            t1 = np.log((cstar - acc) / acc) / g
            m1 = min(1.0, (cstar / 2) * np.min(v / p))
            tau = np.clip(tg - t1, 0, None)
            mt = 1 / (1 + (1 / m1 - 1) * np.exp(-kap * tau))  # overflow-free logistic
            later = tg >= t1
            phase2 = np.max((mt[None, later] * p[:, None]) - Xr[:, later])
            gap = np.max((p.mean() - Xr.mean(0))[later] - p.mean() * (1 - mt[later]))
            # mixed initial state in [acc, 1]^n sandwiched between acc*u and u
            x0 = rng.uniform(acc, 1, n)
            Xm = integrate(sis(B), x0, tg)
            Xu = integrate(sis(B), np.ones(n), tg)
            sandwich = max(np.max(Xr - Xm), np.max(Xm - Xu))
            parts = np.array([phase1, phase2, gap, sandwich])
            worst9 = wmax(worst9, np.max(parts)) if np.all(np.isfinite(parts)) else np.inf
    check("N9_growth_subsolution_and_sandwich", worst9 < TOL, max_excess=worst9)

    # N10. SIR and SIRS: infected mean obeys the SIS envelopes.
    worst10 = -np.inf
    for xi in (0.0, 0.5):
        n = 15
        A = rand_sym(rng, n)
        for rho in (0.8, 1.0, 2.0):
            B = rho * A / perron_sym(A)[0]

            def rhs(_t, zv, B=B, xi=xi):
                S, I = zv[:n], zv[n:]
                Rr = 1 - S - I
                f = S * (B @ I)
                return np.r_[-f + xi * Rr, -I + f]

            I0 = rng.uniform(0.2, 0.95, n)
            a0 = hinv(float(h(I0).mean()))
            Z = integrate(rhs, np.r_[1 - I0, I0], ts)
            worst10 = wmax(worst10, np.max(Z[n:].mean(0) - envelope(ts, rho, 1.0, a0)))
    check("N10_sir_sirs", worst10 < TOL, max_excess=worst10)

    # N11. Sharpness: regular-graph equality; subcritical exact exponential rate.
    Kn = np.ones((6, 6)) - np.eye(6)
    eq_err = 0.0
    for rho in (0.5, 1.0, 2.0):
        X = integrate(sis(rho * Kn / 5), np.ones(6), ts)
        eq_err = wmax(eq_err, np.max(np.abs(X.mean(0) - envelope(ts, rho, 1.0))))
    check("N11a_regular_equality", eq_err < 1e-8, max_abs=eq_err)
    X = integrate(sis(Kn / 5), np.ones(6), ts)
    control("C4_constant_below_one_fails", np.any(X.mean(0) > 0.99 / (1 + ts) + 1e-6))
    A = nx.to_numpy_array(nx.barabasi_albert_graph(60, 2, seed=3))
    rho = 0.7
    B = rho * A / perron_sym(A)[0]
    tsub = np.linspace(0, 80, 161)
    X = integrate(sis(B), np.ones(60), tsub, rtol=1e-11, atol=1e-30)
    scaled = np.exp((1 - rho) * tsub) * X.mean(0)
    check("N11b_subcritical_rate_exact", scaled[-1] > 0 and abs(scaled[-1] / scaled[-21] - 1) < 1e-3,
          ratio=scaled[-1] / scaled[-21])
    faster = np.exp((1 - rho + 1e-3) * tsub) * X.mean(0)
    control("C5_faster_rate_fails", faster[-1] / faster[-21] > 1.005, ratio=faster[-1] / faster[-21])
    kb = asym["BA"]
    control("C7_kappa_one_percent_low_rejected",
            abs(kb["t_y_at_2e4"] - 0.99 * kb["kappa"]) > 5 * abs(kb["t_y_at_2e4"] - kb["kappa"]))

    # N7e. Theorem 7, directed and heterogeneous-recovery limits:
    # t x_i -> v_i (u^T v) / sum_j u_j delta_j v_j^2 with B v = Delta v, u^T B = u^T Delta.
    worst7e = -np.inf
    tq = np.array([0.0, 2e4])
    for kind in ("directed", "heterogeneous"):
        n = 12
        dl = np.ones(n) if kind == "directed" else rng.uniform(0.5, 2.0, n)
        A = rand_dir(rng, n) if kind == "directed" else rand_sym(rng, n)
        lamh = np.max(np.linalg.eigvals(A / dl[:, None]).real)
        B = A / lamh
        _, v, u = perron_pair(B - np.diag(dl))
        X = integrate(sis(B, dl), np.ones(n), tq, rtol=1e-11, atol=1e-16)
        lim = v * (u @ v) / np.sum(u * dl * v**2)
        worst7e = wmax(worst7e, np.max(np.abs(tq[-1] * X[:, -1] - lim) / lim))
    check("N7e_kappa_directed_heterogeneous", worst7e < 5e-3, max_rel_err=worst7e)

    # N12. Theorem 14 end to end: for t >= max(T_D, T_G), |y(t) - y_inf| <= r for u, r u and mixed V(0),
    # across r on both sides of c*/2 and above pbar; regular-graph exactness of T_G.
    def t13(B, R0, r):
        p = endemic(B)
        s_ = 1 - p
        kap = np.min(p / s_)
        sbar, pbar = s_.mean(), p.mean()
        TD = 0.0 if r >= sbar else max(0.0, np.log((1 + sbar * kap / r) / (1 + kap)) / kap)
        _, v = perron_sym(B)
        v = v / v.max()
        cst = (R0 - 1) / R0
        t1 = max(0.0, np.log((cst - r) / r) / (R0 - 1)) if r < cst else 0.0
        m1 = (cst / 2) * np.min(v / p)
        t2 = 0.0 if pbar * (1 - m1) <= r else np.log((1 - m1) * (pbar - r) / (m1 * r)) / kap
        return p, TD, t1 + t2, m1

    worst12 = -np.inf
    m1max = 0.0
    tt = np.linspace(0, 300, 6001)
    for R0 in (1.5, 2.0, 4.0):
        for _ in range(3):
            n = int(rng.integers(3, 9))
            A = rand_sym(rng, n)
            B = R0 * A / perron_sym(A)[0]
            cst = (R0 - 1) / R0
            for frac in (0.01, 0.3, 0.49, 0.51, 0.7, 0.95):
                r = frac * cst
                p, TD, TG, m1 = t13(B, R0, r)
                m1max = max(m1max, m1)
                T = max(TD, TG)
                late = tt >= T
                for x0 in (np.ones(n), r * np.ones(n), rng.uniform(r, 1, n)):
                    Y = integrate(sis(B), x0, tt).mean(0)
                    worst12 = wmax(worst12, np.max(np.abs(Y[late] - p.mean()) - r))
    check("N12a_T14_end_to_end", worst12 < TOL and m1max <= 0.5 + 1e-12, max_excess=worst12, m1_max=m1max)
    # regular graph (cycle C_8, R0 = 2): T_G = 2/(R0-1) log^+((p-r)/r), and it is the exact growth time
    C8 = nx.to_numpy_array(nx.cycle_graph(8))
    reg_err = 0.0
    for r in (1e-3, 0.1, 0.3, 0.5, 0.6):  # p = 1/2: below p/2, between p/2 and p, at p, above p
        B = 2.0 * C8 / 2.0
        p, TD, TG, _ = t13(B, 2.0, r)
        exact = 2 * np.log((p[0] - r) / r) if r < p[0] / 2 else 0.0
        tf = np.linspace(0, max(exact, 1e-9) + 5, 20001)
        Y = integrate(sis(B), r * np.ones(8), tf).mean(0)
        bad = np.nonzero(p.mean() - Y > r)[0]
        obs = tf[bad[-1] + 1] if len(bad) else 0.0
        reg_err = max(reg_err, abs(TG - exact), abs(obs - exact) - (tf[1] - tf[0]))
    check("N12b_T14_regular_exact", reg_err < 1e-6 + 1e-9, max_err=reg_err)

    # C6. Arithmetic initial mean in place of the entropy parameter fails (star K_{1,2} at threshold).
    Bst = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]], float) / np.sqrt(2)
    x0 = np.array([0.5, 0.0, 0.0])
    tsm = np.linspace(0, 0.05, 11)
    X = integrate(sis(Bst), x0, tsm)
    y0 = x0.mean()
    control("C6_arithmetic_mean_fails", np.any(X.mean(0)[1:] > y0 / (1 + y0 * tsm[1:])))

    out = "numerical_checks.json"
    if "--out" in sys.argv:
        out = sys.argv[sys.argv.index("--out") + 1]
    failed = [x["name"] for x in RESULTS if not x["pass"]]
    summary = {
        "script": "numerical_checks.py",
        "numpy": np.__version__,
        "networkx": nx.__version__,
        "python_optimize_flag": sys.flags.optimize,
        "n_positive": sum(x["kind"] == "positive" for x in RESULTS),
        "n_negative_controls": sum(x["kind"] == "negative_control" for x in RESULTS),
        "failed": failed,
        "status": "PASS" if not failed else "FAIL",
        "results": RESULTS,
    }
    with open(out, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")
    print(json.dumps({k: summary[k] for k in ("status", "n_positive", "n_negative_controls", "failed")}))
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
