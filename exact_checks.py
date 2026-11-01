"""Exact symbolic checks for the NIMFA envelope paper (SymPy, rational arithmetic).

Every identity used in a written proof that is purely algebraic is re-derived here.
No `assert` statements are used: each check is recorded by `check()`, and the script
exits with status 1 if any positive check fails or any negative control is NOT
detected. It therefore behaves identically under `python -O`.

Negative controls are deliberately false variants of the identities; the script
requires each of them to be rejected.

Fault injection (for testing the harness): set NIMFA_INJECT_FAULT=<check name>
to perturb that check; the run must then exit with status 1.

Usage: python3 -I -O exact_checks.py [--out exact_checks.json]
"""
import json
import os
import sys

import sympy as sp

INJECT = os.environ.get("NIMFA_INJECT_FAULT", "")
RESULTS = []


def check(name, ok, detail=""):
    if name == INJECT:
        ok = not ok
    RESULTS.append({"name": name, "kind": "positive", "pass": bool(ok), "detail": str(detail)})


def control(name, rejected, detail=""):
    """A negative control passes when the false statement is rejected."""
    if name == INJECT:
        rejected = not rejected
    RESULTS.append({"name": name, "kind": "negative_control", "pass": bool(rejected), "detail": str(detail)})


def is_zero(expr):
    return sp.simplify(sp.together(sp.expand(expr))) == 0


SAMPLE = {sp.Symbol("s"): sp.Rational(1, 3), sp.Symbol("r"): sp.Rational(5, 7), sp.Symbol("d"): sp.Rational(3, 2),
          sp.Symbol("m"): sp.Rational(2, 5), sp.Symbol("p_i"): sp.Rational(1, 4), sp.Symbol("t"): sp.Rational(1, 2)}


def differs(expr):
    """Robust negative-control test: the expression is nonzero at an exact rational sample point."""
    val = sp.nsimplify(expr.subs(SAMPLE))
    return sp.simplify(val) != 0


s, r, d, b, t, a, m, c = sp.symbols("s r d b t a m c", real=True)
delta = sp.symbols("delta", positive=True)
h = -s - sp.log(1 - s)
hp = sp.diff(h, s)
q = d * s**2 / (1 - s) - r * s**2

# E1. h'(s) = s/(1-s); convexity identity of Lemma 3(b): d/ds (q'/h') = d/(1-s)^2 + 2r.
check("E1a_hprime", is_zero(hp - s / (1 - s)))
ratio = sp.simplify(sp.diff(q, s) / hp)
check("E1b_ratio", is_zero(ratio - (d * (2 - s) / (1 - s) - 2 * r * (1 - s))))
check("E1c_convexity", is_zero(sp.diff(ratio, s) - (d / (1 - s) ** 2 + 2 * r)))
control("N1_convexity_wrong_power", differs(sp.diff(ratio, s) - (d / (1 - s) + 2 * r)))
# Right limit of q'/h' at 0 is 2(d - r).
check("E1d_limit_at_zero", is_zero(sp.limit(ratio, s, 0, "+") - 2 * (d - r)))

# E2. Entropy dissipation along one coordinate:
# h'(s) * (-delta s + (1-s) b) = -delta s^2/(1-s) + s b.
check("E2_dissipation", is_zero(hp * (-delta * s + (1 - s) * b) - (-delta * s**2 / (1 - s) + s * b)))

# E3. If z' = (r-d) z - r z^2 then K = h(z) satisfies K' = -q(z).
zdot = (r - d) * s - r * s**2
check("E3_envelope_ode", is_zero(hp * zdot + q))
control("N3_envelope_wrong_sign", differs(hp * zdot - q))

# E4. Closed forms (3.3) solve the Riccati equation with z(0)=a.
g = r - d
z_ne = a * sp.exp(g * t) / (1 + r * a * (sp.exp(g * t) - 1) / g)
check("E4a_closed_form_r_ne_d", is_zero(sp.diff(z_ne, t) - ((r - d) * z_ne - r * z_ne**2)))
check("E4b_initial_value", is_zero(z_ne.subs(t, 0) - a))
z_eq = a / (1 + r * a * t)
check("E4c_closed_form_r_eq_d", is_zero(sp.diff(z_eq, t) - (-r * z_eq**2)))
# Subcritical, a = 1, d = 1, r = rho: z = (1-rho)/(e^{(1-rho)t} - rho).
rho = sp.symbols("rho", positive=True)
z_sub = (1 - rho) / (sp.exp((1 - rho) * t) - rho)
check("E4d_subcritical_form", is_zero(z_sub - z_ne.subs({a: 1, d: 1, r: rho})))
check("E4e_subcritical_to_critical", is_zero(sp.limit(z_sub, rho, 1) - 1 / (1 + t)))
rr = sp.symbols("rr", positive=True)  # accuracy tolerance
T_sub = sp.log(1 + (1 - rho) * (1 / rr - 1)) / (1 - rho)
check("E4f_T_sub_solves", is_zero(sp.simplify(z_sub.subs(t, T_sub)) - rr))
check("E4g_T_sub_limit", is_zero(sp.limit(T_sub, rho, 1) - (1 / rr - 1)))

# E5. Residual transform (4.3): x = p + s_i z with B p = p/(1-p) gives
# zdot = -z/s_i + (1-z)(B S z)_i.  Check componentwise with a symbolic 3x3 symmetric B.
n = 3
Bs = sp.Matrix(n, n, lambda i, j: sp.Symbol(f"B{min(i, j)}{max(i, j)}") if i != j else 0)
pv = sp.Matrix(sp.symbols("p0:3", positive=True))
zv = sp.Matrix(sp.symbols("z0:3"))
sv = sp.Matrix([1 - pv[i] for i in range(n)])
Bp = Bs * pv
eq_sub = {}
# Impose B p = p/(1-p) by solving for the three off-diagonal entries.
sol = sp.solve([Bp[i] - pv[i] / (1 - pv[i]) for i in range(n)], [Bs[0, 1], Bs[0, 2], Bs[1, 2]], dict=True)
check("E5a_equilibrium_solvable", len(sol) == 1)
Bsub = Bs.subs(sol[0])
xv = sp.Matrix([pv[i] + sv[i] * zv[i] for i in range(n)])
Bx = Bsub * xv
BSz = Bsub * sp.Matrix([sv[i] * zv[i] for i in range(n)])
ok = True
for i in range(n):
    lhs = (-xv[i] + (1 - xv[i]) * Bx[i]) / sv[i]
    rhs = -zv[i] / sv[i] + (1 - zv[i]) * BSz[i]
    ok = ok and is_zero(lhs - rhs)
check("E5b_residual_transform", ok)

# E6. Endemic-ray subsolution: with (Bp)_i = p_i/s_i, F_i(m p) = m(1-m) p_i^2/s_i.
pi, si = sp.symbols("p_i s_i", positive=True)
F_mp = -m * pi + (1 - m * pi) * m * (pi / si)
check("E6_endemic_ray", is_zero(F_mp.subs(si, 1 - pi) - (m * (1 - m) * pi**2 / (1 - pi))))
control("N6_endemic_ray_wrong", differs(F_mp.subs(si, 1 - pi) - (m * (1 - m) * pi / (1 - pi))))

# E7. Perron-ray subsolution: with (Bv)_i = R v_i, F_i(c v) = c v_i (R - 1 - R c v_i).
R, vi = sp.symbols("R v_i", positive=True)
F_cv = -c * vi + (1 - c * vi) * c * R * vi
check("E7_perron_ray", is_zero(F_cv - c * vi * (R - 1 - R * c * vi)))

# E8. Diagonal-stability identity: if S is symmetric with S v = 0, then
# x^T S x = -1/2 sum_{i,j} S_ij v_i v_j (x_i/v_i - x_j/v_j)^2.
off = {(0, 1): sp.Symbol("S01"), (0, 2): sp.Symbol("S02"), (1, 2): sp.Symbol("S12")}
vv = sp.Matrix(sp.symbols("v0:3", positive=True))
xx = sp.Matrix(sp.symbols("x0:3"))
S = sp.zeros(n, n)
for (i, j), val in off.items():
    S[i, j] = val
    S[j, i] = val
for i in range(n):
    S[i, i] = -sum(S[i, j] * vv[j] for j in range(n) if j != i) / vv[i]
check("E8a_kernel", all(is_zero(e) for e in (S * vv)))
quad = (xx.T * S * xx)[0]
rhs = -sp.Rational(1, 2) * sum(S[i, j] * vv[i] * vv[j] * (xx[i] / vv[i] - xx[j] / vv[j]) ** 2
                              for i in range(n) for j in range(n) if i != j)
check("E8b_quadratic_identity", is_zero(quad - rhs))

# E9. Star K_{1,m}: Perron vector (sqrt(m), 1, ..., 1); kappa = 2 sqrt(m)/(m+1);
# hub constant v_hub sum v^2 / sum v^3 = 2 sqrt(m)/(sqrt(m)+1).
mm = sp.symbols("mm", positive=True)
v_hub, v_leaf = sp.sqrt(mm), sp.Integer(1)
S1 = v_hub + mm * v_leaf
S2 = v_hub**2 + mm * v_leaf**2
S3 = v_hub**3 + mm * v_leaf**3
check("E9a_star_eigen", is_zero(mm * v_leaf - sp.sqrt(mm) * v_hub) and is_zero(v_hub - sp.sqrt(mm) * v_leaf))
check("E9b_star_kappa", is_zero(S1 * S2 / ((mm + 1) * S3) - 2 * sp.sqrt(mm) / (mm + 1)))
check("E9c_star_hub", is_zero(v_hub * S2 / S3 - 2 * sp.sqrt(mm) / (sp.sqrt(mm) + 1)))

# E10. Hub joined to a d-regular graph on N vertices: Perron vector (1, eps,...,eps),
# lambda = N eps, lambda eps = 1 + d eps, lambda^2 - d lambda - N = 0, lambda - d = N/lambda.
N, dd = sp.symbols("N dd", positive=True)
lam = (dd + sp.sqrt(dd**2 + 4 * N)) / 2
eps = lam / N
check("E10a_hub_eq1", is_zero(lam * eps - (1 + dd * eps)))
check("E10b_quadratic", is_zero(lam**2 - dd * lam - N))
check("E10c_gap", is_zero((lam - dd) - N / lam))
# Numerical instance of the lower bound R >= (d^2/N)/(1 + (d + N/d)^3/N^2) for N=2000, d=158.
Nn, dn = 2000, 158
lam_n = (dn + (dn**2 + 4 * Nn) ** 0.5) / 2
eps_n = lam_n / Nn
R_exact = (1 + Nn * eps_n**2) / (1 + Nn * eps_n**3)
R_lower = (dn**2 / Nn) / (1 + (dn + Nn / dn) ** 3 / Nn**2)
check("E10d_hub_bound_instance", R_exact >= R_lower > 1, f"R={R_exact:.4f} lower={R_lower:.4f}")

# E11. Directed two-node counterexample B=[[0,4],[1/16,0]]: rho = 1/2, mean second derivative 97/32,
# and the Perron weights w = u/v = (1/8, 8) satisfy W B + B^T W <= 2 rho W (i.e. x^T W B x <= rho x^T W x).
Bd = sp.Matrix([[0, 4], [sp.Rational(1, 16), 0]])
eigs = Bd.eigenvals()
check("E11a_rho_half", max(eigs, key=lambda e: sp.re(e)) == sp.Rational(1, 2))
x1, x2 = sp.symbols("x1 x2")
X = sp.Matrix([x1, x2])
F = sp.Matrix([-X[i] + (1 - X[i]) * (Bd * X)[i] for i in range(2)])
J = F.jacobian(X)
Fdd = J * F  # second time derivative of x along the flow
ydd = sp.Rational(1, 2) * sum(Fdd.subs({x1: 1, x2: 1}))
check("E11b_second_derivative", ydd == sp.Rational(97, 32))
control("N11_second_derivative_not_two", ydd != 2)
W = sp.diag(sp.Rational(1, 8), 8)
M = W * Bd + Bd.T * W - W  # should be NSD: 2 x^T W B x <= x^T W x
check("E11c_weighted_nsd", all(sp.re(e) <= 0 for e in M.eigenvals()))

# E12. Arithmetic-mean counterexample on the star K_{1,2} at threshold: y'(0) = (sqrt2 - 1)/6 > 0.
Bst = sp.Matrix([[0, 1, 1], [1, 0, 0], [1, 0, 0]]) / sp.sqrt(2)
x0 = sp.Matrix([sp.Rational(1, 2), 0, 0])
yd = sum(-x0[i] + (1 - x0[i]) * (Bst * x0)[i] for i in range(3)) / 3
check("E12_arith_mean", is_zero(yd - (sp.sqrt(2) - 1) / 6) and yd > 0)

# E13. h bounds: s^2/2 <= h(s) <= s^2/(2(1-s)) on [0,1): derivatives of the differences.
check("E13a_lower", is_zero(sp.diff(h - s**2 / 2, s) - s**2 / (1 - s)))
check("E13b_upper", is_zero(sp.diff(s**2 / (2 * (1 - s)) - h, s) - s**2 / (2 * (1 - s) ** 2)))

# E14. Chebyshev: n sum v^3 - sum v sum v^2 = 1/2 sum_{i,j} (v_i - v_j)(v_i^2 - v_j^2).
v4 = sp.symbols("w0:4", positive=True)
lhs = 4 * sum(x**3 for x in v4) - sum(v4) * sum(x**2 for x in v4)
rhs = sp.Rational(1, 2) * sum((v4[i] - v4[j]) * (v4[i] ** 2 - v4[j] ** 2) for i in range(4) for j in range(4))
check("E14_chebyshev", is_zero(lhs - rhs))

# E15. Regular endemic exactness: on a rho-regular graph the residual satisfies
# z' = -rho z + (1-z) z = (1-rho) z - z^2, i.e. the envelope with r=1, d=rho.
zz = sp.symbols("zz")
check("E15_regular_endemic", is_zero((-rho * zz + (1 - zz) * zz) - ((1 - rho) * zz - zz**2)))

# E16. Exponential form in (4.2): envelope with r=1,d=D,a=1 equals (D-1)/(D e^{(D-1)t} - 1).
Dd = sp.symbols("Dd", positive=True)
check("E16_endemic_exp", is_zero(z_ne.subs({a: 1, r: 1, d: Dd}) - (Dd - 1) / (Dd * sp.exp((Dd - 1) * t) - 1)))

# E17. Logistic solution of m' = k m(1-m), m(0)=m0: 1-m = (1-m0)/(1-m0+m0 e^{kt}).
k, m0 = sp.symbols("k m0", positive=True)
mt = m0 * sp.exp(k * t) / (1 - m0 + m0 * sp.exp(k * t))
check("E17_logistic", is_zero(sp.diff(mt, t) - k * mt * (1 - mt)) and is_zero(1 - mt - (1 - m0) / (1 - m0 + m0 * sp.exp(k * t))))

# E18. Hoelder step of Prop 9(d): with max v = 1, sum v^2 / sum v^3 <= n^{1/3}; exact check on rational vectors.
import random as _random
_rng = _random.Random(20261101)
ok18 = True
for _ in range(200):
    nn = _rng.randint(2, 12)
    vv18 = [sp.Rational(_rng.randint(1, 50), 50) for _ in range(nn - 1)] + [sp.Integer(1)]
    lhs18 = (sum(x**2 for x in vv18) / sum(x**3 for x in vv18)) ** 3
    ok18 = ok18 and bool(lhs18 <= nn)
check("E18_hoelder_node_order", ok18)


def main():
    out = "exact_checks.json"
    if "--out" in sys.argv:
        out = sys.argv[sys.argv.index("--out") + 1]
    failed = [x["name"] for x in RESULTS if not x["pass"]]
    summary = {
        "script": "exact_checks.py",
        "sympy": sp.__version__,
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
