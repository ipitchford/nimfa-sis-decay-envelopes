# SPDX-License-Identifier: MIT
"""Clean replay of the package (read-only on ROOT), run in a temporary copy.

Usage: python -I -B tools/replay_compare.py ROOT [--quick]

Runs every checker with the interpreter running this script and flags -I -O -B
(the same commands as build.sh, without the PDF build), each writing its retained
output file names in the copy. Requires every checker to exit 0, then compares each
regenerated JSON output with the retained one (including JSON-encoded detail strings):
numbers within rel 1e-7 / abs 1e-10,
everything else exactly; package-version fields are ignored. Figures are
regenerated, not byte-compared. Fault injection (NIMFA_INJECT_FAULT) must make
exact_checks.py and numerical_checks.py fail. Fails if a checker writes an
unexpected file or if ROOT changes.

--quick runs only the exact and numerical checkers (seconds) and fault injection;
the full run adds eq20_check.py (about 2 minutes) and illustrations.py (about 5).
"""
import hashlib, json, math, os, shutil, subprocess, sys, tempfile
from pathlib import Path

root = Path(sys.argv[1]).resolve()
quick = '--quick' in sys.argv
VOLATILE = {'numpy', 'scipy', 'sympy', 'networkx', 'mpmath', 'matplotlib'}
FULL = [('exact_checks.py', ['--out', 'exact_checks.json']),
        ('numerical_checks.py', ['--out', 'numerical_checks.json']),
        ('eq20_check.py', ['--out', 'eq20_check.json']),
        ('illustrations.py', ['--outdir', '.'])]
FAULTS = [('exact_checks.py', 'E3_envelope_ode'), ('numerical_checks.py', 'N5a_endemic_algebraic')]


def tree(p):
    return {q.relative_to(p).as_posix(): hashlib.sha256(q.read_bytes()).hexdigest()
            for q in p.rglob('*') if q.is_file() and not q.relative_to(p).as_posix().startswith('.git/')}


def cmp(a, b, path, diffs, stats):
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            diffs.append(f'{path}: keys differ {sorted(set(a) ^ set(b))[:5]}')
        for k in set(a) & set(b):
            if k in VOLATILE:
                stats['volatile'] += 1
                continue
            cmp(a[k], b[k], f'{path}.{k}', diffs, stats)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            diffs.append(f'{path}: length {len(a)} vs {len(b)}')
            return
        for i, (x, y) in enumerate(zip(a, b)):
            cmp(x, y, f'{path}[{i}]', diffs, stats)
    elif isinstance(a, str) and isinstance(b, str) and a != b and a[:1] in '{[' and b[:1] in '{[':
        try:  # JSON-encoded detail strings are compared value by value
            pa, pb = json.loads(a), json.loads(b)
        except ValueError:
            stats['exact'] += 1
            diffs.append(f'{path}: {a[:60]!r} vs {b[:60]!r}')
        else:
            cmp(pa, pb, path, diffs, stats)
    elif isinstance(a, bool) or isinstance(b, bool) or not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
        stats['exact'] += 1
        if a != b:
            diffs.append(f'{path}: {str(a)[:60]!r} vs {str(b)[:60]!r}')
    else:
        stats['numeric'] += 1
        if a == b:
            return
        if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
            return
        if not math.isclose(a, b, rel_tol=1e-7, abs_tol=1e-10):
            diffs.append(f'{path}: {a!r} vs {b!r}')
        else:
            stats['last_place'] += 1


before = tree(root)
with tempfile.TemporaryDirectory(prefix='nimfa-replay-') as td:
    w = Path(td) / 'pkg'
    shutil.copytree(root, w, ignore=shutil.ignore_patterns('.git'))
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg', MPLCONFIGDIR=str(Path(td) / 'mpl'),
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    env.pop('NIMFA_INJECT_FAULT', None)
    plan = FULL[:2] if quick else FULL
    for script, args in plan:
        r = subprocess.run([sys.executable, '-I', '-O', '-B', script] + args, cwd=w, capture_output=True,
                           text=True, env=env, timeout=3600)
        if r.returncode:
            sys.exit(f'FAIL {script} exit {r.returncode}: {r.stderr[-2000:]}')
        print(f'ok {script} {" ".join(args)} (-I -O -B)')
    for script, check in FAULTS:
        r = subprocess.run([sys.executable, '-I', '-O', '-B', script, '--out', str(Path(td) / f'fault-{check}.json')],
                           cwd=w, capture_output=True, text=True, env=dict(env, NIMFA_INJECT_FAULT=check), timeout=3600)
        if r.returncode == 0:
            sys.exit(f'FAIL fault injection {check} in {script} was not detected')
        print(f'ok fault injection {check} in {script} detected (exit {r.returncode}, -I -O -B)')
    after = tree(w)
    extra = sorted(set(after) - set(before))
    if extra:
        sys.exit(f'FAIL replay wrote unexpected files: {extra}')
    changed = sorted(k for k in before if after.get(k) != before[k])
    diffs, stats, compared = [], dict(exact=0, numeric=0, last_place=0, volatile=0), []
    outputs = ['exact_checks.json', 'numerical_checks.json'] + ([] if quick else ['eq20_check.json', 'illustrations.json'])
    for k in outputs:
        cmp(json.loads((root / k).read_text()), json.loads((w / k).read_text()), k, diffs, stats)
        compared.append(k)
    for k in changed:
        if k not in outputs and not k.endswith(('.png', '.pdf')):
            diffs.append(f'{k}: unexpected output changed')
    print(f'outputs compared: {compared}; byte-changed files: {changed}')
    print(f'leaves exact={stats["exact"]} numeric={stats["numeric"]} (within tolerance but not identical: '
          f'{stats["last_place"]}); package-version fields skipped {stats["volatile"]}')
    if diffs:
        print('\n'.join(diffs[:40]))
        sys.exit(f'FAIL {len(diffs)} differences beyond tolerance')
if tree(root) != before:
    sys.exit('FAIL root changed')
mode = 'quick (exact and numerical checkers)' if quick else 'full (all four checkers)'
print(f'PASS clean replay, {mode}: every checker passed under -I -O -B; regenerated outputs equal the retained ones '
      'within rel 1e-7 (package-version fields ignored; figures regenerated, not byte-compared); fault injection '
      'detected in both fast checkers; no unexpected files; root unchanged')
