# Replay receipt

- **Pre-freeze replay (1 November 2026, release preparation):** `python -I -B tools/replay_compare.py <package>` on a copy of the publication working tree taken just before DOI insertion, whose checkers and outputs are byte-identical to the released ones (CPython 3.11.14, pinned packages). All four checkers ran with `-I -O -B` in a temporary copy and passed; `exact_checks.json`, `numerical_checks.json`, `eq20_check.json` and `illustrations.json` were regenerated and compared with the retained ones; fault injection (`NIMFA_INJECT_FAULT=E3_envelope_ode`, `N5a_endemic_algebraic`) made both fast checkers fail as required; no unexpected files.
- **Quick replay:** `tools/replay_compare.py . --quick` (exact and numerical checkers and fault injection, under half a minute) passed with byte-identical outputs.
- **Frozen-archive replay:** the full comparer was run from a fresh extraction of the frozen release ZIP as part of the release gate plan (retained privately with the publication receipts).
- **Earlier replays:** before release preparation, the editor replayed the checkers from a fresh copy (outputs identical), and an independent reviewer re-ran both fast checkers on the revised version, also with identical outputs.
- **CI:** `.github/workflows/replay.yml` reruns the quick comparer on Linux (Python 3.11, pinned packages) as a blocking job on every push, and the full comparer as a non-blocking job.
- Environment: [ENVIRONMENT.txt](ENVIRONMENT.txt).
