# PyWATER Modernization Execplan

Date started: 2026-07-10

Scope: track the ongoing modernization of `pywater.py` (single-file PyMOL
plugin) — what is done, what remains, and important findings so future work
does not re-litigate settled decisions.

Working constraints:

- Do not work directly on `master`; use `hp_*` feature branches.
- Keep the installed plugin copy in sync after each change:
  `cp pywater.py .venv/lib/python3.12/site-packages/pmg_tk/startup/pywater.py`
  (PyMOL loads that copy and must be restarted to pick up changes).
- Preserve unrelated local changes (`.gitignore`, `CLAUDE.md`).

## Done

### Environment / tooling (earlier sessions)
- Open-source PyMOL 3.2.0a0 installed via `uv` into `.venv` on macOS arm64,
  incl. the non-self-contained wheel's dylib/@rpath fix (see memory
  `pymol-uv-install`). GUI PyMOL runs.

### Merged to `master`
- **PR #15** — migrate to modern RCSB APIs (Data/Search/GraphQL), batch
  requests, threaded downloads. Retired the dead `customReport`/XML endpoints.
- **PR #16** (`hp_qt_gui`) — port GUI Tkinter -> Qt (`pymol.Qt`); run the
  search off the main thread so the GUI stays responsive; graceful
  structure-cap guard (`DEFAULT_MAX_STRUCTURES=200`, resolution-sorted with the
  query forced first) to stay under the 50000-water clustering limit; GUI
  structure-selection summary.
- **PR #17 / #18** — README + troubleshooting docs.
- **PR #19** (`hp_modernization`) — drop Python 2 compatibility
  (`sys.version_info`/`urllib2`/`xrange`/`xml.minidom` removed); centralize
  RCSB HTTP helpers (`_get_json`/`_post_json`/`_graphql`, module URL/timeout
  constants); `try/finally` temp-dir cleanup; `ImportError`/`ValueError` in
  place of `sys.exit`/`assert`; add `pyproject.toml` + `requirements.txt`;
  document shell CLI usage; replace `multiprocessing.dummy.Pool` with
  `concurrent.futures.ThreadPoolExecutor`.

### GitHub issue triage
- All 6 open issues (#8, #10, #11, #12, #13, #14) closed as resolved/superseded
  by the RCSB + Qt modernization. See `ISSUE_TRIAGE_EXECPLAN.md`. 0 open issues.

### This branch (`hp_modernization_v2`) — log-noise cleanup
- Per-structure refinement logs moved INFO -> DEBUG in `okMobility` /
  `okBfactor` (mobility/B-factor counts, included/excluded/considered lines).
- Added one INFO **summary** at the caller in `makePDBwithConservedWaters`:
  `<method> refinement filtering: N of M structures retained (K discarded ...)`.
- Downloads: per-structure "Retrieving structure" and per-failure logs moved
  INFO/ERROR -> DEBUG; added a single WARNING summary counting skipped
  (e.g. 404) structures. Normal runs no longer look alarming.
- Status: code done, syntax-checked, installed copy synced. Not yet committed.

### This branch (`hp_modernization_v2`) — file-I/O context managers (#5)
- All manual `open(...)` calls now use `with` (or a list-and-write-once pattern).
  Converted: `okMobility` / `okBfactor` (two reads + rewrite each), Protein
  `calculate_water_coordinates`, and the conserved-waters PDB writer (nested
  `with` over input + output).
- The deeply-nested `clusterPresence.txt` writer (~66 lines) was NOT mass
  re-indented (too risky); instead rows are collected into
  `clusterPresenceLines` and written once via a small `with` block. Same
  output, flat indentation, guaranteed close.
- Only `urllib.urlopen` remains outside a plain `with open`, and it is already
  context-managed. Syntax-checked, installed copy synced. Not yet committed.

### This branch (`hp_modernization_v2`) — automated tests (#1)
- Added `tests/` unit suite (stdlib `unittest`, no pytest) run under a stubbed
  PyMOL: `tests/pymol_stubs.py` installs minimal `pymol` / `pymol.Qt` /
  `pymol.cmd` stubs and imports `pywater`, so tests need only NumPy/SciPy.
- 16 tests over: `_decode_response`, `filterbyResolution` (best-first sort +
  cutoff), `_prioritize_and_cap` (query-first, cap keeps query, no-mutate),
  `toPyWATER` arg coercion/plumbing, `_LogCollector.summary()/found_waters()`.
- Refactor: extracted the inline query-first/cap logic from
  `FindConservedWaters` into pure helper `_prioritize_and_cap(...)` so it is
  unit-testable (behavior preserved; logging unchanged).
- Run: `.venv/bin/python -m unittest discover -s tests`. All green.
- CLAUDE.md updated (removed "no test suite"; documented how to run).

## Remaining

All "worth doing next" items are done. Remaining are lower-priority / later.

### Probably later
2. Split `pywater.py` into modules (GUI / network / clustering / output).
   High-value maintainability but high-risk refactor; wait until behavior is
   stable and tests exist.
3. Real package/plugin install flow beyond current `pyproject.toml`.
4. Type hints / formatting tooling (noisy, low urgency).
5. `pathlib` adoption (cosmetic).

Note: for a large deeply-nested I/O block, prefer collect-into-a-list +
write-once through a `with` rather than re-indenting the whole block into a
context manager — lower risk, same guarantee (used for `clusterPresence.txt`).

## Important findings

- **Do NOT preemptively move `cmd.*` calls back to the main thread.**
  Codex flagged the GUI worker calling PyMOL `cmd.*` off-thread as a crash
  risk. Investigation says otherwise:
  - `_worker` deliberately runs `FindConservedWaters` off the main thread to
    fix the beachball/GUI freeze; results marshal back via the `_finished`
    Qt signal.
  - PyMOL's `cmd.*` layer serializes calls through its own internal API lock,
    so off-thread `cmd.load/super/create/save` is a supported pattern.
  - The macOS SIGABRT seen earlier came from a **test harness** creating a
    Qt/Cocoa dialog on a background thread — Cocoa window creation must be on
    the main thread — NOT from `cmd.*` in the worker.
  - The user tested this threaded build: "pretty fast and not hanging Mac
    window."
  - Fully marshaling every `cmd.*` back to the main thread would likely
    reintroduce the freeze (all `displayInPyMOL` rendering back on the event
    loop). Treat this as fix-when-a-real-crash-is-observed, not preemptive.
- **RCSB downloads 404 routinely** for broad clusters; skipping them is normal,
  not an error (drove the log-noise cleanup above).
- **Water identity string slicing is load-bearing** (`pdb_chain_atomNumber`,
  `[:6]`/`[7:]`) — see CLAUDE.md; don't disturb during refactors.
- **50000-water clustering cap** is a hard memory limit (scipy O(n^2) pairwise);
  the structure-cap guard is the mitigation.

## Log

- 2026-07-10: Merged PR #19. Closed all 6 issues. Started
  `hp_modernization_v2` with log-noise cleanup (#2/#3 from Codex's remaining
  list). This file created.
- 2026-07-10: On `hp_modernization_v2`, converted all file I/O to context
  managers (#5).
- 2026-07-10: Added `tests/` unit suite (#1) + extracted `_prioritize_and_cap`
  helper. 16 tests green. Log-noise + file-I/O + tests bundled for one PR.
