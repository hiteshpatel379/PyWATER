# Implementation Plan — Conserved Waters from Local Files

Branch: `hp_local_files`
Date: 2026-07-10
Status: **Implemented** — `collectLocalStructures`, `_chainHasCA`,
`_writeLocalMap`, `FindConservedWatersLocal`, `toPyWATERLocal`
(`pywater_local` command), and the GUI "Use local files" mode are in
`pywater.py`. 27 unit tests pass; verified end-to-end headless on 3 real PDBs
(40 conserved waters, mapping + friendly-named output written; bad-chain and
missing-reference error paths confirmed).

## Goal

Let a crystallographer find conserved waters across a set of **local,
unpublished structures** (files on their workstation) instead of RCSB
homologs. Same clustering/scoring/visualization; only the *source of
structures* changes.

## Chosen UX (decided)

- **Input:** a **folder path**. Every `*.pdb` in it is a comparison structure.
- **Chain:** a **single chain id** (existing `Chain id` field) applied to all
  files — matches the common case (same construct, chain `A`).
- **Reference/query:** the user **types the reference filename** (option (a)).
  Conserved waters are annotated onto and written out for this structure.

Mental model:

```
Local files folder: /Users/me/my_structures/   (all *.pdb compared)
Chain:              A                           (existing field, all files)
Reference file:     my_apo.pdb                  (query; waters reported on it)
```

## Key constraint (load-bearing)

Water identity strings are `pdbid_chain_atomNumber`, sliced positionally as
`[:6]` (4-char id + `_` + 1-char chain) and `[7:]` throughout clustering /
scoring. The save-superimposed glob `cwm_????_?.pdb` (pywater.py:538) and
`pdbIdFormat` (`^[a-z0-9]{4}$`) / `chainIdFormat` (`^[A-Z0-9]$`) bake in the
same 4-char-id + 1-char-chain shape.

=> Each local file is assigned a **synthetic 4-char id** (`lc00`, `lc01`, …),
with the reference forced to `lc00`. A mapping (synthetic id -> original
filename) is logged and written to the output folder so results stay legible.
Chain must be a single character (fine for normal crystallographic chains).

## Scope of v1

- `.pdb` files only (PDB fixed-column format — `okMobility`/`okBfactor` read
  `line[54:60]`/`line[60:66]`). mmCIF is a follow-up (convert-on-ingest via
  `cmd.load`+`cmd.save`).
- Local mode is **exclusive**: these files and nothing else (no RCSB mixing).
- Refinement filter, clustering method, inconsistency coeff., degree of
  conservation, save-superimposed: all still apply.
- No structure cap needed (user controls the file set), but warn if a very
  large folder risks the 50000-water clustering limit.

## Code changes

### 1. Local-structure ingestion (new, pure + testable)

```python
def collectLocalStructures(local_dir, reference_filename):
    """
    Return (structures, query_id) where structures is an ordered list of
    (synthetic_id, source_path) with the reference first (synthetic_id 'lc00').
    Raises ValueError on: missing/*empty* folder, < 2 .pdb files, or a
    reference filename not present in the folder.
    """
```

- `glob` `*.pdb` in `local_dir` (case-insensitive), sorted for determinism.
- Validate ≥ 2 files and that `reference_filename` (basename match) is present.
- Assign ids: reference -> `lc00`; the rest -> `lc01`, `lc02`, … (cap the
  count at 100 → ids stay 4 chars `lc00`..`lc99`; error if more, pointing to
  the follow-up).
- Return the mapping; caller copies each `source_path` -> `tmp_dir/<id>.pdb`.

A separate tiny helper `_writeLocalMap(outdir, query_chain_dir, mapping)`
writes `local_files_map.txt` and logs the table.

### 2. `FindConservedWaters` — add a local branch

Add trailing keyword params (keeps positional CLI/API order intact):

```python
def FindConservedWaters(..., max_structures=DEFAULT_MAX_STRUCTURES,
                        local_files_dir=None, local_reference=None):
```

When `local_files_dir` is set:
- **Skip** the RCSB reachability check, `isXray`, `chainPresent`,
  `fetchpdbChainsList`, `filterbyResolution`, the cap, and the RCSB download
  loop.
- Still run: `chainIdFormat(chain)` and the numeric range checks.
- `collectLocalStructures(...)` -> build `up` with synthetic ids; the query is
  `Protein(query_id, chain)`; `selectedStruturePDB = query_id`.
- Copy files into `tmp_dir` as `<id>.pdb`; guard each with a **chain-present /
  has-atoms** check (load, verify `chain` exists; skip with a WARNING if not —
  mirrors the failed-download skip, then re-check ≥ 2 remain).
- Converge on the **existing** `makePDBwithConservedWaters(up, tmp_dir, outdir,
  save_sup_files)` — no downstream changes.
- After the run, also copy `<query_id>_<chain>_withConservedWaters.pdb` to a
  friendly `<referencestem>_withConservedWaters.pdb` in the output folder.

Refactor note: extract the current "populate tmp_dir + assemble `up.proteins`
+ drop failures" tail into a small internal path so the RCSB and local branches
share the `makePDBwithConservedWaters` call and the `finally: rmtree` cleanup.

### 3. GUI (`ConservedWaters` dialog)

- New checkbox **"Use local files (skip RCSB)"**.
- New **folder** `QLineEdit` + **"Browse…"** button
  (`QtWidgets.QFileDialog.getExistingDirectory`).
- New **reference filename** `QLineEdit`.
- `varcheck`-style handler: when local mode is on, **disable** PDB id,
  sequence-identity, resolution, and the user-defined-list field; **enable**
  the folder + reference fields. (Chain, refinement, clustering, inconsistency,
  prob, save-superimposed stay active.)
- `run()`: in local mode, call `FindConservedWaters(local_files_dir=…,
  local_reference=…, chain=…, …)` and ignore PDB id / seq_id / resolution.
- Runs on the existing worker thread; `_LogCollector.summary()` gains a
  local-mode-friendly outcome line (e.g. "No structures with chain X").

### 4. CLI / Python API

Avoid overloading the comma-delimited `pywater` command. Add a dedicated one:

```python
def toPyWATERLocal(folder, chain, reference,
                   refinement='Mobility', clustering='complete',
                   inconsistency=2.4, prob=0.7, save_sup_files=False): ...
cmd.extend('pywater_local', toPyWATERLocal)
```

Usage: `cmd.pywater_local('/path/to/folder', 'A', 'my_apo.pdb')`.

## Edge cases & validation

- Folder missing / not a dir / no `.pdb` -> clear error, no crash.
- < 2 usable files (after chain filtering) -> error (need ≥ 2 to superimpose).
- Reference filename not in folder / lacks the chain -> error.
- Chain id not single alphanumeric -> existing `chainIdFormat` error.
- > 100 files -> error pointing to the (future) larger-id follow-up.
- Reference with no waters -> existing "no conserved waters" messaging applies.
- All popup-free (log/`_qt_info`), consistent with current error handling.

## Tests (`tests/test_pywater.py`, stubbed PyMOL)

Pure-Python, no `cmd`/network:
- `collectLocalStructures`: id assignment, reference-first (`lc00`), sorted
  order, ≥2 enforcement, reference-missing error, >100 error, non-`.pdb`
  ignored. Use `tempfile` dirs with empty `.pdb` files.
- `toPyWATERLocal` arg coercion/plumbing (spy on `FindConservedWaters`, assert
  `local_files_dir`/`local_reference`/chain forwarded, numerics coerced).
- `_LogCollector` new local outcome line.

## Docs

- README: new "Local / unpublished structures" section (GUI steps +
  `pywater_local` example + the `.pdb`-only / chain caveats).
- `docs/troubleshooting.md`: chain-not-found and <2-files notes.
- CLAUDE.md: mention the local-files path in the control-flow section.
- Update `MODERNIZATION_EXECPLAN.md` log.

## Out of scope (follow-ups)

- mmCIF input (convert-on-ingest).
- Per-file chain selection / multi-char chains.
- > 100 local files (needs a wider synthetic-id scheme).
- Mixing local files with RCSB homologs.

## Risk assessment

Low–medium. The heavy pipeline (`makePDBwithConservedWaters` onward) is reused
unchanged; new code is confined to structure gathering + a GUI branch. The main
risk is the `FindConservedWaters` refactor to share the converge point — kept
behavior-preserving for the RCSB path and covered by the existing 16 tests plus
a manual RCSB smoke run.
