# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

PyWATER is a single-file PyMOL plugin (`pywater.py`) that finds conserved water molecules in X-ray protein structures. Given a query PDB id + chain, it fetches homologous structures from RCSB PDB, superimposes them, spatially clusters their crystallographic waters, and reports waters that recur across most structures. The whole program is `pywater.py`; `pyproject.toml`/`requirements.txt` document its NumPy/SciPy deps, and `tests/` holds a small unit suite for the non-GUI logic (run under a stubbed PyMOL — see below).

## Running it

PyWATER only runs inside PyMOL (it imports `pymol.cmd` at module load and cannot execute standalone). Requires NumPy and SciPy in PyMOL's Python environment.

- **As a plugin**: install `pywater.py` via PyMOL's Plugin Manager (Plugins → Manage Plugins → Install), restart PyMOL. Adds a "PyWATER" GUI entry under the Plugin menu.
- **From the PyMOL command line**: `pywater 4lyw, A, 95` (registered via `cmd.extend('pywater', toPyWATER)`). Note `cmd.extend` registers a *command keyword* usable as `pywater ...` in the PyMOL command language — it is **not** exposed as a `cmd.pywater(...)` attribute (that older doc claim is inaccurate; `hasattr(cmd,'pywater')` is False).
- **Local / unpublished structures**: `pywater_local /path/to/folder, A, my_ref.pdb` (registered via `cmd.extend('pywater_local', toPyWATERLocal)`) finds conserved waters across local `.pdb` files instead of RCSB homologs. See `LOCAL_FILES_PLAN.md`.

Argument order (all after PDB id and chain are optional): `PDB id, Chain id, seq_id, resolution, refinement, user_def_list, clustering_method, inconsistency_coefficient, prob`. Defaults: seq_id `95`, resolution `2.0`, refinement `Mobility`, clustering `complete`, inconsistency `2.4`, prob `0.7`.

To exercise the full pipeline, launch PyMOL, run a query, and inspect the results. The non-GUI logic (RCSB response parsing, resolution sorting, structure cap / query-first selection, `toPyWATER` argument plumbing, `_LogCollector` summaries) has unit tests: `.venv/bin/python -m unittest discover -s tests` from the repo root. They import `pywater` under a stubbed PyMOL (`tests/pymol_stubs.py`), so they need only NumPy/SciPy, not a running PyMOL. `Benchmark/` and `Gallery/` hold reference outputs (per-structure result folders and PNG renderings), not test fixtures.

Output goes to `~/PyWATER_outdir/` (created at import time): `pywater.log`, and a per-run subfolder `PDBid_CHAINid/` containing the query PDB with conserved waters, a `*_clusterPresence.txt` table, and (optionally) superimposed intermediate PDBs.

## Architecture / control flow

The pipeline lives in `FindConservedWaters()` (the RCSB entry point that the GUI and CLI call through) → `makePDBwithConservedWaters()` (the algorithmic core). Reading these two functions top-to-bottom is the fastest way to understand the program. `FindConservedWatersLocal()` is a parallel entry point for local files: it skips every RCSB step (reachability, `isXray`, `chainPresent`, `fetchpdbChainsList`, `filterbyResolution`, the cap, download), assigns each local file a 4-char id via `collectLocalStructures()` (the real PDB id when the filename starts with one, e.g. `4lyw_a.pdb`→`4lyw`, else synthetic `lc00`…), copies files into the temp dir, guards on chain presence (`_chainHasCA`), then converges on the same `makePDBwithConservedWaters()`. The 4-char-id scheme exists because the water-identity slicing (`[:6]`/`[7:]`) and the `cwm_????_?.pdb` glob require a 4-char id + 1-char chain.

1. **Validation & metadata lookup** (`FindConservedWaters`): validate PDB/chain id format, confirm X-ray via `isXray()`, confirm chain exists via `chainPresent()`, range-check numeric params. All error reporting is dual: `logger.error(...)` plus a Tkinter `tkMessageBox` popup.
2. **Building the structure set**: unless a user-defined list is given, `fetchpdbChainsList()` queries RCSB's sequence-cluster REST API for homologs at the chosen identity cutoff, then `filterbyResolution()` drops low-resolution structures. The query structure is always forced into the list.
3. **Download & superpose** (`makePDBwithConservedWaters`): each PDB is fetched to a temp dir and loaded into PyMOL as `cwm_*` objects; `cmd.super(...////CA...)` aligns each onto the query. Water-only objects are saved per structure.
4. **Refinement filtering**: `okMobility()` / `okBfactor()` remove poorly-refined waters per structure; a structure is discarded entirely if >50% of its waters fail. Selected via the `refinement` param (`Mobility`, `Normalized B-factor`, or `No refinement`). These read B-factor/occupancy from fixed PDB column offsets (`line[54:60]`, `line[60:66]`).
5. **Clustering**: all retained water coordinates are pooled and clustered with `scipy.cluster.hierarchy.fclusterdata` (euclidean, `criterion='distance'`, threshold = inconsistency coefficient). Hard cap of 50000 waters.
6. **Conservation scoring**: for each cluster, duplicate waters from the same PDB are pruned, then degree of conservation = (unique PDBs in cluster) / (total structures). Clusters at/above the `prob` cutoff and containing the query structure become conserved waters.
7. **Output**: writes the query PDB augmented with conserved waters, and `displayInPyMOL()` renders the PyMOL session — surface, H-bond distance objects (protein/ligand/water donors & acceptors within 4.0 Å), and waters colored by conservation via `b`-factor spectrum (`red_blue`).

## Conventions & gotchas

- **Water identity strings**: waters are tracked as strings like `1abc_A_1234` (`pdb_chain_atomNumber`). Slicing `[:6]` extracts `pdb_ch` and `[7:]` the atom number — this positional slicing is load-bearing throughout the clustering/scoring loop. `Protein.__repr__` produces the `pdbid_chain` prefix.
- **PyMOL object namespace**: all working objects are prefixed `cwm_` and bulk-deleted with `cmd.delete('cwm_*')` between phases. Don't introduce unprefixed temp objects into that flow.
- **Python 2/3 compat**: the top of the file branches on `sys.version_info` for `urllib`, Tkinter, and `xrange`. Preserve this dual-support pattern when editing imports.
- **RCSB API is stale**: the REST endpoints (`pdb.org/pdb/rest/...`, `rcsb.org/pdb/rest/...`) are legacy URLs that RCSB has since retired/changed. Network-dependent functions (`isXray`, `chainPresent`, `fetchpdbChainsList`, `filterbyResolution`, PDB file download) may fail against current RCSB services — expect this when debugging "no structures found" issues rather than assuming a logic bug.
- **The `~/PyWATER_outdir` and logger are created at import time**, not inside a function — importing the module has filesystem side effects.
