# PyWATER Troubleshooting

This page collects common installation and runtime issues for PyWATER.

## Plugin Does Not Initialize

Symptoms:

- PyMOL reports `Unable to initialize plugin 'PyWATER'`.
- The `pywater` command is not recognized after installation.
- The plugin appears to install but does not open from the Plugins menu.

Checks:

- Restart PyMOL after installing `pywater.py`.
- Make sure PyWATER was installed into the Python environment used by PyMOL,
  not a separate system Python or virtual environment.
- Install NumPy and SciPy in PyMOL's Python environment.
- Use PyMOL with Python 3; current PyWATER no longer supports Python 2.
- Use a recent PyMOL build with Qt plugin support.
- Use the latest `pywater.py` from this repository.

PyWATER imports NumPy, SciPy, PyMOL's Qt bindings, and `pymol.cmd` when the
plugin is loaded. If one of those imports fails, PyMOL may report a generic
plugin initialization error.

`requirements.txt` and `pyproject.toml` document the Python dependencies for
modern tooling, but PyMOL itself should come from your PyMOL installation.

## RCSB HTTP 404 or XML Parse Errors

Older PyWATER versions used retired RCSB `customReport` XML endpoints. Those
failures often looked like:

- `HTTPError: HTTP Error 404: Not Found`
- `xml.parsers.expat.ExpatError`
- tracebacks inside `isXray()`, `chainPresent()`, or `parseString(...)`

Current PyWATER uses the modern RCSB Data/Search APIs instead:

- `https://data.rcsb.org/rest/v1/core/entry/...`
- `https://data.rcsb.org/graphql`
- `https://search.rcsb.org/rcsbsearch/v2/query`

If you see one of the old `customReport` URLs in a traceback, update your
installed plugin file to the latest `pywater.py`.

## Network or Download Failures

PyWATER needs internet access to:

- confirm the query structure metadata;
- find sequence-cluster neighbors through RCSB;
- fetch structure resolution metadata;
- download PDB files from RCSB.

If structure downloads fail, PyWATER logs the failed PDB ID and skips that
structure. If too few structures remain, the run cannot continue.

Checks:

- Confirm the machine running PyMOL can reach `https://data.rcsb.org`,
  `https://search.rcsb.org`, and `https://files.rcsb.org`.
- Check whether a proxy, firewall, VPN, or institutional network filter is
  blocking PyMOL's Python process.
- Try a small query such as `pywater 4lyw, A` from the PyMOL command line.
- Inspect `~/PyWATER_outdir/pywater.log` for the exact failing URL or PDB ID.

## Very Large Protein Families

Large sequence clusters can return hundreds of structures and tens of
thousands of waters. PyWATER now limits broad sequence-cluster runs with the
`Maximum structures` setting.

Default:

- `Maximum structures`: `200`

Behavior:

- structures are filtered by resolution;
- structures are sorted by best resolution first;
- the query structure is kept at the front;
- the selected list is capped before clustering;
- user-defined structure lists are not capped.

The GUI and log include a line like:

```text
Structure selection summary: 438 candidate chains, 217 passed the 2.00 A resolution cutoff, 200 selected for clustering.
```

If a run still reports too many waters to cluster, lower `Maximum structures`
or provide a smaller user-defined list.

## Where Results Are Written

PyWATER creates:

```text
~/PyWATER_outdir/
```

The main log is:

```text
~/PyWATER_outdir/pywater.log
```

Each run creates a subfolder such as:

```text
~/PyWATER_outdir/4lyw_A/
```

That folder contains the processed query PDB, conserved-water PDB, cluster
presence table, and optionally the superimposed intermediate PDB files.

## Reporting a New Issue

When reporting a new problem, include:

- PyWATER version or commit, if known;
- PyMOL version;
- Python version used by PyMOL;
- operating system;
- exact PyWATER inputs;
- whether the GUI or command line was used;
- the relevant traceback;
- the matching lines from `~/PyWATER_outdir/pywater.log`.
