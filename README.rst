=======
PyWATER
=======

PyWATER is a PyMOL plugin for identifying conserved water molecules in
X-ray protein structures. Given a query PDB ID and chain, PyWATER fetches
homologous structures from RCSB PDB, superimposes them in PyMOL, clusters
their crystallographic waters, and reports waters that recur across the
structure set.

PyWATER source code and documentation are available at:

https://github.com/hiteshpatel379/PyWATER

PyMOL is available at:

https://www.pymol.org/

Copyright 2013 Hitesh Patel and B. Gruening


Installation
============

PyWATER runs inside PyMOL as a plugin. Install the latest ``pywater.py`` file
into the Python environment used by PyMOL.

Requirements
------------

- PyMOL with Python 3 and plugin support
- NumPy installed in PyMOL's Python environment
- SciPy installed in PyMOL's Python environment
- Internet access to query RCSB PDB and download PDB structures

Dependencies are listed in ``requirements.txt`` and ``pyproject.toml`` for
modern Python tooling. PyMOL itself is still provided by your PyMOL
installation, not by PyPI.

Plugin Installation
-------------------

1. Start PyMOL.
2. Open ``Plugins -> Plugin Manager -> Install New Plugin``.
3. Install from the local file ``pywater.py``.
4. Restart PyMOL.
5. Open ``Plugins -> PyWATER``.

PyWATER writes results to ``~/PyWATER_outdir``. Avoid launching PyMOL with
administrator privileges unless necessary, because output files may then be
owned by the administrator account.


Usage
=====

PyWATER can be run from the graphical plugin dialog, the PyMOL command line,
or PyMOL's Python API.

Graphical Interface
-------------------

Open the plugin from ``Plugins -> PyWATER``. Enter a PDB ID and chain ID, then
adjust optional parameters as needed. The GUI is implemented with PyMOL's Qt
plugin API and runs the water search off the main GUI thread so the PyMOL
window remains responsive during longer searches.

PyMOL Command Line
------------------

The shortest form is:

.. code-block:: text

   PyMOL> pywater 4lyw, A

With optional parameters:

.. code-block:: text

   PyMOL> pywater PDB_ID, CHAIN_ID, SEQ_ID, RESOLUTION, REFINEMENT, USER_LIST, LINKAGE, INCONSISTENCY, CONSERVATION, MAX_STRUCTURES

Example:

.. code-block:: text

   PyMOL> pywater 4lyw, A, 95, 2.0, Mobility, , complete, 2.4, 0.7, 200

Shell Command Line
------------------

If you launch PyMOL from a terminal, pass PyWATER through PyMOL's ``-d``
command option. Positional shell arguments such as ``pymol 4lyw A`` are treated
as filenames by PyMOL, not as PyWATER inputs, and will fail with an unsupported
file type error.

Headless/quiet mode:

.. code-block:: bash

   .venv/bin/pymol -cq -d "run .venv/lib/python3.12/site-packages/pmg_tk/startup/pywater.py; pywater 4lyw, A"

Visible PyMOL window:

.. code-block:: bash

   .venv/bin/pymol -q -d "run .venv/lib/python3.12/site-packages/pmg_tk/startup/pywater.py; pywater 4lyw, A"

With explicit optional parameters:

.. code-block:: bash

   .venv/bin/pymol -q -d "run .venv/lib/python3.12/site-packages/pmg_tk/startup/pywater.py; pywater 4lyw, A, 95, 2.0, Mobility, , complete, 2.4, 0.7, 200"

The ``-c`` flag runs PyMOL without the GUI, ``-q`` quiets startup output, and
``-d`` executes PyMOL commands after startup.

Local / Unpublished Structures
------------------------------

If you have your own structures that are not in the PDB (e.g. unpublished
crystal forms on your workstation), PyWATER can find conserved waters across
those local files instead of RCSB homologs.

- **GUI**: tick **"Use local files (skip RCSB)"**, set **Local files folder**
  to a folder containing your ``.pdb`` files (use *Browse...*), set **Chain id**
  to the chain to analyse (applied to every file), and enter the **Reference
  file** whose conserved waters are reported. Sequence identity, resolution and
  maximum-structures are ignored in this mode.

- **Command line** (``pywater_local FOLDER, CHAIN, REFERENCE`` with optional
  ``REFINEMENT, LINKAGE, INCONSISTENCY, CONSERVATION, SAVE_SUPERIMPOSED``):

  .. code-block:: text

     PyMOL> pywater_local /path/to/my_structures, A, my_apo.pdb

  .. code-block:: bash

     .venv/bin/pymol -cq -d "run .venv/lib/python3.12/site-packages/pmg_tk/startup/pywater.py; pywater_local /path/to/my_structures, A, my_apo.pdb"

Notes:

- At least two ``.pdb`` files are required; all files must share the requested
  chain id (single character). Files lacking that chain are skipped with a
  warning; if the reference lacks it, the run stops with an error.
- Each file gets a 4-char id: the real PDB id when the filename starts with a
  valid one (e.g. ``4lyw_a.pdb`` → ``4lyw``), otherwise a synthetic
  ``lc00``/``lc01``/... id. ``local_files_map.txt`` records the id → filename
  mapping. The output folder ``~/PyWATER_outdir/<query_id>_<chain>/`` contains
  that map and the result as both
  ``cwm_<query_id>_<chain>_withConservedWaters.pdb`` and a copy named after your
  reference file (e.g. ``my_apo_withConservedWaters.pdb``).
- Current limits: ``.pdb`` format only, single-character chain ids, up to 100
  files.

PyMOL Python API
----------------

.. code-block:: python

   from pymol import cmd
   cmd.pywater("4lyw", "A")


Input Parameters
================

+----------------------+----------------+-------------------------------------------------------+
| Parameter            | Default        | Description                                           |
+======================+================+=======================================================+
| PDB ID               | required       | PDB identifier for the query protein.                 |
+----------------------+----------------+-------------------------------------------------------+
| Chain ID             | required       | Chain identifier for the query protein.               |
+----------------------+----------------+-------------------------------------------------------+
| Sequence identity    | 95             | RCSB sequence-cluster identity cutoff. Allowed        |
| cutoff               |                | values: 30, 40, 50, 70, 90, 95, 100.                 |
+----------------------+----------------+-------------------------------------------------------+
| Resolution cutoff    | 2.0            | Keep structures with resolution at or below this      |
|                      |                | value. Maximum accepted value is 3.0.                 |
+----------------------+----------------+-------------------------------------------------------+
| Refinement method    | Mobility       | Water filtering method: ``Mobility``,                 |
|                      |                | ``Normalized B-factor``, or ``No refinement``.        |
+----------------------+----------------+-------------------------------------------------------+
| User-defined list    | disabled       | Comma-separated list such as ``1abc_A,2def_B``.       |
|                      |                | When provided, sequence identity and resolution       |
|                      |                | filtering are skipped.                                |
+----------------------+----------------+-------------------------------------------------------+
| Linkage method       | complete       | Hierarchical clustering linkage: ``single``,          |
|                      |                | ``complete``, or ``average``.                         |
+----------------------+----------------+-------------------------------------------------------+
| Inconsistency        | 2.4            | Distance threshold used for water clustering.         |
| coefficient          |                | Maximum accepted value is 2.8.                        |
+----------------------+----------------+-------------------------------------------------------+
| Degree of            | 0.7            | Conservation cutoff. Conserved waters must appear     |
| conservation         |                | in at least this fraction of selected structures.     |
|                      |                | Allowed range: 0.4 to 1.0.                            |
+----------------------+----------------+-------------------------------------------------------+
| Maximum structures   | 200            | Maximum number of sequence-cluster structures used    |
|                      |                | for clustering. The query is kept first, and the      |
|                      |                | remaining structures are prioritized by best          |
|                      |                | resolution. Ignored for user-defined lists.           |
+----------------------+----------------+-------------------------------------------------------+


Structure Selection
===================

For sequence-cluster searches, PyWATER now uses the current RCSB Data and
Search APIs. The structure list is filtered by resolution, sorted from best
to worst resolution, and capped by ``Maximum structures`` to keep clustering
within memory limits for very large protein families.

During a run, the GUI status and ``pywater.log`` include a concise summary:

.. code-block:: text

   Structure selection summary: 438 candidate chains, 217 passed the 2.00 A resolution cutoff, 200 selected for clustering.


Results
=======

PyWATER writes output to ``~/PyWATER_outdir``. Each run creates a subfolder
named for the query structure, for example ``4lyw_A``.

Common output files include:

- ``pywater.log``: input parameters, structure-selection summary, warnings,
  and final status messages.
- ``cwm_<PDB>_<CHAIN>.pdb``: the query structure after PyWATER processing.
- ``cwm_<PDB>_<CHAIN>_withConservedWaters.pdb``: the query structure with
  conserved waters added.
- ``<PDB>_<CHAIN>_clusterPresence.txt``: tabular conservation information for
  water clusters across selected structures.
- Optional superimposed intermediate PDB files when ``Save superimposed pdb
  files`` is enabled.

In PyMOL, conserved waters are displayed as spheres and colored by degree of
conservation. Hydrogen-bond distance objects are created for interactions
between conserved waters and nearby protein, ligand, or water atoms.


Troubleshooting
===============

See ``docs/troubleshooting.md`` for common installation and runtime issues.

Quick checks:

- Install NumPy and SciPy into the same Python environment used by PyMOL.
- Use the latest ``pywater.py`` from this repository.
- Confirm internet access from the machine running PyMOL.
- Check ``~/PyWATER_outdir/pywater.log`` for the full run log.
- If a run fails, include the PyMOL version, Python version, input parameters,
  and the traceback/log excerpt when reporting an issue.


History
=======

- v2.0.0: Major modernization release after the original public release.
  Migrated RCSB access to current Data/Search APIs, modernized the PyMOL Qt
  GUI, moved long searches off the GUI thread, added best-resolution structure
  prioritization with a configurable structure cap, documented shell CLI usage,
  added local-files mode for unpublished/private structures, and refreshed
  package metadata for modern Python tooling. Removed the outdated PDF tutorial;
  this README is now the maintained usage reference.
- v1.0: Initial public release.


Citation
========

Please cite the following article if you use PyWATER:

Patel,H. *et al.* (2014) `PyWATER: a PyMOL plug-in to find conserved water
molecules in proteins by clustering
<https://doi.org/10.1093/bioinformatics/btu424>`_. *Bioinformatics*, **30**,
2978-2980.


Licence (MIT)
=============

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
