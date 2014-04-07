======
PyWATER
======

PyWATER finds conserved water molecules in X-ray protein structure (pdb).


Important or conserved waters are the water molecules which are present in most or all available pdb structures when superimposed.

Copyright 2013 Hitesh Patel and B. Gruening


Installation
------------

PyWATER can be used as PyMOL plugin and you can install from inside PyMOL.  
- Run PyMOL ( on some systems you need to be an administrator.)  
- Install the PyWATER plugin in pymol by following the path: `Plugins -> Manage Plugins -> install`  
- Restart the PyMol  

**Enabling Plugins in MacPyMOL**   
In order to install plug-ins from inside MacPyMOL, you have to rename the "MacPyMOL.app" to "PyMOLX11Hybrid.app" in your Applications folder, which will start PyMOL in the X Window System (X11) mode (requires an installation of [XQuartz](http://xquartz.macosforge.org/landing/)).



Usage:
======


After installation as plugin. It can be run from PyMOL GUI by following the path: `Plugins -> PyWATER`

Enter the required parameters and click on 'Find Conserved Water Molecules'.

OR It can be run from commandline in pymol:
    `pywater [PDB id , Chain id [, sequence identity cutoff [, resolution cutoff [, refinement assessing method [, user defined proteins list [, linkage method [, inconsistency coefficient threshold [, degree of conservation]]]]]]]]`

OR you can use it directly from your python script.

```python
from pymol import cmd

cmd.pywater(PDB id , Chain id [, sequence identity cutoff [, resolution cutoff [, refinement assessing method [, user defined proteins list [, linkage method [, inconsistency coefficient threshold [, degree of conservation]]]]]]])

```


Parameter    | Default Value | Explanation
-------- | -------- | -----------
PDB id | -- | The PDB id of the protein for which you like to find conserved waters.
Chain id | -- | The chain identifier of the protein for which you like to find conserved waters in above mentioned PDB.
sequence identity cutoff | 95 % | All the protein structures, clustered by BlastClust, having sequence identity more than given cutoff will be superimposed to find the conserved water molecules in query protein chain.
resolution cutoff | 2.0 A | All the protein structures to be superimposed will be filtered first according to the structure resolution cutoff. Only structures with better resolution than given cutoff will be used further.
refinement assessing method | Mobility | Choose either 'Mobility' or 'Normalized B-factor' as criteria to assess the refinement quality of crystal structure. Program will filter out the water molecules with bad refinement quality.
user defined proteins list | -- | Give a custom list of protein structures to superimpose. Specifying this list will disable 'sequence identity' and 'resolution cutoff' parameters.
linkage method | complete | Linkage method for hierarchical clustering. Choose one from single, complete, average, weighted, median centroid, ward.
inconsistency coefficient threshold | 2.0 A | Any two clusters of water molecules will not be closer than given inconsistency coefficient threshold. Value ranges from 0 to 2.4.
degree of conservation | 0.7 | Water molecules will be considered CONSERVED if their probability of being conserved is above given cutoff. Value ranges from 0.4 to 1.


=======
History
=======

- v1.0: Initial public release


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

