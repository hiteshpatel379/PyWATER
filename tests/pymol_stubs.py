"""
Test bootstrap: install minimal stand-ins for PyMOL's modules so ``pywater``
can be imported and its non-GUI logic tested under a plain Python interpreter
(no running PyMOL). Importing this module registers the stubs and imports
``pywater``, which it re-exports as ``pywater``.

Only the attributes touched at import time are stubbed:
- ``pymol.cmd.extend`` (called at the bottom of pywater.py),
- ``pymol.Qt.QtWidgets.QDialog`` (base class of ConservedWaters),
- ``pymol.Qt.QtCore.Signal`` (class attribute of ConservedWaters).
GUI methods are not exercised by these tests.
"""

import os
import sys
import types

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)


def _install_pymol_stubs():
    if isinstance(sys.modules.get('pymol'), types.ModuleType) and \
            getattr(sys.modules['pymol'], '_pywater_stub', False):
        return

    pymol = types.ModuleType('pymol')
    pymol._pywater_stub = True

    cmd = types.ModuleType('pymol.cmd')
    cmd.extend = lambda name, fn: fn
    pymol.cmd = cmd

    qt = types.ModuleType('pymol.Qt')

    QtWidgets = types.ModuleType('pymol.Qt.QtWidgets')

    class _QDialog(object):
        def __init__(self, *args, **kwargs):
            pass

    QtWidgets.QDialog = _QDialog
    qt.QtWidgets = QtWidgets

    QtCore = types.ModuleType('pymol.Qt.QtCore')
    QtCore.Signal = lambda *args, **kwargs: object()

    class _Qt(object):
        WaitCursor = 0

    QtCore.Qt = _Qt
    qt.QtCore = QtCore

    sys.modules['pymol'] = pymol
    sys.modules['pymol.cmd'] = cmd
    sys.modules['pymol.Qt'] = qt
    sys.modules['pymol.Qt.QtWidgets'] = QtWidgets
    sys.modules['pymol.Qt.QtCore'] = QtCore


_install_pymol_stubs()

import pywater  # noqa: E402
