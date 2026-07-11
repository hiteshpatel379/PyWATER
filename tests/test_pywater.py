"""
Unit tests for PyWATER's non-GUI logic.

Run from the repo root with:

    .venv/bin/python -m unittest discover -s tests

These use a stubbed PyMOL (see ``tests/pymol_stubs.py``) so they need only
NumPy/SciPy, not a running PyMOL.
"""

import logging
import os
import tempfile
import unittest

from tests.pymol_stubs import pywater as pw


def _touch(path, content=''):
    with open(path, 'w') as f:
        f.write(content)


def _atom_line(atomname, chain):
    line = list(' ' * 80)
    line[0:6] = list('ATOM  ')
    name = atomname.ljust(4)
    line[12:16] = list(name)
    line[21] = chain
    return ''.join(line) + '\n'


class DecodeResponseTests(unittest.TestCase):
    class _Resp(object):
        def __init__(self, data):
            self._data = data

        def read(self):
            return self._data

    def test_decodes_bytes(self):
        self.assertEqual(pw._decode_response(self._Resp(b'{"a": 1}')), '{"a": 1}')

    def test_passes_through_str(self):
        self.assertEqual(pw._decode_response(self._Resp('already text')), 'already text')


class FilterByResolutionTests(unittest.TestCase):
    def setUp(self):
        self._orig_graphql = pw._graphql

    def tearDown(self):
        pw._graphql = self._orig_graphql

    def test_sorts_best_first_and_applies_cutoff(self):
        # 3def is above the cutoff; 4ghi has no resolution reported.
        canned = {
            "data": {
                "entries": [
                    {"rcsb_id": "1ABC", "rcsb_entry_info": {"resolution_combined": [2.0]}},
                    {"rcsb_id": "2XYZ", "rcsb_entry_info": {"resolution_combined": [1.0]}},
                    {"rcsb_id": "3DEF", "rcsb_entry_info": {"resolution_combined": [3.5]}},
                    {"rcsb_id": "4GHI", "rcsb_entry_info": {}},
                ]
            }
        }
        pw._graphql = lambda query, timeout=pw.HTTP_TIMEOUT: canned

        chains = ['1abc:A', '2xyz:A', '3def:A', '4ghi:A']
        result = pw.filterbyResolution(chains, 2.5)

        # Best (lowest) resolution first; over-cutoff and null-resolution dropped.
        self.assertEqual(result, ['2xyz:A', '1abc:A'])

    def test_empty_when_none_pass(self):
        pw._graphql = lambda query, timeout=pw.HTTP_TIMEOUT: {"data": {"entries": []}}
        self.assertEqual(pw.filterbyResolution(['1abc:A'], 2.0), [])


class PrioritizeAndCapTests(unittest.TestCase):
    def test_moves_query_to_front_preserving_order(self):
        chains = ['1abc:A', '2xyz:A', '3def:A']
        result, capped_from = pw._prioritize_and_cap(chains, '3def:A', 200)
        self.assertEqual(result, ['3def:A', '1abc:A', '2xyz:A'])
        self.assertIsNone(capped_from)

    def test_inserts_query_when_absent(self):
        chains = ['1abc:A', '2xyz:A']
        result, capped_from = pw._prioritize_and_cap(chains, '9zzz:A', 200)
        self.assertEqual(result[0], '9zzz:A')
        self.assertIn('1abc:A', result)
        self.assertIsNone(capped_from)

    def test_caps_keeping_query_first(self):
        chains = ['q:A'] + ['s%d:A' % i for i in range(10)]
        result, capped_from = pw._prioritize_and_cap(chains, 'q:A', 3)
        self.assertEqual(len(result), 3)
        self.assertEqual(result[0], 'q:A')
        self.assertEqual(capped_from, 11)

    def test_no_cap_when_disabled(self):
        chains = ['a:A', 'b:A', 'c:A']
        result, capped_from = pw._prioritize_and_cap(chains, 'a:A', 0)
        self.assertEqual(len(result), 3)
        self.assertIsNone(capped_from)

    def test_does_not_mutate_input(self):
        chains = ['1abc:A', '2xyz:A']
        pw._prioritize_and_cap(chains, '2xyz:A', 200)
        self.assertEqual(chains, ['1abc:A', '2xyz:A'])


class ToPyWATERArgTests(unittest.TestCase):
    def setUp(self):
        self._orig = pw.FindConservedWaters
        self.captured = {}

        def _spy(*args, **kwargs):
            self.captured['args'] = args
            self.captured['kwargs'] = kwargs

        pw.FindConservedWaters = _spy

    def tearDown(self):
        pw.FindConservedWaters = self._orig

    def test_defaults_and_type_coercion(self):
        pw.toPyWATER('4LYW', 'a')
        args = self.captured['args']
        kwargs = self.captured['kwargs']
        # pdb lowercased, chain uppercased
        self.assertEqual(args[0], '4lyw')
        self.assertEqual(args[1], 'A')
        # seq_id stays a string; numeric params coerced
        self.assertEqual(args[2], '95')
        self.assertEqual(args[3], 2.0)
        self.assertIsInstance(args[3], float)
        self.assertEqual(args[4], 'Mobility')
        self.assertEqual(args[8], 0.7)
        self.assertIsInstance(args[8], float)
        self.assertEqual(kwargs['max_structures'], pw.DEFAULT_MAX_STRUCTURES)
        self.assertIsInstance(kwargs['max_structures'], int)

    def test_explicit_values_passed_through(self):
        pw.toPyWATER('1abc', 'B', '90', '1.5', 'No refinement', '', 'average', '2.0', '0.5', '50')
        args = self.captured['args']
        kwargs = self.captured['kwargs']
        self.assertEqual(args[2], '90')
        self.assertEqual(args[3], 1.5)
        self.assertEqual(args[4], 'No refinement')
        self.assertEqual(args[6], 'average')
        self.assertEqual(args[7], 2.0)
        self.assertEqual(kwargs['max_structures'], 50)


class LogCollectorTests(unittest.TestCase):
    def _record(self, msg):
        return logging.LogRecord('PyWATER', logging.INFO, __file__, 0, msg, None, None)

    def test_found_waters_true_on_success_line(self):
        c = pw._LogCollector()
        c.emit(self._record('Found 12 conserved water molecules.'))
        self.assertTrue(c.found_waters())

    def test_found_waters_false_otherwise(self):
        c = pw._LogCollector()
        c.emit(self._record('Filtering by resolution ...'))
        self.assertFalse(c.found_waters())

    def test_summary_prefers_outcome_and_includes_selection(self):
        c = pw._LogCollector()
        c.emit(self._record('Structure selection summary: 100 candidate chains, 40 passed, 40 selected.'))
        c.emit(self._record('4lyw_A has too many waters to cluster. Memory is not enough...'))
        summary = c.summary()
        self.assertIn('too many waters to cluster', summary)
        self.assertIn('Structure selection summary:', summary)

    def test_summary_falls_back_to_selection_only(self):
        c = pw._LogCollector()
        c.emit(self._record('Structure selection summary: 3 candidate chains, 3 passed, 3 selected.'))
        self.assertEqual(c.summary(), 'Structure selection summary: 3 candidate chains, 3 passed, 3 selected.')

    def test_summary_default_when_nothing_matches(self):
        c = pw._LogCollector()
        c.emit(self._record('some unrelated debug line'))
        self.assertIn('pywater.log', c.summary())


class CollectLocalStructuresTests(unittest.TestCase):
    def test_reference_first_then_sorted_ids(self):
        with tempfile.TemporaryDirectory() as d:
            for name in ('a.pdb', 'b.pdb', 'c.pdb'):
                _touch(os.path.join(d, name))
            structures, query_id = pw.collectLocalStructures(d, 'b.pdb')
            self.assertEqual(query_id, 'lc00')
            ids = [pid for pid, _ in structures]
            names = [os.path.basename(p) for _, p in structures]
            self.assertEqual(ids, ['lc00', 'lc01', 'lc02'])
            self.assertEqual(names, ['b.pdb', 'a.pdb', 'c.pdb'])  # reference first, rest sorted

    def test_accepts_full_path_reference(self):
        with tempfile.TemporaryDirectory() as d:
            for name in ('a.pdb', 'b.pdb'):
                _touch(os.path.join(d, name))
            structures, query_id = pw.collectLocalStructures(d, os.path.join(d, 'a.pdb'))
            self.assertEqual(os.path.basename(structures[0][1]), 'a.pdb')

    def test_ignores_non_pdb(self):
        with tempfile.TemporaryDirectory() as d:
            for name in ('a.pdb', 'b.pdb', 'notes.txt'):
                _touch(os.path.join(d, name))
            structures, _ = pw.collectLocalStructures(d, 'a.pdb')
            self.assertEqual(len(structures), 2)

    def test_requires_at_least_two(self):
        with tempfile.TemporaryDirectory() as d:
            _touch(os.path.join(d, 'only.pdb'))
            with self.assertRaises(ValueError):
                pw.collectLocalStructures(d, 'only.pdb')

    def test_reference_must_exist(self):
        with tempfile.TemporaryDirectory() as d:
            for name in ('a.pdb', 'b.pdb'):
                _touch(os.path.join(d, name))
            with self.assertRaises(ValueError):
                pw.collectLocalStructures(d, 'missing.pdb')

    def test_missing_folder(self):
        with self.assertRaises(ValueError):
            pw.collectLocalStructures('/no/such/folder/here', 'a.pdb')

    def test_too_many_files(self):
        with tempfile.TemporaryDirectory() as d:
            for i in range(pw.MAX_LOCAL_STRUCTURES + 1):
                _touch(os.path.join(d, 'f%03d.pdb' % i))
            with self.assertRaises(ValueError):
                pw.collectLocalStructures(d, 'f000.pdb')


class ChainHasCATests(unittest.TestCase):
    def test_detects_chain(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 's.pdb')
            _touch(p, _atom_line('CA', 'A') + _atom_line('CB', 'A'))
            self.assertTrue(pw._chainHasCA(p, 'A'))
            self.assertFalse(pw._chainHasCA(p, 'B'))

    def test_missing_file(self):
        self.assertFalse(pw._chainHasCA('/no/such/file.pdb', 'A'))


class ToPyWATERLocalArgTests(unittest.TestCase):
    def setUp(self):
        self._orig = pw.FindConservedWatersLocal
        self.captured = {}

        def _spy(*args, **kwargs):
            self.captured['args'] = args
            self.captured['kwargs'] = kwargs

        pw.FindConservedWatersLocal = _spy

    def tearDown(self):
        pw.FindConservedWatersLocal = self._orig

    def test_defaults_and_coercion(self):
        pw.toPyWATERLocal('/data/mine', 'a', 'ref.pdb')
        args = self.captured['args']
        self.assertEqual(args[0], '/data/mine')
        self.assertEqual(args[1], 'A')            # chain uppercased
        self.assertEqual(args[2], 'ref.pdb')
        self.assertEqual(args[3], 'Mobility')
        self.assertEqual(args[4], 'complete')
        self.assertEqual(args[5], 2.4)
        self.assertIsInstance(args[5], float)
        self.assertEqual(args[6], 0.7)
        self.assertFalse(args[7])                 # save_sup default off

    def test_save_sup_string_truthy(self):
        pw.toPyWATERLocal('/d', 'A', 'r.pdb', 'Mobility', 'complete', 2.4, 0.7, 'true')
        self.assertTrue(self.captured['args'][7])


if __name__ == '__main__':
    unittest.main()
