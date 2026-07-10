"""
Unit tests for PyWATER's non-GUI logic.

Run from the repo root with:

    .venv/bin/python -m unittest discover -s tests

These use a stubbed PyMOL (see ``tests/pymol_stubs.py``) so they need only
NumPy/SciPy, not a running PyMOL.
"""

import logging
import unittest

from tests.pymol_stubs import pywater as pw


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


if __name__ == '__main__':
    unittest.main()
