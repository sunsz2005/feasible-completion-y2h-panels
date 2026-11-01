import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import balanced_assay_initialization_time as experiment
from balanced_assay_pilot import State,matrix


class Clock:
    value=0.
    def read(self): return self.value


class Checks(unittest.TestCase):
    def setUp(self):
        self.a=matrix(4,np.array([[0,2],[1,3]],dtype=np.int32))
        self.state=State(self.a,[0,1],[2,3])

    def run_fake(self,method,allowance):
        clock=Clock(); searches=[]
        def initialize(a,arm,floor):
            clock.value+=3 if arm=='A' else 6
            return self.state.copy(),self.state.copy()
        def solve(initial,floor,seconds,seed,replace):
            searches.append((seconds,seed,replace)); clock.value+=seconds
            return initial.copy(),dict(search_seconds=seconds,starts=1,perturb_accepted=0)
        with patch.object(experiment.time,'perf_counter',clock.read),patch.object(experiment,'initialize',initialize),patch.object(experiment,'multistart',solve):
            state,meta=experiment.pipeline(self.a,method,0.,allowance,123)
        return state,meta,searches

    def test_initialization_is_paid_and_identical_solver_used(self):
        _,a,sa=self.run_fake('A',10.)
        _,b,sb=self.run_fake('B',10.)
        self.assertEqual(sa,[(7.,123,True)])
        self.assertEqual(sb,[(4.,123,True)])
        self.assertEqual(a['total_seconds'],b['total_seconds'])

    def test_initialization_over_budget_fallback_is_visible(self):
        _,m,calls=self.run_fake('B',5.)
        self.assertTrue(m['initialization_fallback'])
        self.assertFalse(m['timing_valid'])
        self.assertEqual(calls,[])

    def test_atomic_initialization_within_tolerance_retained(self):
        _,m,calls=self.run_fake('B',5.9)
        self.assertFalse(m['initialization_fallback'])
        self.assertTrue(m['initialization_over_target'])
        self.assertTrue(m['timing_valid'])
        self.assertEqual(calls,[])

    def test_no_free_search_when_budget_consumed(self):
        _,m,calls=self.run_fake('A',3.)
        self.assertEqual(m['search_allowance'],0.)
        self.assertEqual(calls,[])

    def test_360_unique_outputs(self):
        files=experiment.expected_files()
        self.assertEqual(len(files),360)
        self.assertEqual(len(set(files)),360)

    def test_scoring_requires_complete_seal(self):
        with tempfile.TemporaryDirectory() as directory,patch.object(experiment,'OUT',Path(directory)):
            with self.assertRaises(FileNotFoundError): experiment.verify_seal()


if __name__=='__main__':
    unittest.main()
