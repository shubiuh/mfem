"""Exercise the real ROL shared library, including constraints and failures."""
import math
from pathlib import Path
import tempfile
import unittest

from udar_inversion.core import minimize


class ROLTests(unittest.TestCase):
    def solve(self, predict, initial, observed, **kw):
        with tempfile.TemporaryDirectory() as temporary:
            return minimize(predict, initial, [-2.0]*len(initial), [2.0]*len(initial),
                            observed, [1.0]*len(observed),
                            log_path=Path(temporary)/"rol.log", **kw)

    def test_coupled_least_squares(self):
        result = self.solve(lambda x: [2*x[0]+x[1], x[0]-3*x[1]],
                            [0.0, 0.0], [1.0, 4.0])
        self.assertTrue(result["converged"], result)
        self.assertAlmostEqual(result["log_conductivity"][0], 1, places=5)
        self.assertAlmostEqual(result["log_conductivity"][1], -1, places=5)

    def test_active_bound(self):
        seen = []
        def prediction(x):
            seen.append(x[0])
            return x
        result = self.solve(prediction, [-2.0], [4.0])
        self.assertTrue(result["converged"], result)
        self.assertAlmostEqual(result["log_conductivity"][0], 2.0, places=6)
        self.assertTrue(all(-2.0 <= x <= 2.0 for x in seen))

    def test_regularization(self):
        result = self.solve(lambda x: x, [0.0], [1.0],
                            regularization=3.0, prior=[0.0])
        self.assertTrue(result["converged"], result)
        self.assertAlmostEqual(result["log_conductivity"][0], 0.25, places=5)

    def test_nonbinary_bound_roundoff(self):
        with tempfile.TemporaryDirectory() as temporary:
            seen = []
            def prediction(x):
                seen.append(x[0])
                return [100*x[0]]
            lo, hi = math.log(0.01), math.log(2.0)
            result = minimize(prediction, [math.log(0.1)], [lo], [hi], [400], [1],
                              log_path=Path(temporary)/"rol.log")
            self.assertTrue(result["converged"], result)
            self.assertTrue(all(lo <= x <= hi for x in seen))
            self.assertAlmostEqual(result["log_conductivity"][0], hi, places=12)

    def test_failure_crosses_abi(self):
        def broken(_):
            raise RuntimeError("deliberate forward failure")
        with self.assertRaisesRegex(RuntimeError, "deliberate forward failure"):
            self.solve(broken, [0.0], [1.0])

    def test_nonfinite_and_wrong_size(self):
        with self.assertRaisesRegex(ValueError, "non-finite"):
            self.solve(lambda _: [math.nan], [0.0], [1.0])
        with self.assertRaisesRegex(ValueError, "wrong data size"):
            self.solve(lambda _: [1.0, 2.0], [0.0], [1.0])

    def test_iteration_limit_is_not_convergence(self):
        result = self.solve(lambda x: [math.exp(x[0])], [-1.0], [2.0], max_iterations=1)
        self.assertFalse(result["converged"], result)
        self.assertEqual(result["status"], "iteration_limit")


if __name__ == "__main__":
    unittest.main()
