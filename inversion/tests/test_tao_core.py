"""Exercise the real PETSc/TAO shared library through the stable C ABI."""
import math
from pathlib import Path
import tempfile
import unittest

from udar_inversion.core import default_library, minimize


@unittest.skipUnless(default_library("tao").exists(), "PETSc/TAO core is not built")
class TAOTests(unittest.TestCase):
    def solve(self, predict, initial, observed, **options):
        with tempfile.TemporaryDirectory() as temporary:
            return minimize(
                predict, initial, [-2.0] * len(initial), [2.0] * len(initial),
                observed, [1.0] * len(observed), backend="tao",
                log_path=Path(temporary) / "tao.log", **options)

    def test_coupled_least_squares(self):
        result = self.solve(lambda x: [2*x[0] + x[1], x[0] - 3*x[1]],
                            [0.0, 0.0], [1.0, 4.0])
        self.assertTrue(result["converged"], result)
        self.assertAlmostEqual(result["log_conductivity"][0], 1.0, places=5)
        self.assertAlmostEqual(result["log_conductivity"][1], -1.0, places=5)

    def test_active_bound(self):
        seen = []

        def prediction(x):
            seen.append(x[0])
            return x

        result = self.solve(prediction, [-2.0], [4.0])
        self.assertTrue(result["converged"], result)
        self.assertAlmostEqual(result["log_conductivity"][0], 2.0, places=6)
        self.assertTrue(all(-2.0 <= value <= 2.0 for value in seen))

    def test_regularization(self):
        result = self.solve(lambda x: x, [0.0], [1.0],
                            regularization=3.0, prior=[0.0])
        self.assertTrue(result["converged"], result)
        self.assertAlmostEqual(result["log_conductivity"][0], 0.25, places=5)

    def test_callback_failure(self):
        def broken(_):
            raise RuntimeError("deliberate TAO forward failure")

        with self.assertRaisesRegex(RuntimeError, "deliberate TAO forward failure"):
            self.solve(broken, [0.0], [1.0])

    def test_iteration_limit_is_not_convergence(self):
        result = self.solve(lambda x: [math.exp(x[0])], [-1.0], [2.0],
                            max_iterations=1)
        self.assertFalse(result["converged"], result)
        self.assertEqual(result["status"], "iteration_limit")


if __name__ == "__main__":
    unittest.main()
