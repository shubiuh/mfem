import unittest

from udar_inversion.forward import CHANNELS, GEOMETRY, TensorData


def row(station=0):
    return {"station": station, "frequency_hz": 2000, "spacing_m": 1,
            **{k: 0.0 for k in GEOMETRY}, **{k: float(station+1) for k in CHANNELS}}


class TensorTests(unittest.TestCase):
    def test_reordered_rows_align(self):
        a, b = TensorData([row(0), row(1)]), TensorData([row(1), row(0)])
        self.assertEqual(a.vector(), b.vector(reference=a))

    def test_geometry_mismatch_rejected(self):
        changed = row()
        changed["rx_x"] = 1
        with self.assertRaisesRegex(ValueError, "geometry differs"):
            TensorData([changed]).vector(reference=TensorData([row()]))

    def test_keys_and_duplicates_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            TensorData([row(), row()])
        with self.assertRaisesRegex(ValueError, "keys differ"):
            TensorData([row(1)]).vector(reference=TensorData([row(0)]))

    def test_nonfinite_rejected(self):
        bad = row()
        bad["Hxx_imag"] = float("nan")
        with self.assertRaisesRegex(ValueError, "non-finite"):
            TensorData([bad])


if __name__ == "__main__":
    unittest.main()
