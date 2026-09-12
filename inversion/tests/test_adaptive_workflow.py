"""Validate adaptive-mesh grouping and staged inversion configuration."""
import json
from pathlib import Path
import tempfile
import unittest

from udar_inversion.forward import CHANNELS, GEOMETRY, ForwardModel, TensorData
from udar_inversion.workflow import _prepare_incremental_update


def row(station, frequency):
    return {
        "station": station,
        "frequency_hz": frequency,
        "spacing_m": 1.0,
        **{name: 0.0 for name in GEOMETRY},
        **{name: float(station + frequency) for name in CHANNELS},
    }


class AdaptiveWorkflowTests(unittest.TestCase):
    def test_station_selection_and_frequency_merge(self):
        data = TensorData([row(0, 2000), row(1, 2000), row(0, 48000)])
        selected = data.select_stations([1])
        self.assertEqual(list(selected.by_key)[0][0], 1.0)
        low = TensorData([row(0, 2000)])
        high = TensorData([row(0, 48000)])
        self.assertEqual(len(TensorData.merge_frequency_groups([low, high]).rows), 2)

    def test_frequency_meshes_partition_the_acquisition(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "udar3p.exe").touch()
            (root / "mesh.mesh").touch()
            config = {
                "executable": "udar3p.exe",
                "mesh": "mesh.mesh",
                "solver": "mumps",
                "persistent": False,
                "stations": 2,
                "frequencies": [2000, 48000],
                "frequency_meshes": [
                    {"frequencies": [2000], "octree": True,
                     "octree_levels": 1},
                    {"frequencies": [48000], "octree": True,
                     "octree_levels": 3, "octree_radius": 0.1},
                ],
            }
            forward = ForwardModel(config, root / "forward", root)
            self.assertEqual(len(forward.forward_contexts), 2)
            self.assertEqual(forward.forward_contexts[1]["octree_levels"], 3)
            forward.close()

    def test_incremental_update_freezes_prior_and_selects_new_stations(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            previous = root / "result.json"
            previous.write_text(json.dumps({"conductivity_s_per_m": {
                "sigma_background": 0.2, "sigma_target": 0.05}}),
                                encoding="utf-8")
            config = {
                "forward": {"sigma_background": 0.1, "sigma_target": 0.1,
                            "stations": 3},
                "parameters": [
                    {"name": "sigma_background", "initial": 0.1,
                     "lower": 0.01, "upper": 2.0},
                    {"name": "sigma_target", "initial": 0.1,
                     "lower": 0.01, "upper": 2.0},
                ],
                "incremental_update": {
                    "previous_result": "result.json",
                    "freeze_parameters": ["sigma_background"],
                    "active_stations": [2],
                },
            }
            observations = TensorData([row(0, 2000), row(2, 2000)])
            observed, parameters, forward, frozen, metadata = \
                _prepare_incremental_update(config, root, observations)
            self.assertEqual(len(observed.rows), 1)
            self.assertEqual([parameter["name"] for parameter in parameters],
                             ["sigma_target"])
            self.assertEqual(frozen, {"sigma_background": 0.2})
            self.assertEqual(forward["station_indices"], [2])
            self.assertEqual(metadata["active_stations"], [2])


if __name__ == "__main__":
    unittest.main()
