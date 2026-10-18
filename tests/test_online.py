import unittest
import numpy as np

from validacao_online import OnlineBenchmarks
from visual.modelo import VisualModel, nafems_t4_case


class PublicBenchmarkTests(unittest.TestCase):
    def test_public_references_and_heat_generation(self):
        for result in OnlineBenchmarks().run():
            with self.subTest(result=result.name):
                self.assertTrue(result.passed, result.to_dict())

    def test_t4_mesh_refinement_approaches_published_surface_value(self):
        errors = []
        for nx, ny in ((6, 10), (12, 20), (18, 30)):
            result = VisualModel(nafems_t4_case(nx, ny)).solve()
            errors.append(abs(result["probes"][0]["T_C"] - 18.25))
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])
        self.assertLess(errors[2], .15)

    def test_generated_slab_whole_profile_refines_quadratically(self):
        errors = []
        for cells in (10, 20, 40):
            positions, computed, exact = OnlineBenchmarks.generated_slab(cells)
            errors.append(float(np.sqrt(np.mean((computed - exact) ** 2))))
            self.assertTrue(np.all(np.diff(computed) < 0))
        self.assertAlmostEqual(errors[0] / errors[1], 4, places=5)
        self.assertAlmostEqual(errors[1] / errors[2], 4, places=5)

    def test_boundary_probe_is_surface_temperature_and_rejects_interior(self):
        case = nafems_t4_case(6, 10)
        result = VisualModel(case).solve()
        boundary_cells = [c["T_C"] for c in result["cells"] if c["i"] == result["nx"] - 1 and c["j"] in (1, 2)]
        self.assertLess(result["probes"][0]["T_C"], sum(boundary_cells) / len(boundary_cells))
        case["probes"][0]["x"] = .3
        with self.assertRaisesRegex(ValueError, "face externa"):
            VisualModel(case).solve()
