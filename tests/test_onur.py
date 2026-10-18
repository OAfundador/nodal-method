import unittest
import numpy as np

from validacao_onur import OnurCase, OnurSlab, combined_benchmark, pure_convection_exact


class OnurTests(unittest.TestCase):
    def test_combined_heat_transfer_and_generation(self):
        result = combined_benchmark()
        self.assertTrue(result["passed"], result)
        self.assertGreater(result["theta_insulated"], result["theta_surface"])

    def test_steady_profiles_at_published_parameters(self):
        for ratio in (.25, .5, .75):
            for Bi, Nr in ((1, 1), (5, 1), (10, 1), (1, 5), (1, 10)):
                for u in (0, 1, 5):
                    with self.subTest(ratio=ratio, Bi=Bi, Nr=Nr, u=u):
                        case = OnurCase(Bi, Nr, u, ratio, ratio)
                        slab = OnurSlab(case, intervals=10)
                        result = slab.steady()
                        np.testing.assert_allclose(result.x, case.analytical_steady(slab.x), atol=1e-8, rtol=0)
                        self.assertAlmostEqual(case.outward_flux(result.x[-1]), u, places=8)
                        self.assertTrue(np.all(result.x > 0))

    def test_pure_convection_transient_and_temporal_refinement(self):
        case = OnurCase(Bi=1, Nr=0, u=0)
        errors = []
        for dt in (.04, .02, .01):
            slab = OnurSlab(case, intervals=30)
            result = slab.transient(fo_end=.4, dt=dt)
            exact = pure_convection_exact(slab.x, .4)
            errors.append(float(np.max(np.abs(result.snapshots[-1] - exact))))
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])
        self.assertLess(errors[-1], .002)
        self.assertGreater(errors[0] / errors[1], 1.7)

    def test_transient_energy_and_cooling_bounds(self):
        slab = OnurSlab(OnurCase(u=0), intervals=12)
        result = slab.transient(fo_end=.5, dt=.01)
        values = np.asarray(result.snapshots)
        self.assertLess(slab.transient_balance_error(result), 1e-7)
        self.assertTrue(np.all(values >= .25))
        self.assertTrue(np.all(values <= 1. + 1e-10))
        self.assertTrue(np.all(np.diff(values, axis=0) <= 1e-10))

    def test_plot_normalization_and_validation(self):
        case = OnurCase()
        np.testing.assert_allclose(case.paper_ordinate([1., .25, 1.75]), [0., 1., -1.])
        for values in ({"Bi": -1}, {"Nr": -1}, {"Bi": 0, "Nr": 0}, {"u": -1},
                       {"fluid_ratio": 0}, {"environment_ratio": float("nan")}):
            with self.assertRaises(ValueError):
                OnurCase(**values)


if __name__ == "__main__":
    unittest.main()
