import math
import unittest

import numpy as np

from materiais import Material
from nos import NodalNetwork, NodeKind, TransferKind, LinkDirection
from solver import solve_steady_state, solve_transient


def simple_network():
    net = NodalNetwork()
    hot = net.add_node("hot", NodeKind.DIFFUSION, 0, 0, source=10, temperature=0)
    fixed = net.add_node("fixed", NodeKind.BOUNDARY, 1, 0, fixed_temperature=0)
    net.add_link(hot, fixed, TransferKind.CONDUCTION, conductance=2)
    return net, hot


class SteadyTests(unittest.TestCase):
    def test_analytical_balance_both_backends(self):
        for scipy in (False, True):
            with self.subTest(scipy=scipy):
                net, hot = simple_network()
                result = solve_steady_state(net, prefer_scipy=scipy)
                self.assertTrue(result.success, result.message)
                self.assertAlmostEqual(net.nodes[hot].temperature, 5, places=8)
                self.assertLessEqual(result.residual_norm, 1e-8)

    def test_final_newton_update_counts_as_convergence(self):
        net, hot = simple_network()
        result = solve_steady_state(net, prefer_scipy=False, max_iter=1)
        self.assertTrue(result.success)
        self.assertAlmostEqual(result.x[0], 5, places=8)

    def test_empty_and_boundary_only(self):
        for fixed in (False, True):
            net = NodalNetwork()
            if fixed:
                net.add_node("fixed", NodeKind.BOUNDARY, 0, 0, fixed_temperature=30)
            result = solve_steady_state(net)
            self.assertTrue(result.success)
            self.assertEqual(result.iterations, 0)
            self.assertNotIn("Improper", result.message)

    def test_no_update_has_no_side_effect(self):
        net, hot = simple_network()
        result = solve_steady_state(net, update_network=False)
        self.assertTrue(result.success)
        self.assertEqual(net.nodes[hot].temperature, 0)

    def test_rejects_invalid_solver_arguments(self):
        for kwargs in ({"tol": 0}, {"tol": math.nan}, {"max_iter": 0}, {"max_iter": 1.5},
                       {"z0": [math.nan]}, {"z0": [1, 2]}, {"z0": [[1]]}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                solve_steady_state(simple_network()[0], **kwargs)

    def test_singular_heated_node_does_not_report_success(self):
        net = NodalNetwork()
        net.add_node("isolated", NodeKind.DIFFUSION, 0, 0, source=1)
        result = solve_steady_state(net, max_iter=2)
        self.assertFalse(result.success)

    def test_residual_conserves_undirected_heat_and_transport_direction(self):
        for direction, expected in ((LinkDirection.UNDIRECTED, [-20, 20]),
                                    (LinkDirection.I_TO_J, [0, 20]),
                                    (LinkDirection.J_TO_I, [-20, 0])):
            net = NodalNetwork()
            a = net.add_node("a", NodeKind.FLUID, 0, 0)
            b = net.add_node("b", NodeKind.FLUID, 1, 0)
            net.add_link(a, b, TransferKind.FLUID_TRANSPORT, conductance=2, direction=direction)
            np.testing.assert_allclose(net.residual_steady(np.array([30., 20.])), expected)


class TransientTests(unittest.TestCase):
    def heated_node(self):
        net = NodalNetwork()
        node = net.add_node("mass", NodeKind.DIFFUSION, 0, 0, volume=1,
                            source=1, temperature=0, material=Material("test", rho=1, cp=1))
        return net, node

    def test_reaches_requested_end_time_with_partial_last_step(self):
        for dt in (0.3, 0.6, 2):
            net, node = self.heated_node()
            result = solve_transient(net, t_end=1, dt=dt, max_newton_iter=1)
            self.assertTrue(result.success, result.message)
            self.assertEqual(result.t_values[-1], 1)
            self.assertAlmostEqual(net.nodes[node].temperature, 1, places=8)

    def test_zero_duration_keeps_initial_state(self):
        net, node = self.heated_node()
        result = solve_transient(net, t_end=0, dt=1)
        self.assertEqual(result.t_values, [0])
        self.assertEqual(net.nodes[node].temperature, 0)

    def test_decimal_multiple_does_not_create_near_zero_final_step(self):
        net, node = self.heated_node()
        end = .1 - .01
        result = solve_transient(net, t_end=end, dt=.002)
        self.assertTrue(result.success, result.message)
        self.assertEqual(len(result.t_values), 46)  # Inicial + 45 passos.
        self.assertEqual(result.t_values[-1], end)
        self.assertGreater(min(np.diff(result.t_values)), .0019)
        self.assertAlmostEqual(net.nodes[node].temperature, end, places=9)

    def test_variable_heat_capacity_is_not_frozen(self):
        net, node = self.heated_node()
        net.nodes[node].material = Material("variable", rho=1, cp_func=lambda T, P: 1 + T)
        result = solve_transient(net, t_end=1, dt=0.001, tol=1e-8)
        self.assertTrue(result.success)
        # Integral de (1+T)dT = t; solução contínua independente do integrador.
        self.assertAlmostEqual(net.nodes[node].temperature, math.sqrt(3) - 1, delta=0.0003)

    def test_transient_no_update(self):
        net, node = self.heated_node()
        result = solve_transient(net, t_end=1, dt=0.3, update_network=False)
        self.assertTrue(result.success)
        self.assertEqual(net.nodes[node].temperature, 0)
        self.assertAlmostEqual(result.snapshots[-1][0], 1)

    def test_invalid_transient_arguments(self):
        for kwargs in ({"dt": 0}, {"dt": -1}, {"t_end": -1}, {"save_every": 0},
                       {"max_newton_iter": 0}, {"tol": 0}, {"z0": [math.nan]}):
            options = {"t_end": 1, "dt": 0.1, **kwargs}
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                solve_transient(self.heated_node()[0], **options)
