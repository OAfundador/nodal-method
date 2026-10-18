import unittest

from solver import solve_transient
from validacao_livro import BookBenchmarks


class LiteratureTests(unittest.TestCase):
    def test_four_published_cases(self):
        for result in BookBenchmarks().run():
            with self.subTest(case=result.name):
                self.assertTrue(result.passed, result.to_dict())

    def test_sphere_time_refinement_and_full_history(self):
        errors = []
        for dt in (2, 1, 0.5):
            net, node = BookBenchmarks.sphere_network()
            solution = solve_transient(net, t_end=100, dt=dt, tol=1e-9)
            self.assertTrue(solution.success)
            error = max(abs(snapshot[0] - BookBenchmarks.sphere_exact(t))
                        for t, snapshot in zip(solution.t_values, solution.snapshots))
            errors.append(error)
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])
        self.assertLess(errors[2], 0.032)
