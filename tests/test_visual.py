import json
import math
import threading
import unittest
import urllib.error
import urllib.request

import numpy as np

from visual.modelo import VisualModel, demo_cases
from visual.servidor import create_server
from scipy.optimize import brentq


class VisualModelTests(unittest.TestCase):
    def test_thin_layers_keep_dimensions_and_use_conforming_mesh(self):
        data = demo_cases()["composite"]
        data["regions"][0]["x1"] = .0007
        data["regions"][1]["x0"] = .0007
        data["mesh"].update(nx=3, ny=2)
        result = VisualModel(data).solve()
        self.assertEqual(result["nx"], 3)
        self.assertIn(.0007, result["x_edges"])
        flux = 80 / (.0007 / 10 + (.06 - .0007) / 20)
        for cell in result["cells"]:
            x = (cell["x0"] + cell["x1"]) / 2
            expected = 100 - flux * (x / 10 if x < .0007 else .0007 / 10 + (x - .0007) / 20)
            self.assertAlmostEqual(cell["T_C"], expected, places=8)

    def test_uniform_wall_matches_exact_field_at_every_cell(self):
        data = demo_cases()["uniform"]
        result = VisualModel(data).solve()
        for cell in result["cells"]:
            x = (cell["x0"] + cell["x1"]) / 2
            self.assertAlmostEqual(cell["T_C"], 100 - 80 * x / .06, places=9)

    def test_composite_wall_matches_exact_resistances_and_connects_interface(self):
        result = VisualModel(demo_cases()["composite"]).solve()
        flux = 80 / (.03 / 10 + .03 / 20)
        interface_T = 100 - flux * .03 / 10
        self.assertEqual(len(result["contacts"]), result["ny"])
        for cell in result["cells"]:
            x = (cell["x0"] + cell["x1"]) / 2
            expected = 100 - flux * x / 10 if x < .03 else interface_T - flux * (x - .03) / 20
            self.assertAlmostEqual(cell["T_C"], expected, places=9)

    def test_generated_heat_has_second_order_spatial_convergence(self):
        errors = []
        for nx in (6, 12, 24):
            data = demo_cases()["uniform"]
            data["mesh"].update(nx=nx, ny=2)
            data["regions"][0]["source"] = 10000
            result = VisualModel(data).solve()
            error = max(abs(c["T_C"] - (100 - 80 * ((c["x0"] + c["x1"]) / 2) / .06
                        + 10000 / (2 * 10) * ((c["x0"] + c["x1"]) / 2) * (.06 - (c["x0"] + c["x1"]) / 2)))
                        for c in result["cells"])
            errors.append(error)
            self.assertAlmostEqual(result["source_W"], .0024 * 10, places=10)
        self.assertAlmostEqual(errors[0] / errors[1], 4, places=5)
        self.assertAlmostEqual(errors[1] / errors[2], 4, places=5)

    def test_chip_2d_regression_is_distinct_from_lumped_chip(self):
        result = VisualModel(demo_cases()["chip"]).solve()
        self.assertAlmostEqual(result["maximum_C"], 66.485404, places=5)

    def test_overlap_is_rejected(self):
        data = demo_cases()["composite"]
        data["regions"][1]["x0"] = .025
        with self.assertRaises(ValueError):
            VisualModel(data)

    def test_unanchored_component_is_rejected_even_with_zero_source(self):
        data = demo_cases()["uniform"]
        data["boundaries"] = {side: {"kind": "adiabatic"} for side in ("left", "right", "top", "bottom")}
        with self.assertRaisesRegex(ValueError, "referência térmica"):
            VisualModel(data).solve()

    def test_large_mesh_and_invalid_properties_are_rejected(self):
        for change in (lambda d: d["mesh"].update(nx=100, ny=100),
                       lambda d: d["materials"][0].update(k=-1),
                       lambda d: d["domain"].update(width=math.nan)):
            data = demo_cases()["uniform"]
            change(data)
            with self.assertRaises(ValueError):
                VisualModel(data)

    def test_missing_boundary_value_is_not_silently_zero(self):
        data = demo_cases()["uniform"]
        data["boundaries"]["left"] = {"kind": "temperature", "T": None}
        with self.assertRaises(ValueError):
            VisualModel(data)

    def test_generated_wall_with_shared_convective_and_radiative_surface(self):
        # Caso didático com solução permanente independente: q'''L sai pela direita.
        for kind, h in (("radiation", 0), ("convection_radiation", 20)):
            with self.subTest(kind=kind):
                data = demo_cases()["uniform"]
                data["regions"][0]["source"] = 1e5
                data["boundaries"]["left"] = {"kind": "adiabatic"}
                data["boundaries"]["right"] = {"kind": kind, "h": h, "T_inf": 20,
                    "emissivity": .8, "view_factor": 1, "T_env": 30}
                data["probes"] = [{"name": "superfície", "x": .06, "y": .02}]
                surface = brentq(lambda T: h * (T - 20) + .8 * 5.670374419e-8 *
                                  ((T + 273.15) ** 4 - 303.15 ** 4) - 1e5 * .06, 30, 1500)
                result = VisualModel(data).solve()
                self.assertAlmostEqual(result["probes"][0]["T_C"], surface, places=7)
                self.assertGreater(result["maximum_C"], surface)
                self.assertLess(result["residual_W"], 1e-7)

    def test_radiative_inputs_require_valid_absolute_temperature_and_factors(self):
        for key, value in (("emissivity", 1.1), ("emissivity", None), ("view_factor", -.1),
                           ("T_env", -274), ("T_env", None)):
            data = demo_cases()["uniform"]
            data["boundaries"]["right"] = {"kind": "radiation", "emissivity": .8,
                                           "view_factor": 1, "T_env": 20, key: value}
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                VisualModel(data)

    def test_fluid_requires_explicit_inputs_and_closes_energy_balance(self):
        data = demo_cases()["uniform"]
        data["materials"][0].update(phase="fluid", cp=4000)
        data["regions"][0]["source"] = 100000
        data["boundaries"] = {side: {"kind": "adiabatic"} for side in ("left", "right", "top", "bottom")}
        with self.assertRaisesRegex(ValueError, "vazão"):
            VisualModel(data).solve()
        data["regions"][0]["flow"] = {"direction": "x+", "mdot": .01, "T_in": 20}
        result = VisualModel(data).solve()
        outlet = np.mean([cell["T_C"] for cell in result["cells"] if cell["i"] == result["nx"] - 1])
        self.assertAlmostEqual(.01 * 4000 * (outlet - 20), result["source_W"], places=7)


class LocalApiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = create_server(0)
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.base = f"http://127.0.0.1:{cls.server.server_port}"

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join(timeout=2)

    def test_serves_only_interface_assets(self):
        with urllib.request.urlopen(self.base + "/") as response:
            self.assertEqual(response.status, 200)
            self.assertIn("Editor visual", response.read().decode("utf-8"))
        with self.assertRaises(urllib.error.HTTPError) as error:
            urllib.request.urlopen(self.base + "/solver.py")
        self.assertEqual(error.exception.code, 404)

    def test_solve_endpoint_runs_real_solver(self):
        request = urllib.request.Request(self.base + "/api/solve", data=json.dumps(demo_cases()["uniform"]).encode(), headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(request) as response:
            result = json.load(response)
        self.assertTrue(result["success"])
        self.assertEqual(len(result["cells"]), 96)
        self.assertLessEqual(result["residual_W"], 1e-7)

    def test_local_case_files_are_not_served(self):
        for path in ("/comparacao/relatorio.html", "/comparacao/resumo.json", "/comparacao/saida_placa.csv"):
            with self.subTest(path=path), self.assertRaises(urllib.error.HTTPError) as error:
                urllib.request.urlopen(self.base + path)
            self.assertEqual(error.exception.code, 404)
        request = urllib.request.Request(self.base + "/api/compare", data=b"{}")
        with self.assertRaises(urllib.error.HTTPError) as error:
            urllib.request.urlopen(request)
        self.assertEqual(error.exception.code, 404)

    def test_foreign_origin_and_bad_data_are_rejected(self):
        for headers, payload, expected in (({"Origin": "https://example.invalid"}, {}, 403),
                                            ({}, {"domain": {}}, 400), ({}, {"domain": []}, 400),
                                            ({}, {"materials": [None]}, 400)):
            request = urllib.request.Request(self.base + "/api/solve", data=json.dumps(payload).encode(), headers=headers)
            with self.assertRaises(urllib.error.HTTPError) as error:
                urllib.request.urlopen(request)
            self.assertEqual(error.exception.code, expected)
