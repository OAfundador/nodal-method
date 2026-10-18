from pathlib import Path
import tempfile
import unittest

import numpy as np

from exemplos.nodal_repl import CaseRunner, EstadoREPL, executar
from visual.modelo import VisualModel, demo_cases
from escoamento import water_properties


class InputModesTests(unittest.TestCase):
    def assert_same_field(self, data):
        visual = VisualModel(data)
        result = visual.solve()
        commands = visual.to_repl()
        interactive = EstadoREPL()
        for line in commands.splitlines():
            executar(interactive, line)
        temp_root = Path(__file__).resolve().parents[1] / 'tmp'
        temp_root.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=temp_root) as directory:
            path = Path(directory) / 'modelo.txt'
            path.write_text(commands, encoding='utf-8-sig')
            text = CaseRunner().run_file(path)
        for state in (interactive, text):
            self.assertTrue(state.last_result.success)
            self.assertEqual(len(state.net.nodes), result['nodes'])
            self.assertEqual(len(state.net.links), result['links'])
            self.assertAlmostEqual(sum(n.source for n in state.net.nodes.values()), result['source_W'], places=10)
            self.assertLess(np.linalg.norm(state.net.residual_steady(state.last_result.x)), 1e-7)
            for cell in result['cells']:
                node = state.net.nodes[state.net.cell_to_node[(cell['j'], cell['i'])]]
                self.assertAlmostEqual(node.temperature, cell['T_C'], places=8)

    def test_all_public_visual_examples_roundtrip(self):
        for name, data in demo_cases().items():
            with self.subTest(name=name):
                self.assert_same_field(data)

    def test_generation_radiation_and_nonuniform_interfaces_roundtrip(self):
        data = demo_cases()['composite']
        data['regions'][0].update(x1=.0007, source=100000, conduction_axes='x')
        data['regions'][1]['x0'] = .0007
        data['mesh'].update(nx=3, ny=2)
        data['boundaries']['left'] = {'kind': 'adiabatic'}
        data['boundaries']['right'] = {'kind': 'convection_radiation', 'h': 20, 'T_inf': 20,
                                       'emissivity': .8, 'view_factor': 1, 'T_env': 30}
        self.assert_same_field(data)

    def test_mass_flow_roundtrip(self):
        data = demo_cases()['uniform']
        data['materials'][0].update(phase='fluid', cp=4000)
        data['regions'][0].update(source=100000, flow={'direction': 'x+', 'mdot': .01, 'T_in': 20})
        data['boundaries'] = {s: {'kind': 'adiabatic'} for s in ('left', 'right', 'top', 'bottom')}
        self.assert_same_field(data)

    def test_nonuniform_parallel_flow_weights_close_total_heat_balance(self):
        data = demo_cases()['uniform']
        data['mesh'].update(nx=3, ny=3, h_internal=20)
        data['materials'][0].update(phase='fluid', cp=4000)
        data['materials'].append({'id': 'solid', 'name': 'Sólido', 'phase': 'solid', 'k': 10})
        data['regions'][0].update(y0=.01, source=100000, flow={'direction': 'y-', 'mdot': .01, 'T_in': 20})
        data['regions'].append({'id': 'block', 'material': 'solid', 'x0': .0137, 'x1': .06, 'y0': 0., 'y1': .01})
        data['boundaries'] = {s: {'kind': 'adiabatic'} for s in ('left', 'right', 'top', 'bottom')}
        result = VisualModel(data).solve()
        outlet = [c for c in result['cells'] if c['region'] == 'r1' and abs(c['y0'] - .01) < 1e-12]
        mixed = sum(c['T_C'] * (c['x1'] - c['x0']) / .06 for c in outlet)
        self.assertAlmostEqual(.01 * 4000 * (mixed - 20), result['source_W'], places=8)

    def test_registered_fuel_conductivity_roundtrip(self):
        data = demo_cases()['uniform']
        data['materials'][0]['property_model'] = 'u3si2_al'
        self.assert_same_field(data)

    def test_optional_variable_water_roundtrip_and_independent_axial_march(self):
        try:
            water_properties(20, 2e5)
        except (ImportError, RuntimeError):
            self.skipTest('IAPWS97 é uma dependência opcional.')
        data = demo_cases()['uniform']
        data['mesh'].update(nx=1, ny=6)
        data['materials'][0].update(phase='fluid', property_model='iapws97')
        data['regions'][0].update(source=1e6, flow={'direction': 'y-', 'model': 'iapws97',
            'vdot': 1e-6, 'T_in': 20, 'P_in': 2e5, 'P_out': 1.5e5})
        data['boundaries'] = {s: {'kind': 'adiabatic'} for s in ('left', 'right', 'top', 'bottom')}
        self.assert_same_field(data)
        result = VisualModel(data).solve()
        temperature, pressure = 20., 2e5
        power_per_cell = result['source_W'] / 6
        for index, cell in enumerate(sorted(result['cells'], key=lambda c: c['j'], reverse=True)):
            water = water_properties(temperature, pressure)
            temperature += power_per_cell / (water.rho * 1e-6 * water.cp)
            self.assertAlmostEqual(cell['T_C'], temperature, places=8)
            pressure = 2e5 + (index + 1) / 6 * (1.5e5 - 2e5)
