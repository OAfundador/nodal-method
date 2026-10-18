from copy import deepcopy
import csv
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from associacao import ModelAssembly, component_from_geometry
from escoamento import water_properties
from exemplos.associar_componentes import build_example
from exemplos.nodal_repl import CaseRunner, EstadoREPL, executar
from visual.modelo import VisualModel, demo_cases, heated_channel_case
from visual.relatorio_associacao import export_report

ROOT = Path(__file__).resolve().parents[1]


class AssemblyTests(unittest.TestCase):
    def setUp(self):
        (ROOT / 'tmp').mkdir(exist_ok=True)

    def assert_same_results(self, first, second):
        self.assertEqual(first['count'], second['count'])
        for a, b in zip(first['instances'], second['instances']):
            self.assertEqual((a['id'], a['component']), (b['id'], b['component']))
            self.assertAlmostEqual(a['result']['source_W'], b['result']['source_W'], places=10)
            self.assertEqual(a['result']['nodes'], b['result']['nodes'])
            self.assertEqual([(c['region'], c['material'], c['i'], c['j']) for c in a['result']['cells']],
                             [(c['region'], c['material'], c['i'], c['j']) for c in b['result']['cells']])
            np.testing.assert_allclose([c['T_C'] for c in a['result']['cells']],
                                       [c['T_C'] for c in b['result']['cells']], rtol=0, atol=1e-8)
            self.assertLess(b['result']['residual_W'], 1e-7)

    def test_25_channels_close_analytic_energy_balance_and_preserve_template(self):
        group = build_example()
        original = group.to_dict()
        progress = []
        result = group.solve(on_progress=lambda *args: progress.append(args))
        self.assertEqual(group.to_dict(), original)
        self.assertEqual((result['count'], result['total_cells']), (25, 500))
        self.assertEqual(progress[-1][:2], (25, 25))
        for instance in result['instances']:
            factor = instance['source_factor']
            self.assertAlmostEqual(instance['result']['source_W'], 8 * factor, places=10)
            fluid = next(r for r in instance['result']['region_results'] if r['id'] == 'agua')
            self.assertAlmostEqual(fluid['outlet_C'], 20 + 8 * factor / (.01 * 4200), places=9)
        self.assertTrue(all(abs(c['T_C'] - 20) < 1e-9 for c in result['instances'][0]['result']['cells']))
        self.assertAlmostEqual(result['maximum_C'], 28.647619047619, places=8)
        # Returned definitions and caller inputs cannot mutate the library.
        original['components']['canal']['regions'][0]['source'] = -999
        self.assertNotEqual(group.components['canal']['regions'][0]['source'], -999)

    def test_python_json_portable_txt_and_interactive_have_identical_fields(self):
        group = build_example()
        expected = group.solve()
        for imported in (ModelAssembly(json.loads(json.dumps(group.to_dict()))),
                         ModelAssembly.from_repl(group.to_repl())):
            self.assert_same_results(expected, imported.solve())
        interactive = EstadoREPL()
        for line in group.to_repl().splitlines():
            executar(interactive, line)
        self.assert_same_results(expected, interactive.assembly.result)
        with tempfile.TemporaryDirectory(dir=ROOT / 'tmp') as directory:
            path = Path(directory) / 'conjunto.txt'
            path.write_text(group.to_repl(), encoding='utf-8-sig')
            state = CaseRunner().run_file(path)
            self.assert_same_results(expected, state.assembly.result)

    def test_public_command_file_cli_and_csv(self):
        expected = build_example().solve()
        with tempfile.TemporaryDirectory(dir=ROOT / 'tmp') as directory:
            output = Path(directory) / 'temperaturas.csv'
            process = subprocess.run([sys.executable, '-X', 'utf8', 'exemplos/nodal_repl.py',
                'exemplos/associacao_25_canais.txt', '--quiet', '--csv', str(output)],
                cwd=ROOT, capture_output=True, text=True, encoding='utf-8', timeout=60)
            self.assertEqual(process.returncode, 0, process.stderr)
            self.assertIn('25 instâncias calculadas', process.stdout)
            rows = list(csv.DictReader(output.read_text(encoding='utf-8-sig').splitlines()))
        self.assertEqual(len(rows), expected['total_cells'])
        reference = [c['T_C'] for v in expected['instances'] for c in v['result']['cells']]
        self.assertEqual([row['regiao'] for row in rows], [c['region'] for v in expected['instances'] for c in v['result']['cells']])
        np.testing.assert_allclose([float(row['temperatura_C']) for row in rows], reference, rtol=0, atol=1e-8)
        self.assertAlmostEqual(float(rows[-1]['temperatura_K']) - float(rows[-1]['temperatura_C']), 273.15)

    def test_different_components_keep_their_meshes_materials_and_boundaries(self):
        group = ModelAssembly()
        wall = demo_cases()['uniform']
        wall['materials'][0]['k_poly'] = [10, .01]
        channel = heated_channel_case()
        group.add_component('parede', wall)
        group.add_component('canal', channel)
        group.associate('parede_1', 'parede', source_factor=0)
        group.associate('canal_1', 'canal', source_factor=1.5)
        result = group.solve()
        for spec, instance in zip((wall, channel), result['instances']):
            spec = deepcopy(spec)
            for region in spec['regions']:
                region['source'] = region.get('source', 0) * instance['source_factor']
            single = VisualModel(spec).solve()
            self.assertEqual(single['cells'], instance['result']['cells'])
        self.assertNotEqual(result['instances'][0]['result']['nx'], result['instances'][1]['result']['nx'])
        # A zero source factor does not erase the wall's prescribed temperatures.
        self.assertGreater(result['instances'][0]['result']['maximum_C'], 20)

    def test_capture_geometry_retains_sources_flow_polynomial_and_thickness(self):
        data = heated_channel_case()
        data['materials'][0]['k_poly'] = [15, .002]
        runner = CaseRunner()
        runner.run_lines(VisualModel(data).to_repl().splitlines())
        captured = component_from_geometry(runner.state.geom, **runner.state.build_options)
        self.assertEqual(captured['domain'], data['domain'])
        self.assertEqual(captured['regions'][1]['flow'], data['regions'][1]['flow'])
        self.assertEqual(captured['materials'][0]['k_poly'], data['materials'][0]['k_poly'])
        np.testing.assert_allclose([c['T_C'] for c in VisualModel(captured).solve()['cells']],
                                   [c['T_C'] for c in VisualModel(data).solve()['cells']], rtol=0, atol=1e-8)
        executar(runner.state, 'component_save canal')
        executar(runner.state, 'assembly "Dois canais"')
        executar(runner.state, 'associate primeiro component=canal')
        executar(runner.state, 'associate segundo component=canal source_factor=0')
        executar(runner.state, 'solve_assembly')
        self.assertEqual(runner.state.assembly.result['count'], 2)

    def test_optional_registered_water_and_fuel_survive_component_capture(self):
        try:
            water_properties(20, 2e5)
        except (ImportError, RuntimeError):
            self.skipTest('IAPWS97 é uma dependência opcional.')
        data = heated_channel_case()
        data['materials'][0]['property_model'] = 'u3si2_al'
        data['materials'][1]['property_model'] = 'iapws97'
        data['regions'][1]['flow'] = {'model': 'iapws97', 'direction': 'y-', 'vdot': 1e-5,
            'P_in': 2e5, 'P_out': 1.5e5, 'T_in': 20, 'h_model': 'dittus_boelter', 'hydraulic_diameter': .01}
        runner = CaseRunner()
        runner.run_lines(VisualModel(data).to_repl().splitlines())
        captured = component_from_geometry(runner.state.geom, **runner.state.build_options)
        self.assertEqual([m['property_model'] for m in captured['materials']], ['u3si2_al', 'iapws97'])
        group = ModelAssembly()
        group.add_component('canal', captured)
        group.associate('um', 'canal')
        single = VisualModel(data).solve()
        np.testing.assert_allclose([c['T_C'] for c in group.solve()['instances'][0]['result']['cells']],
                                   [c['T_C'] for c in single['cells']], rtol=0, atol=1e-8)

    def test_invalid_references_factors_and_resource_limits_are_rejected(self):
        group = ModelAssembly()
        group.add_component('canal', heated_channel_case())
        for factor in (-1, float('nan'), float('inf'), True):
            with self.subTest(factor=factor), self.assertRaises(ValueError):
                group.associate('ruim', 'canal', source_factor=factor)
        with self.assertRaises(ValueError):
            group.associate('ruim', 'inexistente')
        group.associate('um', 'canal')
        with self.assertRaises(ValueError):
            group.associate('um', 'canal')
        group.MAX_INSTANCES = 1
        with self.assertRaises(ValueError):
            group.associate('dois', 'canal')
        group.MAX_TOTAL_CELLS = 19
        with self.assertRaisesRegex(ValueError, 'células'):
            group.solve()
        self.assertIsNone(group.result)

    def test_new_or_loaded_unsolved_group_never_exports_an_old_single_solution(self):
        with tempfile.TemporaryDirectory(dir=ROOT / 'tmp') as directory:
            path = Path(directory) / 'conjunto.json'
            path.write_text(json.dumps(build_example().to_dict()), encoding='utf-8')
            state = CaseRunner().run_lines(['node fixo boundary T_fixed=25', 'solve'])
            for command in ('assembly "Novo conjunto"', f'assembly_load file="{path.as_posix()}"'):
                executar(state, command)
                output = Path(directory) / 'nao_exportar.csv'
                with self.assertRaisesRegex(ValueError, 'Resolva'):
                    executar(state, f'export csv="{output.as_posix()}"')
                self.assertFalse(output.exists())

    def test_relative_txt_component_import_preserves_context_and_skips_exports(self):
        with tempfile.TemporaryDirectory(dir=ROOT / 'tmp') as directory:
            folder = Path(directory)
            forbidden = folder / 'exportacao_omitida.csv'
            (folder / 'componente.txt').write_text('reset\n' +
                (ROOT / 'exemplos/componentes/canal_aquecido.txt').read_text(encoding='utf-8') +
                f'export csv="{forbidden.as_posix()}"\n', encoding='utf-8')
            path = folder / 'grupo.txt'
            path.write_text('reset\ncomponent_load canal file=componente.txt\nassembly "Grupo"\n'
                            'associate um component=canal\nsolve_assembly\n', encoding='utf-8')
            state = CaseRunner().run_file(path)
            self.assertFalse(forbidden.exists())
            self.assertAlmostEqual(state.assembly.result['maximum_C'], 24.32380952381, places=8)
            with self.assertRaisesRegex(ValueError, 'pressão global'):
                component_from_geometry(VisualModel(heated_channel_case()).geom, pressure=1e5)

    def test_failed_solve_identifies_instance_and_does_not_retain_old_result(self):
        group = ModelAssembly()
        group.add_component('canal', heated_channel_case())
        group.associate('canal_ruim', 'canal')
        group.solve()
        group.components['canal']['regions'][1]['flow']['mdot'] = 0
        with self.assertRaisesRegex(ValueError, 'canal_ruim'):
            group.solve()
        self.assertIsNone(group.result)

    def test_portable_names_are_literal_and_general_commands_are_not_imported(self):
        group = ModelAssembly(name='Conjunto $nome $(1/0)\nsegunda linha')
        data = heated_channel_case()
        data['name'] = 'Canal $nome $(1/0)'
        group.add_component('canal', data)
        group.associate('um', 'canal', name='Canal "um" $nome $(1/0)')
        state = CaseRunner().run_lines(group.to_repl().splitlines())
        imported = ModelAssembly.from_repl(group.to_repl())
        self.assertEqual(state.assembly.to_dict(), group.to_dict())
        self.assertEqual(imported.to_dict(), group.to_dict())
        for command in ('let x = 1/0', 'component_load x file=segredo.json', 'export csv=saida.csv'):
            with self.subTest(command=command), self.assertRaises(ValueError):
                ModelAssembly.from_repl(command)

    def test_portable_report_contains_escaped_data_and_all_assets(self):
        group = ModelAssembly(name='Teste </script><div>literal</div>')
        group.add_component('canal', heated_channel_case())
        group.associate('um', 'canal')
        group.solve()
        with tempfile.TemporaryDirectory(dir=ROOT / 'tmp') as directory:
            path = Path(directory) / 'relatorio.html'
            export_report(group, path)
            html = path.read_text(encoding='utf-8')
        tags = []
        class Assets(HTMLParser):
            def handle_starttag(self, tag, attrs):
                tags.append((tag, dict(attrs)))
        Assets().feed(html)
        self.assertFalse(any(tag == 'script' and 'src' in attrs for tag, attrs in tags))
        self.assertFalse(any(tag == 'link' and attrs.get('rel') == 'stylesheet' for tag, attrs in tags))
        payload = re.search(r'<script id="saved-assembly" type="application/json">(.*?)</script>', html, re.S)[1]
        saved = json.loads(payload)
        self.assertNotIn('<', payload)
        self.assertEqual(saved['assembly']['name'], group.name)
        self.assert_same_results(saved['result'], group.result)


if __name__ == '__main__':
    unittest.main()
