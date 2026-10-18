"""Componentes reutilizáveis e associações de modelos térmicos independentes.

Uma instância referencia um componente e escala suas fontes volumétricas.
Materiais, malha, contornos e escoamento continuam definidos no componente.
Não são criadas ligações térmicas ou hidráulicas entre instâncias.
"""
from __future__ import annotations

from copy import deepcopy
import csv
import io
import json
from pathlib import Path
import re
import shlex

from visual.modelo import VisualModel, number


def identifier(value, label='Identificador'):
    if not isinstance(value, str) or not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_-]{0,63}', value):
        raise ValueError(f'{label}: use até 64 letras ASCII, números, _ ou -, começando por letra ou _.')
    return value


class ModelAssembly:
    MAX_COMPONENTS = 50
    MAX_INSTANCES = 100
    MAX_TOTAL_CELLS = 20000

    def __init__(self, data=None, *, name='Associação de componentes'):
        data = deepcopy(data) if data is not None else {'kind': 'assembly', 'name': name, 'components': {}, 'instances': []}
        if not isinstance(data, dict) or data.get('kind') != 'assembly':
            raise ValueError('O arquivo deve conter uma associação com kind="assembly".')
        if data.get('coupling', 'independent') != 'independent':
            raise ValueError('Esta associação calcula componentes independentes. Defina contatos físicos dentro de cada componente.')
        self.name = str(data.get('name', name))[:120]
        self.components = {}
        self.instances = []
        self.result = None
        view = data.get('view', {})
        if not isinstance(view, dict) or view.get('axis', 'y') not in ('x', 'y'):
            raise ValueError('O eixo do mapa deve ser x ou y.')
        self.view = {'axis': view.get('axis', 'y'), 'offset_m': number(view.get('offset_m', 0), 'Origem do eixo'),
                     'region': str(view.get('region', ''))[:80], 'unit': view.get('unit', 'C')}
        if self.view['unit'] not in ('C', 'K'):
            raise ValueError('Unidade do mapa deve ser C ou K.')
        components = data.get('components', {})
        if not isinstance(components, dict) or len(components) > self.MAX_COMPONENTS:
            raise ValueError(f'Use uma biblioteca com até {self.MAX_COMPONENTS} componentes.')
        for key, model in components.items():
            self.add_component(key, model)
        instances = data.get('instances', [])
        if not isinstance(instances, list):
            raise ValueError('As instâncias devem ser uma lista.')
        for instance in instances:
            if not isinstance(instance, dict):
                raise ValueError('Cada instância deve ser um objeto.')
            self.associate(instance.get('id'), instance.get('component'),
                           source_factor=instance.get('source_factor', 1), name=instance.get('name'))

    def add_component(self, key, model, *, replace=False):
        key = identifier(key, 'Componente')
        if key in self.components and not replace:
            raise ValueError(f'O componente {key} já existe.')
        if key not in self.components and len(self.components) >= self.MAX_COMPONENTS:
            raise ValueError('A biblioteca de componentes atingiu o limite.')
        VisualModel(model)  # Mesma validação nos três modos; não executa comandos.
        self.components[key] = deepcopy(model)
        self.result = None

    def associate(self, key, component, *, source_factor=1., name=None):
        key = identifier(key, 'Instância')
        if component not in self.components:
            raise ValueError(f'Componente {component!r} não encontrado. Cadastre ou importe o componente primeiro.')
        if any(instance['id'] == key for instance in self.instances):
            raise ValueError(f'A instância {key} já existe.')
        if len(self.instances) >= self.MAX_INSTANCES:
            raise ValueError(f'Use até {self.MAX_INSTANCES} instâncias por associação.')
        factor = number(source_factor, 'Fator de geração', minimum=0)
        self.instances.append({'id': key, 'name': str(name or key)[:120], 'component': component, 'source_factor': factor})
        self.result = None

    def to_dict(self):
        return deepcopy({'version': 1, 'kind': 'assembly', 'name': self.name, 'coupling': 'independent',
                         'components': self.components, 'instances': self.instances, 'view': self.view})

    def solve(self, *, on_progress=None):
        self.result = None
        if not self.instances:
            raise ValueError('Associe pelo menos uma instância antes de resolver.')
        models = []
        cell_count = 0
        for instance in self.instances:
            data = deepcopy(self.components[instance['component']])
            for region in data['regions']:
                region['source'] = number(region.get('source', 0), 'Geração volumétrica') * instance['source_factor']
            model = VisualModel(data)
            cell_count += model.geom.nx * model.geom.ny
            if cell_count > self.MAX_TOTAL_CELLS:
                raise ValueError(f'A associação excede {self.MAX_TOTAL_CELLS} células. Reduza a malha ou o número de instâncias.')
            models.append(model)
        results = []
        for index, (instance, model) in enumerate(zip(self.instances, models)):
            try:
                result = model.solve()
            except (ValueError, RuntimeError) as error:
                raise ValueError(f'Instância {instance["id"]}: {error}') from error
            results.append({**instance, 'domain': deepcopy(model.data['domain']),
                            'regions': [{'id': r['id'], 'name': r.get('name', r['id'])} for r in model.data['regions']],
                            'result': result})
            if on_progress:
                on_progress(index + 1, len(models), instance['id'])
        self.result = {'kind': 'assembly_result', 'success': True, 'name': self.name, 'coupling': 'independent',
                       'instances': results, 'count': len(results),
                       'maximum_C': max(r['result']['maximum_C'] for r in results),
                       'minimum_C': min(r['result']['minimum_C'] for r in results),
                       'source_W': sum(r['result']['source_W'] for r in results),
                       'max_residual_W': max(r['result']['residual_W'] for r in results),
                       'total_cells': sum(len(r['result']['cells']) for r in results)}
        return self.result

    def to_repl(self):
        # O literal protege nomes com $, aspas ou quebras de linha da interpolação do REPL.
        literal = shlex.quote(json.dumps(self.to_dict(), ensure_ascii=False, separators=(',', ':')))
        return ('# Associação portátil: componentes incorporados; comprimentos em m, temperaturas em °C.\n'
                'assembly_json ' + literal + '\nsolve_assembly\n')

    @classmethod
    def from_repl(cls, text):
        """Lê somente o formato portátil de associação, sem eval nem acesso a arquivos."""
        assembly = cls()
        if not isinstance(text, str) or len(text.encode('utf-8')) > 1048576:
            raise ValueError('Use um arquivo de até 1 MB.')
        for line_number, line in enumerate(text.splitlines(), 1):
            try:
                tokens = shlex.split(line, comments=True)
                if not tokens:
                    continue
                command, *args = tokens
                if command == 'assembly_json' and len(args) == 1:
                    assembly = cls(json.loads(args[0]))
                elif command == 'assembly_reset' and not args:
                    assembly = cls()
                elif command == 'component_json' and len(args) == 2:
                    assembly.add_component(args[0], json.loads(args[1]))
                elif command == 'assembly' and len(args) == 1:
                    assembly.name = args[0][:120]
                elif command == 'assembly_view':
                    options = dict(arg.split('=', 1) for arg in args)
                    assembly.set_view(options)
                elif command == 'associate' and args:
                    options = dict(arg.split('=', 1) for arg in args[1:])
                    if set(options) - {'component', 'source_factor', 'name'}:
                        raise ValueError('Opção de instância desconhecida.')
                    assembly.associate(args[0], options.get('component'), source_factor=options.get('source_factor', 1), name=options.get('name'))
                elif command == 'solve_assembly' and not args:
                    continue  # Importação nunca dispara o cálculo.
                else:
                    raise ValueError('Este importador aceita o .txt portátil exportado pela associação. Arquivos de comandos gerais são executados pelo REPL.')
            except (ValueError, TypeError) as error:
                raise ValueError(f'Linha {line_number}: {error}') from error
        if not assembly.components:
            raise ValueError('O arquivo não contém componentes.')
        return assembly

    def set_view(self, options):
        if set(options) - {'axis', 'offset_m', 'region', 'unit'}:
            raise ValueError('Opção de visualização desconhecida.')
        view = {**self.view, **options}
        if view['axis'] not in ('x', 'y') or view['unit'] not in ('C', 'K'):
            raise ValueError('Revise eixo e unidade do mapa.')
        view['offset_m'] = number(view['offset_m'], 'Origem do eixo')
        view['region'] = str(view['region'])[:80]
        self.view = view

    @staticmethod
    def result_csv(result):
        stream = io.StringIO(newline='')
        writer = csv.writer(stream)
        writer.writerow(['instancia', 'nome', 'componente', 'fator_geracao', 'regiao', 'material',
                         'i', 'j', 'x_m', 'y_m', 'temperatura_C', 'temperatura_K'])
        for instance in result['instances']:
            for cell in instance['result']['cells']:
                writer.writerow([instance['id'], instance['name'], instance['component'], instance['source_factor'],
                                 cell['region'], cell['material'], cell['i'], cell['j'],
                                 (cell['x0'] + cell['x1']) / 2, (cell['y0'] + cell['y1']) / 2,
                                 cell['T_C'], cell['T_C'] + 273.15])
        return '\ufeff' + stream.getvalue()

    def export_csv(self, path):
        if self.result is None:
            raise ValueError('Resolva a associação antes de exportar os resultados.')
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.result_csv(self.result), encoding='utf-8')


def component_from_geometry(geom, *, thickness=1., h_internal=None, pressure=None, name='Componente'):
    """Captura a geometria atual do REPL em um componente do editor.

Recusa construções que o formato geométrico não representa, sem descartá-las.
"""
    if geom is None or geom.nx is None:
        raise ValueError('Defina a geometria e sua malha antes de component_save.')
    if pressure is not None:
        raise ValueError('A pressão global de build_from_geom não é representada neste formato. Para água variável, defina P_in e P_out em flow_region.')
    if geom.interface_rules or any(r.thermal_mode is not None for r in geom.regions):
        raise ValueError('O componente geométrico usa as interações padrão entre as fases dos materiais.')
    materials, regions = [], []
    for index, region in enumerate(geom.regions):
        material = region.material
        model = getattr(material, 'property_model', 'constant')
        if model == 'constant' and any(getattr(material, key + '_func') is not None or getattr(material, key + '_table') is not None for key in ('k', 'rho', 'cp', 'mu')):
            raise ValueError('Para component_save, use propriedades constantes, k polinomial ou correlações registradas. Funções/tabelas personalizadas continuam disponíveis no núcleo.')
        if any(getattr(material, key + '_coeffs') is not None for key in ('rho', 'cp', 'mu')):
            raise ValueError('Este formato de componente aceita polinômio em k e correlações registradas para água.')
        spec = {'id': f'mat{index}', 'name': material.name, 'phase': material.phase.value,
                'k': material.k if material.k is not None else material.prop('k', 30.), 'property_model': model}
        for key in ('rho', 'cp', 'mu'):
            if getattr(material, key) is not None:
                spec[key] = getattr(material, key)
        if material.k_coeffs is not None:
            spec['k_poly'] = list(material.k_coeffs)
        materials.append(spec)
        sources = geom.source_for_region(region.name)
        if any(source.function is not None or source.kind != 'volumetric' for source in sources):
            raise ValueError('O componente geométrico requer fontes volumétricas constantes por região.')
        regions.append({'id': region.name, 'name': region.name, 'material': spec['id'],
                        'x0': region.x0, 'x1': region.x1, 'y0': region.y0, 'y1': region.y1,
                        'source': sum(source.value for source in sources),
                        'conduction_axes': region.metadata.get('conduction_axes', 'xy'),
                        **({'flow': deepcopy(region.metadata['flow'])} if region.metadata.get('flow') else {})})
    boundaries = {side: {'kind': 'adiabatic'} for side in ('left', 'right', 'top', 'bottom')}
    for boundary in geom.boundary_conditions:
        boundaries[boundary.side] = {'kind': boundary.kind, **deepcopy(boundary.data)}
    data = {'version': 1, 'name': name, 'domain': {'width': geom.width, 'height': geom.height, 'thickness': thickness},
            'mesh': {'nx': geom.nx, 'ny': geom.ny, 'align': True, 'h_internal': h_internal},
            'materials': materials, 'regions': regions, 'boundaries': boundaries}
    rebuilt = VisualModel(data)
    import numpy as np
    if any(len(a) != len(b) or not np.allclose(a, b, rtol=0, atol=1e-12)
           for a, b in ((rebuilt.geom.x_edges, geom.x_edges), (rebuilt.geom.y_edges, geom.y_edges))):
        raise ValueError('A malha atual não pode ser representada por este componente; use fix=interfaces antes de salvar.')
    return data
