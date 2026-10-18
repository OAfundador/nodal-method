"""Leis de propriedades e escoamento usadas pelo desenho e pelo REPL.

Unidades: °C, Pa, m³/s, kg/s e SI. O modelo IAPWS usa água líquida,
vazão volumétrica constante e cp avaliado a montante de cada intervalo.
Trata-se da marcha em temperatura; não substitui uma formulação em entalpia.
"""
from dataclasses import dataclass
from functools import lru_cache
import importlib
import math
from pathlib import Path
import sys

from materiais import Material
from nos import (NodeKind, TransferKind, LinkDirection, build_network_from_geometry,
                 fluid_node_chains_by_region)


@lru_cache(maxsize=1)
def _iapws_class():
    try:
        return importlib.import_module('iapws').IAPWS97
    except ModuleNotFoundError as error:
        if error.name != 'iapws':
            raise
        local = Path(__file__).resolve().parent / '.deps'
        if local.is_dir() and str(local) not in sys.path:
            sys.path.insert(0, str(local))
            try:
                return importlib.import_module('iapws').IAPWS97
            except ModuleNotFoundError as local_error:
                if local_error.name != 'iapws':
                    raise
        raise RuntimeError('Para água variável, instale a dependência opcional: python -m pip install -r requirements-agua.txt') from error


@dataclass(frozen=True)
class WaterProperties:
    rho: float
    cp: float
    mu: float
    k: float


@lru_cache(maxsize=16384)
def water_properties(T_C, P_Pa):
    if P_Pa is None:
        raise ValueError('IAPWS97 exige pressão em Pa.')
    water = _iapws_class()(T=float(T_C) + 273.15, P=float(P_Pa) / 1e6)
    if water.x != 0:
        raise ValueError('Este modelo de escoamento exige água líquida monofásica.')
    return WaterProperties(water.rho, water.cp * 1000, water.mu, water.k)


def property_material(name, phase, property_model='constant', **properties):
    if property_model == 'u3si2_al':
        if phase != 'solid':
            raise ValueError('A correlação U3Si2-Al exige fase sólida.')
        # Mesma correlação já disponível na biblioteca do REPL; T em °C.
        properties['k_func'] = lambda T, P=None: 1.73073 * (
            3978.1 / (724.61 + 1.8 * T) + 6.02366e-12 * (1.8 * T + 492) ** 3)
    elif property_model == 'iapws97':
        if phase != 'fluid':
            raise ValueError('IAPWS97 exige fase fluida.')
        for key in ('k', 'rho', 'cp', 'mu'):
            properties[key + '_func'] = lambda T, P, key=key: getattr(water_properties(T, P), key)
    elif property_model != 'constant':
        raise ValueError('Modelo de propriedades desconhecido.')
    return Material(name=name, phase=phase, **properties)


def _positive(value, name):
    if isinstance(value, bool):
        raise ValueError(f'{name}: informe um número positivo.')
    try:
        value = float(value)
    except (ValueError, TypeError) as error:
        raise ValueError(f'{name}: informe um número positivo.') from error
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f'{name}: informe um número positivo.')
    return value


class RegionFlow:
    """Uma vazão por região, dividida proporcionalmente entre as cadeias."""
    def __init__(self, region, thickness, data):
        if not isinstance(data, dict) or not data:
            raise ValueError(f'Fluido em {region.name}: informe vazão, direção e temperatura de entrada.')
        self.region, self.thickness = region, thickness
        self.direction = data.get('direction', 'y-')
        if self.direction not in ('x+', 'x-', 'y+', 'y-'):
            raise ValueError('Direção do fluido deve ser x+, x-, y+ ou y-.')
        try:
            self.T_in = float(data['T_in'])
        except (KeyError, ValueError, TypeError) as error:
            raise ValueError('Informe a temperatura de entrada do fluido.') from error
        if not math.isfinite(self.T_in) or self.T_in <= -273.15:
            raise ValueError('Temperatura de entrada inválida.')
        self.model = data.get('model', 'constant')
        self.h_model = data.get('h_model', 'constant')
        if self.h_model not in ('constant', 'dittus_boelter'):
            raise ValueError('Modelo de h desconhecido.')
        self.axial_x = self.direction.startswith('x')
        self.cross_width = (region.y1 - region.y0) if self.axial_x else (region.x1 - region.x0)
        self.area = self.cross_width * thickness
        if self.model == 'iapws97':
            self.vdot = _positive(data.get('vdot'), 'Vazão volumétrica')
            self.P_in = _positive(data.get('P_in'), 'Pressão de entrada')
            self.P_out = _positive(data.get('P_out'), 'Pressão de saída')
            water_properties(self.T_in, self.P_in)  # Valida também a dependência opcional.
        elif self.model == 'constant':
            self.mdot = _positive(data.get('mdot'), 'Vazão mássica')
            self.cp = _positive(data.get('cp', region.material.cp), 'cp do fluido')
            if self.h_model != 'constant':
                raise ValueError('Dittus–Boelter neste editor exige o modelo de água IAPWS97.')
        else:
            raise ValueError('Modelo de escoamento desconhecido.')
        if self.h_model == 'dittus_boelter':
            self.diameter = _positive(data.get('hydraulic_diameter'), 'Diâmetro hidráulico')

    def pressure(self, node=None):
        if node is None or node.bounds is None:
            fraction = 0.
        else:
            x0, x1, y0, y1 = node.bounds
            r = self.region
            fraction = {'x+': (x1 - r.x0) / (r.x1 - r.x0),
                        'x-': (r.x1 - x0) / (r.x1 - r.x0),
                        'y+': (y1 - r.y0) / (r.y1 - r.y0),
                        'y-': (r.y1 - y0) / (r.y1 - r.y0)}[self.direction]
        return self.P_in + fraction * (self.P_out - self.P_in)

    def film(self, temperatures, fluid_node):
        if self.h_model != 'dittus_boelter':
            raise ValueError('Informe h interno ou selecione uma correlação para este fluido.')
        water = water_properties(temperatures[fluid_node.id], self.pressure(fluid_node))
        reynolds = water.rho * self.vdot / self.area * self.diameter / water.mu
        prandtl = water.cp * water.mu / water.k
        return .023 * water.k / self.diameter * reynolds ** .8 * prandtl ** .4

    def connect(self, net):
        chains = fluid_node_chains_by_region(net, self.region.name, self.direction)
        for chain_index, chain in enumerate(chains):
            first = net.nodes[chain[0]]
            x0, x1, y0, y1 = first.bounds
            weight = ((y1 - y0) if self.axial_x else (x1 - x0)) / self.cross_width
            inlet = net.add_node(name=f'inlet_{self.region.name}_{chain_index}', kind=NodeKind.BOUNDARY,
                                 x=first.x, y=first.y, fixed_temperature=self.T_in,
                                 temperature=self.T_in, material=first.material)
            for up, down in zip([inlet] + chain[:-1], chain):
                if self.model == 'constant':
                    options = {'conductance': weight * self.mdot * self.cp}
                else:
                    def conductance(T, upstream=up, weight=weight):
                        water = water_properties(T[upstream], self.pressure(net.nodes[upstream]))
                        return weight * water.rho * self.vdot * water.cp
                    options = {'conductance_func': conductance}
                net.add_link(up, down, TransferKind.FLUID_TRANSPORT, direction=LinkDirection.I_TO_J,
                             name=f'flow_{up}_{down}', metadata={'region_flow': self.region.name}, **options)
        return chains


def build_with_flows(geom, *, thickness_z=1., internal_convection_h=None, pressure=None):
    flows = {r.name: RegionFlow(r, thickness_z, r.metadata.get('flow')) for r in geom.regions if r.material.is_fluid()}
    def film(T, node):
        flow = flows[node.region.name]
        if flow.h_model == 'dittus_boelter':
            return flow.film(T, node)
        if internal_convection_h is None:
            raise ValueError('Informe h entre sólido e fluido ou selecione Dittus–Boelter na região fluida.')
        return _positive(internal_convection_h, 'h interno')
    net = build_network_from_geometry(geom, thickness_z=thickness_z,
                                      internal_convection_h=film if flows else internal_convection_h,
                                      pressure=pressure)
    for flow in flows.values():
        flow.connect(net)
    net.region_flows = flows
    return net
