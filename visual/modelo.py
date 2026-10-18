"""Tradução validada do editor de retângulos para Geometry2D/NodalNetwork."""

from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np
from scipy.optimize import brentq

from geometria import Geometry2D
from escoamento import property_material, build_with_flows
from solver import solve_steady_state


def number(value, label, *, minimum=None, positive=False):
    if isinstance(value, bool):
        raise ValueError(f"{label}: informe um número.")
    try:
        value = float(value)
    except (ValueError, TypeError) as error:
        raise ValueError(f"{label}: informe um número.") from error
    if not math.isfinite(value) or (positive and value <= 0) or (minimum is not None and value < minimum):
        raise ValueError(f"{label}: valor inválido.")
    return value


@dataclass(frozen=True)
class MaterialSpec:
    id: str
    name: str
    phase: str
    k: float
    rho: float | None = None
    cp: float | None = None
    mu: float | None = None
    k_poly: tuple | None = None
    property_model: str = "constant"

    @classmethod
    def from_dict(cls, data):
        name = str(data.get("name", "")).strip()
        identifier = str(data.get("id", "")).strip()
        if not name or not identifier or len(name) > 80 or len(identifier) > 80:
            raise ValueError("Cada material precisa de nome e identificador (até 80 caracteres).")
        phase = data.get("phase", "solid")
        if phase not in ("solid", "fluid"):
            raise ValueError("Fase do material deve ser solid ou fluid.")
        optional = {key: number(data[key], key, positive=True) if data.get(key) not in (None, "") else None
                    for key in ("rho", "cp", "mu")}
        coefficients = data.get("k_poly")
        if coefficients is not None:
            if not isinstance(coefficients, list) or not 1 <= len(coefficients) <= 6:
                raise ValueError("k(T) aceita de 1 a 6 coeficientes, do termo constante em diante.")
            coefficients = tuple(number(c, "coeficiente de k(T)") for c in coefficients)
        property_model = data.get("property_model", "constant")
        if property_model not in ("constant", "u3si2_al", "iapws97"):
            raise ValueError("Modelo de propriedades desconhecido.")
        return cls(identifier, name, phase, number(data.get("k"), "k", positive=True),
                   **optional, k_poly=coefficients, property_model=property_model)

    def material(self):
        return property_material(self.name, self.phase, self.property_model, k=self.k,
                                 rho=self.rho, cp=self.cp, mu=self.mu, k_coeffs=self.k_poly)


class VisualModel:
    MAX_CELLS = 900
    MAX_UNKNOWNS = 1600

    def __init__(self, data):
        if not isinstance(data, dict):
            raise ValueError("O modelo deve ser um objeto JSON.")
        if any(not isinstance(data.get(key, {}), dict) for key in ("domain", "mesh", "boundaries")):
            raise ValueError("Domínio, malha e contornos devem ser objetos JSON.")
        if not isinstance(data.get("materials", []), list) or len(data.get("materials", [])) > 100:
            raise ValueError("A biblioteca deve conter uma lista de até 100 materiais.")
        if any(not isinstance(raw, dict) for raw in data.get("materials", [])):
            raise ValueError("Cada material deve ser um objeto JSON.")
        self.data = data
        domain = data.get("domain", {})
        self.width = number(domain.get("width"), "Largura do domínio", positive=True)
        self.height = number(domain.get("height"), "Altura do domínio", positive=True)
        self.thickness = number(domain.get("thickness"), "Espessura em z", positive=True)
        self.geom = Geometry2D(width=self.width, height=self.height)
        self.materials = {}
        for raw in data.get("materials", []):
            spec = MaterialSpec.from_dict(raw)
            if spec.id in self.materials:
                raise ValueError("Identificador de material repetido.")
            self.materials[spec.id] = spec
        regions = data.get("regions", [])
        if not isinstance(regions, list) or not 1 <= len(regions) <= 60:
            raise ValueError("Desenhe de 1 a 60 regiões para gerar a malha.")
        if any(not isinstance(raw, dict) for raw in regions):
            raise ValueError("Cada região deve ser um objeto JSON.")
        self.regions = {}
        for raw in regions:
            identifier = str(raw.get("id", "")).strip()
            if not identifier or identifier in self.regions:
                raise ValueError("Identificador de região ausente ou repetido.")
            material_id = raw.get("material")
            if material_id not in self.materials:
                raise ValueError("Escolha um material para cada retângulo.")
            coordinates = {key: number(raw.get(key), key, minimum=0) for key in ("x0", "x1", "y0", "y1")}
            self.geom.material(identifier, self.materials[material_id].material(), **coordinates,
                               conduction_axes=raw.get("conduction_axes", "xy"), flow=raw.get("flow"))
            source = number(raw.get("source", 0), "Geração volumétrica")
            if source:
                self.geom.source(identifier, "volumetric", value=source)
            self.regions[identifier] = raw
        mesh = data.get("mesh", {})
        nx = number(mesh.get("nx", 12), "nx", positive=True)
        ny = number(mesh.get("ny", 8), "ny", positive=True)
        if int(nx) != nx or int(ny) != ny or nx * ny > self.MAX_CELLS:
            raise ValueError(f"Use nx e ny inteiros, com até {self.MAX_CELLS} células no domínio.")
        self.geom.mesh(int(nx), int(ny), fix="interfaces" if mesh.get("align", True) else "raise")
        if self.geom.nx * self.geom.ny > self.MAX_CELLS:
            raise ValueError("O alinhamento gerou uma malha grande. Encaixe as bordas na grade ou reduza nx/ny.")
        for side in ("left", "right", "top", "bottom"):
            boundary = data.get("boundaries", {}).get(side, {"kind": "adiabatic"})
            if not isinstance(boundary, dict):
                raise ValueError("Cada condição de contorno deve ser um objeto JSON.")
            kind = boundary.get("kind", "adiabatic")
            if kind == "temperature":
                self.geom.bc(side, kind, T=number(boundary.get("T"), "Temperatura de contorno", minimum=-273.15))
            elif kind in ("convection", "radiation", "convection_radiation"):
                values = {}
                if kind in ("convection", "convection_radiation"):
                    values.update(h=number(boundary.get("h"), "h", positive=True),
                                  T_inf=number(boundary.get("T_inf"), "Temperatura do fluido", minimum=-273.15))
                if kind in ("radiation", "convection_radiation"):
                    epsilon = number(boundary.get("emissivity"), "Emissividade", minimum=0)
                    factor = number(boundary.get("view_factor"), "Fator de vista", minimum=0)
                    if epsilon > 1 or factor > 1:
                        raise ValueError("Emissividade e fator de vista devem estar entre 0 e 1.")
                    values.update(emissivity=epsilon, view_factor=factor,
                                  T_env=number(boundary.get("T_env"), "Temperatura radiativa", minimum=-273.15))
                self.geom.bc(side, kind, **values)
            elif kind == "heat_flux":
                self.geom.bc(side, kind, q=number(boundary.get("q"), "Fluxo de calor"))
            elif kind == "adiabatic":
                self.geom.bc(side, kind)
            else:
                raise ValueError("Tipo de condição de contorno desconhecido.")

    def mesh_description(self):
        grid = self.geom.region_grid()
        cells = []
        contacts = []
        for j in range(self.geom.ny):
            for i in range(self.geom.nx):
                identifier = grid[j, i]
                if identifier == "VOID":
                    continue
                x0, x1, y0, y1 = self.geom.cell_bounds(i, j)
                cells.append({"i": i, "j": j, "region": str(identifier),
                              "x0": x0, "x1": x1, "y0": y0, "y1": y1})
                for di, dj in ((1, 0), (0, 1)):
                    ni, nj = i + di, j + dj
                    if ni < self.geom.nx and nj < self.geom.ny:
                        neighbor = grid[nj, ni]
                        if neighbor != "VOID" and neighbor != identifier:
                            contacts.append({"a": str(identifier), "b": str(neighbor),
                                             "x0": x1 if di else x0, "y0": y0 if di else y1,
                                             "x1": x1, "y1": y1})
        if not cells:
            raise ValueError("A malha não contém células com material.")
        return {"nx": self.geom.nx, "ny": self.geom.ny,
                "x_edges": self.geom.x_edges.tolist(), "y_edges": self.geom.y_edges.tolist(),
                "cells": cells, "contacts": contacts,
                "void_cells": self.geom.nx * self.geom.ny - len(cells)}

    @staticmethod
    def _check_anchors(net):
        adjacency = {node: [] for node in net.nodes}
        temperatures = net.temperature_map_from_vector(net.initial_guess())
        for link in net.links:
            if link.G(temperatures) > 0:
                adjacency[link.node_i].append(link.node_j)
                adjacency[link.node_j].append(link.node_i)
        remaining = set(net.nodes)
        while remaining:
            stack = [remaining.pop()]
            component = []
            while stack:
                node = stack.pop()
                component.append(node)
                for neighbor in adjacency[node]:
                    if neighbor in remaining:
                        remaining.remove(neighbor)
                        stack.append(neighbor)
            if any(net.nodes[n].is_unknown() for n in component) and not any(net.nodes[n].is_boundary() for n in component):
                region = next((getattr(net.nodes[n].region, "name", "região") for n in component
                               if net.nodes[n].is_unknown()), "região")
                raise ValueError(f"A região {region} está sem referência térmica. Defina temperatura, convecção ou entrada de fluido conectada a ela.")

    def solve(self):
        description = self.mesh_description()
        mesh_options = self.data.get("mesh", {})
        h = mesh_options.get("h_internal")
        h = number(h, "h interno", positive=True) if h not in (None, "") else None
        net = build_with_flows(self.geom, thickness_z=self.thickness, internal_convection_h=h)
        if len(net.unknown_node_ids()) > self.MAX_UNKNOWNS:
            raise ValueError("A rede tem muitos nós de interface. Reduza a malha.")
        self._check_anchors(net)
        result = solve_steady_state(net, tol=1e-7, max_iter=100)
        if not result.success:
            raise ValueError(f"A solução não convergiu. Resíduo {result.residual_norm:.3g} W. Revise materiais e contornos.")
        self.network, self.solution = net, result
        temperatures = []
        for cell in description["cells"]:
            node = net.nodes[net.cell_to_node[(cell["j"], cell["i"])]]
            cell["T_C"] = node.temperature
            cell["material"] = node.material.name
            temperatures.append(node.temperature)
        description.update({"success": True, "minimum_C": min(temperatures), "maximum_C": max(temperatures),
                            "residual_W": result.residual_norm, "nodes": len(net.nodes), "links": len(net.links),
                            "source_W": sum(node.source for node in net.nodes.values()),
                            "iterations": result.iterations})
        probes = self.data.get("probes", [])
        if not isinstance(probes, list) or len(probes) > 20:
            raise ValueError("Use uma lista de até 20 pontos de avaliação.")
        description["probes"] = [self.boundary_probe(description, probe) for probe in probes]
        description["region_results"] = []
        temperatures_map = net.temperatures_dict()
        for identifier, region in self.regions.items():
            cells = [c for c in description["cells"] if c["region"] == identifier]
            summary = {"id": identifier, "name": region.get("name", identifier),
                       "minimum_C": min(c["T_C"] for c in cells), "maximum_C": max(c["T_C"] for c in cells)}
            flow = net.region_flows.get(identifier)
            if flow:
                r = self.geom.get_region(identifier)
                direction = flow.direction
                axis, edge = ("x", r.x1 if direction == "x+" else r.x0) if direction.startswith("x") else ("y", r.y1 if direction == "y+" else r.y0)
                end_key = axis + ("1" if direction.endswith("+") else "0")
                outlet = [c for c in cells if abs(c[end_key] - edge) < 1e-12]
                weights = [(c["y1"] - c["y0"]) if axis == "x" else (c["x1"] - c["x0"]) for c in outlet]
                summary["outlet_C"] = sum(c["T_C"] * w for c, w in zip(outlet, weights)) / sum(weights)
                if flow.h_model == "dittus_boelter":
                    hs = [flow.film(temperatures_map, net.nodes[net.cell_to_node[(c["j"], c["i"])]]) for c in cells]
                    summary.update(h_min=min(hs), h_max=max(hs), sampling="Temperatura do fluido na saída de cada intervalo axial.")
            description["region_results"].append(summary)
        description["surfaces"] = [{"type": n.metadata["type"], "x": n.x, "y": n.y, "T_C": n.temperature}
                                   for n in net.nodes.values() if n.metadata.get("type") in ("conductive_interface", "convective_surface")]
        return description

    def to_repl(self):
        """Exporta as entradas físicas; os três modos usam a mesma montagem."""
        lines = ["# Caso exportado pelo editor visual. Unidades SI; temperatura em °C.",
                 f"domain W={self.width!r} H={self.height!r}"]
        materials = {}
        for index, spec in enumerate(self.materials.values()):
            alias = f"mat{index}"
            materials[spec.id] = alias
            values = [f"phase={spec.phase}", f"k={spec.k!r}", f"property_model={spec.property_model}"]
            for key in ("rho", "cp", "mu"):
                if getattr(spec, key) is not None:
                    values.append(f"{key}={getattr(spec, key)!r}")
            if spec.k_poly:
                values.append("k_poly=" + ",".join(repr(c) for c in spec.k_poly))
            lines.append(f"material {alias} " + " ".join(values))
        for index, region in enumerate(self.regions.values()):
            alias = f"reg{index}"
            values = " ".join(f"{key}={float(region[key])!r}" for key in ("x0", "x1", "y0", "y1"))
            lines.append(f"region {alias} mat={materials[region['material']]} {values} conduction_axes={region.get('conduction_axes', 'xy')}")
            lines.append(f"source {alias} q={float(region.get('source', 0))!r}")
            if self.materials[region["material"]].phase == "fluid":
                flow = region.get("flow")
                if not isinstance(flow, dict):
                    raise ValueError("Defina o escoamento antes de exportar o caso.")
                # Somente parâmetros conhecidos; nomes do usuário não viram comandos.
                values = []
                for key in ("direction", "model", "h_model", "T_in", "mdot", "cp", "vdot", "P_in", "P_out", "hydraulic_diameter"):
                    if flow.get(key) is None:
                        continue
                    value = flow[key]
                    if key in ("direction", "model", "h_model"):
                        allowed = {"direction": ("x+", "x-", "y+", "y-"), "model": ("constant", "iapws97"), "h_model": ("constant", "dittus_boelter")}
                        if value not in allowed[key]:
                            raise ValueError("Opção de escoamento inválida.")
                    else:
                        value = repr(number(value, key))
                    values.append(f"{key}={value}")
                lines.append(f"flow_region {alias} " + " ".join(values))
        for side in ("left", "right", "top", "bottom"):
            boundary = self.data.get("boundaries", {}).get(side, {"kind": "adiabatic"})
            values = []
            for key in ("T", "h", "T_inf", "emissivity", "view_factor", "T_env", "q"):
                if boundary.get(key) is not None:
                    values.append(f"{key}={number(boundary[key], key)!r}")
            lines.append(f"bc {side} {boundary['kind']} " + " ".join(values))
        lines.append(f"mesh nx={self.geom.nx} ny={self.geom.ny} fix={'interfaces' if self.data['mesh'].get('align', True) else 'raise'}")
        h = self.data['mesh'].get('h_internal')
        lines.append(f"build_from_geom tz={self.thickness!r}" + (f" h_int={number(h, 'h', positive=True)!r}" if h not in (None, '') else ''))
        lines.append("solve tol=1e-7 max_iter=100 linhas=0")
        return "\n".join(lines) + "\n"

    def boundary_probe(self, result, probe):
        """Temperatura de superfície, interpolada entre centros das faces.

        O valor no centro da célula não é a temperatura da superfície.
        Usa a mesma resistência de meia célula + convecção montada no núcleo.
        """
        if not isinstance(probe, dict):
            raise ValueError("Ponto de avaliação inválido.")
        x, y = number(probe.get("x"), "x do ponto"), number(probe.get("y"), "y do ponto")
        if not 0 <= x <= self.width or not 0 <= y <= self.height:
            raise ValueError("O ponto deve ficar no domínio.")
        sides = [side for side, on_edge in (("left", abs(x) < 1e-12),
                 ("right", abs(x - self.width) < 1e-12), ("bottom", abs(y) < 1e-12),
                 ("top", abs(y - self.height) < 1e-12)) if on_edge]
        if len(sides) != 1:
            raise ValueError("Nesta versão, o ponto deve estar em uma face externa, fora dos cantos.")
        side = sides[0]
        boundary = self.data["boundaries"].get(side, {"kind": "adiabatic"})
        if boundary["kind"] == "temperature":
            temperature = number(boundary["T"], "Temperatura")
        elif boundary["kind"] in ("convection", "radiation", "convection_radiation"):
            samples = []
            vertical = side in ("left", "right")
            for cell in result["cells"]:
                on_edge = {"left": cell["i"] == 0, "right": cell["i"] == result["nx"] - 1,
                           "bottom": cell["j"] == 0, "top": cell["j"] == result["ny"] - 1}[side]
                if not on_edge:
                    continue
                material_id = self.regions[cell["region"]]["material"]
                k = self.materials[material_id].material().prop("k", cell["T_C"])
                distance = (cell["x1"] - cell["x0"] if vertical else cell["y1"] - cell["y0"]) / 2
                has_convection = boundary["kind"] in ("convection", "convection_radiation")
                has_radiation = boundary["kind"] in ("radiation", "convection_radiation")
                h = boundary["h"] if has_convection else 0
                ambient = boundary["T_inf"] if has_convection else cell["T_C"]
                if not has_radiation:
                    surface_T = (k / distance * cell["T_C"] + h * ambient) / (k / distance + h)
                else:
                    environment = boundary["T_env"]
                    coefficient = boundary["emissivity"] * boundary["view_factor"] * 5.670374419e-8
                    def balance(T):
                        return (k / distance * (cell["T_C"] - T) - h * (T - ambient) -
                                coefficient * ((T + 273.15) ** 4 - (environment + 273.15) ** 4))
                    surface_T = brentq(balance, min(cell["T_C"], ambient, environment),
                                       max(cell["T_C"], ambient, environment))
                coordinate = (cell["y0"] + cell["y1"]) / 2 if vertical else (cell["x0"] + cell["x1"]) / 2
                samples.append((coordinate, surface_T))
            samples.sort()
            target = y if vertical else x
            if not samples or not samples[0][0] <= target <= samples[-1][0]:
                raise ValueError("O ponto não está entre centros das faces com material; refine a malha.")
            coordinates, values = zip(*samples)
            temperature = float(np.interp(target, coordinates, values))
        else:
            raise ValueError("A avaliação na superfície exige temperatura, convecção ou radiação.")
        return {"name": str(probe.get("name", "Ponto")), "x": x, "y": y, "T_C": temperature,
                "method": "surface_temperature_with_linear_face_interpolation"}


def demo_cases():
    """Casos didáticos: valores de teste explícitos, sem dados de projeto."""
    boundaries = {side: {"kind": "adiabatic"} for side in ("left", "right", "top", "bottom")}
    boundaries["left"] = {"kind": "temperature", "T": 100}
    boundaries["right"] = {"kind": "temperature", "T": 20}
    common = {"version": 1, "name": "Parede composta", "domain": {"width": .06, "height": .04, "thickness": .001},
              "mesh": {"nx": 12, "ny": 8, "align": True, "h_internal": None}, "boundaries": boundaries,
              "materials": [{"id": "a", "name": "Material A", "phase": "solid", "k": 10, "color": "#7190a7"},
                            {"id": "b", "name": "Material B", "phase": "solid", "k": 20, "color": "#60a897"}],
              "regions": [{"id": "r1", "name": "Parede A", "material": "a", "x0": 0, "x1": .03, "y0": 0, "y1": .04, "source": 0},
                          {"id": "r2", "name": "Parede B", "material": "b", "x0": .03, "x1": .06, "y0": 0, "y1": .04, "source": 0}],
              "benchmark": {"kind": "composite", "source": "Bastos & Andrade, cap. I, condução em placas justapostas; valores didáticos de teste.",
                            "description": "Solução exata por resistências em série. T da interface = 46,6667 °C."}}
    import copy
    uniform = copy.deepcopy(common)
    uniform.update(name="Parede uniforme")
    uniform["regions"] = [{"id": "r1", "name": "Parede", "material": "a", "x0": 0, "x1": .06, "y0": 0, "y1": .04, "source": 0}]
    uniform["benchmark"] = {"kind": "uniform", "source": "Bastos & Andrade, cap. II: condução permanente em parede plana; valores didáticos de teste.",
                             "description": "Solução exata linear: T(x) = 100 − 80 x / 0,06. O campo nos centros das células deve coincidir."}
    chip = {"version": 1, "name": "Chip do livro · geometria 2D", "domain": {"width": .03, "height": .03, "thickness": .03},
            "mesh": {"nx": 12, "ny": 12, "align": True, "h_internal": None},
            "materials": [{"id": "legs", "name": "Pernas de cobre (equivalente)", "phase": "solid", "k": 3.979, "color": "#be9873"},
                          {"id": "silicon", "name": "Silício", "phase": "solid", "k": 150, "color": "#7190a7"}],
            "regions": [{"id": "r1", "name": "Pernas equivalentes", "material": "legs", "x0": 0, "x1": .03, "y0": 0, "y1": .02, "source": 0},
                        {"id": "r2", "name": "Componente", "material": "silicon", "x0": 0, "x1": .03, "y0": .02, "y1": .03, "source": 4.444e5}],
            "boundaries": {"bottom": {"kind": "temperature", "T": 44}, "left": {"kind": "adiabatic"}, "right": {"kind": "adiabatic"}, "top": {"kind": "adiabatic"}},
            "benchmark": {"kind": "chip2d", "source": "Adaptação 2D do exemplo existente interativo_chip.txt; livro, exercício 1(a), pp. 58-60.",
                          "description": "A rede geométrica distribuída resulta em aproximadamente 66,4854 °C. O modelo concentrado do livro resulta em 66,3375 °C; são discretizações diferentes."}}
    return {"composite": common, "uniform": uniform, "chip": chip, "nafems_t4": nafems_t4_case()}


def nafems_t4_case(nx=18, ny=30):
    """Dados públicos da documentação COMSOL do benchmark NAFEMS T4."""
    return {"version": 1, "name": "NAFEMS T4 · placa com convecção",
            "domain": {"width": .6, "height": 1., "thickness": 1.},
            "materials": [{"id": "nafems", "name": "Condutor NAFEMS · k = 52", "phase": "solid", "k": 52, "color": "#7190a7"}],
            "regions": [{"id": "plate", "name": "Placa", "material": "nafems", "x0": 0, "x1": .6, "y0": 0, "y1": 1., "source": 0}],
            "mesh": {"nx": nx, "ny": ny, "align": True, "h_internal": None},
            "boundaries": {"left": {"kind": "adiabatic"}, "bottom": {"kind": "temperature", "T": 100},
                           "top": {"kind": "convection", "h": 750, "T_inf": 0},
                           "right": {"kind": "convection", "h": 750, "T_inf": 0}},
            "probes": [{"name": "E", "x": .6, "y": .2}],
            "benchmark": {"kind": "nafems_t4", "target_C": 18.25,
                          "source": "https://www.comsol.com/model/download/1329071/models.mph.heat_convection_2d.pdf",
                          "description": "NAFEMS T4: temperatura de referência de 18,25 °C na superfície, em E = (600, 200) mm. O ponto usa temperatura de superfície interpolada entre faces; refine a malha para conferir a convergência."}}
