"""
nos.py

Rede nodal térmica construída a partir de Geometry2D.

Esta versão faz a ponte:

    Geometry2D  ->  NodalNetwork

Ela cria:
- um nó de célula para cada célula física da malha;
- fonte integrada Q_i em cada nó;
- ligações internas por condução;
- nós aritméticos em interfaces entre materiais diferentes;
- nós aritméticos de superfície em convecção/radiação;
- nós de contorno para temperatura prescrita e ambiente;
- residual de regime permanente.

Convenções:
- Temperatura em °C.
- Comprimentos em m.
- Área em m².
- Volume em m³.
- Fonte em W.
- Condutância em W/K.
- Fluxo positivo = calor entrando no nó balanceado.

Regra de sinal:
    fluxo entrando no nó atual = G * (T_vizinho - T_atual)

Regime permanente:
    R_i = Q_i + soma(fluxos entrando em i) = 0
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Optional
import math
import numpy as np

from condutancias import (
    conduction_G,
    convection_G,
    fluid_transport_G,
    radiation_linearized_G,
    equivalent_series_G,
)


class NodeKind(Enum):
    DIFFUSION = "diffusion"
    FLUID = "fluid"
    ARITHMETIC = "arithmetic"
    BOUNDARY = "boundary"


class TransferKind(Enum):
    CONDUCTION = "conduction"
    CONVECTION = "convection"
    FLUID_TRANSPORT = "fluid_transport"
    RADIATION = "radiation"
    EQUIVALENT = "equivalent"


class LinkDirection(Enum):
    UNDIRECTED = "undirected"
    I_TO_J = "i_to_j"
    J_TO_I = "j_to_i"


TemperatureMap = dict[int, float]
ConductanceFunction = Callable[[TemperatureMap], float]


def material_name(material: Any) -> str:
    return getattr(material, "name", str(material))


def region_name(region: Any) -> str:
    return getattr(region, "name", str(region))


def is_fluid_material(material: Any) -> bool:
    if hasattr(material, "is_fluid"):
        return bool(material.is_fluid())
    phase = getattr(material, "phase", None)
    if phase is None:
        return False
    return str(getattr(phase, "value", phase)).lower() in ("fluid", "fluido")


def conductivity(material: Any, T: float, P: Optional[float] = None) -> float:
    if material is None:
        raise ValueError("Material ausente.")

    if hasattr(material, "conductivity"):
        value = material.conductivity(T, P)
    elif hasattr(material, "k"):
        value = material.k
    else:
        raise ValueError(f"Material {material!r} não possui condutividade.")

    value = float(value)

    if value <= 0.0 or not math.isfinite(value):
        raise ValueError(f"Condutividade inválida para {material_name(material)}: {value}")

    return value


def density(material: Any, T: float, P: Optional[float] = None) -> float:
    if hasattr(material, "density"):
        return float(material.density(T, P))
    if hasattr(material, "rho") and material.rho is not None:
        return float(material.rho)
    raise ValueError(f"Material {material_name(material)} não possui rho.")


def specific_heat(material: Any, T: float, P: Optional[float] = None) -> float:
    if hasattr(material, "specific_heat"):
        return float(material.specific_heat(T, P))
    if hasattr(material, "cp") and material.cp is not None:
        return float(material.cp)
    raise ValueError(f"Material {material_name(material)} não possui cp.")


def conduction_func_single_material(
    material: Any,
    node_temperature_id: int,
    area: float,
    distance: float,
    pressure: Optional[float] = None,
) -> ConductanceFunction:
    def func(T: TemperatureMap) -> float:
        k = conductivity(material, T[node_temperature_id], pressure)
        return conduction_G(k, area, distance)

    return func


def conduction_func_between_nodes(
    material: Any,
    node_i: int,
    node_j: int,
    area: float,
    distance: float,
    pressure: Optional[float] = None,
) -> ConductanceFunction:
    def func(T: TemperatureMap) -> float:
        Tm = 0.5 * (T[node_i] + T[node_j])
        k = conductivity(material, Tm, pressure)
        return conduction_G(k, area, distance)

    return func


def radiation_func_between_nodes(
    node_i: int,
    node_j: int,
    emissivity: float,
    area: float,
    view_factor: float = 1.0,
) -> ConductanceFunction:
    def func(T: TemperatureMap) -> float:
        return radiation_linearized_G(
            emissivity=emissivity,
            area=area,
            T_i_C=T[node_i],
            T_j_C=T[node_j],
            view_factor=view_factor,
        )

    return func


@dataclass
class ThermalNode:
    id: int
    name: str
    kind: NodeKind

    x: float
    y: float
    z: float = 0.0

    material: Optional[Any] = None
    region: Optional[Any] = None

    volume: float = 0.0
    source: float = 0.0

    temperature: float = 30.0
    fixed_temperature: Optional[float] = None

    i: Optional[int] = None
    j: Optional[int] = None
    bounds: Optional[tuple[float, float, float, float]] = None

    metadata: dict[str, Any] = field(default_factory=dict)

    def is_unknown(self) -> bool:
        return self.kind != NodeKind.BOUNDARY

    def is_boundary(self) -> bool:
        return self.kind == NodeKind.BOUNDARY

    def is_arithmetic(self) -> bool:
        return self.kind == NodeKind.ARITHMETIC

    def has_capacitance(self) -> bool:
        return self.kind in (NodeKind.DIFFUSION, NodeKind.FLUID) and self.volume > 0.0

    def add_source(self, Q: float) -> None:
        self.source += float(Q)

    def heat_capacity(self, pressure: Optional[float] = None) -> float:
        if not self.has_capacitance():
            return 0.0

        if self.material is None:
            raise ValueError(f"Nó {self.name!r} não possui material.")

        rho = density(self.material, self.temperature, pressure)
        cp = specific_heat(self.material, self.temperature, pressure)
        return rho * cp * self.volume


@dataclass
class ThermalLink:
    node_i: int
    node_j: int
    kind: TransferKind
    direction: LinkDirection = LinkDirection.UNDIRECTED

    name: str = ""
    conductance: Optional[float] = None
    conductance_func: Optional[ConductanceFunction] = None

    metadata: dict[str, Any] = field(default_factory=dict)

    def connects(self, node_id: int) -> bool:
        return node_id == self.node_i or node_id == self.node_j

    def G(self, temperatures: TemperatureMap) -> float:
        if self.conductance_func is not None:
            value = self.conductance_func(temperatures)
        elif self.conductance is not None:
            value = self.conductance
        else:
            raise ValueError(f"Ligação {self.name!r} não possui condutância.")

        value = float(value)

        if not math.isfinite(value):
            raise ValueError(f"Ligação {self.name!r} retornou G não finito.")
        if value < 0.0:
            raise ValueError(f"Ligação {self.name!r} retornou G negativo: {value}")

        return value

    def flux_into(self, node_id: int, temperatures: TemperatureMap) -> float:
        if not self.connects(node_id):
            raise ValueError(f"Nó {node_id} não pertence à ligação {self.name!r}.")

        Ti = temperatures[self.node_i]
        Tj = temperatures[self.node_j]
        G = self.G(temperatures)

        if self.direction == LinkDirection.UNDIRECTED:
            if node_id == self.node_i:
                return G * (Tj - Ti)
            return G * (Ti - Tj)

        if self.direction == LinkDirection.I_TO_J:
            if node_id == self.node_j:
                return G * (Ti - Tj)
            return 0.0

        if self.direction == LinkDirection.J_TO_I:
            if node_id == self.node_i:
                return G * (Tj - Ti)
            return 0.0

        raise ValueError(f"Direção inválida na ligação {self.name!r}.")


@dataclass
class NodalNetwork:
    nodes: dict[int, ThermalNode] = field(default_factory=dict)
    links: list[ThermalLink] = field(default_factory=list)
    cell_to_node: dict[tuple[int, int], int] = field(default_factory=dict)

    _next_node_id: int = 0

    def add_node(
        self,
        name: str,
        kind: NodeKind,
        x: float,
        y: float,
        *,
        z: float = 0.0,
        material: Optional[Any] = None,
        region: Optional[Any] = None,
        volume: float = 0.0,
        source: float = 0.0,
        temperature: float = 30.0,
        fixed_temperature: Optional[float] = None,
        i: Optional[int] = None,
        j: Optional[int] = None,
        bounds: Optional[tuple[float, float, float, float]] = None,
        metadata: Optional[dict[str, Any]] = None,
    ) -> int:
        if kind == NodeKind.ARITHMETIC and abs(volume) > 0.0:
            raise ValueError("Nó aritmético deve ter volume zero.")

        if kind == NodeKind.BOUNDARY and fixed_temperature is None:
            raise ValueError("Nó de contorno precisa de fixed_temperature.")

        node_id = self._next_node_id
        self._next_node_id += 1

        node = ThermalNode(
            id=node_id,
            name=name,
            kind=kind,
            x=float(x),
            y=float(y),
            z=float(z),
            material=material,
            region=region,
            volume=float(volume),
            source=float(source),
            temperature=float(temperature),
            fixed_temperature=fixed_temperature,
            i=i,
            j=j,
            bounds=bounds,
            metadata=metadata or {},
        )

        self.nodes[node_id] = node
        return node_id

    def add_link(
        self,
        node_i: int,
        node_j: int,
        kind: TransferKind,
        *,
        direction: LinkDirection = LinkDirection.UNDIRECTED,
        name: str = "",
        conductance: Optional[float] = None,
        conductance_func: Optional[ConductanceFunction] = None,
        metadata: Optional[dict[str, Any]] = None,
    ) -> None:
        if node_i not in self.nodes:
            raise KeyError(f"node_i={node_i} não existe.")
        if node_j not in self.nodes:
            raise KeyError(f"node_j={node_j} não existe.")
        if node_i == node_j:
            raise ValueError("Uma ligação não pode conectar um nó a ele mesmo.")

        self.links.append(
            ThermalLink(
                node_i=node_i,
                node_j=node_j,
                kind=kind,
                direction=direction,
                name=name,
                conductance=conductance,
                conductance_func=conductance_func,
                metadata=metadata or {},
            )
        )

    def links_of(self, node_id: int) -> list[ThermalLink]:
        return [link for link in self.links if link.connects(node_id)]

    def unknown_node_ids(self) -> list[int]:
        return [node_id for node_id, node in self.nodes.items() if node.is_unknown()]

    def initial_guess(self) -> np.ndarray:
        return np.array(
            [self.nodes[node_id].temperature for node_id in self.unknown_node_ids()],
            dtype=float,
        )

    def temperature_map_from_vector(self, z: np.ndarray) -> TemperatureMap:
        unknowns = self.unknown_node_ids()

        if len(z) != len(unknowns):
            raise ValueError(f"Esperado vetor com {len(unknowns)} temperaturas.")

        T: TemperatureMap = {}

        for node_id, value in zip(unknowns, z):
            T[node_id] = float(value)

        for node_id, node in self.nodes.items():
            if node.is_boundary():
                T[node_id] = float(node.fixed_temperature)

        return T

    def residual_steady(self, z: np.ndarray) -> np.ndarray:
        T = self.temperature_map_from_vector(z)
        residuals = []

        for node_id in self.unknown_node_ids():
            node = self.nodes[node_id]
            balance = node.source

            for link in self.links_of(node_id):
                balance += link.flux_into(node_id, T)

            residuals.append(balance)

        return np.array(residuals, dtype=float)

    def update_temperatures(self, z: np.ndarray) -> None:
        """
        Atualiza as temperaturas dos nós desconhecidos a partir do vetor solução.

        Nós BOUNDARY mantêm fixed_temperature.
        """
        unknowns = self.unknown_node_ids()

        if len(z) != len(unknowns):
            raise ValueError(f"Esperado vetor com {len(unknowns)} temperaturas.")

        for node_id, value in zip(unknowns, z):
            self.nodes[node_id].temperature = float(value)

        for node in self.nodes.values():
            if node.is_boundary():
                node.temperature = float(node.fixed_temperature)

    def temperatures_dict(self) -> dict[int, float]:
        """
        Retorna {node_id: temperatura_atual}.
        """
        return {node_id: node.temperature for node_id, node in self.nodes.items()}

    def nodes_in_region(self, region_name_value: str, *, only_cell_nodes: bool = True) -> list[int]:
        """
        Retorna ids de nós associados a uma região.

        Por padrão, retorna apenas nós vindos de células físicas da geometria,
        evitando nós aritméticos e de contorno.
        """
        found: list[int] = []

        for node_id, node in self.nodes.items():
            if node.region is None:
                continue

            if region_name(node.region) != region_name_value:
                continue

            if only_cell_nodes and not node.metadata.get("from_geometry_cell", False):
                continue

            found.append(node_id)

        return found

    def count_nodes_by_kind(self) -> dict[str, int]:
        counts = {kind.value: 0 for kind in NodeKind}
        for node in self.nodes.values():
            counts[node.kind.value] += 1
        return counts

    def count_links_by_kind(self) -> dict[str, int]:
        counts = {kind.value: 0 for kind in TransferKind}
        for link in self.links:
            counts[link.kind.value] += 1
        return counts

    def summary(self, max_nodes: int = 30, max_links: int = 40) -> str:
        lines = [
            "=== NodalNetwork ===",
            f"Nós totais: {len(self.nodes)}",
            f"Nós desconhecidos: {len(self.unknown_node_ids())}",
            f"Ligações totais: {len(self.links)}",
            "",
            "Nós por tipo:",
        ]

        for key, value in self.count_nodes_by_kind().items():
            lines.append(f"  - {key}: {value}")

        lines.append("")
        lines.append("Ligações por tipo:")

        for key, value in self.count_links_by_kind().items():
            lines.append(f"  - {key}: {value}")

        lines.append("")
        lines.append("Amostra de nós:")

        for idx, (node_id, node) in enumerate(self.nodes.items()):
            if idx >= max_nodes:
                lines.append(f"  ... {len(self.nodes) - max_nodes} nó(s) omitido(s)")
                break

            mat = "-" if node.material is None else material_name(node.material)
            reg = "-" if node.region is None else region_name(node.region)
            fixed = "" if node.fixed_temperature is None else f", Tfix={node.fixed_temperature:g}"

            lines.append(
                f"  {node_id}: {node.name}, kind={node.kind.value}, "
                f"region={reg}, mat={mat}, x={node.x:g}, y={node.y:g}, "
                f"V={node.volume:g}, Q={node.source:g}{fixed}"
            )

        lines.append("")
        lines.append("Amostra de ligações:")

        for idx, link in enumerate(self.links):
            if idx >= max_links:
                lines.append(f"  ... {len(self.links) - max_links} ligação(ões) omitida(s)")
                break

            lines.append(
                f"  {link.name}: {link.node_i} - {link.node_j}, "
                f"kind={link.kind.value}, direction={link.direction.value}"
            )

        return "\n".join(lines)


def build_network_from_geometry(
    geom: Any,
    *,
    thickness_z: float = 1.0,
    default_temperature: float = 30.0,
    pressure: Optional[float] = None,
    internal_convection_h: Optional[float] = None,
    add_arithmetic_interfaces: bool = True,
    apply_boundary_conditions: bool = True,
) -> NodalNetwork:
    if thickness_z <= 0.0:
        raise ValueError("thickness_z deve ser positiva.")

    if getattr(geom, "nx", None) is None:
        raise ValueError("A geometria precisa estar discretizada: chame geom.mesh(...).")

    net = NodalNetwork()
    region_grid = geom.region_object_grid()

    for j in range(geom.ny):
        for i in range(geom.nx):
            region = region_grid[j, i]

            if region is None or region == "VOID":
                continue

            x0, x1, y0, y1 = geom.cell_bounds(i, j)
            x, y = geom.cell_center(i, j)
            volume = (x1 - x0) * (y1 - y0) * thickness_z

            material = region.material
            kind = NodeKind.FLUID if is_fluid_material(material) else NodeKind.DIFFUSION

            Q = integrated_source_for_cell(
                geom=geom,
                region=region,
                x=x,
                y=y,
                x0=x0,
                x1=x1,
                y0=y0,
                y1=y1,
                thickness_z=thickness_z,
            )

            node_id = net.add_node(
                name=f"N_{j}_{i}_{region.name}",
                kind=kind,
                x=x,
                y=y,
                material=material,
                region=region,
                volume=volume,
                source=Q,
                temperature=default_temperature,
                i=i,
                j=j,
                bounds=(x0, x1, y0, y1),
                metadata={"from_geometry_cell": True, "cell": (j, i)},
            )

            net.cell_to_node[(j, i)] = node_id

    for j in range(geom.ny):
        for i in range(geom.nx):
            if (j, i) not in net.cell_to_node:
                continue

            if i + 1 < geom.nx and (j, i + 1) in net.cell_to_node:
                add_internal_link_between_cells(
                    geom=geom,
                    net=net,
                    cell_a=(j, i),
                    cell_b=(j, i + 1),
                    orientation="vertical",
                    thickness_z=thickness_z,
                    pressure=pressure,
                    internal_convection_h=internal_convection_h,
                    add_arithmetic_interfaces=add_arithmetic_interfaces,
                )

            if j + 1 < geom.ny and (j + 1, i) in net.cell_to_node:
                add_internal_link_between_cells(
                    geom=geom,
                    net=net,
                    cell_a=(j, i),
                    cell_b=(j + 1, i),
                    orientation="horizontal",
                    thickness_z=thickness_z,
                    pressure=pressure,
                    internal_convection_h=internal_convection_h,
                    add_arithmetic_interfaces=add_arithmetic_interfaces,
                )

    if apply_boundary_conditions:
        apply_external_boundary_conditions(
            geom=geom,
            net=net,
            thickness_z=thickness_z,
            pressure=pressure,
        )

    return net


def integrated_source_for_cell(
    geom: Any,
    region: Any,
    x: float,
    y: float,
    x0: float,
    x1: float,
    y0: float,
    y1: float,
    thickness_z: float,
) -> float:
    if not hasattr(geom, "source_for_region"):
        return 0.0

    Q = 0.0
    dx = x1 - x0
    dy = y1 - y0
    volume = dx * dy * thickness_z

    for source in geom.source_for_region(region.name):
        q = source.evaluate(x, y)

        if source.kind == "volumetric":
            Q += q * volume
        elif source.kind == "surface":
            Q += q * dx * dy
        else:
            raise ValueError(f"Tipo de fonte não suportado: {source.kind}")

    return Q


def face_data_between_cells(
    geom: Any,
    cell_a: tuple[int, int],
    cell_b: tuple[int, int],
    orientation: str,
    thickness_z: float,
) -> tuple[float, float, float, float, float]:
    ja, ia = cell_a
    jb, ib = cell_b

    ax0, ax1, ay0, ay1 = geom.cell_bounds(ia, ja)
    bx0, bx1, by0, by1 = geom.cell_bounds(ib, jb)

    if orientation == "vertical":
        x_face = ax1
        y_face = 0.5 * (max(ay0, by0) + min(ay1, by1))
        overlap = min(ay1, by1) - max(ay0, by0)
        area = overlap * thickness_z
        distance_a = 0.5 * (ax1 - ax0)
        distance_b = 0.5 * (bx1 - bx0)
        return x_face, y_face, area, distance_a, distance_b

    if orientation == "horizontal":
        y_face = ay1
        x_face = 0.5 * (max(ax0, bx0) + min(ax1, bx1))
        overlap = min(ax1, bx1) - max(ax0, bx0)
        area = overlap * thickness_z
        distance_a = 0.5 * (ay1 - ay0)
        distance_b = 0.5 * (by1 - by0)
        return x_face, y_face, area, distance_a, distance_b

    raise ValueError("orientation deve ser 'vertical' ou 'horizontal'.")


def add_internal_link_between_cells(
    geom: Any,
    net: NodalNetwork,
    cell_a: tuple[int, int],
    cell_b: tuple[int, int],
    orientation: str,
    thickness_z: float,
    pressure: Optional[float],
    internal_convection_h: Optional[float],
    add_arithmetic_interfaces: bool,
) -> None:
    node_a = net.cell_to_node[cell_a]
    node_b = net.cell_to_node[cell_b]

    A = net.nodes[node_a]
    B = net.nodes[node_b]

    region_a = A.region
    region_b = B.region

    if region_a is None or region_b is None:
        return

    mode = geom.transfer_mode_between_regions(region_a, region_b)

    x_face, y_face, area, dist_a, dist_b = face_data_between_cells(
        geom=geom,
        cell_a=cell_a,
        cell_b=cell_b,
        orientation=orientation,
        thickness_z=thickness_z,
    )

    if area <= 0.0:
        return

    if mode == "conduction":
        add_conduction_between_cells(
            net=net,
            node_a=node_a,
            node_b=node_b,
            x_face=x_face,
            y_face=y_face,
            area=area,
            dist_a=dist_a,
            dist_b=dist_b,
            pressure=pressure,
            add_arithmetic_interface=add_arithmetic_interfaces,
        )
        return

    if mode == "convection":
        add_internal_convection_between_cells(
            net=net,
            node_a=node_a,
            node_b=node_b,
            x_face=x_face,
            y_face=y_face,
            area=area,
            dist_a=dist_a,
            dist_b=dist_b,
            h=internal_convection_h,
            pressure=pressure,
        )
        return

    if mode == "radiation":
        add_internal_radiation_between_cells(
            net=net,
            node_a=node_a,
            node_b=node_b,
            x_face=x_face,
            y_face=y_face,
            area=area,
            dist_a=dist_a,
            dist_b=dist_b,
            pressure=pressure,
            emissivity=1.0,
            view_factor=1.0,
        )
        return

    if mode == "fluid_transport":
        # Em Geometry2D, dois materiais fluidos sugerem "fluid_transport".
        # Porém o transporte entálpico exige direção de escoamento e vazão mássica.
        # Nesta etapa transversal da rede, não inferimos isso automaticamente.
        # As ligações fluido-fluido direcionais devem ser adicionadas depois com:
        #     add_fluid_transport_link(...)
        # ou
        #     add_fluid_transport_chain(...)
        return

    raise ValueError(f"Modo térmico não suportado: {mode}")


def add_conduction_between_cells(
    net: NodalNetwork,
    node_a: int,
    node_b: int,
    x_face: float,
    y_face: float,
    area: float,
    dist_a: float,
    dist_b: float,
    pressure: Optional[float],
    add_arithmetic_interface: bool,
) -> None:
    A = net.nodes[node_a]
    B = net.nodes[node_b]

    same_region = A.region is B.region
    same_material = material_name(A.material) == material_name(B.material)

    if same_region and same_material:
        net.add_link(
            node_i=node_a,
            node_j=node_b,
            kind=TransferKind.CONDUCTION,
            name=f"cond_{A.name}_to_{B.name}",
            conductance_func=conduction_func_between_nodes(
                material=A.material,
                node_i=node_a,
                node_j=node_b,
                area=area,
                distance=dist_a + dist_b,
                pressure=pressure,
            ),
            metadata={"area": area, "distance": dist_a + dist_b},
        )
        return

    if add_arithmetic_interface:
        interface_id = net.add_node(
            name=f"I_{A.name}_to_{B.name}",
            kind=NodeKind.ARITHMETIC,
            x=x_face,
            y=y_face,
            volume=0.0,
            source=0.0,
            temperature=0.5 * (A.temperature + B.temperature),
            metadata={
                "type": "conductive_interface",
                "between": (node_a, node_b),
                "area": area,
            },
        )

        net.add_link(
            node_i=node_a,
            node_j=interface_id,
            kind=TransferKind.CONDUCTION,
            name=f"cond_{A.name}_to_I{interface_id}",
            conductance_func=conduction_func_single_material(
                material=A.material,
                node_temperature_id=node_a,
                area=area,
                distance=dist_a,
                pressure=pressure,
            ),
            metadata={"area": area, "distance": dist_a},
        )

        net.add_link(
            node_i=interface_id,
            node_j=node_b,
            kind=TransferKind.CONDUCTION,
            name=f"cond_I{interface_id}_to_{B.name}",
            conductance_func=conduction_func_single_material(
                material=B.material,
                node_temperature_id=node_b,
                area=area,
                distance=dist_b,
                pressure=pressure,
            ),
            metadata={"area": area, "distance": dist_b},
        )
        return

    def geq_func(T: TemperatureMap) -> float:
        k_a = conductivity(A.material, T[node_a], pressure)
        k_b = conductivity(B.material, T[node_b], pressure)
        G_a = conduction_G(k_a, area, dist_a)
        G_b = conduction_G(k_b, area, dist_b)
        return equivalent_series_G(G_a, G_b)

    net.add_link(
        node_i=node_a,
        node_j=node_b,
        kind=TransferKind.EQUIVALENT,
        name=f"geq_{A.name}_to_{B.name}",
        conductance_func=geq_func,
        metadata={"area": area, "dist_a": dist_a, "dist_b": dist_b},
    )


def add_internal_convection_between_cells(
    net: NodalNetwork,
    node_a: int,
    node_b: int,
    x_face: float,
    y_face: float,
    area: float,
    dist_a: float,
    dist_b: float,
    h: Optional[float],
    pressure: Optional[float],
) -> None:
    if h is None:
        raise ValueError(
            "Interface interna por convecção encontrada, mas internal_convection_h=None. "
            "Forneça internal_convection_h ao chamar build_network_from_geometry(...)."
        )

    A = net.nodes[node_a]
    B = net.nodes[node_b]

    a_fluid = is_fluid_material(A.material)
    b_fluid = is_fluid_material(B.material)

    if a_fluid and not b_fluid:
        fluid_node = node_a
        solid_node = node_b
        solid_distance = dist_b
        solid_material = B.material
    else:
        fluid_node = node_b
        solid_node = node_a
        solid_distance = dist_a
        solid_material = A.material

    solid = net.nodes[solid_node]
    fluid = net.nodes[fluid_node]

    surface_id = net.add_node(
        name=f"S_{solid.name}_to_{fluid.name}",
        kind=NodeKind.ARITHMETIC,
        x=x_face,
        y=y_face,
        volume=0.0,
        source=0.0,
        temperature=solid.temperature,
        metadata={
            "type": "convective_surface",
            "solid_node": solid_node,
            "fluid_node": fluid_node,
            "area": area,
        },
    )

    net.add_link(
        node_i=solid_node,
        node_j=surface_id,
        kind=TransferKind.CONDUCTION,
        name=f"cond_{solid.name}_to_S{surface_id}",
        conductance_func=conduction_func_single_material(
            material=solid_material,
            node_temperature_id=solid_node,
            area=area,
            distance=solid_distance,
            pressure=pressure,
        ),
        metadata={"area": area, "distance": solid_distance},
    )

    net.add_link(
        node_i=surface_id,
        node_j=fluid_node,
        kind=TransferKind.CONVECTION,
        name=f"conv_S{surface_id}_to_{fluid.name}",
        conductance=convection_G(h, area),
        metadata={"h": h, "area": area},
    )


def add_internal_radiation_between_cells(
    net: NodalNetwork,
    node_a: int,
    node_b: int,
    x_face: float,
    y_face: float,
    area: float,
    dist_a: float,
    dist_b: float,
    pressure: Optional[float],
    emissivity: float,
    view_factor: float,
) -> None:
    A = net.nodes[node_a]
    B = net.nodes[node_b]

    surf_a = net.add_node(
        name=f"Srad_{A.name}",
        kind=NodeKind.ARITHMETIC,
        x=x_face,
        y=y_face,
        volume=0.0,
        temperature=A.temperature,
        metadata={"type": "radiation_surface", "side": "A"},
    )

    surf_b = net.add_node(
        name=f"Srad_{B.name}",
        kind=NodeKind.ARITHMETIC,
        x=x_face,
        y=y_face,
        volume=0.0,
        temperature=B.temperature,
        metadata={"type": "radiation_surface", "side": "B"},
    )

    net.add_link(
        node_i=node_a,
        node_j=surf_a,
        kind=TransferKind.CONDUCTION,
        name=f"cond_{A.name}_to_Srad{surf_a}",
        conductance_func=conduction_func_single_material(
            material=A.material,
            node_temperature_id=node_a,
            area=area,
            distance=dist_a,
            pressure=pressure,
        ),
    )

    net.add_link(
        node_i=surf_b,
        node_j=node_b,
        kind=TransferKind.CONDUCTION,
        name=f"cond_Srad{surf_b}_to_{B.name}",
        conductance_func=conduction_func_single_material(
            material=B.material,
            node_temperature_id=node_b,
            area=area,
            distance=dist_b,
            pressure=pressure,
        ),
    )

    net.add_link(
        node_i=surf_a,
        node_j=surf_b,
        kind=TransferKind.RADIATION,
        name=f"rad_S{surf_a}_to_S{surf_b}",
        conductance_func=radiation_func_between_nodes(
            node_i=surf_a,
            node_j=surf_b,
            emissivity=emissivity,
            area=area,
            view_factor=view_factor,
        ),
        metadata={"area": area, "emissivity": emissivity, "view_factor": view_factor},
    )


def apply_external_boundary_conditions(
    geom: Any,
    net: NodalNetwork,
    thickness_z: float,
    pressure: Optional[float],
) -> None:
    for bc in getattr(geom, "boundary_conditions", []):
        side = bc.side
        kind = bc.kind
        data = dict(bc.data or {})

        for j, i in cells_on_side(geom, net, side):
            node_id = net.cell_to_node[(j, i)]
            node = net.nodes[node_id]

            x0, x1, y0, y1 = geom.cell_bounds(i, j)
            face_x, face_y, area, distance = boundary_face_data(
                side=side,
                x0=x0,
                x1=x1,
                y0=y0,
                y1=y1,
                thickness_z=thickness_z,
            )

            if kind in ("symmetry", "adiabatic"):
                continue

            if kind in ("temperature", "dirichlet"):
                T_bc = get_temperature_from_data(data)

                bnode = net.add_node(
                    name=f"BC_{side}_{j}_{i}",
                    kind=NodeKind.BOUNDARY,
                    x=face_x,
                    y=face_y,
                    volume=0.0,
                    fixed_temperature=T_bc,
                    temperature=T_bc,
                    metadata={"bc": kind, "side": side, "data": data},
                )

                net.add_link(
                    node_i=node_id,
                    node_j=bnode,
                    kind=TransferKind.CONDUCTION,
                    name=f"bc_temp_{side}_{j}_{i}",
                    conductance_func=conduction_func_single_material(
                        material=node.material,
                        node_temperature_id=node_id,
                        area=area,
                        distance=distance,
                        pressure=pressure,
                    ),
                    metadata={"side": side, "area": area, "distance": distance},
                )
                continue

            if kind in ("heat_flux", "flux", "neumann"):
                q_in = get_heat_flux_from_data(data)
                node.add_source(q_in * area)
                continue

            if kind in ("convection", "convective"):
                h = get_h_from_data(data)
                T_inf = get_temperature_from_data(data, keys=("T_inf", "t_inf", "T", "temperature", "value"))

                surface = net.add_node(
                    name=f"Sbc_{side}_{j}_{i}",
                    kind=NodeKind.ARITHMETIC,
                    x=face_x,
                    y=face_y,
                    volume=0.0,
                    temperature=node.temperature,
                    metadata={"type": "boundary_surface", "side": side, "bc": data},
                )

                ambient = net.add_node(
                    name=f"AMB_{side}_{j}_{i}",
                    kind=NodeKind.BOUNDARY,
                    x=face_x,
                    y=face_y,
                    volume=0.0,
                    fixed_temperature=T_inf,
                    temperature=T_inf,
                    metadata={"type": "ambient", "side": side, "bc": data},
                )

                net.add_link(
                    node_i=node_id,
                    node_j=surface,
                    kind=TransferKind.CONDUCTION,
                    name=f"cond_{node.name}_to_Sbc{surface}",
                    conductance_func=conduction_func_single_material(
                        material=node.material,
                        node_temperature_id=node_id,
                        area=area,
                        distance=distance,
                        pressure=pressure,
                    ),
                    metadata={"side": side, "area": area, "distance": distance},
                )

                net.add_link(
                    node_i=surface,
                    node_j=ambient,
                    kind=TransferKind.CONVECTION,
                    name=f"conv_Sbc{surface}_to_AMB{ambient}",
                    conductance=convection_G(h, area),
                    metadata={"side": side, "h": h, "area": area},
                )
                continue

            if kind in ("radiation", "rad"):
                T_env = get_temperature_from_data(data, keys=("T_env", "T_inf", "T", "temperature", "value"))
                emissivity = float(data.get("emissivity", data.get("epsilon", 1.0)))
                view_factor = float(data.get("view_factor", data.get("F", 1.0)))

                surface = net.add_node(
                    name=f"Sradbc_{side}_{j}_{i}",
                    kind=NodeKind.ARITHMETIC,
                    x=face_x,
                    y=face_y,
                    volume=0.0,
                    temperature=node.temperature,
                    metadata={"type": "boundary_radiation_surface", "side": side, "bc": data},
                )

                env = net.add_node(
                    name=f"RADENV_{side}_{j}_{i}",
                    kind=NodeKind.BOUNDARY,
                    x=face_x,
                    y=face_y,
                    volume=0.0,
                    fixed_temperature=T_env,
                    temperature=T_env,
                    metadata={"type": "radiation_environment", "side": side, "bc": data},
                )

                net.add_link(
                    node_i=node_id,
                    node_j=surface,
                    kind=TransferKind.CONDUCTION,
                    name=f"cond_{node.name}_to_Sradbc{surface}",
                    conductance_func=conduction_func_single_material(
                        material=node.material,
                        node_temperature_id=node_id,
                        area=area,
                        distance=distance,
                        pressure=pressure,
                    ),
                )

                net.add_link(
                    node_i=surface,
                    node_j=env,
                    kind=TransferKind.RADIATION,
                    name=f"rad_Sradbc{surface}_to_ENV{env}",
                    conductance_func=radiation_func_between_nodes(
                        node_i=surface,
                        node_j=env,
                        emissivity=emissivity,
                        area=area,
                        view_factor=view_factor,
                    ),
                    metadata={
                        "side": side,
                        "area": area,
                        "emissivity": emissivity,
                        "view_factor": view_factor,
                    },
                )
                continue

            if kind in ("inlet", "outlet"):
                continue

            raise ValueError(f"Condição de contorno não suportada: {kind}")


def cells_on_side(geom: Any, net: NodalNetwork, side: str) -> list[tuple[int, int]]:
    cells: list[tuple[int, int]] = []

    if side == "left":
        i = 0
        for j in range(geom.ny):
            if (j, i) in net.cell_to_node:
                cells.append((j, i))
        return cells

    if side == "right":
        i = geom.nx - 1
        for j in range(geom.ny):
            if (j, i) in net.cell_to_node:
                cells.append((j, i))
        return cells

    if side == "bottom":
        j = 0
        for i in range(geom.nx):
            if (j, i) in net.cell_to_node:
                cells.append((j, i))
        return cells

    if side == "top":
        j = geom.ny - 1
        for i in range(geom.nx):
            if (j, i) in net.cell_to_node:
                cells.append((j, i))
        return cells

    raise ValueError(f"Lado inválido: {side}")


def boundary_face_data(
    side: str,
    x0: float,
    x1: float,
    y0: float,
    y1: float,
    thickness_z: float,
) -> tuple[float, float, float, float]:
    if side == "left":
        return x0, 0.5 * (y0 + y1), (y1 - y0) * thickness_z, 0.5 * (x1 - x0)

    if side == "right":
        return x1, 0.5 * (y0 + y1), (y1 - y0) * thickness_z, 0.5 * (x1 - x0)

    if side == "bottom":
        return 0.5 * (x0 + x1), y0, (x1 - x0) * thickness_z, 0.5 * (y1 - y0)

    if side == "top":
        return 0.5 * (x0 + x1), y1, (x1 - x0) * thickness_z, 0.5 * (y1 - y0)

    raise ValueError(f"Lado inválido: {side}")


def get_temperature_from_data(
    data: dict[str, Any],
    keys: tuple[str, ...] = ("T", "temperature", "T_inf", "t_inf", "T_env", "value"),
) -> float:
    for key in keys:
        if key in data and data[key] is not None:
            return float(data[key])
    raise ValueError(f"Dados de contorno sem temperatura: {data}")


def get_h_from_data(data: dict[str, Any]) -> float:
    for key in ("h", "coef", "coefficient", "value"):
        if key in data and data[key] is not None:
            return float(data[key])
    raise ValueError(f"Dados de convecção sem h: {data}")


def get_heat_flux_from_data(data: dict[str, Any]) -> float:
    for key in ("q", "q_flux", "heat_flux", "value"):
        if key in data and data[key] is not None:
            return float(data[key])
    raise ValueError(f"Dados de fluxo sem q: {data}")


def add_fluid_transport_link(
    net: NodalNetwork,
    upstream: int,
    downstream: int,
    m_dot: float,
    cp: float,
    *,
    name: Optional[str] = None,
) -> None:
    net.add_link(
        node_i=upstream,
        node_j=downstream,
        kind=TransferKind.FLUID_TRANSPORT,
        direction=LinkDirection.I_TO_J,
        name=name or f"fluid_{upstream}_to_{downstream}",
        conductance=fluid_transport_G(m_dot, cp),
        metadata={"m_dot": m_dot, "cp": cp},
    )


def add_fluid_transport_chain(
    net: NodalNetwork,
    node_ids: list[int],
    m_dot: float,
    cp: float,
    *,
    name_prefix: str = "fluid_chain",
) -> None:
    if len(node_ids) < 2:
        return

    for n, (upstream, downstream) in enumerate(zip(node_ids[:-1], node_ids[1:])):
        add_fluid_transport_link(
            net=net,
            upstream=upstream,
            downstream=downstream,
            m_dot=m_dot,
            cp=cp,
            name=f"{name_prefix}_{n}",
        )


# ============================================================
# Transporte fluido automático por região
# ============================================================

def fluid_node_chains_by_region(
    net: NodalNetwork,
    region_name_value: str,
    direction: str,
) -> list[list[int]]:
    """
    Agrupa e ordena os nós de uma região fluida em cadeias direcionais.

    direction:
        "x+" : escoamento da esquerda para a direita.
        "x-" : escoamento da direita para a esquerda.
        "y+" : escoamento de baixo para cima.
        "y-" : escoamento de cima para baixo.

    Para direction="x+", cada linha j vira uma cadeia:
        (j, i0) -> (j, i1) -> ...

    Para direction="y-", cada coluna i vira uma cadeia:
        (j_top, i) -> ... -> (j_bottom, i)
    """
    direction = str(direction).strip().lower()

    allowed = {"x+", "x-", "y+", "y-"}
    if direction not in allowed:
        raise ValueError(f"direction deve ser um de {sorted(allowed)}.")

    node_ids = net.nodes_in_region(region_name_value, only_cell_nodes=True)

    if not node_ids:
        raise ValueError(f"Nenhum nó de célula encontrado para a região {region_name_value!r}.")

    groups: dict[int, list[int]] = {}

    if direction in ("x+", "x-"):
        for node_id in node_ids:
            node = net.nodes[node_id]
            if node.j is None or node.i is None:
                continue
            groups.setdefault(node.j, []).append(node_id)

        chains = []
        reverse = direction == "x-"
        for key in sorted(groups):
            chain = sorted(groups[key], key=lambda nid: net.nodes[nid].i, reverse=reverse)
            if len(chain) >= 1:
                chains.append(chain)

        return chains

    # direction in ("y+", "y-")
    for node_id in node_ids:
        node = net.nodes[node_id]
        if node.j is None or node.i is None:
            continue
        groups.setdefault(node.i, []).append(node_id)

    chains = []
    reverse = direction == "y-"
    for key in sorted(groups):
        chain = sorted(groups[key], key=lambda nid: net.nodes[nid].j, reverse=reverse)
        if len(chain) >= 1:
            chains.append(chain)

    return chains


def add_fluid_transport_link_auto_cp(
    net: NodalNetwork,
    upstream: int,
    downstream: int,
    m_dot: float,
    *,
    cp: Optional[float] = None,
    pressure: Optional[float] = None,
    name: Optional[str] = None,
) -> None:
    """
    Adiciona ligação direcional upstream -> downstream.

    Se cp for dado:
        Gf = m_dot * cp

    Se cp=None:
        cp é calculado com o material do nó downstream na temperatura média
        entre upstream e downstream.
    """
    if cp is not None:
        add_fluid_transport_link(
            net=net,
            upstream=upstream,
            downstream=downstream,
            m_dot=m_dot,
            cp=float(cp),
            name=name,
        )
        return

    def gf_func(T: TemperatureMap) -> float:
        node = net.nodes[downstream]
        Tm = 0.5 * (T[upstream] + T[downstream])
        cp_value = specific_heat(node.material, Tm, pressure)
        return fluid_transport_G(m_dot, cp_value)

    net.add_link(
        node_i=upstream,
        node_j=downstream,
        kind=TransferKind.FLUID_TRANSPORT,
        direction=LinkDirection.I_TO_J,
        name=name or f"fluid_{upstream}_to_{downstream}",
        conductance_func=gf_func,
        metadata={"m_dot": m_dot, "cp": cp, "auto_cp": cp is None},
    )


def add_fluid_transport_for_region(
    net: NodalNetwork,
    region_name_value: str,
    *,
    direction: str,
    m_dot_total: float,
    cp: Optional[float] = None,
    pressure: Optional[float] = None,
    inlet_temperature: Optional[float] = None,
    distribute_m_dot: bool = True,
    name_prefix: Optional[str] = None,
) -> list[list[int]]:
    """
    Adiciona transporte entálpico direcional em uma região fluida.

    Parâmetros
    ----------
    region_name_value:
        Nome da região fluida, por exemplo "canal_agua".

    direction:
        "x+", "x-", "y+" ou "y-".

    m_dot_total:
        Vazão mássica total associada à região.

    distribute_m_dot:
        Se True, divide m_dot_total pelo número de cadeias paralelas.
        Exemplo: em direction="x+", cada linha j é uma cadeia paralela.

    inlet_temperature:
        Se fornecida, cria um nó BOUNDARY na entrada de cada cadeia e liga:
            entrada_fixa -> primeiro_nó_fluido
        Isso fecha o balanço entálpico do primeiro nó.

    Retorna
    -------
    chains:
        Lista das cadeias de nós fluido usadas.
    """
    if m_dot_total < 0.0:
        raise ValueError("m_dot_total não pode ser negativo.")

    chains = fluid_node_chains_by_region(
        net=net,
        region_name_value=region_name_value,
        direction=direction,
    )

    if not chains:
        raise ValueError(f"Nenhuma cadeia fluida encontrada para {region_name_value!r}.")

    m_dot_chain = m_dot_total / len(chains) if distribute_m_dot else m_dot_total
    prefix = name_prefix or f"fluid_{region_name_value}_{direction}"

    for chain_index, chain in enumerate(chains):
        if inlet_temperature is not None:
            first = chain[0]
            first_node = net.nodes[first]

            inlet = net.add_node(
                name=f"INLET_{region_name_value}_{chain_index}",
                kind=NodeKind.BOUNDARY,
                x=first_node.x,
                y=first_node.y,
                volume=0.0,
                fixed_temperature=float(inlet_temperature),
                temperature=float(inlet_temperature),
                metadata={
                    "type": "fluid_inlet",
                    "region": region_name_value,
                    "chain_index": chain_index,
                    "direction": direction,
                },
            )

            add_fluid_transport_link_auto_cp(
                net=net,
                upstream=inlet,
                downstream=first,
                m_dot=m_dot_chain,
                cp=cp,
                pressure=pressure,
                name=f"{prefix}_inlet_{chain_index}",
            )

        for link_index, (upstream, downstream) in enumerate(zip(chain[:-1], chain[1:])):
            add_fluid_transport_link_auto_cp(
                net=net,
                upstream=upstream,
                downstream=downstream,
                m_dot=m_dot_chain,
                cp=cp,
                pressure=pressure,
                name=f"{prefix}_{chain_index}_{link_index}",
            )

    return chains


def extract_region_temperatures(net: NodalNetwork, region_name_value: str) -> list[dict[str, float]]:
    """
    Extrai temperaturas dos nós físicos de uma região como lista de dicionários.
    """
    rows = []

    for node_id in net.nodes_in_region(region_name_value, only_cell_nodes=True):
        node = net.nodes[node_id]
        rows.append(
            {
                "node_id": node_id,
                "i": node.i,
                "j": node.j,
                "x": node.x,
                "y": node.y,
                "temperature": node.temperature,
                "source": node.source,
                "volume": node.volume,
            }
        )

    rows.sort(key=lambda row: (row["j"], row["i"]))
    return rows



if __name__ == "__main__":
    from materiais import Material, MaterialPhase
    from geometria import Geometry2D

    aluminio = Material(name="Aluminio", phase=MaterialPhase.SOLID, k=180.0)
    ar = Material(name="Ar", phase=MaterialPhase.FLUID, k=0.026, rho=1.2, cp=1007.0, mu=1.8e-5)
    agua = Material(name="Agua", phase=MaterialPhase.FLUID, k=0.60, rho=997.0, cp=4180.0, mu=8.9e-4)

    geom = Geometry2D(width=0.010, height=0.004)

    geom.material("aluminio_esq", aluminio, 0.000, 0.003, 0.000, 0.004)
    geom.material("gap_ar", ar, 0.003, 0.004, 0.000, 0.004, thermal_mode="conduction")
    geom.material("aluminio_dir", aluminio, 0.004, 0.007, 0.000, 0.004)
    geom.material("canal_agua", agua, 0.007, 0.010, 0.000, 0.004)

    geom.bc("left", "temperature", T=100.0)
    geom.bc("right", "convection", h=1000.0, T_inf=30.0)
    geom.bc("top", "symmetry")
    geom.bc("bottom", "symmetry")

    geom.source("aluminio_esq", "volumetric", value=2.5e6, name="fonte_aluminio_esq")
    geom.interface("gap_ar", "aluminio_dir", mode="conduction", name="interface_ar_aluminio_conducao")

    geom.mesh(nx=7, ny=3, fix="auto")

    net = build_network_from_geometry(
        geom,
        thickness_z=1.0,
        default_temperature=30.0,
        internal_convection_h=1000.0,
    )

    print(net.summary())
    print()
    z0 = net.initial_guess()
    print("Tamanho de z0:", len(z0))
    print("Norma do resíduo inicial:", float(np.linalg.norm(net.residual_steady(z0))))
