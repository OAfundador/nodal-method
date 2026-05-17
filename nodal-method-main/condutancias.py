"""
condutancias.py

Cálculos fundamentais de condutâncias térmicas para o método nodal.

Este módulo concentra as fórmulas físicas elementares. O solver não deve
conhecer detalhes geométricos; ele apenas resolve R(T)=0. As condutâncias
entram na montagem da rede nodal.

Convenções:
- comprimento: m
- área: m²
- temperatura: °C ou K apenas para diferenças; radiação usa conversão para K
- condutância: W/K
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Iterable


def _positive(name: str, value: float) -> float:
    value = float(value)
    if not isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} deve ser positivo. Recebido: {value!r}")
    return value


def _non_negative(name: str, value: float) -> float:
    value = float(value)
    if not isfinite(value) or value < 0.0:
        raise ValueError(f"{name} não pode ser negativo. Recebido: {value!r}")
    return value


def conduction_G(k: float, area: float, distance: float) -> float:
    """
    Condutância condutiva cartesiana:

        G = k A / L

    onde:
        k        condutividade térmica [W/(m K)]
        area     área perpendicular ao fluxo [m²]
        distance distância entre nós ou meia-distância apropriada [m]
    """
    k = _positive("k", k)
    area = _positive("area", area)
    distance = _positive("distance", distance)
    return k * area / distance


def resistance_from_conduction(k: float, area: float, distance: float) -> float:
    """
    Resistência condutiva cartesiana:

        R = L / (k A)
    """
    return 1.0 / conduction_G(k, area, distance)


def convection_G(h: float, area: float) -> float:
    """
    Condutância convectiva:

        G = h A
    """
    h = _non_negative("h", h)
    area = _positive("area", area)
    return h * area


def resistance_from_convection(h: float, area: float) -> float:
    """
    Resistência convectiva:

        R = 1 / (h A)
    """
    G = convection_G(h, area)
    if G <= 0.0:
        raise ValueError("A resistência convectiva não é definida para h=0.")
    return 1.0 / G


def fluid_transport_G(m_dot: float, cp: float) -> float:
    """
    Condutância entálpica de transporte de massa:

        G_f = m_dot cp

    Esta ligação é direcional: montante -> jusante.
    """
    m_dot = _non_negative("m_dot", m_dot)
    cp = _positive("cp", cp)
    return m_dot * cp


def equivalent_series_G(*conductances: float) -> float:
    """
    Condutância equivalente em série:

        1/G_eq = Σ 1/G_i

    Para duas condutâncias:
        G_eq = G1 G2 / (G1 + G2)
    """
    if len(conductances) == 1 and not isinstance(conductances[0], (int, float)):
        conductances = tuple(conductances[0])

    if len(conductances) == 0:
        raise ValueError("Informe pelo menos uma condutância.")

    inv_sum = 0.0
    for idx, G in enumerate(conductances, start=1):
        G = _positive(f"G{idx}", G)
        inv_sum += 1.0 / G

    return 1.0 / inv_sum


def equivalent_parallel_G(*conductances: float) -> float:
    """
    Condutância equivalente em paralelo:

        G_eq = Σ G_i
    """
    if len(conductances) == 1 and not isinstance(conductances[0], (int, float)):
        conductances = tuple(conductances[0])

    if len(conductances) == 0:
        raise ValueError("Informe pelo menos uma condutância.")

    return sum(_positive(f"G{idx}", G) for idx, G in enumerate(conductances, start=1))


def interface_conduction_G(
    k_left: float,
    k_right: float,
    area: float,
    distance_left: float,
    distance_right: float,
) -> float:
    """
    Condutância equivalente entre dois materiais diferentes em coordenadas cartesianas.

    Modelo:
        nó esquerdo -- interface -- nó direito

        G_left  = k_left  A / L_left
        G_right = k_right A / L_right
        G_eq    = série(G_left, G_right)
    """
    G_left = conduction_G(k_left, area, distance_left)
    G_right = conduction_G(k_right, area, distance_right)
    return equivalent_series_G(G_left, G_right)


def wall_with_convection_G(
    k: float,
    thickness: float,
    area: float,
    h: float,
) -> float:
    """
    Condutância equivalente de parede plana + convecção:

        R_total = L/(kA) + 1/(hA)
        G_eq = 1/R_total
    """
    G_cond = conduction_G(k, area, thickness)
    G_conv = convection_G(h, area)
    return equivalent_series_G(G_cond, G_conv)


def radiation_linearized_G(
    emissivity: float,
    area: float,
    T_i_C: float,
    T_j_C: float,
    view_factor: float = 1.0,
    sigma: float = 5.670374419e-8,
) -> float:
    """
    Condutância radiativa linearizada:

        q = ε σ F A (Ti^4 - Tj^4)
        q = G_rad (Ti - Tj)

        G_rad = ε σ F A (Ti + Tj)(Ti² + Tj²)

    Temperaturas de entrada em °C; cálculo interno em K.
    """
    emissivity = float(emissivity)
    view_factor = float(view_factor)

    if not (0.0 <= emissivity <= 1.0):
        raise ValueError("emissivity deve estar entre 0 e 1.")
    if not (0.0 <= view_factor <= 1.0):
        raise ValueError("view_factor deve estar entre 0 e 1.")

    area = _positive("area", area)
    Ti = float(T_i_C) + 273.15
    Tj = float(T_j_C) + 273.15

    return emissivity * sigma * view_factor * area * (Ti + Tj) * (Ti**2 + Tj**2)


@dataclass(frozen=True)
class CartesianLinkData:
    """
    Registro simples para documentar a origem geométrica de uma ligação cartesiana.
    """
    name: str
    kind: str
    area: float
    distance: float | None = None
    k: float | None = None
    h: float | None = None
    conductance: float | None = None


def print_conductance_table(rows: list[tuple[str, float]]) -> None:
    """
    Impressão padronizada de condutâncias.
    """
    print("nome                 G [W/K]")
    print("-------------------  ----------------")
    for name, value in rows:
        print(f"{name:19s}  {value:16.10e}")
