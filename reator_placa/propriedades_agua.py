"""
propriedades_agua.py

Propriedades termofísicas da água líquida subresfriada para uso didático
em problemas de troca convectiva forçada.

Faixa válida: 20-100 °C, P em torno de 100-200 kPa (dependência fraca de P
para a fase líquida nessa janela). Correlações simples, ajustadas a tabelas
padrão (Incropera Tab. A.6 e similares).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from materiais import Material, MaterialPhase


P_REF_PA = 150e3


def rho_agua(T_C: float, P_Pa: Optional[float] = None) -> float:
    """Densidade [kg/m³]."""
    T = float(T_C)
    return (
        999.842594
        + 6.793952e-2 * T
        - 9.095290e-3 * T ** 2
        + 1.001685e-4 * T ** 3
        - 1.120083e-6 * T ** 4
        + 6.536332e-9 * T ** 5
    )


def mu_agua(T_C: float, P_Pa: Optional[float] = None) -> float:
    """Viscosidade dinâmica [Pa·s] (correlação de Reichardt)."""
    T_K = float(T_C) + 273.15
    return 2.414e-5 * 10.0 ** (247.8 / (T_K - 140.0))


def cp_agua(T_C: float, P_Pa: Optional[float] = None) -> float:
    """Calor específico isobárico [J/(kg·K)]."""
    T = float(T_C)
    return 4217.0 - 3.358 * T + 0.04148 * T ** 2 - 1.6e-4 * T ** 3


def k_agua(T_C: float, P_Pa: Optional[float] = None) -> float:
    """Condutividade térmica [W/(m·K)]."""
    T = float(T_C)
    return 0.5615 + 1.939e-3 * T - 7.51e-6 * T ** 2


def Pr_agua(T_C: float, P_Pa: Optional[float] = None) -> float:
    return cp_agua(T_C, P_Pa) * mu_agua(T_C, P_Pa) / k_agua(T_C, P_Pa)


def T_sat_agua(P_Pa: float) -> float:
    """Temperatura de saturação [°C] - equação de Antoine."""
    P_mmHg = float(P_Pa) / 133.322
    A, B, C = 8.07131, 1730.63, 233.426
    if P_mmHg <= 0.0:
        raise ValueError("P deve ser positiva.")
    return B / (A - math.log10(P_mmHg)) - C


def criar_agua_didatica() -> Material:
    return Material(
        name="Água",
        phase=MaterialPhase.FLUID,
        k_func=lambda T, P=None: k_agua(T, P),
        rho_func=lambda T, P=None: rho_agua(T, P),
        cp_func=lambda T, P=None: cp_agua(T, P),
        mu_func=lambda T, P=None: mu_agua(T, P),
    )


def coeficiente_pelicula(
    T_C: float, P_Pa: float,
    m_dot_canal: float, area_flow: float, Dh: float,
) -> dict:
    """
    Coeficiente convectivo h por Dittus-Boelter:

        Nu = 0.023 · Re^0.8 · Pr^0.4
        h  = Nu · k_fluido / Dh
    """
    rho = rho_agua(T_C, P_Pa)
    mu = mu_agua(T_C, P_Pa)
    cp = cp_agua(T_C, P_Pa)
    k = k_agua(T_C, P_Pa)
    v = m_dot_canal / (rho * area_flow)
    Re = rho * v * Dh / mu
    Pr = cp * mu / k
    Nu = 0.023 * (Re ** 0.8) * (Pr ** 0.4)
    h = Nu * k / Dh
    return {"rho": rho, "mu": mu, "cp": cp, "k": k,
            "v": v, "Re": Re, "Pr": Pr, "Nu": Nu, "h": h}
