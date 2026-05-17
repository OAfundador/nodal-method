"""
modelo_nodal_reator.py

Construção da rede nodal de um EC tipo placa pelo método nodal,
seguindo os 5 passos:

    Passo 1: Discretização (N camadas axiais, 9 nós cada)
    Passo 2: Tipos de nó (DIFFUSION, ARITHMETIC, FLUID, BOUNDARY)
    Passo 3: Condutâncias G (conduction_func / convection_func / fluid_transport)
    Passo 4: Equações de balanço (auto via NodalNetwork.residual_steady)
    Passo 5: Resolução R(T)=0 (via solver.solve_steady_state)

Topologia por camada axial j ∈ {0, ..., N-1}, com 9 nós:

    Tf_t —Gfi— Ti_t —Gcli— Tcl_t —Gcls— Ts_t —Gchs— Tch
                                                       |
    Tf_b —Gfi— Ti_b —Gcli— Tcl_b —Gcls— Ts_b —Gchs ————┘

Mais 1 nó BOUNDARY (Tch_in) com a temperatura de entrada do fluido.

A fonte axial é Q_face[j] aplicada em CADA nó de combustível (top e bot),
com renormalização para que 2·Σ Q_face = P_placa exatamente.

Distribuição axial possível:
    cos(π/2 · x/L)    - perfil cossenoidal típico de reator com refletor
    constante         - para comparação
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from condutancias import conduction_G
from materiais import criar_aluminio
from nos import (
    NodalNetwork, NodeKind, TransferKind, LinkDirection,
    conduction_func_between_nodes,
)

from reator_placa.geometria_reator import GeometriaReator
from reator_placa.propriedades_agua import (
    cp_agua, rho_agua, coeficiente_pelicula, T_sat_agua,
)
from reator_placa.propriedades_combustivel import (
    criar_combustivel_didatico,
)


MODO_COS = "cos"
MODO_CONST = "constante"


@dataclass(frozen=True)
class ConfigCaso:
    """Parâmetros de um caso de cálculo."""

    n_axial: int = 40
    vazao_canal_m3_s: float = 3.0e-4    # vazão por canal
    P_placa_W: float = 5000.0           # potência da placa modelada
    modo_fluxo: str = MODO_COS
    T_in_C: float = 30.0
    P_in_Pa: float = 150e3
    dP_canal_Pa: float = 10e3
    tol: float = 1e-7
    max_iter: int = 200


def distribuicao_axial_potencia(
    cfg: ConfigCaso, g: GeometriaReator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calcula x[j], z[j], Q_face[j] (potência por meia placa).

    x ∈ [-L/2, +L/2]; z = x + L/2 (desde a entrada).
    Renormalização: 2·Σ Q_face[j] == P_placa exatamente.
    """
    L = g.Lx
    Lz = g.Ly
    dx = L / cfg.n_axial
    j_idx = np.arange(cfg.n_axial)
    x = -L / 2.0 + (j_idx + 0.5) * dx
    z = x + L / 2.0

    if cfg.modo_fluxo == MODO_CONST:
        fator = np.ones_like(x)
    elif cfg.modo_fluxo == MODO_COS:
        fator = np.cos(math.pi / 2.0 * x / L)
    else:
        raise ValueError(f"modo_fluxo inválido: {cfg.modo_fluxo}")

    Ax = dx * Lz
    soma_total = 2.0 * float(np.sum(fator * Ax))   # 2 meias placas
    escala = cfg.P_placa_W / soma_total
    Q_face = fator * escala * Ax
    return x, z, Q_face


def pressao_local(cfg: ConfigCaso, j: int) -> float:
    """Queda linear de P entre entrada e saída."""
    if cfg.n_axial <= 1:
        return cfg.P_in_Pa
    frac = j / (cfg.n_axial - 1)
    return cfg.P_in_Pa - frac * cfg.dP_canal_Pa


@dataclass
class IdsCamada:
    Tf_t: int; Ti_t: int; Tcl_t: int; Ts_t: int
    Tch: int
    Ts_b: int; Tcl_b: int; Ti_b: int; Tf_b: int


@dataclass
class RedeReator:
    net: NodalNetwork
    camadas: list[IdsCamada]
    Tch_in: int
    cfg: ConfigCaso
    geom: GeometriaReator
    Q_face: np.ndarray
    x_axial: np.ndarray
    z_axial: np.ndarray
    m_dot_canal: float
    pressoes: np.ndarray = field(default_factory=lambda: np.array([]))


def _make_h_func(ch_id, P_local, m_dot, area_flow, Dh, area_face):
    """h·A para Ts↔Tch: h avaliado na T do fluido (Dittus-Boelter)."""
    def func(T_map):
        T_ch = T_map[ch_id]
        coef = coeficiente_pelicula(T_ch, P_local, m_dot, area_flow, Dh)
        return coef["h"] * area_face
    return func


def _make_gch_func(up_id, down_id, m_dot, P_local):
    """m_dot·cp para o transporte entálpico entre camadas."""
    def func(T_map):
        Tm = 0.5 * (T_map[up_id] + T_map[down_id])
        return m_dot * cp_agua(Tm, P_local)
    return func


def construir_rede(cfg: ConfigCaso, g: GeometriaReator) -> RedeReator:
    """
    Constrói a NodalNetwork seguindo os passos 1-4.

    A não-linearidade k_f(T) e h(T,P) entra via conductance_func, e o
    Newton em solve_steady_state lida com ela automaticamente.
    """
    combustivel = criar_combustivel_didatico()
    revestimento = criar_aluminio()
    net = NodalNetwork()

    L = g.Lx
    dx = L / cfg.n_axial
    Lz = g.Ly
    Ax = dx * Lz
    Lfi = g.df / 2.0
    Lcli = g.dcl / 2.0

    x_axial, z_axial, Q_face = distribuicao_axial_potencia(cfg, g)

    rho_in = rho_agua(cfg.T_in_C, cfg.P_in_Pa)
    m_dot = cfg.vazao_canal_m3_s * rho_in

    pressoes = np.array([pressao_local(cfg, j) for j in range(cfg.n_axial)])

    G_cli = conduction_G(revestimento.k, Ax, Lcli)
    G_cls = conduction_G(revestimento.k, Ax, Lcli)

    camadas: list[IdsCamada] = []
    T0 = cfg.T_in_C

    for j in range(cfg.n_axial):
        z_j = float(z_axial[j])
        Tf_t = net.add_node(f"Tf_top_{j+1:03d}", NodeKind.DIFFUSION,
                            x=z_j, y=0.0, material=combustivel,
                            volume=(g.df/2.0)*dx*Lz, source=float(Q_face[j]),
                            temperature=T0)
        Ti_t = net.add_node(f"Ti_top_{j+1:03d}", NodeKind.ARITHMETIC,
                            x=z_j, y=1.0, temperature=T0)
        Tcl_t = net.add_node(f"Tcl_top_{j+1:03d}", NodeKind.DIFFUSION,
                             x=z_j, y=2.0, material=revestimento,
                             volume=g.dcl*dx*Lz, temperature=T0)
        Ts_t = net.add_node(f"Ts_top_{j+1:03d}", NodeKind.ARITHMETIC,
                            x=z_j, y=3.0, temperature=T0)
        Tch = net.add_node(f"Tch_{j+1:03d}", NodeKind.FLUID,
                           x=z_j, y=4.0, volume=g.area_flow*dx,
                           temperature=T0)
        Ts_b = net.add_node(f"Ts_bot_{j+1:03d}", NodeKind.ARITHMETIC,
                            x=z_j, y=5.0, temperature=T0)
        Tcl_b = net.add_node(f"Tcl_bot_{j+1:03d}", NodeKind.DIFFUSION,
                             x=z_j, y=6.0, material=revestimento,
                             volume=g.dcl*dx*Lz, temperature=T0)
        Ti_b = net.add_node(f"Ti_bot_{j+1:03d}", NodeKind.ARITHMETIC,
                            x=z_j, y=7.0, temperature=T0)
        Tf_b = net.add_node(f"Tf_bot_{j+1:03d}", NodeKind.DIFFUSION,
                            x=z_j, y=8.0, material=combustivel,
                            volume=(g.df/2.0)*dx*Lz, source=float(Q_face[j]),
                            temperature=T0)
        camadas.append(IdsCamada(Tf_t, Ti_t, Tcl_t, Ts_t, Tch,
                                 Ts_b, Tcl_b, Ti_b, Tf_b))

    Tch_in = net.add_node("Tch_inlet", NodeKind.BOUNDARY,
                          x=-1.0, y=4.0,
                          fixed_temperature=cfg.T_in_C,
                          temperature=cfg.T_in_C)

    for j, c in enumerate(camadas):
        P_j = float(pressoes[j])

        # ramo TOP
        net.add_link(c.Tf_t, c.Ti_t, TransferKind.CONDUCTION,
                     conductance_func=conduction_func_between_nodes(
                         combustivel, c.Tf_t, c.Ti_t,
                         area=Ax, distance=Lfi),
                     name=f"Gfi_top_{j+1:03d}")
        net.add_link(c.Ti_t, c.Tcl_t, TransferKind.CONDUCTION,
                     conductance=G_cli, name=f"Gcli_top_{j+1:03d}")
        net.add_link(c.Tcl_t, c.Ts_t, TransferKind.CONDUCTION,
                     conductance=G_cls, name=f"Gcls_top_{j+1:03d}")
        net.add_link(c.Ts_t, c.Tch, TransferKind.CONVECTION,
                     conductance_func=_make_h_func(
                         c.Tch, P_j, m_dot,
                         area_flow=g.area_flow, Dh=g.Dh, area_face=Ax),
                     name=f"Gchs_top_{j+1:03d}")

        # ramo BOTTOM (espelho)
        net.add_link(c.Tf_b, c.Ti_b, TransferKind.CONDUCTION,
                     conductance_func=conduction_func_between_nodes(
                         combustivel, c.Tf_b, c.Ti_b,
                         area=Ax, distance=Lfi),
                     name=f"Gfi_bot_{j+1:03d}")
        net.add_link(c.Ti_b, c.Tcl_b, TransferKind.CONDUCTION,
                     conductance=G_cli, name=f"Gcli_bot_{j+1:03d}")
        net.add_link(c.Tcl_b, c.Ts_b, TransferKind.CONDUCTION,
                     conductance=G_cls, name=f"Gcls_bot_{j+1:03d}")
        net.add_link(c.Ts_b, c.Tch, TransferKind.CONVECTION,
                     conductance_func=_make_h_func(
                         c.Tch, P_j, m_dot,
                         area_flow=g.area_flow, Dh=g.Dh, area_face=Ax),
                     name=f"Gchs_bot_{j+1:03d}")

    # transporte do fluido (direcional, montante -> jusante)
    P0 = float(pressoes[0])
    net.add_link(Tch_in, camadas[0].Tch, TransferKind.FLUID_TRANSPORT,
                 direction=LinkDirection.I_TO_J,
                 conductance_func=_make_gch_func(
                     Tch_in, camadas[0].Tch, m_dot, P0),
                 name="Gch_inlet")
    for j in range(cfg.n_axial - 1):
        P_j = float(pressoes[j])
        net.add_link(camadas[j].Tch, camadas[j+1].Tch,
                     TransferKind.FLUID_TRANSPORT,
                     direction=LinkDirection.I_TO_J,
                     conductance_func=_make_gch_func(
                         camadas[j].Tch, camadas[j+1].Tch, m_dot, P_j),
                     name=f"Gch_{j+1:03d}_{j+2:03d}")

    return RedeReator(net=net, camadas=camadas, Tch_in=Tch_in,
                      cfg=cfg, geom=g,
                      Q_face=Q_face, x_axial=x_axial, z_axial=z_axial,
                      m_dot_canal=m_dot, pressoes=pressoes)


def gerar_chute_inicial(rede: RedeReator) -> np.ndarray:
    """Marcha axial simplificada para acelerar o Newton."""
    cfg = rede.cfg
    g = rede.geom
    n = cfg.n_axial
    m_dot = rede.m_dot_canal

    Q_layer_total = 2.0 * rede.Q_face
    cp_ref = cp_agua(cfg.T_in_C, cfg.P_in_Pa)
    Tch_estim = np.zeros(n)
    Tup = cfg.T_in_C
    for j in range(n):
        dT = Q_layer_total[j] / (m_dot * cp_ref)
        Tch_estim[j] = Tup + dT
        Tup = Tch_estim[j]

    coef_in = coeficiente_pelicula(
        cfg.T_in_C, cfg.P_in_Pa, m_dot, g.area_flow, g.Dh)
    h_ref = coef_in["h"]
    Ax = (g.Lx / n) * g.Ly
    G_chs = h_ref * Ax
    G_cls = conduction_G(180.0, Ax, g.dcl / 2.0)
    G_cli = conduction_G(180.0, Ax, g.dcl / 2.0)
    G_fi = conduction_G(85.0, Ax, g.df / 2.0)

    z_init = []
    for j in range(n):
        Qh = float(rede.Q_face[j])
        Tch = Tch_estim[j]
        Ts = Tch + Qh / G_chs if G_chs > 0 else Tch
        Tcl = Ts + Qh / G_cls
        Ti = Tcl + Qh / G_cli
        Tf = Ti + Qh / G_fi
        # Ordem dos nós: Tf_t, Ti_t, Tcl_t, Ts_t, Tch, Ts_b, Tcl_b, Ti_b, Tf_b
        z_init.extend([Tf, Ti, Tcl, Ts, Tch, Ts, Tcl, Ti, Tf])
    return np.array(z_init, dtype=float)


def extrair_resultados(rede: RedeReator) -> dict:
    """Empacota os vetores axiais Tf, Ti, Tcl, Ts, Tch + diagnósticos."""
    cfg = rede.cfg
    g = rede.geom
    n = cfg.n_axial

    Tf = np.zeros(n); Ti = np.zeros(n); Tcl = np.zeros(n)
    Ts = np.zeros(n); Tch = np.zeros(n)
    for j, c in enumerate(rede.camadas):
        Tf[j] = rede.net.nodes[c.Tf_t].temperature
        Ti[j] = rede.net.nodes[c.Ti_t].temperature
        Tcl[j] = rede.net.nodes[c.Tcl_t].temperature
        Ts[j] = rede.net.nodes[c.Ts_t].temperature
        Tch[j] = rede.net.nodes[c.Tch].temperature

    h = np.zeros(n); Re = np.zeros(n); Pr = np.zeros(n)
    Tsat = np.zeros(n)
    for j in range(n):
        coef = coeficiente_pelicula(float(Tch[j]), float(rede.pressoes[j]),
                                    rede.m_dot_canal, g.area_flow, g.Dh)
        h[j] = coef["h"]; Re[j] = coef["Re"]; Pr[j] = coef["Pr"]
        Tsat[j] = T_sat_agua(float(rede.pressoes[j]))

    return {
        "cfg": cfg, "geom": g,
        "x": rede.x_axial, "z": rede.z_axial,
        "Q_face": rede.Q_face, "Q_total": 2.0 * rede.Q_face,
        "q_flux_face": rede.Q_face / ((g.Lx / n) * g.Ly),
        "Tf": Tf, "Ti": Ti, "Tcl": Tcl, "Ts": Ts, "Tch": Tch,
        "h": h, "Re": Re, "Pr": Pr, "Tsat": Tsat,
        "m_dot": rede.m_dot_canal, "pressoes": rede.pressoes,
    }
