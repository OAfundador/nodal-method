"""
nodal_repl.py

REPL / runner interativo do Metodo Nodal.
Suporta problemas genericos (chip, placa 2D) e EC tipo placa nuclear.

Toda a fisica do reator e embutida aqui — sem dependencia de reator_placa/.
Usa apenas os 5 modulos core: condutancias, geometria, materiais, nos, solver.

Uso:
  python exemplos/nodal_repl.py                    # modo interativo
  python exemplos/nodal_repl.py --demo             # demo rapido
  python exemplos/nodal_repl.py exemplos/interativo_reator.txt
"""

from __future__ import annotations
import csv as _csv
import math as _math
import shlex
import sys
from dataclasses import dataclass as _dataclass, field as _field
from pathlib import Path

import numpy as _np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from condutancias import conduction_G, convection_G
from geometria import Geometry2D
from materiais import Material, MaterialPhase
from nos import (
    NodalNetwork, NodeKind, TransferKind, LinkDirection,
    build_network_from_geometry, conduction_func_between_nodes,
)
from solver import solve_steady_state


# ===========================================================================
# Fisica do reator tipo placa — inline (sem reator_placa/)
# Usa: condutancias, materiais, nos (modulos core)
# ===========================================================================

# --- Constantes ---
MODO_COS   = "cos"
MODO_CONST = "constante"

# --- Propriedades do refrigerante ---
#
# Toda a fisica do refrigerante (rho, mu, cp, k, T_sat) flui agora pelo
# objeto Material instanciado a partir do dict de propriedades guardado
# no REPL (state.materials) ou, como fallback, do _MATERIAIS_BUILTIN["agua"]
# que contem as correlacoes de Incropera (Tab. A.6) e a equacao de Antoine
# embutidas como *_expr. Veja _refrigerante_da_sessao() e _material_from_props().

def _coef_pelicula(T_C, P_Pa, m_dot, area_flow, Dh, refrigerante):
    """Coeficiente convectivo h por Dittus-Boelter: Nu=0.023*Re^0.8*Pr^0.4.

    Todas as propriedades vem do Material refrigerante (avaliacao em T_C, P_Pa).
    """
    rho = refrigerante.prop("rho", T_C, P_Pa)
    mu  = refrigerante.prop("mu",  T_C, P_Pa)
    cp  = refrigerante.prop("cp",  T_C, P_Pa)
    k   = refrigerante.prop("k",   T_C, P_Pa)
    v   = m_dot / (rho * area_flow)
    Re  = rho * v * Dh / mu
    Pr  = cp * mu / k
    Nu  = 0.023 * Re**0.8 * Pr**0.4
    h   = Nu * k / Dh
    return {"rho": rho, "mu": mu, "cp": cp, "k": k,
            "v": v, "Re": Re, "Pr": Pr, "Nu": Nu, "h": h}

# --- Combustivel U3Si2-Al: k(T) pela correlacao do PDF TNR5703 ---

def _k_u3si2al(T_C, P=None):
    """k [W/(m.K)], T em C.  Original em BTU/h.ft.F, convertido para SI."""
    T_F = 1.8 * float(T_C) + 32.0
    k_BTU = 3978.1 / (692.61 + T_F) + 6.02366e-12 * (T_F + 460.0)**3
    return k_BTU * 1.73073   # 1 BTU/(h.ft.F) = 1.73073 W/(m.K)

def _criar_combustivel_padrao():
    """Material U3Si2-Al com k(T) da correlacao do PDF TNR5703."""
    return Material(name="U3Si2-Al", phase=MaterialPhase.SOLID,
                    k_func=_k_u3si2al, rho=4300.0, cp=836.0)

def _criar_aluminio_revestimento():
    """Al-6061 para revestimento."""
    return Material(name="Al-6061", phase=MaterialPhase.SOLID,
                    k=180.0, rho=2700.0, cp=896.0)

# --- Geometria do EC tipo placa ---

@_dataclass(frozen=True)
class GeometriaReator:
    """Dimensoes transversais de um EC tipo placa."""
    Lx:     float = 0.500    # comprimento ativo [m]
    Ly:     float = 0.060    # largura do cerne [m]
    df:     float = 0.001    # espessura do combustivel [m]
    dcl:    float = 0.0005   # espessura do revestimento [m]
    dch:    float = 0.003    # abertura do canal de resfriamento [m]
    Lcanal: float = 0.065    # comprimento transversal do canal [m]

    @property
    def area_flow(self):
        return self.dch * self.Lcanal

    @property
    def perimetro_molhado(self):
        return 2.0 * (self.dch + self.Lcanal)

    @property
    def Dh(self):
        return 4.0 * self.area_flow / self.perimetro_molhado

# --- Configuracao do caso ---

@_dataclass(frozen=True)
class ConfigCaso:
    """Parametros operacionais do EC."""
    n_axial:          int   = 40
    vazao_canal_m3_s: float = 3.0e-4
    P_placa_W:        float = 5000.0
    modo_fluxo:       str   = MODO_COS
    T_in_C:           float = 30.0
    P_in_Pa:          float = 150e3
    dP_canal_Pa:      float = 10e3
    tol:              float = 1e-7
    max_iter:         int   = 200

# --- Distribuicao axial de potencia ---

def _dist_potencia(cfg, g):
    L = g.Lx; Lz = g.Ly; dx = L / cfg.n_axial
    j_idx = _np.arange(cfg.n_axial)
    x = -L / 2.0 + (j_idx + 0.5) * dx
    z = x + L / 2.0
    if cfg.modo_fluxo == MODO_COS:
        fator = _np.cos(_math.pi / 2.0 * x / L)
    elif cfg.modo_fluxo == MODO_CONST:
        fator = _np.ones_like(x)
    else:
        raise ValueError(f"modo_fluxo invalido: {cfg.modo_fluxo!r}")
    Ax = dx * Lz
    escala = cfg.P_placa_W / (2.0 * float(_np.sum(fator * Ax)))
    return x, z, fator * escala * Ax

def _P_local(cfg, j):
    if cfg.n_axial <= 1:
        return cfg.P_in_Pa
    return cfg.P_in_Pa - j / (cfg.n_axial - 1) * cfg.dP_canal_Pa

# --- Estruturas de dados da rede ---

@_dataclass
class _IdsCamada:
    Tf_t: int; Ti_t: int; Tcl_t: int; Ts_t: int
    Tch:  int
    Ts_b: int; Tcl_b: int; Ti_b: int; Tf_b: int

@_dataclass
class RedeReator:
    net:           NodalNetwork
    camadas:       list
    Tch_in:        int
    cfg:           ConfigCaso
    geom:          GeometriaReator
    Q_face:        object   # np.ndarray
    x_axial:       object
    z_axial:       object
    m_dot_canal:   float
    pressoes:      object = _field(default_factory=lambda: _np.array([]))
    combustivel:   object = None
    refrigerante:  object = None

# --- Funcoes de condutancia dependente de temperatura ---

def _make_h_func(ch_id, P_local, m_dot, area_flow, Dh, area_face, refrigerante):
    def func(T_map):
        coef = _coef_pelicula(T_map[ch_id], P_local, m_dot, area_flow,
                              Dh, refrigerante)
        return coef["h"] * area_face
    return func

def _make_gch_func(up_id, dn_id, m_dot, P_local, refrigerante):
    def func(T_map):
        T_avg = 0.5 * (T_map[up_id] + T_map[dn_id])
        return m_dot * refrigerante.prop("cp", T_avg, P_local)
    return func

# --- Construcao da rede nodal ---

def construir_rede(cfg: ConfigCaso, g: GeometriaReator,
                   combustivel=None, refrigerante=None) -> RedeReator:
    """
    Constroi NodalNetwork de 9 nos/camada axial.

    combustivel:  Material opcional; usa U3Si2-Al (PDF TNR5703) se None.
    refrigerante: Material opcional; usa builtin "agua" com correlacoes
                  de Incropera embutidas se None.
    """
    if combustivel is None:
        combustivel = _criar_combustivel_padrao()
    if refrigerante is None:
        refrigerante = _material_from_props("agua", _MATERIAIS_BUILTIN["agua"])
    revestimento = _criar_aluminio_revestimento()
    net = NodalNetwork()

    dx = g.Lx / cfg.n_axial
    Ax = dx * g.Ly
    Lfi  = g.df  / 2.0
    Lcli = g.dcl / 2.0

    x_axial, z_axial, Q_face = _dist_potencia(cfg, g)
    rho_in  = refrigerante.prop("rho", cfg.T_in_C, cfg.P_in_Pa)
    m_dot   = cfg.vazao_canal_m3_s * rho_in
    pressoes = _np.array([_P_local(cfg, j) for j in range(cfg.n_axial)])

    G_cli = conduction_G(revestimento.k, Ax, Lcli)
    G_cls = conduction_G(revestimento.k, Ax, Lcli)

    camadas = []
    T0 = cfg.T_in_C

    for j in range(cfg.n_axial):
        z_j = float(z_axial[j])
        Tf_t  = net.add_node(f"Tf_top_{j+1:03d}", NodeKind.DIFFUSION,
                             x=z_j, y=0.0, material=combustivel,
                             volume=(g.df/2.0)*dx*g.Ly,
                             source=float(Q_face[j]), temperature=T0)
        Ti_t  = net.add_node(f"Ti_top_{j+1:03d}", NodeKind.ARITHMETIC,
                             x=z_j, y=1.0, temperature=T0)
        Tcl_t = net.add_node(f"Tcl_top_{j+1:03d}", NodeKind.DIFFUSION,
                             x=z_j, y=2.0, material=revestimento,
                             volume=g.dcl*dx*g.Ly, temperature=T0)
        Ts_t  = net.add_node(f"Ts_top_{j+1:03d}", NodeKind.ARITHMETIC,
                             x=z_j, y=3.0, temperature=T0)
        Tch   = net.add_node(f"Tch_{j+1:03d}", NodeKind.FLUID,
                             x=z_j, y=4.0, volume=g.area_flow*dx, temperature=T0)
        Ts_b  = net.add_node(f"Ts_bot_{j+1:03d}", NodeKind.ARITHMETIC,
                             x=z_j, y=5.0, temperature=T0)
        Tcl_b = net.add_node(f"Tcl_bot_{j+1:03d}", NodeKind.DIFFUSION,
                             x=z_j, y=6.0, material=revestimento,
                             volume=g.dcl*dx*g.Ly, temperature=T0)
        Ti_b  = net.add_node(f"Ti_bot_{j+1:03d}", NodeKind.ARITHMETIC,
                             x=z_j, y=7.0, temperature=T0)
        Tf_b  = net.add_node(f"Tf_bot_{j+1:03d}", NodeKind.DIFFUSION,
                             x=z_j, y=8.0, material=combustivel,
                             volume=(g.df/2.0)*dx*g.Ly,
                             source=float(Q_face[j]), temperature=T0)
        camadas.append(_IdsCamada(Tf_t, Ti_t, Tcl_t, Ts_t, Tch,
                                  Ts_b, Tcl_b, Ti_b, Tf_b))

    Tch_in = net.add_node("Tch_inlet", NodeKind.BOUNDARY,
                          x=-1.0, y=4.0,
                          fixed_temperature=cfg.T_in_C,
                          temperature=cfg.T_in_C)

    for j, c in enumerate(camadas):
        P_j = float(pressoes[j])
        # Lado superior
        net.add_link(c.Tf_t,  c.Ti_t,  TransferKind.CONDUCTION,
                     conductance_func=conduction_func_between_nodes(
                         combustivel, c.Tf_t, c.Ti_t, area=Ax, distance=Lfi),
                     name=f"Gfi_top_{j+1:03d}")
        net.add_link(c.Ti_t,  c.Tcl_t, TransferKind.CONDUCTION,
                     conductance=G_cli, name=f"Gcli_top_{j+1:03d}")
        net.add_link(c.Tcl_t, c.Ts_t,  TransferKind.CONDUCTION,
                     conductance=G_cls, name=f"Gcls_top_{j+1:03d}")
        net.add_link(c.Ts_t,  c.Tch,   TransferKind.CONVECTION,
                     conductance_func=_make_h_func(
                         c.Tch, P_j, m_dot, g.area_flow, g.Dh, Ax,
                         refrigerante),
                     name=f"Gchs_top_{j+1:03d}")
        # Lado inferior
        net.add_link(c.Tf_b,  c.Ti_b,  TransferKind.CONDUCTION,
                     conductance_func=conduction_func_between_nodes(
                         combustivel, c.Tf_b, c.Ti_b, area=Ax, distance=Lfi),
                     name=f"Gfi_bot_{j+1:03d}")
        net.add_link(c.Ti_b,  c.Tcl_b, TransferKind.CONDUCTION,
                     conductance=G_cli, name=f"Gcli_bot_{j+1:03d}")
        net.add_link(c.Tcl_b, c.Ts_b,  TransferKind.CONDUCTION,
                     conductance=G_cls, name=f"Gcls_bot_{j+1:03d}")
        net.add_link(c.Ts_b,  c.Tch,   TransferKind.CONVECTION,
                     conductance_func=_make_h_func(
                         c.Tch, P_j, m_dot, g.area_flow, g.Dh, Ax,
                         refrigerante),
                     name=f"Gchs_bot_{j+1:03d}")

    # Transporte de entalpia axial
    P0 = float(pressoes[0])
    net.add_link(Tch_in, camadas[0].Tch, TransferKind.FLUID_TRANSPORT,
                 direction=LinkDirection.I_TO_J,
                 conductance_func=_make_gch_func(Tch_in, camadas[0].Tch,
                                                 m_dot, P0, refrigerante),
                 name="Gch_inlet")
    for j in range(cfg.n_axial - 1):
        net.add_link(camadas[j].Tch, camadas[j+1].Tch,
                     TransferKind.FLUID_TRANSPORT,
                     direction=LinkDirection.I_TO_J,
                     conductance_func=_make_gch_func(
                         camadas[j].Tch, camadas[j+1].Tch,
                         m_dot, float(pressoes[j]), refrigerante),
                     name=f"Gch_{j+1:03d}_{j+2:03d}")

    return RedeReator(net=net, camadas=camadas, Tch_in=Tch_in,
                      cfg=cfg, geom=g, Q_face=Q_face,
                      x_axial=x_axial, z_axial=z_axial,
                      m_dot_canal=m_dot, pressoes=pressoes,
                      combustivel=combustivel, refrigerante=refrigerante)


def gerar_chute_inicial(rede: RedeReator) -> _np.ndarray:
    cfg = rede.cfg; g = rede.geom; n = cfg.n_axial; m_dot = rede.m_dot_canal
    refrigerante = rede.refrigerante
    combustivel  = rede.combustivel
    Q2 = 2.0 * rede.Q_face
    cp_ref = refrigerante.prop("cp", cfg.T_in_C, cfg.P_in_Pa)
    Tch_e  = _np.zeros(n); Tup = cfg.T_in_C
    for j in range(n):
        dT = Q2[j] / (m_dot * cp_ref)
        Tch_e[j] = Tup + dT; Tup = Tch_e[j]
    coef_in = _coef_pelicula(cfg.T_in_C, cfg.P_in_Pa, m_dot,
                             g.area_flow, g.Dh, refrigerante)
    h_ref = coef_in["h"]; Ax = (g.Lx / n) * g.Ly
    G_chs = h_ref * Ax
    G_cls = conduction_G(180.0, Ax, g.dcl / 2.0)
    G_cli = conduction_G(180.0, Ax, g.dcl / 2.0)
    k_fuel_ref = combustivel.prop("k", cfg.T_in_C)
    G_fi  = conduction_G(k_fuel_ref, Ax, g.df / 2.0)
    z_init = []
    for j in range(n):
        Qh  = float(rede.Q_face[j]); Tch = Tch_e[j]
        Ts  = Tch + Qh / G_chs if G_chs > 0 else Tch
        Tcl = Ts  + Qh / G_cls
        Ti  = Tcl + Qh / G_cli
        Tf  = Ti  + Qh / G_fi
        z_init.extend([Tf, Ti, Tcl, Ts, Tch, Ts, Tcl, Ti, Tf])
    return _np.array(z_init, dtype=float)


def extrair_resultados(rede: RedeReator) -> dict:
    cfg = rede.cfg; g = rede.geom; n = cfg.n_axial
    Tf  = _np.zeros(n); Ti  = _np.zeros(n)
    Tcl = _np.zeros(n); Ts  = _np.zeros(n); Tch = _np.zeros(n)
    for j, c in enumerate(rede.camadas):
        Tf[j]  = rede.net.nodes[c.Tf_t].temperature
        Ti[j]  = rede.net.nodes[c.Ti_t].temperature
        Tcl[j] = rede.net.nodes[c.Tcl_t].temperature
        Ts[j]  = rede.net.nodes[c.Ts_t].temperature
        Tch[j] = rede.net.nodes[c.Tch].temperature
    h   = _np.zeros(n); Re  = _np.zeros(n)
    Pr  = _np.zeros(n); Tsat = _np.zeros(n)
    refrigerante = rede.refrigerante
    t_sat_func   = getattr(refrigerante, "t_sat_func", None)
    for j in range(n):
        coef = _coef_pelicula(float(Tch[j]), float(rede.pressoes[j]),
                              rede.m_dot_canal, g.area_flow, g.Dh,
                              refrigerante)
        h[j]    = coef["h"]; Re[j]   = coef["Re"]; Pr[j] = coef["Pr"]
        if t_sat_func is not None:
            Tsat[j] = t_sat_func(float(rede.pressoes[j]))
        else:
            Tsat[j] = float("nan")
    return {
        "cfg": cfg, "geom": g,
        "x": rede.x_axial, "z": rede.z_axial,
        "Q_face": rede.Q_face, "Q_total": 2.0 * rede.Q_face,
        "q_flux_face": rede.Q_face / ((g.Lx / n) * g.Ly),
        "Tf": Tf, "Ti": Ti, "Tcl": Tcl, "Ts": Ts, "Tch": Tch,
        "h": h, "Re": Re, "Pr": Pr, "Tsat": Tsat,
        "m_dot": rede.m_dot_canal, "pressoes": rede.pressoes,
    }


def texto_resumo(res: dict, caso: str) -> str:
    cfg = res["cfg"]; g = res["geom"]
    margem = res["Tsat"] - res["Tch"]
    return (
        f"EC tipo placa — Metodo Nodal\n\n"
        f"Caso = {caso}\n"
        f"Modo de fluxo = {cfg.modo_fluxo}\n"
        f"N axial = {cfg.n_axial}\n"
        f"Vazao por canal = {cfg.vazao_canal_m3_s:.4e} m3/s\n"
        f"m_dot por canal = {res['m_dot']:.4e} kg/s\n\n"
        f"Geometria:\n"
        f"  comprimento ativo Lx = {g.Lx*1e3:.1f} mm\n"
        f"  largura cerne Ly     = {g.Ly*1e3:.1f} mm\n"
        f"  espessura combust.   = {g.df*1e3:.3f} mm\n"
        f"  espessura clad       = {g.dcl*1e3:.3f} mm\n"
        f"  canal interno        = {g.dch*1e3:.2f} x {g.Lcanal*1e3:.1f} mm\n"
        f"  Dh                   = {g.Dh*1e3:.3f} mm\n\n"
        f"Rede nodal:\n"
        f"  n nos   = {res['n_nos']}\n"
        f"  n links = {res['n_links']}\n\n"
        f"Solver:\n"
        f"  success        = {res['solver_success']}\n"
        f"  residual_norm  = {res['solver_residual_norm']:.3e}\n"
        f"  iteracoes      = {res['solver_iterations']}\n\n"
        f"Potencia:\n"
        f"  P_placa configurada = {cfg.P_placa_W:.2f} W\n"
        f"  Potencia integrada  = {float(_np.sum(res['Q_total'])):.2f} W\n\n"
        f"Temperaturas:\n"
        f"  T_fluido_entrada    = {cfg.T_in_C:.2f} C\n"
        f"  T_fluido_saida      = {float(res['Tch'][-1]):.2f} C\n"
        f"  T_superficie_max    = {float(_np.max(res['Ts'])):.2f} C\n"
        f"  T_revestimento_max  = {float(_np.max(res['Tcl'])):.2f} C\n"
        f"  T_interface_max     = {float(_np.max(res['Ti'])):.2f} C\n"
        f"  T_combustivel_max   = {float(_np.max(res['Tf'])):.2f} C\n"
        f"  margem T_sat-T_ch (min) = {float(_np.min(margem)):.2f} C\n\n"
        f"Hidraulica:\n"
        f"  Re_min, Re_max = {float(_np.min(res['Re'])):.1f}, {float(_np.max(res['Re'])):.1f}\n"
        f"  h_min, h_max   = {float(_np.min(res['h'])):.1f}, {float(_np.max(res['h'])):.1f} W/m2K\n"
    )


def salvar_csv(res: dict, path: Path) -> None:
    cfg = res["cfg"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = _csv.writer(f, delimiter=";")
        w.writerow([
            "j", "x_m", "z_m", "q_flux_W_m2", "Q_total_W",
            "T_fluido_C", "T_superficie_C", "T_revestimento_C",
            "T_interface_C", "T_combustivel_max_C",
            "h_W_m2K", "Re", "Pr", "T_sat_C", "margem_sat_C",
        ])
        for j in range(cfg.n_axial):
            w.writerow([
                j+1, float(res["x"][j]), float(res["z"][j]),
                float(res["q_flux_face"][j]), float(res["Q_total"][j]),
                float(res["Tch"][j]), float(res["Ts"][j]),
                float(res["Tcl"][j]), float(res["Ti"][j]), float(res["Tf"][j]),
                float(res["h"][j]), float(res["Re"][j]), float(res["Pr"][j]),
                float(res["Tsat"][j]), float(res["Tsat"][j] - res["Tch"][j]),
            ])


def plot_temperaturas(res: dict, path: Path, titulo: str) -> None:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    z_cm = res["z"] * 100.0
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    ax.plot(z_cm, res["Tch"], label="fluido")
    ax.plot(z_cm, res["Ts"],  label="superficie")
    ax.plot(z_cm, res["Tcl"], label="centro do revestimento")
    ax.plot(z_cm, res["Ti"],  label="interface combustivel/revestimento")
    ax.plot(z_cm, res["Tf"],  label="centro do combustivel (T max.)")
    ax.set_xlabel("posicao axial desde a entrada [cm]")
    ax.set_ylabel("Temperatura [C]"); ax.set_title(titulo)
    ax.grid(True, alpha=0.3); ax.legend(loc="upper left")
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)


def plot_fluxo(res: dict, path: Path, titulo: str) -> None:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4.7))
    ax.plot(res["x"] * 100.0, res["q_flux_face"] / 1e6)
    ax.set_xlabel("posicao axial x (centro=0) [cm]")
    ax.set_ylabel("q'' por face [MW/m2]"); ax.set_title(titulo)
    ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)


def _make_expr_func(expr_str):
    """Compila uma expressao em T (e opcionalmente P) num callable f(T,P=None)."""
    _ns = {"math": _math, **vars(_math)}
    return lambda T, P=None, _e=expr_str, _n=_ns: float(
        eval(_e, {"T": T, "P": P, **_n}))


def _make_t_sat_func(expr_str):
    """Compila uma expressao em P (Pa) num callable f(P_Pa) -> T_sat (C)."""
    _ns = {"math": _math, **vars(_math)}
    return lambda P, _e=expr_str, _n=_ns: float(
        eval(_e, {"P": P, **_n}))


def _material_from_props(nome, props):
    """Cria Material a partir do dict de propriedades armazenado no REPL.

    Para cada propriedade (k, rho, cp, mu), aceita:
      - <prop>_expr : string de expressao em T (e opcionalmente P) -> *_func
      - <prop>_poly : lista de coeficientes [a0, a1, ...]           -> *_coeffs
      - <prop>      : constante                                      -> escalar

    Se "t_sat_expr" estiver presente (expressao em P em Pa), e atachado
    como mat.t_sat_func — usado em extrair_resultados para a margem ate
    a saturacao.
    """
    ph = (MaterialPhase.FLUID if props.get("phase", "solid") == "fluid"
          else MaterialPhase.SOLID)

    kwargs = {"name": nome, "phase": ph}
    for prop in ("k", "rho", "cp", "mu"):
        expr = props.get(prop + "_expr")
        poly = props.get(prop + "_poly")
        const = props.get(prop)
        if expr:
            kwargs[prop + "_func"] = _make_expr_func(expr)
        elif poly:
            kwargs[prop + "_coeffs"] = poly
        elif const is not None:
            kwargs[prop] = const

    mat = Material(**kwargs)

    t_sat_expr = props.get("t_sat_expr")
    if t_sat_expr:
        mat.t_sat_func = _make_t_sat_func(t_sat_expr)
    else:
        mat.t_sat_func = None

    return mat


# ===========================================================================
# REPL
# ===========================================================================

NODE_KIND_MAP = {
    "diffusion": NodeKind.DIFFUSION,
    "arithmetic": NodeKind.ARITHMETIC,
    "fluid": NodeKind.FLUID,
    "boundary": NodeKind.BOUNDARY,
}
LINK_KIND_MAP = {
    "cond": TransferKind.CONDUCTION,
    "conv": TransferKind.CONVECTION,
    "fluid": TransferKind.FLUID_TRANSPORT,
    "rad": TransferKind.RADIATION,
}

HELP = """
============================================================
  METODO NODAL - REPL INTERATIVO
  Comandos disponiveis
============================================================

--- MATERIAIS ---
  material <nome> [phase=solid|fluid] [k=v] [rho=v] [cp=v] [mu=v]
                  [k_poly=a0,a1,...] [k_expr="expressao em T"]
       k pode ser constante, polinomio ou expressao de T. Exemplos:
         material cobre  k=380 rho=8900 cp=385
         material agua   phase=fluid k=0.6 rho=997 cp=4182 mu=8.9e-4
         material uo2    k_poly=10.5,-0.012,4e-6 rho=10970 cp=330
         material u3si2  k_expr="1.73073*(3978.1/(724.61+1.8*T)+6.02366e-12*(1.8*T+492)**3)" rho=4300 cp=836
  materiais        Lista os materiais disponiveis (builtin + sessao)

--- GEOMETRIA 2D ---
  domain W=v H=v
       Cria dominio fisico largura W x altura H [m]
  region <nome> [mat=<material> | k=v] x0=v x1=v y0=v y1=v
       Define regiao retangular com material
  source <regiao> q=v
       Fonte volumetrica q_triplo [W/m3] na regiao
  bc <face> <tipo> [T=v] [h=v] [T_inf=v]
       face: bottom | top | left | right
       tipo: temperature | convection | adiabatic
  mesh nx=v ny=v [fix=auto|none]
  show_geom        Resume dominio, regioes, BCs e malha
  build_from_geom  [tz=v]   Gera rede nodal automaticamente

--- REDE NODAL DIRETA ---
  node <nome> <tipo> [Q=v] [T_fixed=v] [T_init=v] [V=v]
       tipos: diffusion | fluid | arithmetic | boundary
  link <a> <b> <tipo> G=v [dir=fwd|bwd|undirected]
       tipos: cond | conv | fluid | rad | equiv
  fluid_chain n1 n2 [n3 ...] mdot=v cp=v [T_in=v]
       Cadeia FLUID_TRANSPORT direcional n1->n2->...
       T_in cria no BOUNDARY de entrada automaticamente.
  show             Resumo da rede (nos e ligacoes)
  solve            Resolve estado permanente (Newton-Raphson)
  reset            Limpa rede, geometria e materiais

--- CONDUTANCIAS ---
  g_cond k=v A=v L=v   -> G = k*A/L  [W/K]
  g_conv h=v A=v        -> G = h*A    [W/K]

--- VISUALIZACAO ---
  viz  [png=arquivo.png] [titulo=Texto]
       PNG da rede: nos coloridos por T, links com setas
  gif  [gif=arquivo.gif] [fps=10] [steps=40] [dpi=100] [titulo=Texto]
       GIF pseudo-transiente: frio -> regime permanente
  nucleo_init [rows=5] [cols=5]
       Inicializa matriz vazia do nucleo para entrada manual
  nucleo_linha row=N vals=v1,v2,... [tipos=EC,CR,...]
       Define uma linha da matriz (row=0 e linha de topo)
  nucleo_mapa [png=mapa.png] [titulo=Texto]
       Gera PNG do mapa de potencia radial (usa matriz inserida ou padrao 5x5)
  nucleo_solve [P_total=5e6] [n_canais=19] [png=mapa_T.png] [plots=no] [dir=saidas_nucleo]
       Resolve todos os ECs e gera mapa 5x5 de temperatura maxima do combustivel
  ec_associacao [n_int=17] [n_ext=2] [dch_int=0.00289] [dch_ext=0.00452]
                [potencia=P] [vazao=V] [png=ec_assoc.png]
       Distribui vazao entre canais internos e externos (mesmo dP, Blasius)

--- REATOR TIPO PLACA ---
  reator_config [n_axial=40] [vazao=3e-4] [potencia=5000]
                [modo=cos|constante] [T_in=30] [P_in=150000] [dP=10000]
  reator_geom   [Lx=0.5] [Ly=0.06] [df=0.001] [dcl=0.0005]
                [dch=0.003] [Lcanal=0.065]
  reator_solve / reator_show / reator_plot [dir=.]
  reator_csv [dir=.] / reator_gif [dir=.] [fps=12] [dpi=100]

  Nota: o combustivel padrao e U3Si2-Al com k(T) do PDF TNR5703.
  Para usar outro combustivel, defina via "material" e use o nome
  "u3si2", "combustivel" ou "fuel" antes de reator_solve.

--- CONTROLE ---
  help / quit / exit
============================================================
"""


class EstadoREPL:
    def __init__(self):
        self.net = NodalNetwork()
        self.ids = {}
        self.x_counter = 0.0
        self.geom = None
        self.materials = {}
        self.reator_cfg = None
        self.reator_geom = None
        self.reator_result = None
        self.nucleo_data = None

    def reset(self):
        self.__init__()


def parse_kwargs(tokens, lower_keys=True):
    out = {}
    for t in tokens:
        if "=" not in t:
            raise ValueError(f"opcao invalida: {t!r}")
        k, v = t.split("=", 1)
        out[k.strip().lower() if lower_keys else k.strip()] = v.strip()
    return out

def fnum(s): return float(str(s).replace(",", "."))

def cmd_node(state, args):
    nome, tipo = args[0], args[1].lower()
    if tipo not in NODE_KIND_MAP: raise ValueError(f"tipo invalido: {tipo!r}")
    if nome in state.ids: raise ValueError(f"no {nome!r} ja existe.")
    kw = parse_kwargs(args[2:])
    kind = NODE_KIND_MAP[tipo]
    if kind is NodeKind.BOUNDARY:
        Tf = fnum(kw["t_fixed"])
        nid = state.net.add_node(name=nome, kind=kind, x=state.x_counter, y=0.0,
                                  fixed_temperature=Tf, temperature=Tf)
    else:
        nid = state.net.add_node(name=nome, kind=kind, x=state.x_counter, y=0.0,
                                  volume=fnum(kw.get("v", "1")),
                                  source=fnum(kw.get("q", "0")),
                                  temperature=fnum(kw.get("t_init", "30")))
    state.ids[nome] = nid; state.x_counter += 1.0
    return f"  -> no {nome} criado (id={nid}, kind={kind.value})"

def cmd_link(state, args):
    a, b, tipo = args[0], args[1], args[2].lower()
    kw = parse_kwargs(args[3:])
    dir_map = {
        "fwd": LinkDirection.I_TO_J, "i_to_j": LinkDirection.I_TO_J,
        "bwd": LinkDirection.J_TO_I, "j_to_i": LinkDirection.J_TO_I,
        "undirected": LinkDirection.UNDIRECTED,
    }
    dir_kw    = kw.get("dir", "undirected").lower()
    direction = dir_map.get(dir_kw, LinkDirection.UNDIRECTED)
    G_val     = fnum(kw["g"])
    state.net.add_link(state.ids[a], state.ids[b], LINK_KIND_MAP[tipo],
                       conductance=G_val, direction=direction,
                       name=a + "_" + b)
    return "  -> ligacao " + a + "-" + b + f" ({tipo}) G={G_val:.6g} W/K  dir=" + dir_kw

def cmd_show(state): return state.net.summary(max_nodes=20, max_links=30)

def cmd_solve(state):
    if len(state.net.nodes) == 0: return "  rede vazia."
    r = solve_steady_state(state.net, tol=1e-9, max_iter=500)
    out = [f"  convergiu={r.success}, |R|={r.residual_norm:.2e}, iter={r.iterations}"]
    if state.ids:
        out += ["", "  Temperaturas:"] + [
            f"    {n:8s} = {state.net.nodes[i].temperature:9.4f} C"
            for n, i in state.ids.items()
        ]
    return "\n".join(out)

def _req(state):
    if not state.geom: raise ValueError("Use domain W= H= primeiro.")
    return state.geom


# ---------------------------------------------------------------------------
# Materiais
# ---------------------------------------------------------------------------
# Builtin "agua" carrega as correlacoes de Incropera (Tab. A.6) e a equacao
# de Antoine como expressoes em T (C) e P (Pa). Quando o usuario faz
#   material agua phase=fluid k=0.620 rho=994 cp=4182 mu=7.7e-4
# cmd_material derruba os respectivos *_expr e usa as constantes do usuario,
# mas mantem t_sat_expr (saturacao continua disponivel para a margem termica).
_MATERIAIS_BUILTIN = {
    "cobre":    dict(phase="solid",  k=380.0,   rho=8900.0,  cp=385.0),
    "aluminio": dict(phase="solid",  k=205.0,   rho=2700.0,  cp=900.0),
    "aco":      dict(phase="solid",  k=50.0,    rho=7800.0,  cp=500.0),
    "silicio":  dict(phase="solid",  k=150.0,   rho=2330.0,  cp=712.0),
    "uo2":      dict(phase="solid",  k=3.6,     rho=10970.0, cp=247.0),
    "agua":     dict(
        phase="fluid",
        # Incropera Tab. A.6 (T em C, k em W/(m K))
        k_expr="0.5615 + 1.939e-3*T - 7.51e-6*T**2",
        # Polinomio de Kell para densidade da agua liquida (T em C, rho em kg/m3)
        rho_expr=("999.842594 + 6.793952e-2*T - 9.095290e-3*T**2 "
                  "+ 1.001685e-4*T**3 - 1.120083e-6*T**4 + 6.536332e-9*T**5"),
        # Calor especifico (T em C, cp em J/(kg K))
        cp_expr="4217.0 - 3.358*T + 0.04148*T**2 - 1.6e-4*T**3",
        # Viscosidade dinamica (Reichardt) (T em C, mu em Pa.s)
        mu_expr="2.414e-5 * 10.0**(247.8 / (T + 273.15 - 140.0))",
        # Saturacao por Antoine (P em Pa, T_sat em C)
        t_sat_expr="1730.63 / (8.07131 - log10(P/133.322)) - 233.426",
    ),
    "ar":       dict(phase="fluid",  k=0.026,   rho=1.18,    cp=1007.0, mu=1.85e-5),
}

def cmd_material(state, args):
    if not args:
        raise ValueError("material precisa de nome. Ex: material cobre k=380")
    nome = args[0].lower()
    kw   = parse_kwargs(args[1:]) if len(args) > 1 else {}
    # Merge sobre o builtin: usuario sobrepoe campos individuais.
    props = dict(_MATERIAIS_BUILTIN.get(nome, {}))

    # Constantes explicitas derrubam expr/poly herdados para a mesma propriedade.
    for campo in ("k", "rho", "cp", "mu"):
        if campo in kw:
            props[campo] = fnum(kw[campo])
            props.pop(campo + "_expr", None)
            props.pop(campo + "_poly", None)

    # Polinomios derrubam constante e expr da mesma propriedade.
    for campo in ("k_poly", "rho_poly", "cp_poly", "mu_poly"):
        if campo in kw:
            props[campo] = [fnum(v) for v in kw[campo].split(",")]
            base = campo[:-5]
            props.pop(base, None)
            props.pop(base + "_expr", None)

    # Expressoes derrubam constante e poly da mesma propriedade.
    for campo in ("k_expr", "rho_expr", "cp_expr", "mu_expr"):
        if campo in kw:
            props[campo] = kw[campo]
            base = campo[:-5]
            props.pop(base, None)
            props.pop(base + "_poly", None)

    # Saturacao (fluidos): expressao em P (Pa).
    if "t_sat_expr" in kw:
        props["t_sat_expr"] = kw["t_sat_expr"]

    if "phase" in kw:
        props["phase"] = kw["phase"].lower()
    if not props.get("phase"):
        props["phase"] = "solid"

    state.materials[nome] = props

    partes = []
    if "k_expr" in props:
        partes.append("k=f(T)[" + props["k_expr"][:40] + "]")
    elif "k_poly" in props:
        partes.append("k=poly" + str(props["k_poly"]))
    elif "k" in props:
        partes.append(f"k={props['k']} W/(mK)")
    for campo, unid in [("rho", "kg/m3"), ("cp", "J/(kgK)"), ("mu", "Pa.s")]:
        if campo + "_expr" in props:
            partes.append(campo + "=f(T)")
        elif campo + "_poly" in props:
            partes.append(campo + "=poly")
        elif campo in props:
            partes.append(campo + "=" + str(props[campo]) + " " + unid)
    if "t_sat_expr" in props:
        partes.append("T_sat=f(P)")
    return "  -> material " + repr(nome) + " [" + props["phase"] + "]  " + "  ".join(partes)

def _resumo_propriedade_k(d):
    if "k_expr" in d:
        return "k=f(T)"
    if "k_poly" in d:
        return "k=poly"
    if "k" in d:
        return "k=" + str(d["k"])
    return "k=-"

def cmd_materiais(state, args):
    linhas = ["  --- Builtin ---"]
    for n, d in _MATERIAIS_BUILTIN.items():
        linhas.append("    " + n.ljust(10) + " [" + d["phase"] + "]  "
                      + _resumo_propriedade_k(d))
    if state.materials:
        linhas.append("  --- Definidos na sessao ---")
        for n, d in state.materials.items():
            linhas.append("    " + n.ljust(10) + " [" + d.get("phase", "solid") + "]  "
                          + _resumo_propriedade_k(d))
    return "\n".join(linhas)

def cmd_domain(state, args):
    kw = parse_kwargs(args)
    state.geom = Geometry2D(width=fnum(kw["w"]), height=fnum(kw["h"]))
    return f"  -> Geometry2D {kw['w']} x {kw['h']} m"

def cmd_region(state, args):
    geom = _req(state); kw = parse_kwargs(args[1:]); nome = args[0]
    if "mat" in kw:
        mat_nome = kw["mat"].lower()
        # cmd_material ja merge builtin+usuario; preferir o dict da sessao
        # quando existe (preserva a derrubada de *_expr feita por overrides).
        if mat_nome in state.materials:
            props = dict(state.materials[mat_nome])
        else:
            props = dict(_MATERIAIS_BUILTIN.get(mat_nome, {}))
        if not props:
            raise ValueError("material " + repr(mat_nome) + " nao definido.")
        mat = _material_from_props(mat_nome, props)
    else:
        k_val = fnum(kw["k"])
        mat   = Material(name="mat_" + nome, phase=MaterialPhase.SOLID, k=k_val)
    geom.material(nome, mat,
                  x0=fnum(kw["x0"]), x1=fnum(kw["x1"]),
                  y0=fnum(kw["y0"]), y1=fnum(kw["y1"]))
    try:
        k_show = mat.prop("k", 30.0)
        k_str  = "{:.6g}".format(k_show)
    except Exception:
        k_str = "f(T)" if mat.k_func is not None else "?"
    return "  -> regiao " + repr(nome) + " k=" + k_str

def cmd_source(state, args):
    geom = _req(state); kw = parse_kwargs(args[1:])
    geom.source(region=args[0], kind="volumetric", value=fnum(kw["q"]))
    return f"  -> fonte {kw['q']} W/m3 em {args[0]!r}"

def cmd_bc(state, args):
    geom = _req(state); kw = parse_kwargs(args[2:], lower_keys=False)
    geom.bc(args[0].lower(), args[1].lower(), **{k: fnum(v) for k, v in kw.items()})
    return f"  -> BC {args[1]!r} em {args[0]!r}"

def cmd_mesh(state, args):
    geom = _req(state); kw = parse_kwargs(args)
    geom.mesh(nx=int(fnum(kw["nx"])), ny=int(fnum(kw["ny"])), fix=kw.get("fix", "auto"))
    return f"  -> malha nx={geom.nx} ny={geom.ny}"

def cmd_build_from_geom(state, args):
    geom = _req(state); kw = parse_kwargs(args) if args else {}
    state.net = build_network_from_geometry(geom, thickness_z=fnum(kw.get("tz", "1")))
    state.ids.clear()
    return f"  -> NodalNetwork: {len(state.net.nodes)} nos, {len(state.net.links)} ligacoes"


# ---------------------------------------------------------------------------
# Comandos do reator tipo placa
# ---------------------------------------------------------------------------

def _combustivel_da_sessao(state):
    """Retorna Material de combustivel se definido na sessao, ou None."""
    for nome_comb in ("u3si2", "u3si2al", "combustivel", "fuel", "uo2"):
        if nome_comb in state.materials:
            return _material_from_props(nome_comb, state.materials[nome_comb])
    return None   # usa padrao (k(T) do PDF)


def _refrigerante_da_sessao(state):
    """Retorna Material do refrigerante.

    Procura por 'agua', 'refrigerante' ou 'coolant' em state.materials.
    Se nao houver definicao na sessao, usa o builtin 'agua' — que carrega
    as correlacoes de Incropera (k, rho, cp, mu) e a equacao de Antoine
    (T_sat) embutidas como *_expr. Nunca retorna None: construir_rede sempre
    recebe um Material valido.
    """
    for nome in ("agua", "refrigerante", "coolant"):
        if nome in state.materials:
            return _material_from_props(nome, state.materials[nome])
    return _material_from_props("agua", _MATERIAIS_BUILTIN["agua"])

def cmd_reator_config(state, args):
    kw = parse_kwargs(args); b = state.reator_cfg or ConfigCaso()
    n  = int(fnum(kw["n_axial"])) if "n_axial"  in kw else b.n_axial
    v  = fnum(kw["vazao"])        if "vazao"     in kw else b.vazao_canal_m3_s
    p  = fnum(kw["potencia"])     if "potencia"  in kw else b.P_placa_W
    Ti = fnum(kw["t_in"])         if "t_in"      in kw else b.T_in_C
    Pi = fnum(kw["p_in"])         if "p_in"      in kw else b.P_in_Pa
    dP = fnum(kw["dp"])           if "dp"        in kw else b.dP_canal_Pa
    tl = fnum(kw["tol"])          if "tol"       in kw else b.tol
    it = int(fnum(kw["iter"]))    if "iter"      in kw else b.max_iter
    mr = kw.get("modo", b.modo_fluxo).lower()
    if mr in ("cos", "cossenoidal"):
        modo = MODO_COS
    elif mr in ("constante", "const"):
        modo = MODO_CONST
    else:
        raise ValueError(f"modo invalido: {mr!r}")
    state.reator_cfg = ConfigCaso(
        n_axial=n, vazao_canal_m3_s=v, P_placa_W=p,
        modo_fluxo=modo, T_in_C=Ti, P_in_Pa=Pi,
        dP_canal_Pa=dP, tol=tl, max_iter=it)
    state.reator_result = None
    return (f"  -> ConfigCaso: n_axial={n}, vazao={v:.3e} m3/s, "
            f"potencia={p:.1f} W, modo={modo}, T_in={Ti} C")

def cmd_reator_geom(state, args):
    kw = parse_kwargs(args); b = state.reator_geom or GeometriaReator()
    state.reator_geom = GeometriaReator(
        Lx    =fnum(kw["lx"])     if "lx"     in kw else b.Lx,
        Ly    =fnum(kw["ly"])     if "ly"     in kw else b.Ly,
        df    =fnum(kw["df"])     if "df"     in kw else b.df,
        dcl   =fnum(kw["dcl"])    if "dcl"    in kw else b.dcl,
        dch   =fnum(kw["dch"])    if "dch"    in kw else b.dch,
        Lcanal=fnum(kw["lcanal"]) if "lcanal" in kw else b.Lcanal)
    g = state.reator_geom; state.reator_result = None
    return (f"  -> GeometriaReator: Lx={g.Lx*1e3:.1f}mm df={g.df*1e3:.2f}mm "
            f"dcl={g.dcl*1e3:.2f}mm dch={g.dch*1e3:.2f}mm Dh={g.Dh*1e3:.3f}mm")

def cmd_reator_solve(state):
    cfg          = state.reator_cfg  or ConfigCaso()
    g            = state.reator_geom or GeometriaReator()
    combustivel  = _combustivel_da_sessao(state)   # None = usa padrao k(T)
    refrigerante = _refrigerante_da_sessao(state)  # nunca None
    rede = construir_rede(cfg, g, combustivel=combustivel,
                          refrigerante=refrigerante)
    z0   = gerar_chute_inicial(rede)
    sol  = solve_steady_state(rede.net, z0=z0, tol=cfg.tol, max_iter=cfg.max_iter,
                              update_network=True, prefer_scipy=True)
    res  = extrair_resultados(rede)
    res["solver_success"]       = bool(sol.success)
    res["solver_residual_norm"] = float(sol.residual_norm)
    res["solver_iterations"]    = sol.iterations
    res["n_nos"]                = len(rede.net.nodes)
    res["n_links"]              = len(rede.net.links)
    state.reator_result = res
    Tch_out = float(res["Tch"][-1])
    Tf_max  = float(_np.max(res["Tf"]))
    comb_label = (combustivel.name if combustivel else "U3Si2-Al (padrao)")
    refr_label = refrigerante.name
    return (f"  -> convergiu={sol.success}, |R|={sol.residual_norm:.2e}, "
            f"iter={sol.iterations}\n"
            f"  -> combustivel={comb_label}\n"
            f"  -> refrigerante={refr_label}\n"
            f"  -> Tch_saida={Tch_out:.2f} C, Tf_max={Tf_max:.2f} C")

def cmd_reator_show(state):
    if state.reator_result is None: return "  sem resultado. Execute reator_solve."
    return texto_resumo(state.reator_result, state.reator_result["cfg"].modo_fluxo)

def cmd_reator_plot(state, args):
    if state.reator_result is None: return "  sem resultado. Execute reator_solve."
    kw     = parse_kwargs(args) if args else {}
    outdir = Path(kw.get("dir", "saidas_reator_placa"))
    outdir.mkdir(parents=True, exist_ok=True)
    res  = state.reator_result; modo = res["cfg"].modo_fluxo
    plot_temperaturas(res, outdir / f"temperaturas_{modo}.png",
                      f"EC placa q'' {modo}")
    plot_fluxo(res,        outdir / f"fluxo_{modo}.png",
               f"Fluxo axial q'' {modo}")
    return f"  -> salvo em {outdir}: temperaturas_{modo}.png, fluxo_{modo}.png"

def cmd_reator_csv(state, args):
    if state.reator_result is None: return "  sem resultado. Execute reator_solve."
    kw     = parse_kwargs(args) if args else {}
    outdir = Path(kw.get("dir", "saidas_reator_placa"))
    outdir.mkdir(parents=True, exist_ok=True)
    res  = state.reator_result; modo = res["cfg"].modo_fluxo
    path = outdir / f"resultado_{modo}.csv"
    salvar_csv(res, path)
    return f"  -> CSV salvo em {path}"


def cmd_reator_gif(state, args):
    """GIF animado da rede nodal: campo 2D de temperatura com varredura axial."""
    if state.reator_result is None:
        return "  sem resultado. Execute reator_solve primeiro."
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from matplotlib.colors import Normalize

    kw     = parse_kwargs(args) if args else {}
    outdir = Path(kw.get("dir", "saidas_reator_placa"))
    outdir.mkdir(parents=True, exist_ok=True)
    fps  = int(fnum(kw.get("fps", "12")))
    dpi  = int(fnum(kw.get("dpi", "100")))
    res  = state.reator_result
    cfg  = res["cfg"]; modo = cfg.modo_fluxo; N = cfg.n_axial
    z    = res["z"] * 100.0

    y_labels = ["comb.(inf)", "interf.(inf)", "clad(inf)", "sup.(inf)",
                "fluido", "sup.(sup)", "clad(sup)", "interf.(sup)", "comb.(sup)"]
    T2d = _np.array([res["Tf"], res["Ti"], res["Tcl"], res["Ts"], res["Tch"],
                     res["Ts"], res["Tcl"], res["Ti"],  res["Tf"]])

    T_lo, T_hi = T2d.min(), T2d.max()
    norm = Normalize(vmin=T_lo, vmax=T_hi)
    cmap = plt.cm.hot

    fig   = plt.figure(figsize=(13, 5))
    ax_map  = fig.add_axes([0.06, 0.12, 0.58, 0.78])
    ax_prof = fig.add_axes([0.70, 0.12, 0.18, 0.78])
    ax_cb   = fig.add_axes([0.90, 0.12, 0.025, 0.78])

    img = ax_map.imshow(T2d, aspect="auto", origin="lower", cmap=cmap, norm=norm,
                        extent=[z[0], z[-1], -0.5, 8.5], interpolation="nearest")
    ax_map.set_xlabel("Posicao axial desde a entrada [cm]", fontsize=9)
    ax_map.set_yticks(range(9)); ax_map.set_yticklabels(y_labels, fontsize=7)
    ax_map.set_title(f"Campo nodal  |  q'' {modo}  |  {N} camadas", fontsize=9)
    plt.colorbar(img, cax=ax_cb, label="T [C]")

    vline = ax_map.axvline(x=z[0], color="cyan", lw=1.8, ls="--", alpha=0.85)
    ax_prof.set_xlim(T_lo - 0.5, T_hi + 0.5); ax_prof.set_ylim(-0.5, 8.5)
    ax_prof.set_yticks(range(9)); ax_prof.set_yticklabels(y_labels, fontsize=7)
    ax_prof.set_xlabel("T [C]", fontsize=8); ax_prof.set_title("Perfil\nz=--", fontsize=8)
    ax_prof.grid(True, alpha=0.3)
    pts_prof  = ax_prof.scatter(T2d[:, 0], range(9), c=T2d[:, 0],
                                cmap=cmap, norm=norm, s=60, zorder=5)
    line_prof, = ax_prof.plot(T2d[:, 0], range(9), color="gray", lw=1, alpha=0.6)
    txt_tch = ax_map.text(0.02, 0.97, "", transform=ax_map.transAxes,
                          fontsize=8, va="top", color="cyan",
                          bbox=dict(boxstyle="round,pad=0.2", fc="black", alpha=0.5))

    n_frames  = min(N, 60)
    frame_idx = _np.linspace(0, N-1, n_frames, dtype=int)
    frame_idx = _np.concatenate([frame_idx, _np.full(6, N-1, dtype=int)])

    def update(fi):
        j   = frame_idx[fi]
        vline.set_xdata([z[j]])
        col = T2d[:, j]
        pts_prof.set_offsets(_np.column_stack([col, range(9)]))
        pts_prof.set_array(col)
        line_prof.set_xdata(col)
        ax_prof.set_title("Perfil\nz=" + f"{z[j]:.1f}" + " cm", fontsize=8)
        txt_tch.set_text("Tch=" + f"{float(res['Tch'][j]):.2f}" + " C")
        return vline, pts_prof, line_prof, txt_tch

    ani = animation.FuncAnimation(fig, update, frames=len(frame_idx),
                                  interval=1000 // fps, blit=True)
    gif_path = outdir / f"animacao_{modo}.gif"
    ani.save(str(gif_path), writer="pillow", dpi=dpi)
    plt.close(fig)
    return f"  -> GIF salvo em {gif_path}  ({len(frame_idx)} frames, {fps} fps)"


# ===========================================================================
# Transporte entalpico generico
# ===========================================================================

def cmd_fluid_chain(state, args):
    """fluid_chain n1 n2 ... mdot=v cp=v [T_in=v]"""
    nos_names = [a for a in args if "=" not in a]
    kw        = parse_kwargs([a for a in args if "=" in a])
    if len(nos_names) < 2:
        raise ValueError("fluid_chain precisa de pelo menos 2 nos")
    if "mdot" not in kw or "cp" not in kw:
        raise ValueError("fluid_chain precisa de mdot=v e cp=v")
    mdot = fnum(kw["mdot"]); cp = fnum(kw["cp"]); G = mdot * cp
    T_in_val = fnum(kw["t_in"]) if "t_in" in kw else None

    def resolve(name):
        if name in state.ids: return state.ids[name]
        try:
            nid = int(name)
            if nid in state.net.nodes: return nid
        except ValueError:
            pass
        raise ValueError("No nao encontrado: " + repr(name))

    chain = [resolve(n) for n in nos_names]; msgs = []
    if T_in_val is not None:
        fn = state.net.nodes[chain[0]]
        inlet_id = state.net.add_node(
            name="inlet_" + nos_names[0], kind=NodeKind.BOUNDARY,
            x=fn.x - 1.0, y=fn.y,
            fixed_temperature=T_in_val, temperature=T_in_val)
        state.net.add_link(inlet_id, chain[0], TransferKind.FLUID_TRANSPORT,
                           direction=LinkDirection.I_TO_J, conductance=G,
                           name="fluid_inlet_" + nos_names[0])
        msgs.append("  -> inlet BOUNDARY T=" + str(T_in_val) + " C -> " + nos_names[0])

    for k in range(len(chain) - 1):
        na, nb = chain[k], chain[k+1]
        state.net.add_link(na, nb, TransferKind.FLUID_TRANSPORT,
                           direction=LinkDirection.I_TO_J, conductance=G,
                           name="fluid_" + nos_names[k] + "_" + nos_names[k+1])
        msgs.append("  -> " + nos_names[k] + " -> " + nos_names[k+1]
                    + "  G=" + str(round(G, 4)) + " W/K")
    return "\n".join(msgs)


# ===========================================================================
# Visualizacao generica
# ===========================================================================

def _render_network(net, title="Rede Nodal", fig=None, ax=None, alpha_links=0.55):
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    ids   = list(net.nodes.keys())
    xs    = _np.array([net.nodes[i].x for i in ids])
    ys    = _np.array([net.nodes[i].y for i in ids])
    Ts    = _np.array([net.nodes[i].temperature for i in ids])
    kinds = [net.nodes[i].kind for i in ids]

    T_lo = float(Ts.min()); T_hi = float(Ts.max())
    if abs(T_hi - T_lo) < 0.01: T_hi = T_lo + 1.0
    norm = Normalize(vmin=T_lo, vmax=T_hi); cmap = plt.cm.hot

    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(11, 6))

    link_style = {
        TransferKind.CONDUCTION:      ("#888888", "-",  0.9),
        TransferKind.CONVECTION:      ("#4499ff", "--", 1.0),
        TransferKind.FLUID_TRANSPORT: ("#ff4444", "-",  1.4),
        TransferKind.RADIATION:       ("#cc44ff", ":",  1.0),
        TransferKind.EQUIVALENT:      ("#aaaaaa", "-.", 0.8),
    }
    id_pos = {i: (net.nodes[i].x, net.nodes[i].y) for i in ids}
    for lk in net.links:
        xi, yi = id_pos[lk.node_i]; xj, yj = id_pos[lk.node_j]
        col, ls, lw = link_style.get(lk.kind, ("#888888", "-", 0.9))
        ax.plot([xi, xj], [yi, yj], color=col, ls=ls, lw=lw,
                alpha=alpha_links, zorder=1)
        if lk.direction != LinkDirection.UNDIRECTED:
            mx = 0.5*(xi+xj); my = 0.5*(yi+yj)
            dx = xj-xi; dy = yj-yi
            if lk.direction == LinkDirection.J_TO_I: dx=-dx; dy=-dy
            nd = (dx**2+dy**2)**0.5
            if nd > 1e-9:
                ax.annotate("",
                    xy=(mx+dx/nd*0.12, my+dy/nd*0.12), xytext=(mx, my),
                    arrowprops=dict(arrowstyle="->", color=col, lw=1.3), zorder=2)

    mk_map = {"diffusion":("o",90),"fluid":("s",90),"arithmetic":("^",55),"boundary":("D",65)}
    for kv in set(k.value for k in kinds):
        mask = [k.value == kv for k in kinds]; idx = [i for i,m in enumerate(mask) if m]
        xm = xs[idx]; ym = ys[idx]; Tm = Ts[idx]
        mk, sz = mk_map.get(kv, ("o", 80))
        ax.scatter(xm, ym, c=Tm, cmap=cmap, norm=norm,
                   marker=mk, s=sz, edgecolors="black", lw=0.6, zorder=4, label=kv)

    sm = ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    plt.colorbar(sm, ax=ax, label="T [degC]", shrink=0.85)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title(title + "  |  " + str(len(ids)) + " nos  " + str(len(net.links)) + " links")
    ax.legend(loc="upper right", fontsize=7, framealpha=0.6)
    ax.grid(True, alpha=0.18)
    return fig, ax, norm, cmap


def cmd_viz(state, args):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if len(state.net.nodes) == 0: return "  rede vazia."
    kw = parse_kwargs(args) if args else {}
    outfile = kw.get("png", "rede_nodal.png"); titulo = kw.get("titulo", "Rede Nodal")
    fig, ax, _, _ = _render_network(state.net, title=titulo)
    outpath = Path(outfile); outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(outpath), dpi=130, bbox_inches="tight"); plt.close(fig)
    return "  -> viz salvo em " + str(outpath)


def cmd_gif_generico(state, args):
    """GIF pseudo-transiente: frio -> regime permanente."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    if len(state.net.nodes) == 0: return "  rede vazia."
    kw      = parse_kwargs(args) if args else {}
    outfile = kw.get("gif", "rede_nodal.gif")
    fps     = int(fnum(kw.get("fps", "10")))
    dpi     = int(fnum(kw.get("dpi", "100")))
    n_steps = int(fnum(kw.get("steps", "40")))
    titulo  = kw.get("titulo", "Rede Nodal")

    net = state.net; ids = list(net.nodes.keys())
    xs  = _np.array([net.nodes[i].x for i in ids])
    ys  = _np.array([net.nodes[i].y for i in ids])
    T_final = _np.array([net.nodes[i].temperature for i in ids])
    kinds   = [net.nodes[i].kind for i in ids]
    T_lo = float(T_final.min()); T_hi = float(T_final.max())
    if abs(T_hi - T_lo) < 0.01: T_hi = T_lo + 1.0
    norm = Normalize(vmin=T_lo, vmax=T_hi); cmap = plt.cm.hot

    is_boundary = _np.array([net.nodes[i].is_boundary() for i in ids])
    T_start = _np.where(is_boundary, T_final, T_lo)
    alphas  = _np.concatenate([_np.linspace(0.0, 1.0, n_steps), _np.ones(6)])

    fig, ax = plt.subplots(figsize=(11, 6))
    link_style = {
        TransferKind.CONDUCTION:      ("#888888", "-",  0.8),
        TransferKind.CONVECTION:      ("#4499ff", "--", 0.9),
        TransferKind.FLUID_TRANSPORT: ("#ff4444", "-",  1.3),
        TransferKind.RADIATION:       ("#cc44ff", ":",  0.9),
        TransferKind.EQUIVALENT:      ("#aaaaaa", "-.", 0.7),
    }
    id_pos = {i: (net.nodes[i].x, net.nodes[i].y) for i in ids}
    for lk in net.links:
        xi, yi = id_pos[lk.node_i]; xj, yj = id_pos[lk.node_j]
        col, ls, lw = link_style.get(lk.kind, ("#888888", "-", 0.8))
        ax.plot([xi, xj], [yi, yj], color=col, ls=ls, lw=lw, alpha=0.45, zorder=1)
        if lk.direction != LinkDirection.UNDIRECTED:
            mx=0.5*(xi+xj); my=0.5*(yi+yj); dx=xj-xi; dy=yj-yi
            if lk.direction == LinkDirection.J_TO_I: dx=-dx; dy=-dy
            nd=(dx**2+dy**2)**0.5
            if nd > 1e-9:
                ax.annotate("", xy=(mx+dx/nd*0.12, my+dy/nd*0.12), xytext=(mx,my),
                    arrowprops=dict(arrowstyle="->", color=col, lw=1.2), zorder=2)

    mk_map = {"diffusion":"o","fluid":"s","arithmetic":"^","boundary":"D"}
    sz_map = {"diffusion":90,"fluid":90,"arithmetic":55,"boundary":65}
    scs    = {}
    for kv in set(k.value for k in kinds):
        mask = _np.array([k.value == kv for k in kinds])
        scs[kv] = ax.scatter(xs[mask], ys[mask], c=T_start[mask], cmap=cmap, norm=norm,
                              marker=mk_map.get(kv,"o"), s=sz_map.get(kv,80),
                              edgecolors="black", lw=0.6, zorder=4, label=kv)

    sm = ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    plt.colorbar(sm, ax=ax, label="T [degC]", shrink=0.85)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_title(titulo)
    ax.legend(loc="upper right", fontsize=7, framealpha=0.6)
    ax.grid(True, alpha=0.18)
    pct_txt = ax.text(0.02, 0.97, "0%", transform=ax.transAxes, fontsize=9, va="top",
                      bbox=dict(boxstyle="round", fc="white", alpha=0.7))

    def update(fi):
        a = float(alphas[fi]); T_cur = T_start + a*(T_final - T_start)
        for kv in scs:
            mask = _np.array([k.value == kv for k in kinds])
            scs[kv].set_array(T_cur[mask])
        pct_txt.set_text(str(int(round(a*100))) + "%   Tmax=" + f"{float(T_cur.max()):.1f}" + " C")
        return list(scs.values()) + [pct_txt]

    ani = animation.FuncAnimation(fig, update, frames=len(alphas),
                                  interval=1000//fps, blit=True)
    outpath = Path(outfile); outpath.parent.mkdir(parents=True, exist_ok=True)
    ani.save(str(outpath), writer="pillow", dpi=dpi); plt.close(fig)
    return "  -> GIF salvo em " + str(outpath) + "  (" + str(len(alphas)) + " frames, " + str(fps) + " fps)"


# ---------------------------------------------------------------------------
# Mapa de potencia do nucleo
# ---------------------------------------------------------------------------
_MAPA_PADRAO = [
    [1.321, 1.563, 0.981, 1.628, 1.030],
    [0.857, 0.515, 1.129, 0.402, 0.826],
    [1.050, 1.877, 0.000, 1.914, 0.979],
    [0.860, 0.411, 1.151, 0.519, 0.822],
    [0.906, 1.028, 0.878, 1.044, 0.867],
]
_TIPOS_PADRAO = [
    ["EC","EC","EC","EC","EC"],
    ["EC","CR","EC","CR","EC"],
    ["EC","EC","CR","EC","EC"],
    ["EC","CR","EC","CR","EC"],
    ["EC","EC","EC","EC","EC"],
]


def _solve_one_canal(geom, cfg_template, P_canal, combustivel=None,
                     refrigerante=None):
    cfg = ConfigCaso(
        n_axial          = cfg_template.n_axial,
        vazao_canal_m3_s = cfg_template.vazao_canal_m3_s,
        P_placa_W        = P_canal,
        modo_fluxo       = cfg_template.modo_fluxo,
        T_in_C           = cfg_template.T_in_C,
        P_in_Pa          = cfg_template.P_in_Pa,
        dP_canal_Pa      = cfg_template.dP_canal_Pa,
        tol              = cfg_template.tol,
        max_iter         = cfg_template.max_iter,
    )
    rede = construir_rede(cfg, geom, combustivel=combustivel,
                          refrigerante=refrigerante)
    z0   = gerar_chute_inicial(rede)
    solve_steady_state(rede.net, z0=z0, tol=cfg.tol, max_iter=cfg.max_iter,
                       update_network=True, prefer_scipy=True)
    return extrair_resultados(rede)


def _gerar_mapa_temperaturas(mapa_fat, tipos, T_fuel, nrows, ncols, outfile, titulo):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    vals = [T_fuel[i][j] for i in range(nrows) for j in range(ncols)
            if T_fuel[i][j] == T_fuel[i][j]]
    T_min_v = min(vals) if vals else 30.0; T_max_v = max(vals) if vals else 100.0
    pos_max = (0,0); T_abs_max = 0.0
    for i in range(nrows):
        for j in range(ncols):
            v = T_fuel[i][j]
            if v == v and v > T_abs_max: T_abs_max = v; pos_max = (i,j)
    fig, ax = plt.subplots(figsize=(7,7))
    ax.set_xlim(0,ncols); ax.set_ylim(0,nrows); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(titulo, fontsize=12, fontweight="bold", pad=14)
    cmap = plt.cm.plasma
    for i in range(nrows):
        for j in range(ncols):
            row_plot = nrows-1-i; tipo=tipos[i][j]; fat=mapa_fat[i][j]; T=T_fuel[i][j]
            if fat == 0.0:
                face="#c0c0c0"; txt="---"; tcol="#333"
            elif tipo == "CR":
                face="#8ab4d4"; txt="CR"; tcol="white"
            elif T != T:
                face="#eeeeee"; txt="ERR"; tcol="#c00"
            else:
                norm = max(0.0, min(1.0,(T-T_min_v)/(T_max_v-T_min_v+1e-12)))
                rgba = cmap(norm)
                face = "#{:02x}{:02x}{:02x}".format(int(rgba[0]*255),int(rgba[1]*255),int(rgba[2]*255))
                txt = f"{T:.1f}"; tcol = "white" if norm>0.5 else "#222"
            ax.add_patch(Rectangle((j,row_plot),1,1,facecolor=face,edgecolor="white",linewidth=2.5))
            if (i,j) == pos_max:
                ax.add_patch(Rectangle((j+0.04,row_plot+0.04),0.92,0.92,
                    fill=False,edgecolor="#cc0000",linewidth=3.0))
            ax.text(j+0.5,row_plot+0.5,txt,ha="center",va="center",
                    fontsize=11,fontweight="bold",color=tcol)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=T_min_v,vmax=T_max_v))
    sm.set_array([]); cbar=fig.colorbar(sm,ax=ax,fraction=0.035,pad=0.04)
    cbar.set_label("T combustivel max [C]",fontsize=9)
    fig.tight_layout()
    outpath = Path(outfile); outpath.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(str(outpath),dpi=150,bbox_inches="tight"); plt.close(fig)


def cmd_nucleo_solve(state, args):
    if state.nucleo_data is None:
        raise ValueError("use nucleo_init + nucleo_linha para definir o mapa.")
    if state.reator_geom is None:
        raise ValueError("use reator_geom para definir a geometria.")
    if state.reator_cfg is None:
        raise ValueError("use reator_config para definir a configuracao.")
    kw       = parse_kwargs(args) if args else {}
    P_total  = fnum(kw.get("p_total",  "5e6"))
    n_canais = int(fnum(kw.get("n_canais", "19")))
    outfile  = kw.get("png", "saidas/mapa_temperaturas.png")
    do_plots = kw.get("plots", "no").lower() in ("yes","sim","1","true")
    plotdir  = Path(kw.get("dir", "saidas_nucleo"))
    mapa  = state.nucleo_data["vals"]; tipos = state.nucleo_data["tipos"]
    nrows = state.nucleo_data["rows"]; ncols = state.nucleo_data["cols"]
    n_ec  = sum(1 for i in range(nrows) for j in range(ncols)
                if tipos[i][j]=="EC" and mapa[i][j]>0.0)
    if n_ec == 0: raise ValueError("nenhum EC ativo no mapa.")
    P_base_canal = P_total / n_ec / n_canais
    geom = state.reator_geom; cfg = state.reator_cfg
    combustivel  = _combustivel_da_sessao(state)
    refrigerante = _refrigerante_da_sessao(state)
    T_fuel = [[float("nan")]*ncols for _ in range(nrows)]
    if do_plots:
        plotdir.mkdir(parents=True, exist_ok=True)
    log = []
    for i in range(nrows):
        for j in range(ncols):
            fator = mapa[i][j]; tipo = tipos[i][j]
            if tipo != "EC" or fator <= 0.0: continue
            try:
                res  = _solve_one_canal(geom, cfg, fator*P_base_canal,
                                        combustivel, refrigerante)
                Tf   = float(_np.max(res["Tf"]))
                Ts   = float(_np.max(res["Ts"]))
                Tch  = float(res["Tch"][-1])
                T_fuel[i][j] = Tf
                log.append(f"  EC({i+1},{j+1}) f={fator:.3f}  Tf_max={Tf:.1f}C  "
                            f"Ts_max={Ts:.1f}C  Tch_out={Tch:.1f}C")
                if do_plots:
                    plot_temperaturas(res, plotdir/f"ec_{i+1}_{j+1}.png",
                                      f"EC ({i+1},{j+1}) fator={fator:.3f}")
            except Exception as e:
                log.append(f"  EC({i+1},{j+1}) ERRO: {e}")
    _gerar_mapa_temperaturas(mapa, tipos, T_fuel, nrows, ncols, outfile,
                             "T max combustivel por EC [C]")
    log.append(f"  -> mapa salvo em {outfile}")
    if do_plots: log.append(f"  -> plots em {plotdir}/")
    gif_file = kw.get("gif", "")
    if gif_file:
        r = _gif_nucleo_scan(mapa, tipos, T_fuel, nrows, ncols, gif_file,
                             fps=int(fnum(kw.get("fps","6"))))
        log.append(r)
    return "\n".join(log)


def _plot_ec_associacao(res_int, res_ext, n_int, n_ext, outfile, titulo):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12,5))
    fig.suptitle(titulo, fontsize=12, fontweight="bold")
    for ax, res, label in [
        (axes[0], res_int, f"Canal interno  (x{n_int})"),
        (axes[1], res_ext, f"Canal externo  (x{n_ext})"),
    ]:
        z = res["z"] * 100
        ax.plot(z, res["Tf"],  "r-",  lw=2,   label="T combustivel")
        ax.plot(z, res["Ti"],  "m--", lw=1.5, label="T interface")
        ax.plot(z, res["Tcl"], "b--", lw=1.5, label="T revestimento")
        ax.plot(z, res["Ts"],  "g-",  lw=1.5, label="T superficie")
        ax.plot(z, res["Tch"], "c-",  lw=2,   label="T fluido")
        ax.set_xlabel("Posicao axial z [cm]"); ax.set_ylabel("Temperatura [C]")
        ax.set_title(label); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    outpath = Path(outfile); outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(outpath), dpi=130, bbox_inches="tight"); plt.close(fig)


def cmd_ec_associacao(state, args):
    if state.reator_geom is None:
        raise ValueError("use reator_geom primeiro.")
    if state.reator_cfg is None:
        raise ValueError("use reator_config primeiro.")
    kw      = parse_kwargs(args) if args else {}
    n_int   = int(fnum(kw.get("n_int",   "17")))
    n_ext   = int(fnum(kw.get("n_ext",    "2")))
    dch_int = fnum(kw.get("dch_int", str(state.reator_geom.dch)))
    dch_ext = fnum(kw.get("dch_ext", "0.00452"))
    Lcanal  = fnum(kw.get("lcanal",  str(state.reator_geom.Lcanal)))
    P_EC    = fnum(kw.get("potencia", str(state.reator_cfg.P_placa_W)))
    V_EC    = fnum(kw.get("vazao",    str(state.reator_cfg.vazao_canal_m3_s)))
    outfile = kw.get("png", "saidas/ec_associacao.png")
    geom    = state.reator_geom; cfg = state.reator_cfg
    combustivel  = _combustivel_da_sessao(state)
    refrigerante = _refrigerante_da_sessao(state)

    def dh_canal(dch, Lc):
        A = dch*Lc; P = 2*(dch+Lc); return 4*A/P, A

    Dh_int, A_int = dh_canal(dch_int, Lcanal)
    Dh_ext, A_ext = dh_canal(dch_ext, Lcanal)
    w_int = A_int * Dh_int**(5/7)
    w_ext = A_ext * Dh_ext**(5/7)
    W_tot = n_int*w_int + n_ext*w_ext
    V_int = V_EC * (w_int/W_tot); V_ext = V_EC * (w_ext/W_tot)
    n_pl  = n_int + 1
    P_int = P_EC / n_pl; P_ext = P_EC / (2*n_pl)

    linhas = [
        f"  Canais internos: {n_int}  dch={dch_int*1e3:.2f}mm  Dh={Dh_int*1e3:.3f}mm",
        f"  Canais externos: {n_ext}  dch={dch_ext*1e3:.2f}mm  Dh={Dh_ext*1e3:.3f}mm",
        f"  Vazao  V_int={V_int*1e6:.2f}cm3/s ({V_int/V_EC*100:.1f}%)"
        f"  V_ext={V_ext*1e6:.2f}cm3/s ({V_ext/V_EC*100:.1f}%)",
        f"  Potencia  P_int={P_int:.1f}W/canal  P_ext={P_ext:.1f}W/canal",
        f"  Resolvendo...",
    ]
    geom_int = GeometriaReator(Lx=geom.Lx, Ly=geom.Ly, df=geom.df, dcl=geom.dcl,
                                dch=dch_int, Lcanal=Lcanal)
    geom_ext = GeometriaReator(Lx=geom.Lx, Ly=geom.Ly, df=geom.df, dcl=geom.dcl,
                                dch=dch_ext, Lcanal=Lcanal)
    cfg_int = ConfigCaso(n_axial=cfg.n_axial, vazao_canal_m3_s=V_int, P_placa_W=P_int,
                         modo_fluxo=cfg.modo_fluxo, T_in_C=cfg.T_in_C, P_in_Pa=cfg.P_in_Pa,
                         dP_canal_Pa=cfg.dP_canal_Pa, tol=cfg.tol, max_iter=cfg.max_iter)
    cfg_ext = ConfigCaso(n_axial=cfg.n_axial, vazao_canal_m3_s=V_ext, P_placa_W=P_ext,
                         modo_fluxo=cfg.modo_fluxo, T_in_C=cfg.T_in_C, P_in_Pa=cfg.P_in_Pa,
                         dP_canal_Pa=cfg.dP_canal_Pa, tol=cfg.tol, max_iter=cfg.max_iter)
    res_int = _solve_one_canal(geom_int, cfg_int, P_int, combustivel, refrigerante)
    res_ext = _solve_one_canal(geom_ext, cfg_ext, P_ext, combustivel, refrigerante)
    Tf_i  = float(_np.max(res_int["Tf"])); Ts_i = float(_np.max(res_int["Ts"]))
    Tch_i = float(res_int["Tch"][-1])
    Tf_e  = float(_np.max(res_ext["Tf"])); Ts_e = float(_np.max(res_ext["Ts"]))
    Tch_e = float(res_ext["Tch"][-1])
    linhas += [
        f"  Canal INTERNO: Tf_max={Tf_i:.2f}C   Ts_max={Ts_i:.2f}C   Tch_saida={Tch_i:.2f}C",
        f"  Canal EXTERNO: Tf_max={Tf_e:.2f}C   Ts_max={Ts_e:.2f}C   Tch_saida={Tch_e:.2f}C",
    ]
    _plot_ec_associacao(res_int, res_ext, n_int, n_ext, outfile,
                        f"Associacao de canais  P={P_EC:.0f}W  V={V_EC*1e6:.0f}cm3/s")
    linhas.append(f"  -> figura salva em {outfile}")
    gif_file = kw.get("gif", "")
    if gif_file:
        r = _gif_ec_associacao(res_int, res_ext, n_int, n_ext, gif_file,
                               fps=int(fnum(kw.get("fps","10"))),
                               steps=int(fnum(kw.get("steps","40"))))
        linhas.append(r)
    state.reator_result = res_int
    return "\n".join(linhas)


def _gif_nucleo_scan(mapa_fat, tipos, T_fuel, nrows, ncols, outfile, fps=6):
    """GIF: revela o mapa de temperaturas EC por EC."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from matplotlib.patches import Rectangle

    vals = [T_fuel[i][j] for i in range(nrows) for j in range(ncols)
            if T_fuel[i][j] == T_fuel[i][j]]
    T_min_v = min(vals) if vals else 30.0; T_max_v = max(vals) if vals else 100.0
    pos_max=(0,0); T_abs=0.0
    for i in range(nrows):
        for j in range(ncols):
            v=T_fuel[i][j]
            if v==v and v>T_abs: T_abs=v; pos_max=(i,j)
    ordem = [(i,j) for i in range(nrows) for j in range(ncols)
             if tipos[i][j]=="EC" and mapa_fat[i][j]>0.0]
    cmap = plt.cm.plasma

    def draw_frame(ax, revelados):
        ax.cla(); ax.set_xlim(0,ncols); ax.set_ylim(0,nrows)
        ax.set_aspect("equal"); ax.axis("off")
        for i in range(nrows):
            for j in range(ncols):
                row_plot=nrows-1-i; fat=mapa_fat[i][j]; tipo=tipos[i][j]; T=T_fuel[i][j]
                if fat==0.0: face="#c0c0c0"; txt="---"; tcol="#333"
                elif tipo=="CR": face="#8ab4d4"; txt="CR"; tcol="white"
                elif (i,j) not in revelados: face="#f0f0f0"; txt=f"{fat:.3f}"; tcol="#aaa"
                elif T!=T: face="#eeeeee"; txt="ERR"; tcol="#c00"
                else:
                    norm=max(0.0,min(1.0,(T-T_min_v)/(T_max_v-T_min_v+1e-12)))
                    rgba=cmap(norm)
                    face="#{:02x}{:02x}{:02x}".format(int(rgba[0]*255),int(rgba[1]*255),int(rgba[2]*255))
                    txt=f"{T:.1f}"; tcol="white" if norm>0.5 else "#222"
                ax.add_patch(Rectangle((j,row_plot),1,1,facecolor=face,edgecolor="white",linewidth=2))
                if (i,j)==pos_max and (i,j) in revelados:
                    ax.add_patch(Rectangle((j+0.04,row_plot+0.04),0.92,0.92,
                        fill=False,edgecolor="#cc0000",linewidth=3))
                ax.text(j+0.5,row_plot+0.5,txt,ha="center",va="center",
                        fontsize=10,fontweight="bold",color=tcol)

    fig, ax = plt.subplots(figsize=(6,6))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=T_min_v,vmax=T_max_v))
    sm.set_array([]); fig.colorbar(sm,ax=ax,fraction=0.03,pad=0.03).set_label("Tf max [C]",fontsize=8)
    frames=[set()]; rev=set()
    for pos in ordem: rev=rev|{pos}; frames.append(frozenset(rev))
    for _ in range(3): frames.append(frozenset(rev))

    def animate(k):
        ax.set_title(f"Mapa de T combustivel — EC {k}/{len(ordem)}", fontsize=11, fontweight="bold")
        draw_frame(ax, frames[k]); return []

    ani = animation.FuncAnimation(fig, animate, frames=len(frames),
                                  interval=1000//fps, blit=False)
    outpath = Path(outfile); outpath.parent.mkdir(parents=True, exist_ok=True)
    ani.save(str(outpath), writer="pillow", fps=fps); plt.close(fig)
    return f"  -> GIF nucleo salvo em {outpath}  ({len(frames)} frames, {fps} fps)"


def _gif_ec_associacao(res_int, res_ext, n_int, n_ext, outfile, fps=10, steps=40):
    """GIF pseudo-transiente dos dois tipos de canal."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    T_in   = res_int["Tch"][0]
    alphas = _np.concatenate([_np.linspace(0.0, 1.0, steps), _np.ones(6)])
    fig, axes = plt.subplots(1,2,figsize=(12,5))
    fig.suptitle("Associacao de canais — evolucao termica", fontsize=12, fontweight="bold")
    curves = [("Tf","r-",2.0,"T combustivel"),("Ti","m--",1.5,"T interface"),
              ("Tcl","b--",1.5,"T revestimento"),("Ts","g-",1.5,"T superficie"),
              ("Tch","c-",2.0,"T fluido")]
    all_T = [v for res in (res_int,res_ext) for k,_,_,_ in curves for v in res[k].tolist()]
    T_lo = min(all_T)-1; T_hi = max(all_T)+1
    lines_int={}; lines_ext={}
    for ax, res, lns, label in [
        (axes[0],res_int,lines_int,f"Canal interno  (x{n_int})"),
        (axes[1],res_ext,lines_ext,f"Canal externo  (x{n_ext})"),
    ]:
        z=res["z"]*100; ax.set_xlim(z[0],z[-1]); ax.set_ylim(T_lo,T_hi)
        ax.set_xlabel("z [cm]"); ax.set_ylabel("T [C]")
        ax.set_title(label); ax.grid(True,alpha=0.3)
        for k,ls,lw,lab in curves:
            ln,=ax.plot([],[],ls,lw=lw,label=lab); lns[k]=ln
        ax.legend(fontsize=8)

    def animate(frame_idx):
        alpha=alphas[frame_idx]
        for res,lns in [(res_int,lines_int),(res_ext,lines_ext)]:
            z=res["z"]*100
            for k,_,_,_ in curves:
                lns[k].set_data(z, T_in+alpha*(res[k]-T_in))
        fig.suptitle(f"Associacao de canais — evolucao termica  [{int(alpha*100)}%]",
                     fontsize=12,fontweight="bold")
        return list(lines_int.values())+list(lines_ext.values())

    ani=animation.FuncAnimation(fig,animate,frames=len(alphas),interval=1000//fps,blit=True)
    outpath=Path(outfile); outpath.parent.mkdir(parents=True,exist_ok=True)
    ani.save(str(outpath),writer="pillow",fps=fps); plt.close(fig)
    return f"  -> GIF associacao salvo em {outpath}  ({len(alphas)} frames, {fps} fps)"


def cmd_nucleo_init(state, args):
    kw=parse_kwargs(args) if args else {}
    rows=int(fnum(kw.get("rows","5"))); cols=int(fnum(kw.get("cols","5")))
    state.nucleo_data = {"rows":rows,"cols":cols,
                         "vals":[[0.0]*cols for _ in range(rows)],
                         "tipos":[["EC"]*cols for _ in range(rows)]}
    return f"  nucleo {rows}x{cols} inicializado — use nucleo_linha para preencher."

def cmd_nucleo_linha(state, args):
    if state.nucleo_data is None: raise ValueError("use nucleo_init antes de nucleo_linha.")
    kw=parse_kwargs(args); row=int(fnum(kw["row"]))
    vals=[fnum(v) for v in kw["vals"].split(",")]; cols=state.nucleo_data["cols"]
    if len(vals)!=cols: raise ValueError(f"esperados {cols} valores, recebidos {len(vals)}.")
    tipos=kw["tipos"].split(",") if "tipos" in kw else ["EC"]*cols
    if len(tipos)!=cols: raise ValueError(f"esperados {cols} tipos, recebidos {len(tipos)}.")
    state.nucleo_data["vals"][row]=vals; state.nucleo_data["tipos"][row]=tipos
    return "  linha " + str(row) + ": " + "  ".join(f"{v:.3f}" for v in vals)

def cmd_nucleo_mapa(state, args):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Patch

    kw=parse_kwargs(args) if args else {}
    outfile=kw.get("png","saidas/mapa_nucleo.png")
    titulo=kw.get("titulo","Distribuicao de Potencia Radial do Nucleo")
    if state.nucleo_data is not None:
        mapa=state.nucleo_data["vals"]; tipos=state.nucleo_data["tipos"]
    else:
        mapa=_MAPA_PADRAO; tipos=_TIPOS_PADRAO
    nrows,ncols=len(mapa),len(mapa[0])
    val_max=0.0; pos_max=(0,0)
    for i in range(nrows):
        for j in range(ncols):
            if mapa[i][j]>val_max: val_max=mapa[i][j]; pos_max=(i,j)
    fig,ax=plt.subplots(figsize=(7,7))
    ax.set_xlim(0,ncols); ax.set_ylim(0,nrows); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(titulo,fontsize=13,fontweight="bold",pad=14)
    for i in range(nrows):
        for j in range(ncols):
            val=mapa[i][j]; tipo=tipos[i][j]; row_plot=nrows-1-i
            if val==0.0: face="#c0c0c0"; txt_col="#333333"
            elif tipo=="CR": face="#8ab4d4"; txt_col="white"
            else:
                intensity=0.55+0.45*(val/val_max)
                r=int(255*intensity); g2=int(160*(1-0.4*(val/val_max))); b=int(130*(1-0.5*(val/val_max)))
                face="#{:02x}{:02x}{:02x}".format(min(r,255),max(g2,0),max(b,0)); txt_col="#222222"
            ax.add_patch(Rectangle((j,row_plot),1,1,facecolor=face,edgecolor="white",linewidth=2.5))
            if (i,j)==pos_max:
                ax.add_patch(Rectangle((j+0.04,row_plot+0.04),0.92,0.92,
                    fill=False,edgecolor="#cc0000",linewidth=3.0))
            ax.text(j+0.5,row_plot+0.5,f"{val:.3f}" if val>0 else "---",
                    ha="center",va="center",fontsize=13,fontweight="bold",color=txt_col)
    legenda=[Patch(facecolor="#d97b5a",edgecolor="white",label="EC combustivel"),
             Patch(facecolor="#8ab4d4",edgecolor="white",label="Controle/Refletor"),
             Patch(facecolor="#c0c0c0",edgecolor="white",label="Barra inserida (q=0)"),
             Patch(facecolor="white",edgecolor="#cc0000",linewidth=2,label="EC mais quente")]
    ax.legend(handles=legenda,loc="lower center",bbox_to_anchor=(0.5,-0.06),
              ncol=2,fontsize=9,framealpha=0.8)
    fig.tight_layout()
    outpath=Path(outfile); outpath.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(str(outpath),dpi=150,bbox_inches="tight"); plt.close(fig)
    return ("  -> mapa salvo em " + str(outpath)
            + "  | EC mais quente: (" + str(pos_max[0]+1) + "," + str(pos_max[1]+1)
            + ") fator=" + str(val_max))


# ---------------------------------------------------------------------------
# Dispatcher principal
# ---------------------------------------------------------------------------

def executar(state, linha):
    raw = linha.strip()
    if not raw or raw.startswith("#"): return "", False
    tokens = shlex.split(raw); cmd = tokens[0].lower(); args = tokens[1:]
    if cmd in ("quit","exit","q"): return "Saindo.", True
    cmds = {
        "help":            lambda: HELP,
        "reset":           lambda: state.reset() or "  rede zerada.",
        "show":            lambda: cmd_show(state),
        "solve":           lambda: cmd_solve(state),
        "node":            lambda: cmd_node(state, args),
        "link":            lambda: cmd_link(state, args),
        "g_cond":          lambda: cmd_gcond(args),
        "g_conv":          lambda: cmd_gconv(args),
        "material":        lambda: cmd_material(state, args),
        "materiais":       lambda: cmd_materiais(state, args),
        "domain":          lambda: cmd_domain(state, args),
        "region":          lambda: cmd_region(state, args),
        "source":          lambda: cmd_source(state, args),
        "bc":              lambda: cmd_bc(state, args),
        "mesh":            lambda: cmd_mesh(state, args),
        "show_geom":       lambda: (
            state.geom.summary() if state.geom
            else (
                "GeometriaReator:\n"
                + "  Lx="     + str(round(state.reator_geom.Lx*1e3,2))     + " mm"
                + "  Ly="     + str(round(state.reator_geom.Ly*1e3,2))     + " mm"
                + "  df="     + str(round(state.reator_geom.df*1e3,3))     + " mm"
                + "  dcl="    + str(round(state.reator_geom.dcl*1e3,3))    + " mm"
                + "  dch="    + str(round(state.reator_geom.dch*1e3,3))    + " mm"
                + "  Lcanal=" + str(round(state.reator_geom.Lcanal*1e3,2)) + " mm"
                + "  Dh="     + str(round(state.reator_geom.Dh*1e3,4))     + " mm"
                + "  A_flow=" + str(round(state.reator_geom.area_flow*1e6,4))+ " mm2"
            ) if state.reator_geom
            else "sem geometria."
        ),
        "build_from_geom": lambda: cmd_build_from_geom(state, args),
        "reator_config":   lambda: cmd_reator_config(state, args),
        "reator_geom":     lambda: cmd_reator_geom(state, args),
        "reator_solve":    lambda: cmd_reator_solve(state),
        "reator_show":     lambda: cmd_reator_show(state),
        "reator_plot":     lambda: cmd_reator_plot(state, args),
        "reator_csv":      lambda: cmd_reator_csv(state, args),
        "reator_gif":      lambda: cmd_reator_gif(state, args),
        "fluid_chain":     lambda: cmd_fluid_chain(state, args),
        "viz":             lambda: cmd_viz(state, args),
        "gif":             lambda: cmd_gif_generico(state, args),
        "nucleo_init":     lambda: cmd_nucleo_init(state, args),
        "nucleo_linha":    lambda: cmd_nucleo_linha(state, args),
        "nucleo_mapa":     lambda: cmd_nucleo_mapa(state, args),
        "nucleo_solve":    lambda: cmd_nucleo_solve(state, args),
        "ec_associacao":   lambda: cmd_ec_associacao(state, args),
    }
    if cmd not in cmds:
        raise ValueError(f"comando desconhecido: {cmd!r}. Digite help.")
    return cmds[cmd](), False


def cmd_gcond(args):
    kw=parse_kwargs(args); G=conduction_G(fnum(kw["k"]),fnum(kw["a"]),fnum(kw["l"]))
    return f"  G_cond = {G:.6g} W/K"

def cmd_gconv(args):
    kw=parse_kwargs(args); G=convection_G(fnum(kw["h"]),fnum(kw["a"]))
    return f"  G_conv = {G:.6g} W/K"


def rodar_demo():
    state = EstadoREPL()
    demo = [
        "# Exemplo chip - item (a)",
        "g_cond k=380 A=7.854e-7 L=0.02",
        "node Tc diffusion Q=4 T_init=30",
        "node Tp boundary T_fixed=44",
        "link Tc Tp cond G=0.17908",
        "solve",
        "# item (b): conveccao",
        "g_conv h=30 A=9e-4",
        "node Tar boundary T_fixed=20",
        "link Tc Tar conv G=0.027",
        "solve",
    ]
    print("="*72, "\nExemplo chip - demo\n" + "="*72)
    for l in demo:
        if l.strip(): print(f"\n>>> {l}")
        try:
            r, _ = executar(state, l)
            if r: print(r)
        except Exception as e:
            print(f"  ERRO: {e}")


def rodar_arquivo(caminho):
    state = EstadoREPL()
    print("="*72 + f"\nExecutando: {caminho}\n" + "="*72)
    with open(caminho) as f:
        linhas = f.readlines()
    for linha in linhas:
        raw = linha.rstrip()
        if not raw or raw.lstrip().startswith("#"):
            if raw: print(f"  {raw}")
            continue
        print(f"\n>>> {raw}")
        try:
            r, stop = executar(state, raw)
            if r: print(r)
            if stop: break
        except Exception as e:
            print(f"  ERRO: {e}")


def rodar_repl():
    state = EstadoREPL()
    print("="*72 + "\nMetodo Nodal - interativo\nDigite help ou quit\n" + "="*72)
    while True:
        try:
            l = input(">>> ")
        except (EOFError, KeyboardInterrupt):
            print("\nSaindo."); break
        try:
            r, stop = executar(state, l)
            if r: print(r)
            if stop: break
        except Exception as e:
            print(f"  ERRO: {e}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg == "--demo":
            rodar_demo()
        else:
            rodar_arquivo(arg)
    else:
        rodar_repl()
