"""
exemplo_reator_placa.py

Demonstração: aplica o método nodal a um elemento combustível (EC) tipo
placa, modelando 1 canal de refrigeração + 2 meias placas com 9 nós por
camada axial.

DADOS USADOS - TODOS DIDÁTICOS / FICTÍCIOS
-------------------------------------------
Dimensões redondas, fator de pico unitário, vazão razoável para a
geometria. NÃO REPRODUZ nenhum reator real. Use apenas como exemplo
de aplicação do framework.

Casos rodados:
    1. q'' cossenoidal (perfil típico de reator com refletor)
    2. q'' constante (para comparação)

Saídas em saidas_reator_placa/:
    geometria_secao_transversal.png
    rede_nodal_resumo.txt
    resumo_<caso>.txt
    resultado_<caso>.csv
    temperaturas_<caso>.png
    fluxo_<caso>.png
    comparacao_cos_vs_constante.png
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from solver import solve_steady_state

from reator_placa.geometria_reator import (
    GeometriaReator, construir_geometry2d, desenhar_secao_transversal,
)
from reator_placa.modelo_nodal_reator import (
    ConfigCaso, MODO_COS, MODO_CONST,
    construir_rede, extrair_resultados, gerar_chute_inicial,
)


# Potência didática da placa modelada
P_PLACA_W = 5000.0          # W (mockado — escolhido para gerar T_max razoável)
VAZAO_CANAL_M3_S = 3.0e-4   # m³/s ≈ 1.08 m³/h por canal


def resolver_caso(cfg: ConfigCaso) -> dict:
    g = GeometriaReator()
    rede = construir_rede(cfg, g)
    z0 = gerar_chute_inicial(rede)
    result = solve_steady_state(
        rede.net, z0=z0, tol=cfg.tol, max_iter=cfg.max_iter,
        update_network=True, prefer_scipy=True,
    )
    out = extrair_resultados(rede)
    out["solver_success"] = bool(result.success)
    out["solver_residual_norm"] = float(result.residual_norm)
    out["solver_iterations"] = result.iterations
    out["n_nos"] = len(rede.net.nodes)
    out["n_links"] = len(rede.net.links)
    return out


def texto_resumo(res: dict, caso: str) -> str:
    cfg = res["cfg"]
    g = res["geom"]
    margem = res["Tsat"] - res["Tch"]
    return f"""EC tipo placa — Método Nodal (DIDÁTICO, dados mockados)

Caso = {caso}
Modo de fluxo = {cfg.modo_fluxo}
N axial = {cfg.n_axial}
Vazão por canal = {cfg.vazao_canal_m3_s:.4e} m³/s
m_dot por canal = {res['m_dot']:.4e} kg/s

Geometria (dados didáticos):
  comprimento ativo Lx = {g.Lx*1e3:.1f} mm
  largura cerne Ly     = {g.Ly*1e3:.1f} mm
  espessura combust.   = {g.df*1e3:.2f} mm
  espessura clad       = {g.dcl*1e3:.2f} mm
  canal interno        = {g.dch*1e3:.2f} × {g.Lcanal*1e3:.1f} mm
  Dh                   = {g.Dh*1e3:.3f} mm

Rede nodal:
  n nós   = {res['n_nos']}
  n links = {res['n_links']}

Solver:
  success        = {res['solver_success']}
  residual_norm  = {res['solver_residual_norm']:.3e}
  iterações      = {res['solver_iterations']}

Potência:
  P_placa configurada = {cfg.P_placa_W:.2f} W
  Potência integrada  = {float(np.sum(res['Q_total'])):.2f} W

Temperaturas:
  T_fluido_entrada    = {cfg.T_in_C:.2f} °C
  T_fluido_saida      = {float(res['Tch'][-1]):.2f} °C
  T_superficie_max    = {float(np.max(res['Ts'])):.2f} °C
  T_revestimento_max  = {float(np.max(res['Tcl'])):.2f} °C
  T_interface_max     = {float(np.max(res['Ti'])):.2f} °C
  T_combustivel_max   = {float(np.max(res['Tf'])):.2f} °C
  margem T_sat-T_ch (min) = {float(np.min(margem)):.2f} °C

Hidráulica:
  Re_min, Re_max = {float(np.min(res['Re'])):.1f}, {float(np.max(res['Re'])):.1f}
  h_min, h_max   = {float(np.min(res['h'])):.1f}, {float(np.max(res['h'])):.1f} W/m²K
"""


def salvar_csv(res: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cfg = res["cfg"]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
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
                float(res["Tsat"][j]),
                float(res["Tsat"][j] - res["Tch"][j]),
            ])


def plot_temperaturas(res: dict, path: Path, titulo: str) -> None:
    z_cm = res["z"] * 100.0
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    ax.plot(z_cm, res["Tch"], label="fluido")
    ax.plot(z_cm, res["Ts"], label="superfície")
    ax.plot(z_cm, res["Tcl"], label="centro do revestimento")
    ax.plot(z_cm, res["Ti"], label="interface combustível/revestimento")
    ax.plot(z_cm, res["Tf"], label="centro do combustível (T máx.)")
    ax.set_xlabel("posição axial desde a entrada [cm]")
    ax.set_ylabel("Temperatura [°C]")
    ax.set_title(titulo)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_fluxo(res: dict, path: Path, titulo: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.7))
    ax.plot(res["x"] * 100.0, res["q_flux_face"] / 1e6)
    ax.set_xlabel("posição axial x (centro = 0) [cm]")
    ax.set_ylabel("q'' por face [MW/m²]")
    ax.set_title(titulo)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    outdir = ROOT / "saidas_reator_placa"
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("EC tipo placa - Método Nodal")
    print("=" * 78)
    print()
    print("AVISO: este é um exemplo DIDÁTICO. Todos os dados numéricos")
    print("(dimensões, potência, vazão, propriedades) são FICTÍCIOS.")
    print("Para detalhes, leia reator_placa/LEIA-ME.md.")
    print("=" * 78)

    print("\nGeometria (Geometry2D)...")
    geom2d, _ = construir_geometry2d()
    desenhar_secao_transversal(outdir / "geometria_secao_transversal.png")
    print(f"  -> {outdir/'geometria_secao_transversal.png'}")

    casos = [
        ("cos", ConfigCaso(n_axial=40, vazao_canal_m3_s=VAZAO_CANAL_M3_S,
                           P_placa_W=P_PLACA_W, modo_fluxo=MODO_COS)),
        ("constante", ConfigCaso(n_axial=40, vazao_canal_m3_s=VAZAO_CANAL_M3_S,
                                 P_placa_W=P_PLACA_W, modo_fluxo=MODO_CONST)),
    ]

    resultados = {}
    for tag, cfg in casos:
        print(f"\nResolvendo caso '{tag}' ...")
        res = resolver_caso(cfg)
        resultados[tag] = res

        salvar_csv(res, outdir / f"resultado_{tag}.csv")
        (outdir / f"resumo_{tag}.txt").write_text(
            texto_resumo(res, tag), encoding="utf-8")
        plot_temperaturas(res, outdir / f"temperaturas_{tag}.png",
                          f"EC placa — q'' {tag}")
        plot_fluxo(res, outdir / f"fluxo_{tag}.png",
                   f"Fluxo axial — q'' {tag}")

        print(f"  Tch_out = {res['Tch'][-1]:.2f} °C, "
              f"Tf_max = {float(np.max(res['Tf'])):.2f} °C, "
              f"|R| = {res['solver_residual_norm']:.1e}")

    # Comparação
    z = resultados["cos"]["z"] * 100.0
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(z, resultados["cos"]["Tf"], label="Tf - cos")
    ax.plot(z, resultados["constante"]["Tf"], label="Tf - constante")
    ax.plot(z, resultados["cos"]["Tch"], "--", label="fluido - cos")
    ax.plot(z, resultados["constante"]["Tch"], "--", label="fluido - constante")
    ax.set_xlabel("posição axial desde a entrada [cm]")
    ax.set_ylabel("Temperatura [°C]")
    ax.set_title("Comparação q'' cos(π/2 x/L) vs q'' constante")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "comparacao_cos_vs_constante.png", dpi=180)
    plt.close(fig)
    print(f"\nComparação salva em {outdir/'comparacao_cos_vs_constante.png'}")

    # Resumo da rede para documentação
    g = GeometriaReator()
    rede_demo = construir_rede(
        ConfigCaso(n_axial=10, P_placa_W=P_PLACA_W,
                   vazao_canal_m3_s=VAZAO_CANAL_M3_S, modo_fluxo=MODO_COS),
        g,
    )
    (outdir / "rede_nodal_resumo.txt").write_text(
        rede_demo.net.summary(max_nodes=30, max_links=40), encoding="utf-8")

    print(f"\nTodas as saídas em: {outdir}")


if __name__ == "__main__":
    main()
