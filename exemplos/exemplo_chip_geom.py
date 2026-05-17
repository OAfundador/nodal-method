"""
exemplo_chip_geom.py  —  ABORDAGEM GEOMETRY2D (malha automática)
================================================================

Resolve o mesmo problema do chip eletrônico (Q = 4 W, pernas de cobre),
mas pelo caminho automático: Geometry2D → regiões → malha → build_network_from_geometry.

Compare com exemplo_chip.py, que monta a rede manualmente (1 nó por item).
Aqui a malha cria dezenas de nós e a condutividade das pernas é agregada
em k_eff para compatibilidade com a representação 2D.

Útil como referência para geometrias mais complexas, onde o caminho
manual seria impraticável.

Reproduz o Exemplo 1 do Capítulo V do livro (componente eletrônico) usando
o caminho Geometry2D + build_network_from_geometry + solve_steady_state,
em vez do caminho lumped de `exemplo_chip.py`.

ESTRATÉGIA: CONDUTIVIDADE EFETIVA NA CAMADA DAS PERNAS
------------------------------------------------------

O problema 3D real tem 12 pernas cilíndricas (A_cil = π·d²/4) discretas
imersas em ar, com um componente cúbico em cima.

Em uma representação 2D, cada "perna" vira uma barra retangular cuja área
transversal é a largura × espessura_z, em geral muito maior que a área do
cilindro original. Para que a CONDUTÂNCIA térmica do modelo 2D bata com o
problema 3D, agregamos as 12 pernas + o ar circundante em uma ÚNICA
camada com condutividade efetiva k_eff calibrada para satisfazer

    k_eff · A_camada_2D / L_camada  =  12 · k_cobre · A_cilindro / L_perna

Isso é uma técnica padrão para tratar materiais compostos em 2D.
Cálculo: k_eff = G_pernas_real · L / (W · t_z) = 0.179·0.020 / (0.030·0.030)
                ≈ 3.98 W/(m·K)

GEOMETRIA MODELADA
------------------
                                  topo (item a: adiabatic;
                                        item b: convection h, T_inf)
              ___________________
       y=30  |                   |   COMPONENTE (silício, k=150)
             |        Q          |   q''' = Q_total / V_componente
       y=20  |___________________|
             |                   |
             |   camada pernas   |   k_eff = 3.98 W/(m·K)
             |   (12 pernas+ar)  |
       y=0   |___________________|
             x=0               x=30  bottom (placa, T = 44 °C)

  left, right: adiabatic
  thickness_z = 30 mm (profundidade do componente em z)

O componente em silício (k=150) tem temperatura quase uniforme — o livro
o trata como UM nó isotérmico (Tc). Aqui ele aparece como muitos nós da
discretização 2D e a "Tc" é a temperatura média do componente.

RESULTADO ESPERADO
------------------
  Item (a) sem convecção  : Tc ≈ 66.34 °C  (livro: 66.3 °C)
  Item (b) com convecção  : Tc ≈ 60.27 °C
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from materiais import Material, MaterialPhase
from geometria import Geometry2D
from nos import build_network_from_geometry
from solver import solve_steady_state


# =============================================================================
# Dados
# =============================================================================

LADO_COMPONENTE = 30e-3
ALTURA_COMPONENTE = 10e-3
ALTURA_PERNAS = 20e-3
ESPESSURA_Z = 30e-3

N_PERNAS = 12
D_PERNA = 1.0e-3
K_COBRE = 380.0

A_CILINDRO = math.pi * D_PERNA ** 2 / 4.0
G_PERNAS_3D = N_PERNAS * K_COBRE * A_CILINDRO / ALTURA_PERNAS
A_CAMADA_2D = LADO_COMPONENTE * ESPESSURA_Z
K_EFF_PERNAS = G_PERNAS_3D * ALTURA_PERNAS / A_CAMADA_2D

K_SILICIO = 150.0
Q_TOTAL = 4.0
V_COMPONENTE = LADO_COMPONENTE * ALTURA_COMPONENTE * ESPESSURA_Z
Q_VOLUMETRICO = Q_TOTAL / V_COMPONENTE

T_PLACA = 44.0
T_AR_AMB = 20.0
H_CONVECCAO = 30.0


def construir_geometria(item: str = "a", nx: int = 12, ny: int = 12) -> Geometry2D:
    silicio = Material(name="Silício", phase=MaterialPhase.SOLID, k=K_SILICIO)
    pernas_eff = Material(
        name=f"Pernas+ar (k_eff={K_EFF_PERNAS:.2f})",
        phase=MaterialPhase.SOLID, k=K_EFF_PERNAS,
    )

    geom = Geometry2D(width=LADO_COMPONENTE, height=ALTURA_PERNAS + ALTURA_COMPONENTE)

    geom.material(
        "camada_pernas", pernas_eff,
        x0=0.0, x1=LADO_COMPONENTE, y0=0.0, y1=ALTURA_PERNAS,
    )
    geom.material(
        "componente", silicio,
        x0=0.0, x1=LADO_COMPONENTE,
        y0=ALTURA_PERNAS, y1=ALTURA_PERNAS + ALTURA_COMPONENTE,
    )

    geom.source(region="componente", kind="volumetric", value=Q_VOLUMETRICO)

    geom.bc("bottom", "temperature", T=T_PLACA)
    geom.bc("left", "adiabatic")
    geom.bc("right", "adiabatic")
    if item == "a":
        geom.bc("top", "adiabatic")
    elif item == "b":
        geom.bc("top", "convection", h=H_CONVECCAO, T_inf=T_AR_AMB)
    else:
        raise ValueError("item deve ser 'a' ou 'b'.")

    geom.mesh(nx=nx, ny=ny, fix="auto")
    return geom


def resolver(geom: Geometry2D) -> dict:
    net = build_network_from_geometry(
        geom, thickness_z=ESPESSURA_Z, default_temperature=T_PLACA,
    )
    result = solve_steady_state(
        net, tol=1e-9, max_iter=500,
        prefer_scipy=True, update_network=True,
    )

    Ts_componente = [
        n.temperature for n in net.nodes.values()
        if n.kind.value != "boundary"
        and getattr(n.region, "name", "") == "componente"
    ]
    Ts_pernas = [
        n.temperature for n in net.nodes.values()
        if n.kind.value != "boundary"
        and getattr(n.region, "name", "") == "camada_pernas"
    ]

    return {
        "n_nos": len(net.nodes),
        "n_unknowns": len(net.unknown_node_ids()),
        "n_links": len(net.links),
        "Tc_max": max(Ts_componente),
        "Tc_min": min(Ts_componente),
        "Tc_med": sum(Ts_componente) / len(Ts_componente),
        "Tp_med_camada": sum(Ts_pernas) / len(Ts_pernas),
        "convergiu": result.success,
        "residuo": result.residual_norm,
        "iter": result.iterations,
    }


def main() -> None:
    print("=" * 72)
    print("Exemplo do chip via Geometry2D + condutividade efetiva")
    print("Reproduz o Exemplo 1 do Cap. V do livro (Bastos & Andrade)")
    print("=" * 72)

    print("\nDados / calibração:")
    print(f"  Componente  : {LADO_COMPONENTE*1e3:.0f} × {LADO_COMPONENTE*1e3:.0f}"
          f" × {ALTURA_COMPONENTE*1e3:.0f} mm  (k_silício = {K_SILICIO} W/m·K)")
    print(f"  Camada pernas: 30 × 20 mm × {ESPESSURA_Z*1e3:.0f} mm,"
          f" k_eff = {K_EFF_PERNAS:.3f} W/m·K")
    print(f"  G(12 pernas) = {G_PERNAS_3D:.5f} W/K  (= G_2D por construção)")
    print(f"  Q_total = {Q_TOTAL} W,"
          f" q''' = {Q_VOLUMETRICO:.3e} W/m³,"
          f" Tp = {T_PLACA} °C")

    for item in ["a", "b"]:
        descricao = "(a) sem convecção (top adiabatic)" if item == "a" \
                    else "(b) com convecção no topo (h=30, Tar=20)"
        print(f"\n--- Item {descricao} ---")
        geom = construir_geometria(item=item)
        out = resolver(geom)
        print(f"  Geometria: {geom.nx}×{geom.ny} cells,"
              f" dx={geom.width/geom.nx*1e3:.3f} mm,"
              f" dy={geom.height/geom.ny*1e3:.3f} mm")
        print(f"  Rede     : {out['n_nos']} nós ({out['n_unknowns']} incógnitas),"
              f" {out['n_links']} ligações")
        print(f"  Solver   : convergiu={out['convergiu']},"
              f" |R|={out['residuo']:.2e}, iter={out['iter']}")
        print(f"  T_componente  : min={out['Tc_min']:6.3f},"
              f" max={out['Tc_max']:6.3f}, média={out['Tc_med']:6.3f} °C")
        print(f"  T_camada_pernas (média) : {out['Tp_med_camada']:.3f} °C")

    print("\nReferência (livro, lumped):")
    print("  item (a) - Listagem 5.1.1 + eq. TC = TP + Q/G : 66.3 °C")
    print("  item (b) - Listagem 5.1.3 com h=30 W/m²K      : 60.266 °C")
    print()
    print("Reprodução:")
    print("  exemplo_chip.py (lumped, 2-3 nós):  item (a) 66.34 °C / item (b) 60.27 °C")
    print("  este script (Geometry2D + k_eff):  item (a) 66.44 °C / item (b) 60.32 °C")
    print("  (a pequena diferença do Geometry2D mostra a queda de T dentro")
    print("   do componente, que o modelo lumped do livro despreza)")


if __name__ == "__main__":
    main()
