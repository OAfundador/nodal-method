"""
exemplo_placa.py  —  PROBLEMA DIFERENTE: placa retangular 2D
============================================================

ESTE ARQUIVO NÃO É O CHIP. É um problema independente:
uma placa retangular de alumínio com fonte de calor central.

Demonstra o pipeline Geometry2D completo com um caso mais rico que o chip:
múltiplas regiões, fonte volumétrica em sub-região, BC de Dirichlet (bottom)
e BC de convecção (top). Bom ponto de partida para adaptar a novos problemas.

PROBLEMA: placa retangular de alumínio (100 × 60 mm)
    - fonte uniforme q''' = 5×10⁵ W/m³ em região central (40 × 30 mm)
    - bottom: T = 30 °C (Dirichlet)
    - top: convecção h = 50 W/(m²·K), T_ar = 20 °C
    - left e right: adiabáticas

Pipeline completo:
    Geometry2D  ->  regiões com material  ->  fronteiras (BC)
                ->  fonte volumétrica
                ->  discretização (manual ou automática)
                ->  build_network_from_geometry  ->  NodalNetwork
                ->  solve_steady_state
"""

from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from materiais import Material, MaterialPhase
from geometria import Geometry2D
from nos import build_network_from_geometry
from solver import solve_steady_state


# =============================================================================
# Dados do problema
# =============================================================================

W_PLACA = 0.100        # m
H_PLACA = 0.060        # m
ESPESSURA = 0.010      # m (espessura na direção z, fictícia para 2D)

K_ALUMINIO = 200.0     # W/(m·K)

# Região da fonte (centrada)
X0_FONTE = (W_PLACA - 0.040) / 2     # 30 mm
X1_FONTE = X0_FONTE + 0.040           # 70 mm
Y0_FONTE = (H_PLACA - 0.030) / 2     # 15 mm
Y1_FONTE = Y0_FONTE + 0.030           # 45 mm
Q_VOLUME = 5e5         # W/m³

T_BASE = 30.0          # °C, fronteira inferior
T_AR = 20.0            # °C
H_CONVECCAO = 50.0     # W/(m²·K)


# =============================================================================
# Construção da geometria
# =============================================================================

def construir_geometria(nx: int = 20, ny: int = 12) -> Geometry2D:
    """
    Monta uma Geometry2D com:
        - 1 domínio retangular (alumínio)
        - 1 região "fonte" (mesmo material, mas com q''' adicional)
        - 4 fronteiras
        - malha nx × ny
    """
    aluminio = Material(
        name="Alumínio puro",
        phase=MaterialPhase.SOLID,
        k=K_ALUMINIO,
    )

    geom = Geometry2D(width=W_PLACA, height=H_PLACA)

    # Definição de áreas (passo "geometria"):
    # uma região grande (a placa) + uma região menor sobreposta logicamente
    # como "etiqueta da fonte". A região da fonte é definida como uma das
    # áreas; o resto fica coberto pela região "placa". Para evitar
    # sobreposição, descrevemos a placa em 4 retângulos ao redor da fonte.

    # placa - faixa esquerda
    geom.material(
        "placa_esquerda", aluminio,
        x0=0.0, x1=X0_FONTE, y0=0.0, y1=H_PLACA,
    )
    # placa - faixa direita
    geom.material(
        "placa_direita", aluminio,
        x0=X1_FONTE, x1=W_PLACA, y0=0.0, y1=H_PLACA,
    )
    # placa - faixa inferior central
    geom.material(
        "placa_central_inf", aluminio,
        x0=X0_FONTE, x1=X1_FONTE, y0=0.0, y1=Y0_FONTE,
    )
    # placa - faixa superior central
    geom.material(
        "placa_central_sup", aluminio,
        x0=X0_FONTE, x1=X1_FONTE, y0=Y1_FONTE, y1=H_PLACA,
    )
    # região da fonte (mesmo alumínio, mas marcada para receber q''')
    geom.material(
        "fonte", aluminio,
        x0=X0_FONTE, x1=X1_FONTE, y0=Y0_FONTE, y1=Y1_FONTE,
    )

    # Fonte volumétrica
    geom.source(
        region="fonte", kind="volumetric",
        value=Q_VOLUME, name="aquecedor_central",
    )

    # Condições de contorno (passo "fronteiras"):
    geom.bc("bottom", "temperature", T=T_BASE)
    geom.bc("top", "convection", h=H_CONVECCAO, T_inf=T_AR)
    geom.bc("left", "adiabatic")
    geom.bc("right", "adiabatic")

    # Discretização (passo "malha"): tenta nx × ny manual.
    # Caso a malha não esteja alinhada com as bordas das regiões,
    # fix='auto' pede ao framework para sugerir a malha mínima alinhada.
    geom.mesh(nx=nx, ny=ny, fix="auto")
    return geom


# =============================================================================
# Solução
# =============================================================================

def resolver(geom: Geometry2D) -> dict:
    # Conversão automática Geometry2D -> NodalNetwork
    net = build_network_from_geometry(
        geom,
        thickness_z=ESPESSURA,
        default_temperature=T_BASE,
    )

    result = solve_steady_state(
        net,
        tol=1e-9, max_iter=500,
        prefer_scipy=True, update_network=True,
    )

    # Estatísticas básicas
    Ts = [n.temperature for n in net.nodes.values() if n.kind.value != "boundary"]
    Tmin, Tmax = min(Ts), max(Ts)
    Tmedio = sum(Ts) / len(Ts)
    return {
        "n_nos": len(net.nodes),
        "n_unknowns": len(net.unknown_node_ids()),
        "n_links": len(net.links),
        "Tmin": Tmin, "Tmax": Tmax, "Tmedio": Tmedio,
        "convergiu": result.success,
        "residuo": result.residual_norm,
        "iter": result.iterations,
        "net": net,
    }


# =============================================================================
# Saída
# =============================================================================

def main() -> None:
    print("=" * 72)
    print("Exemplo 2D: placa de alumínio com fonte central")
    print("Geometry2D + build_network_from_geometry + solve_steady_state")
    print("=" * 72)

    print("\n--- PASSO 1+2: GEOMETRIA E DISCRETIZAÇÃO ---")
    geom = construir_geometria(nx=20, ny=12)
    print(geom.summary())

    print("\n--- PASSO 3-5: REDE NODAL E RESOLUÇÃO ---")
    out = resolver(geom)
    print(f"  Rede: {out['n_nos']} nós ({out['n_unknowns']} incógnitas), "
          f"{out['n_links']} ligações")
    print(f"  Convergiu: {out['convergiu']}, |R| = {out['residuo']:.2e}, "
          f"iter = {out['iter']}")
    print()
    print(f"  T_min   = {out['Tmin']:7.3f} °C")
    print(f"  T_max   = {out['Tmax']:7.3f} °C")
    print(f"  T_média = {out['Tmedio']:7.3f} °C")

    print("\nPara visualizar a geometria com matplotlib (precisa de matplotlib):")
    print("  >>> from exemplos.exemplo_placa import construir_geometria")
    print("  >>> geom = construir_geometria()")
    print("  >>> geom.show(grid=True, labels=True); import matplotlib.pyplot as plt; plt.show()")


if __name__ == "__main__":
    main()
