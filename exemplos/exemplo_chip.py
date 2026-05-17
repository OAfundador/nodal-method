"""
exemplo_chip.py

Exemplo 1 do Capítulo V do livro
    "Transferência de Calor Computacional - Método Nodal"
    J. L. Ferraz Bastos e D. A. de Andrade, IPEN/CNEN-SP

CÁLCULO DA TEMPERATURA DE UM COMPONENTE ELETRÔNICO

Geometria
    componente: 30 mm × 30 mm soldado a uma placa por 12 conexões de cobre
    conexões  : L = 20 mm, d = 1 mm, k = 380 W/(m·K)

Hipóteses
    Tp (temperatura da placa) = 44 °C   (condição de contorno)
    Q (potência dissipada)    = 4 W
    item (a) sem convecção
    item (b) com convecção para o ar:  Tar = 20 °C, h = 30 W/(m²·K),
             área convectiva = 30 mm × 30 mm

Este script demonstra os 5 passos do método nodal usando a infraestrutura
do framework:

    Passo 1 - Discretização espacial      -> 1 nó (Tc) por modelo
    Passo 2 - Identificação dos tipos      -> Tc DIFFUSION, Tp/Tar BOUNDARY
    Passo 3 - Cálculo das condutâncias     -> conduction_G, convection_G
    Passo 4 - Equações de balanço          -> NodalNetwork.residual_steady
    Passo 5 - Resolução do sistema R(T)=0  -> solve_steady_state
"""

from __future__ import annotations

import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from condutancias import conduction_G, convection_G
from nos import NodalNetwork, NodeKind, TransferKind
from solver import solve_steady_state


# =============================================================================
# Dados do livro (Exemplo 1, Capítulo V)
# =============================================================================

N_CONEXOES = 12
L_CONEXAO = 20e-3        # m
D_CONEXAO = 1e-3         # m
K_COBRE = 380.0          # W/(m·K)

TP_PLACA = 44.0          # °C, condição de contorno
Q_DISSIPADA = 4.0        # W

# Item (b): convecção para o ar
TAR_AMBIENTE = 20.0      # °C
H_CONVECCAO = 30.0       # W/(m²·K)
LADO_COMPONENTE = 30e-3  # m
AREA_CONVECCAO = LADO_COMPONENTE ** 2


# =============================================================================
# Item (a) — somente condução pelas pernas
# =============================================================================

def resolver_item_a() -> dict:
    """
    Apenas 2 nós:
        Tc (DIFFUSION) com fonte Q
        Tp (BOUNDARY)  com Tp fixa em 44 °C

    Ligados por uma única condutância de condução, equivalente às 12 pernas
    de cobre em paralelo:

        G_cond = N · k · A / L,    A = π d² / 4
    """
    # Passo 3: condutância equivalente das 12 pernas
    area_perna = math.pi * D_CONEXAO ** 2 / 4.0
    G_cond_uma_perna = conduction_G(K_COBRE, area_perna, L_CONEXAO)
    G_cond_total = N_CONEXOES * G_cond_uma_perna

    # Passos 1 e 2: rede com 2 nós
    net = NodalNetwork()
    Tc = net.add_node(
        name="Tc", kind=NodeKind.DIFFUSION,
        x=0.0, y=0.0, volume=1.0, source=Q_DISSIPADA,
        temperature=30.0,
    )
    Tp = net.add_node(
        name="Tp", kind=NodeKind.BOUNDARY,
        x=0.0, y=-L_CONEXAO,
        fixed_temperature=TP_PLACA, temperature=TP_PLACA,
    )

    # Passo 4: ligação de condução (G é constante neste caso)
    net.add_link(
        Tc, Tp, TransferKind.CONDUCTION,
        conductance=G_cond_total, name="conducao_pernas",
    )

    # Passo 5: resolução do sistema
    result = solve_steady_state(net, tol=1e-10, max_iter=50)

    Tc_val = net.nodes[Tc].temperature
    Tc_analitico = TP_PLACA + Q_DISSIPADA / G_cond_total

    return {
        "G_cond_uma_perna": G_cond_uma_perna,
        "G_cond_total": G_cond_total,
        "Tc": Tc_val,
        "Tc_analitico": Tc_analitico,
        "residuo": result.residual_norm,
        "convergiu": result.success,
    }


# =============================================================================
# Item (b) — condução pelas pernas + convecção para o ar
# =============================================================================

def resolver_item_b() -> dict:
    """
    Mesma rede do item (a), agora com uma ligação adicional de CONVECTION
    entre Tc e Tar (BOUNDARY).
    """
    # Passo 3: mesma condutância de condução das pernas
    area_perna = math.pi * D_CONEXAO ** 2 / 4.0
    G_cond_total = N_CONEXOES * conduction_G(K_COBRE, area_perna, L_CONEXAO)
    G_conv = convection_G(H_CONVECCAO, AREA_CONVECCAO)

    # Passos 1 e 2: rede com 3 nós
    net = NodalNetwork()
    Tc = net.add_node(
        name="Tc", kind=NodeKind.DIFFUSION,
        x=0.0, y=0.0, volume=1.0, source=Q_DISSIPADA,
        temperature=30.0,
    )
    Tp = net.add_node(
        name="Tp", kind=NodeKind.BOUNDARY,
        x=0.0, y=-L_CONEXAO,
        fixed_temperature=TP_PLACA, temperature=TP_PLACA,
    )
    Tar = net.add_node(
        name="Tar", kind=NodeKind.BOUNDARY,
        x=LADO_COMPONENTE, y=0.0,
        fixed_temperature=TAR_AMBIENTE, temperature=TAR_AMBIENTE,
    )

    # Passo 4: duas ligações
    net.add_link(
        Tc, Tp, TransferKind.CONDUCTION,
        conductance=G_cond_total, name="conducao_pernas",
    )
    net.add_link(
        Tc, Tar, TransferKind.CONVECTION,
        conductance=G_conv, name="conveccao_ar",
    )

    # Passo 5: resolução
    result = solve_steady_state(net, tol=1e-10, max_iter=50)

    Tc_val = net.nodes[Tc].temperature
    Tc_analitico = (
        G_cond_total * TP_PLACA + G_conv * TAR_AMBIENTE + Q_DISSIPADA
    ) / (G_cond_total + G_conv)

    return {
        "G_cond_total": G_cond_total,
        "G_conv": G_conv,
        "Tc": Tc_val,
        "Tc_analitico": Tc_analitico,
        "residuo": result.residual_norm,
        "convergiu": result.success,
    }


# =============================================================================
# Saída
# =============================================================================

def main() -> None:
    print("=" * 72)
    print("Exemplo 1 do livro — componente eletrônico")
    print("Transferência de Calor Computacional, Método Nodal")
    print("J. L. Ferraz Bastos e D. A. de Andrade, IPEN/CNEN-SP")
    print("=" * 72)

    print("\nDados:")
    print(f"  componente {LADO_COMPONENTE*1e3:.0f} × {LADO_COMPONENTE*1e3:.0f} mm,"
          f" {N_CONEXOES} pernas de cobre")
    print(f"  L = {L_CONEXAO*1e3:.0f} mm, d = {D_CONEXAO*1e3:.1f} mm,"
          f" k = {K_COBRE:.0f} W/(m·K)")
    print(f"  Tp = {TP_PLACA:.1f} °C, Q = {Q_DISSIPADA:.1f} W")

    print("\nItem (a) — somente condução pelas pernas")
    a = resolver_item_a()
    print(f"  G_cond (12 pernas)   = {a['G_cond_total']:.4f} W/K")
    print(f"  Tc (NodalNetwork)    = {a['Tc']:.4f} °C")
    print(f"  Tc (analítico)       = {a['Tc_analitico']:.4f} °C")
    print(f"  Tc (livro)           ≈ 66.3 °C   (Tp + Q/G = 44 + 4/0.1791)")
    print(f"  Convergiu: {a['convergiu']}, |R| = {a['residuo']:.2e}")

    print(f"\nItem (b) — condução + convecção (h = {H_CONVECCAO} W/(m²·K),"
          f" Tar = {TAR_AMBIENTE} °C)")
    b = resolver_item_b()
    print(f"  G_cond               = {b['G_cond_total']:.4f} W/K")
    print(f"  G_conv               = {b['G_conv']:.4f} W/K")
    print(f"  Tc (NodalNetwork)    = {b['Tc']:.4f} °C")
    print(f"  Tc (analítico)       = {b['Tc_analitico']:.4f} °C")
    print(f"  Convergiu: {b['convergiu']}, |R| = {b['residuo']:.2e}")


if __name__ == "__main__":
    main()
