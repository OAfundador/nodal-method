"""
geometria_reator.py

Representação geométrica de um elemento combustível tipo placa para uso
didático. As dimensões aqui são FICTÍCIAS (valores redondos escolhidos
por simplicidade); a estrutura segue um EC tipo MTR genérico.

Seção transversal modelada (simetria "1 canal + 2 meias placas"):

    y_topo  ╔═══════════════════╗  ← simetria (centro do combustível superior)
            ║   1/2 combustível ║
            ╠═══════════════════╣
            ║   revestimento    ║
            ╠═══════════════════╣
            ║   canal (água)    ║
            ╠═══════════════════╣
            ║   revestimento    ║
            ╠═══════════════════╣
            ║   1/2 combustível ║
    y=0     ╚═══════════════════╝  ← simetria (centro do combustível inferior)

A largura do cerne é Lz; a profundidade do canal é Lcanal. O eixo X
representa a posição AXIAL ao longo do comprimento ativo Lx.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from materiais import Material, MaterialPhase, criar_aluminio
from geometria import Geometry2D
from reator_placa.propriedades_agua import criar_agua_didatica
from reator_placa.propriedades_combustivel import criar_combustivel_didatico


# Dimensões DIDÁTICAS (mockadas — ajuste ao seu problema real)
COMPRIMENTO_ATIVO_M = 0.500              # 500 mm de comprimento ativo
LARGURA_CERNE_M = 0.060                  # 60 mm de largura do cerne
ESPESSURA_COMBUSTIVEL_M = 0.001          # 1.0 mm
ESPESSURA_REVESTIMENTO_M = 0.0005        # 0.5 mm
ESPESSURA_CANAL_M = 0.003                # 3.0 mm
LARGURA_CANAL_M = 0.065                  # 65 mm


@dataclass(frozen=True)
class GeometriaReator:
    """Snapshot das dimensões e materiais da seção."""

    Lx: float = COMPRIMENTO_ATIVO_M
    Ly: float = LARGURA_CERNE_M
    df: float = ESPESSURA_COMBUSTIVEL_M
    dcl: float = ESPESSURA_REVESTIMENTO_M
    dch: float = ESPESSURA_CANAL_M
    Lcanal: float = LARGURA_CANAL_M

    @property
    def altura_secao(self) -> float:
        return self.df + 2.0 * self.dcl + self.dch

    @property
    def area_flow(self) -> float:
        return self.dch * self.Lcanal

    @property
    def perimetro_molhado(self) -> float:
        return 2.0 * (self.dch + self.Lcanal)

    @property
    def Dh(self) -> float:
        return 4.0 * self.area_flow / self.perimetro_molhado


def construir_geometry2d() -> tuple[Geometry2D, dict[str, Material]]:
    """
    Monta a Geometry2D documentando a seção transversal modelada.

    Útil para visualização e checagem; a rede nodal em si é construída
    em modelo_nodal_reator.py com a topologia clássica EC placa.
    """
    g = GeometriaReator()
    combustivel = criar_combustivel_didatico()
    revestimento = criar_aluminio()
    agua = criar_agua_didatica()

    geom = Geometry2D(width=g.Lx, height=g.altura_secao)

    y0 = 0.0
    y1 = y0 + g.df / 2.0
    y2 = y1 + g.dcl
    y3 = y2 + g.dch
    y4 = y3 + g.dcl
    y5 = y4 + g.df / 2.0

    geom.material("meia_placa_inf", combustivel, 0.0, g.Lx, y0, y1)
    geom.material("revestimento_inf", revestimento, 0.0, g.Lx, y1, y2)
    geom.material("canal", agua, 0.0, g.Lx, y2, y3,
                  thermal_mode="convection")
    geom.material("revestimento_sup", revestimento, 0.0, g.Lx, y3, y4)
    geom.material("meia_placa_sup", combustivel, 0.0, g.Lx, y4, y5)

    geom.bc("bottom", "symmetry")
    geom.bc("top", "symmetry")
    geom.bc("left", "fluid_inlet", T=30.0)
    geom.bc("right", "outlet")

    return geom, {"combustivel": combustivel,
                  "revestimento": revestimento,
                  "agua": agua}


def desenhar_secao_transversal(outpath: Path) -> None:
    """Diagrama em escala da seção transversal (espessura × largura do canal)."""
    g = GeometriaReator()
    df_mm = g.df * 1000.0
    dcl_mm = g.dcl * 1000.0
    dch_mm = g.dch * 1000.0
    H_mm = g.Lcanal * 1000.0

    parts = [
        ("1/2 combustível", df_mm / 2.0, "#e6b35a"),
        ("revestimento", dcl_mm, "#b5b5b5"),
        ("canal de água", dch_mm, "#9bd3f5"),
        ("revestimento", dcl_mm, "#b5b5b5"),
        ("1/2 combustível", df_mm / 2.0, "#e6b35a"),
    ]

    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    x = 0.0
    for label, w_mm, color in parts:
        ax.add_patch(Rectangle((x, 0), w_mm, H_mm,
                               facecolor=color, edgecolor="black", alpha=0.85))
        ax.text(x + w_mm / 2.0, H_mm * 0.5, label,
                ha="center", va="center", fontsize=9,
                rotation=90 if w_mm < 0.7 else 0)
        x += w_mm

    ax.set_xlim(0, x)
    ax.set_ylim(-H_mm * 0.05, H_mm * 1.05)
    ax.set_aspect("equal")
    ax.set_xlabel("espessura [mm]")
    ax.set_ylabel("largura do canal [mm]")
    ax.set_title("Seção transversal modelada (1 canal + 2 meias placas)")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)
