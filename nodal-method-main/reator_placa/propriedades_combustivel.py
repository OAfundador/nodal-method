"""
propriedades_combustivel.py

Material genérico do tipo "combustível em dispersão metálica" para uso
didático nos exemplos deste módulo.

Os valores aqui são fictícios e servem apenas para exercitar a estrutura
do método nodal (k(T) tabelado, condução acoplada à temperatura média via
conductance_func). Para um cálculo real, substitua a tabela pelos dados
do combustível específico do seu projeto.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from materiais import Material, MaterialPhase


# Tabela genérica: k decresce suavemente com T (típico de dispersão metálica)
T_TABLE_CELSIUS = (20.0, 50.0, 100.0, 150.0, 200.0, 300.0, 400.0)
K_TABLE_W_MK = (90.0, 88.0, 85.0, 82.0, 80.0, 76.0, 72.0)


def criar_combustivel_didatico() -> Material:
    """Material com k(T) tabelado, interpolado linearmente."""
    return Material(
        name="Combustível didático",
        phase=MaterialPhase.SOLID,
        k_table=(T_TABLE_CELSIUS, K_TABLE_W_MK),
    )
