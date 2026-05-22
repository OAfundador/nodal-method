"""
materiais.py

Classe simples para representar materiais e suas propriedades termofísicas
no modelo nodal.

Convenções:
- Temperatura T em °C
- Pressão P em Pa
- k em W/(m K)
- rho em kg/m³
- cp em J/(kg K)
- mu em Pa s

Cada propriedade pode ser definida por:
1. valor constante;
2. polinômio em T;
3. tabela com interpolação linear;
4. função customizada f(T, P).

Atualização principal:
- O material agora possui uma classificação física:
    MaterialPhase.SOLID
    MaterialPhase.FLUID

Essa classificação serve para sugerir o tipo padrão de troca térmica:
- sólido-sólido  -> condução
- sólido-fluido  -> convecção
- fluido-fluido  -> transporte fluido, se houver direção de escoamento

Atenção:
- Fonte de calor NÃO é propriedade do material.
- Fonte de calor deve entrar depois na geometria/região ou diretamente no nó.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Callable, Optional, Sequence
import numpy as np


# ============================================================
# Classificação física do material
# ============================================================

class MaterialPhase(Enum):
    """
    Classificação física básica do material.

    SOLID:
        Material sólido. Normalmente participa de condução.

    FLUID:
        Material fluido. Pode participar de convecção em interface sólido-fluido
        e de transporte entálpico quando existir direção de escoamento.
    """

    SOLID = "solid"
    FLUID = "fluid"


# ============================================================
# Classe Material
# ============================================================

@dataclass
class Material:
    """
    Material com propriedades termofísicas.

    Parâmetros
    ----------
    name:
        Nome do material.

    phase:
        MaterialPhase.SOLID ou MaterialPhase.FLUID.
        Também aceita strings:
            "solid", "solido", "sólido", "s"
            "fluid", "fluido", "f"

    Propriedades
    ------------
    Para cada propriedade k, rho, cp e mu, a prioridade é:

        função customizada > tabela > polinômio > constante

    Exemplo sólido:
        aluminio = Material(
            name="Alumínio",
            phase="solid",
            k=180.0,
        )

    Exemplo fluido:
        agua = Material(
            name="Água",
            phase="fluid",
            k_func=...,
            rho_func=...,
            cp_func=...,
            mu_func=...,
        )
    """

    name: str
    phase: MaterialPhase | str = MaterialPhase.SOLID

    # Valores constantes
    k: Optional[float] = None
    rho: Optional[float] = None
    cp: Optional[float] = None
    mu: Optional[float] = None

    # Coeficientes polinomiais em T:
    # prop(T) = a0 + a1*T + a2*T² + ...
    k_coeffs: Optional[Sequence[float]] = None
    rho_coeffs: Optional[Sequence[float]] = None
    cp_coeffs: Optional[Sequence[float]] = None
    mu_coeffs: Optional[Sequence[float]] = None

    # Tabelas para interpolação linear:
    # prop(T) = interp(T, T_table, prop_table)
    k_table: Optional[tuple[Sequence[float], Sequence[float]]] = None
    rho_table: Optional[tuple[Sequence[float], Sequence[float]]] = None
    cp_table: Optional[tuple[Sequence[float], Sequence[float]]] = None
    mu_table: Optional[tuple[Sequence[float], Sequence[float]]] = None

    # Funções customizadas f(T, P)
    k_func: Optional[Callable[[float, Optional[float]], float]] = None
    rho_func: Optional[Callable[[float, Optional[float]], float]] = None
    cp_func: Optional[Callable[[float, Optional[float]], float]] = None
    mu_func: Optional[Callable[[float, Optional[float]], float]] = None

    def __post_init__(self) -> None:
        self.phase = self._normalize_phase(self.phase)

    @staticmethod
    def _normalize_phase(phase: MaterialPhase | str) -> MaterialPhase:
        if isinstance(phase, MaterialPhase):
            return phase

        text = str(phase).strip().lower()

        solid_aliases = {
            "solid",
            "solido",
            "sólido",
            "s",
            "solid_phase",
        }

        fluid_aliases = {
            "fluid",
            "fluido",
            "f",
            "fluid_phase",
        }

        if text in solid_aliases:
            return MaterialPhase.SOLID

        if text in fluid_aliases:
            return MaterialPhase.FLUID

        raise ValueError(
            f"Fase de material inválida: {phase!r}. "
            "Use 'solid'/'sólido' ou 'fluid'/'fluido'."
        )

    # ========================================================
    # Classificação física
    # ========================================================

    def is_solid(self) -> bool:
        return self.phase == MaterialPhase.SOLID

    def is_fluid(self) -> bool:
        return self.phase == MaterialPhase.FLUID

    def phase_name(self) -> str:
        return self.phase.value

    def default_interaction_with(self, other: "Material") -> str:
        """
        Sugere o tipo padrão de troca térmica entre dois materiais.

        Retornos possíveis:
            "conduction"
            "convection"
            "fluid_transport"

        Importante:
        Esta função só sugere o padrão. Exceções físicas devem ser dadas
        explicitamente pelo usuário na geometria/região/interface.
        """

        if self.is_solid() and other.is_solid():
            return "conduction"

        if self.is_fluid() and other.is_fluid():
            return "fluid_transport"

        return "convection"

    # ========================================================
    # Verificações de propriedade
    # ========================================================

    def has_property_definition(self, prop_name: str) -> bool:
        """
        Verifica se a propriedade foi definida de alguma forma.
        Não calcula o valor.
        """

        prop_name = prop_name.lower().strip()

        if prop_name == "k":
            return any(
                item is not None
                for item in [self.k, self.k_coeffs, self.k_table, self.k_func]
            )

        if prop_name == "rho":
            return any(
                item is not None
                for item in [self.rho, self.rho_coeffs, self.rho_table, self.rho_func]
            )

        if prop_name == "cp":
            return any(
                item is not None
                for item in [self.cp, self.cp_coeffs, self.cp_table, self.cp_func]
            )

        if prop_name == "mu":
            return any(
                item is not None
                for item in [self.mu, self.mu_coeffs, self.mu_table, self.mu_func]
            )

        raise ValueError("Use 'k', 'rho', 'cp' ou 'mu'.")

    def missing_core_properties(self) -> list[str]:
        """
        Lista propriedades essenciais ausentes segundo a fase.

        Para sólido:
            k é essencial.

        Para fluido:
            k, rho, cp e mu são essenciais para cálculo típico de Re, Pr e h.
        """

        if self.is_solid():
            required = ["k"]
        else:
            required = ["k", "rho", "cp", "mu"]

        return [
            prop
            for prop in required
            if not self.has_property_definition(prop)
        ]

    def validate_core_properties(self) -> None:
        """
        Gera erro se faltar propriedade essencial para a fase do material.
        """

        missing = self.missing_core_properties()

        if missing:
            raise ValueError(
                f"Material '{self.name}' está incompleto para phase={self.phase.value}. "
                f"Propriedades ausentes: {', '.join(missing)}."
            )

    # ========================================================
    # Avaliação de propriedades
    # ========================================================

    def prop(self, prop_name: str, T: float, P: Optional[float] = None) -> float:
        """
        Calcula uma propriedade do material.

        Parâmetros
        ----------
        prop_name:
            "k", "rho", "cp" ou "mu".

        T:
            Temperatura em °C.

        P:
            Pressão em Pa. Pode ser None se a propriedade não depender de pressão.
        """

        prop_name = prop_name.lower().strip()

        if prop_name == "k":
            return self._eval_property(
                prop_name="k",
                constant_value=self.k,
                coeffs=self.k_coeffs,
                table=self.k_table,
                func=self.k_func,
                T=T,
                P=P,
            )

        if prop_name == "rho":
            return self._eval_property(
                prop_name="rho",
                constant_value=self.rho,
                coeffs=self.rho_coeffs,
                table=self.rho_table,
                func=self.rho_func,
                T=T,
                P=P,
            )

        if prop_name == "cp":
            return self._eval_property(
                prop_name="cp",
                constant_value=self.cp,
                coeffs=self.cp_coeffs,
                table=self.cp_table,
                func=self.cp_func,
                T=T,
                P=P,
            )

        if prop_name == "mu":
            return self._eval_property(
                prop_name="mu",
                constant_value=self.mu,
                coeffs=self.mu_coeffs,
                table=self.mu_table,
                func=self.mu_func,
                T=T,
                P=P,
            )

        raise ValueError(
            f"Propriedade desconhecida '{prop_name}'. "
            "Use 'k', 'rho', 'cp' ou 'mu'."
        )

    def _eval_property(
        self,
        prop_name: str,
        constant_value: Optional[float],
        coeffs: Optional[Sequence[float]],
        table: Optional[tuple[Sequence[float], Sequence[float]]],
        func: Optional[Callable[[float, Optional[float]], float]],
        T: float,
        P: Optional[float],
    ) -> float:
        """
        Avalia uma propriedade seguindo a prioridade:
            função > tabela > polinômio > constante
        """

        if func is not None:
            value = func(T, P)

        elif table is not None:
            T_table, y_table = table
            value = self._interp(T, T_table, y_table, prop_name)

        elif coeffs is not None:
            value = self._poly(T, coeffs)

        elif constant_value is not None:
            value = constant_value

        else:
            raise ValueError(
                f"Material '{self.name}' não possui definição para "
                f"a propriedade '{prop_name}'."
            )

        value = float(value)

        if not np.isfinite(value):
            raise ValueError(
                f"Material '{self.name}', propriedade '{prop_name}' "
                "retornou valor não finito."
            )

        if value <= 0.0:
            raise ValueError(
                f"Material '{self.name}', propriedade '{prop_name}' "
                f"retornou valor não físico: {value}."
            )

        return value

    @staticmethod
    def _poly(T: float, coeffs: Sequence[float]) -> float:
        """
        Avalia:
            a0 + a1*T + a2*T² + ...
        """
        return sum(a * T**i for i, a in enumerate(coeffs))

    def _interp(
        self,
        T: float,
        T_table: Sequence[float],
        y_table: Sequence[float],
        prop_name: str,
    ) -> float:
        """
        Interpolação linear em temperatura.

        Por segurança, esta função não extrapola.
        Se T estiver fora da tabela, gera erro.
        """

        T_array = np.asarray(T_table, dtype=float)
        y_array = np.asarray(y_table, dtype=float)

        if T_array.ndim != 1 or y_array.ndim != 1:
            raise ValueError(
                f"Material '{self.name}', propriedade '{prop_name}': "
                "tabelas precisam ser vetores 1D."
            )

        if len(T_array) != len(y_array):
            raise ValueError(
                f"Material '{self.name}', propriedade '{prop_name}': "
                "T_table e prop_table precisam ter o mesmo tamanho."
            )

        if len(T_array) < 2:
            raise ValueError(
                f"Material '{self.name}', propriedade '{prop_name}': "
                "interpolação exige pelo menos dois pontos."
            )

        if np.any(np.diff(T_array) <= 0.0):
            raise ValueError(
                f"Material '{self.name}', propriedade '{prop_name}': "
                "T_table precisa estar em ordem crescente e sem repetição."
            )

        Tmin = T_array[0]
        Tmax = T_array[-1]

        if T < Tmin or T > Tmax:
            raise ValueError(
                f"Material '{self.name}', propriedade '{prop_name}': "
                f"T={T:.3f} °C fora da tabela [{Tmin:.3f}, {Tmax:.3f}] °C."
            )

        return float(np.interp(T, T_array, y_array))

    # ========================================================
    # Atalhos
    # ========================================================

    def conductivity(self, T: float, P: Optional[float] = None) -> float:
        return self.prop("k", T, P)

    def density(self, T: float, P: Optional[float] = None) -> float:
        return self.prop("rho", T, P)

    def specific_heat(self, T: float, P: Optional[float] = None) -> float:
        return self.prop("cp", T, P)

    def viscosity(self, T: float, P: Optional[float] = None) -> float:
        return self.prop("mu", T, P)

    def describe(self) -> str:
        missing = self.missing_core_properties()
        missing_text = "nenhuma" if not missing else ", ".join(missing)

        return (
            f"Material(name={self.name!r}, phase={self.phase.value!r}, "
            f"missing_core_properties={missing_text})"
        )


# ============================================================
# Materiais padrão do projeto
# ============================================================

def criar_aluminio() -> Material:
    """
    Revestimento de alumínio:
        k_r = 180 W/(m K), constante.
    """

    return Material(
        name="Alumínio",
        phase=MaterialPhase.SOLID,
        k=180.0,
    )


def water_rho_exemplo(T: float, P: Optional[float] = None) -> float:
    return 997.0


def water_cp_exemplo(T: float, P: Optional[float] = None) -> float:
    return 4180.0


def water_mu_exemplo(T: float, P: Optional[float] = None) -> float:
    return 8.9e-4


def water_k_exemplo(T: float, P: Optional[float] = None) -> float:
    return 0.60


def criar_agua_exemplo() -> Material:
    """
    Água desmineralizada.

    No projeto:
        rho = f(T, P)
        cp  = f(T, P)
        mu  = f(T, P)
        k   = f(T, P)

    Aqui usamos funções placeholder.
    """

    return Material(
        name="Água desmineralizada",
        phase=MaterialPhase.FLUID,
        rho_func=water_rho_exemplo,
        cp_func=water_cp_exemplo,
        mu_func=water_mu_exemplo,
        k_func=water_k_exemplo,
    )


# ============================================================
# Teste rápido
# ============================================================

if __name__ == "__main__":
    aluminio = criar_aluminio()
    agua = criar_agua_exemplo()

    T = 30.0
    P = 160e3

    print("=== Teste da classe Material ===")
    print()

    for material in [aluminio, agua]:
        print(material.describe())

    print()
    print(f"Material: {aluminio.name}")
    print(f"k = {aluminio.conductivity(T):.6g} W/(m K)")

    print()
    print(f"Material: {agua.name}")
    print(f"rho = {agua.density(T, P):.6g} kg/m³")
    print(f"cp  = {agua.specific_heat(T, P):.6g} J/(kg K)")
    print(f"mu  = {agua.viscosity(T, P):.6g} Pa s")
    print(f"k   = {agua.conductivity(T, P):.6g} W/(m K)")
