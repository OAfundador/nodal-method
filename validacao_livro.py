"""Benchmarks de Bastos & Andrade (1996).

Fontes: capítulo V, exercício 1(a,b,e), pp. 58-64, e exercício 2, pp. 64-65.
Nenhuma consulta externa é necessária.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import argparse
import csv
import json
import math

from condutancias import conduction_G, convection_G, radiation_linearized_G
from materiais import Material
from nos import NodalNetwork, NodeKind, TransferKind
from solver import solve_steady_state, solve_transient


@dataclass(frozen=True)
class BenchmarkResult:
    name: str
    source: str
    reference_C: float
    computed_C: float
    tolerance_C: float

    @property
    def error_C(self):
        return self.computed_C - self.reference_C

    @property
    def passed(self):
        return abs(self.error_C) <= self.tolerance_C

    def to_dict(self):
        return {**asdict(self), "error_C": self.error_C, "passed": self.passed}


class BookBenchmarks:
    """Modelos pequenos cuja solução de referência vem do livro."""

    @staticmethod
    def chip_network(*, convection=False, radiation=False):
        net = NodalNetwork()
        chip = net.add_node("chip", NodeKind.DIFFUSION, 0, 0, source=4, temperature=50)
        # Preserva os 44 C/20 C do enunciado. A listagem radiativa arredonda
        # a conversão para K usando +273 e volta a C com -273; aqui usamos
        # +273.15 em radiation_linearized_G. A diferença é registrada na revisão.
        plate_C = 44
        ambient_C = 20
        plate = net.add_node("placa", NodeKind.BOUNDARY, 0, -1, fixed_temperature=plate_C)
        conductance = 12 * conduction_G(380, math.pi * 0.001 ** 2 / 4, 0.020)
        net.add_link(chip, plate, TransferKind.CONDUCTION, conductance=conductance)
        if convection or radiation:
            air = net.add_node("ar", NodeKind.BOUNDARY, 1, 0, fixed_temperature=ambient_C)
            net.add_link(chip, air, TransferKind.CONVECTION, conductance=convection_G(30, 0.030 ** 2))
        if radiation:
            wall = net.add_node("parede", NodeKind.BOUNDARY, -1, 0, fixed_temperature=ambient_C)
            def radiative_G(T):
                return radiation_linearized_G(1, 0.030 ** 2, T[chip], T[wall], sigma=5.67e-8)
            net.add_link(chip, wall, TransferKind.RADIATION, conductance_func=radiative_G)
        return net, chip

    @staticmethod
    def sphere_network():
        net = NodalNetwork()
        material = Material("esfera_livro", rho=3000, cp=1000)
        radius = 0.005
        sphere = net.add_node("esfera", NodeKind.DIFFUSION, 0, 0, material=material,
                              volume=4 / 3 * math.pi * radius ** 3, temperature=400)
        fluid = net.add_node("fluido", NodeKind.BOUNDARY, 1, 0, fixed_temperature=20)
        net.add_link(sphere, fluid, TransferKind.CONVECTION,
                     conductance=convection_G(10, 4 * math.pi * radius ** 2))
        return net, sphere

    @staticmethod
    def sphere_exact(time_s):
        # h*A/(rho*cp*V) = 3*h/(rho*cp*r) = 1/500 s^-1.
        return 20 + (400 - 20) * math.exp(-time_s / 500)

    def run(self):
        results = []
        for name, conv, rad, expected, pages in (
            ("Chip: condução", False, False, 66.3, "Ex. 1(a), pp. 58-60"),
            ("Chip: condução + convecção", True, False, 60.3, "Ex. 1(b), pp. 60-61"),
            ("Chip: condução + convecção + radiação", True, True, 59.1, "Ex. 1(e), pp. 62-64"),
        ):
            net, node = self.chip_network(convection=conv, radiation=rad)
            solution = solve_steady_state(net, tol=1e-9)
            if not solution.success:
                raise RuntimeError(f"{name}: {solution.message}")
            results.append(BenchmarkResult(name, pages, expected, net.nodes[node].temperature, 0.05))
        net, sphere = self.sphere_network()
        solution = solve_transient(net, t_end=100, dt=0.5, tol=1e-9)
        if not solution.success:
            raise RuntimeError(solution.message)
        results.append(BenchmarkResult("Esfera: resfriamento em 100 s", "Ex. 2, pp. 64-65",
                                       self.sphere_exact(100), net.nodes[sphere].temperature, 0.032))
        return results


def main(argv=None):
    parser = argparse.ArgumentParser(description="Valida o núcleo numérico com quatro casos de Bastos & Andrade.")
    parser.add_argument("--saida", type=Path, default=Path("saidas_comparacao"))
    args = parser.parse_args(argv)
    results = BookBenchmarks().run()
    args.saida.mkdir(parents=True, exist_ok=True)
    dictionaries = [result.to_dict() for result in results]
    (args.saida / "literatura.json").write_text(json.dumps(dictionaries, ensure_ascii=False, indent=2), encoding="utf-8")
    with (args.saida / "literatura.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dictionaries[0]))
        writer.writeheader()
        writer.writerows(dictionaries)
    for result in results:
        print(f"{'OK' if result.passed else 'FALHOU':7s} {result.name}: {result.computed_C:.4f} C "
              f"(referência {result.reference_C:.4f} C; erro {result.error_C:+.4f} C)")
    return 0 if all(result.passed for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
