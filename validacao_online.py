"""Benchmarks públicos NAFEMS; execução inteiramente local.

Os URLs registram a origem dos dados. Este programa não acessa a internet.
T4: COMSOL, heat_convection_2d; T2: documentação de benchmarks do Abaqus.
"""
from pathlib import Path
import json
import numpy as np

from condutancias import conduction_G, radiation_linearized_G
from nos import NodalNetwork, NodeKind, TransferKind
from solver import solve_steady_state
from validacao_livro import BenchmarkResult
from visual.modelo import VisualModel, nafems_t4_case

T4_SOURCE = "https://www.comsol.com/model/download/1329071/models.mph.heat_convection_2d.pdf"
T2_SOURCE = "https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-t2.htm"
GENERATION_SOURCE = "https://tmap8.inl.gov/verification_and_validation/ver-1fa.html"


class OnlineBenchmarks:
    @staticmethod
    def generated_slab(cells=20):
        """TMAP8 ver-1fa, limite permanente; seção unitária cancela no perfil."""
        length, k, generation, wall_C = 1.6, 10., 1e4, 300 - 273.15
        dx = length / cells
        net = NodalNetwork()
        ids = [net.add_node(f"slab_{i}", NodeKind.DIFFUSION, (i + .5) * dx, 0,
                            source=generation * dx, temperature=wall_C) for i in range(cells)]
        for a, b in zip(ids[:-1], ids[1:]):
            net.add_link(a, b, TransferKind.CONDUCTION, conductance=k / dx)
        wall = net.add_node("parede", NodeKind.BOUNDARY, length, 0, fixed_temperature=wall_C)
        net.add_link(ids[-1], wall, TransferKind.CONDUCTION, conductance=2 * k / dx)
        solution = solve_steady_state(net, tol=1e-7)
        if not solution.success:
            raise RuntimeError(solution.message)
        positions = np.array([net.nodes[node].x for node in ids])
        computed = np.array([net.nodes[node].temperature for node in ids])
        exact = wall_C + generation / (2 * k) * (length ** 2 - positions ** 2)
        return positions, computed, exact

    @staticmethod
    def radiation_bar():
        net = NodalNetwork()
        start = net.add_node("A", NodeKind.BOUNDARY, 0, 0, fixed_temperature=1000 - 273.15)
        end = net.add_node("B", NodeKind.ARITHMETIC, .1, 0, temperature=650)
        ambient = net.add_node("ambiente", NodeKind.BOUNDARY, .1, 0, fixed_temperature=300 - 273.15)
        area = .01 * .01
        net.add_link(start, end, TransferKind.CONDUCTION, conductance=conduction_G(55.6, area, .1))
        net.add_link(end, ambient, TransferKind.RADIATION, conductance_func=lambda T:
                     radiation_linearized_G(.98, area, T[end], T[ambient], sigma=5.67e-8))
        result = solve_steady_state(net, tol=1e-9)
        if not result.success:
            raise RuntimeError(result.message)
        return net.nodes[end].temperature

    def run(self):
        t4 = VisualModel(nafems_t4_case()).solve()["probes"][0]["T_C"]
        t2 = self.radiation_bar()
        _, generated, exact = self.generated_slab(40)
        rms = float(np.sqrt(np.mean((generated - exact) ** 2)))
        return [BenchmarkResult("NAFEMS T4 · ponto E", T4_SOURCE, 18.25, t4, .15),
                BenchmarkResult("NAFEMS T2 · ponta radiativa", T2_SOURCE, 927 - 273.15, t2, .1),
                BenchmarkResult("TMAP8 ver-1fa · erro RMS do perfil", GENERATION_SOURCE, 0., rms, .21)]


def main():
    results = OnlineBenchmarks().run()
    path = Path(__file__).resolve().parent / "saidas_comparacao" / "benchmarks_online.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps([r.to_dict() for r in results], ensure_ascii=False, indent=2), encoding="utf-8")
    for result in results:
        print(f"{'OK' if result.passed else 'REVER'} {result.name}: {result.computed_C:.5f} C; "
              f"referência {result.reference_C:.5f} C; diferença {result.error_C:+.5f} C")
    return 0 if all(r.passed for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
