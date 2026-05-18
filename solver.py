"""
solver.py

Solver estacionário para NodalNetwork.

Resolve:
    R(T) = 0

onde:
    R_i = Q_i + soma(fluxos entrando no nó i)

Preferência:
- usa scipy.optimize.root se SciPy estiver disponível;
- caso contrário, usa Newton amortecido com diferenças finitas.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class SolverResult:
    success: bool
    message: str
    x: np.ndarray
    residual: np.ndarray
    residual_norm: float
    iterations: Optional[int] = None


def solve_steady_state(
    net,
    *,
    z0: Optional[np.ndarray] = None,
    tol: float = 1e-8,
    max_iter: int = 100,
    update_network: bool = True,
    prefer_scipy: bool = True,
) -> SolverResult:
    """
    Resolve a rede nodal em regime permanente.

    Parâmetros
    ----------
    net:
        NodalNetwork.

    z0:
        Chute inicial. Se None, usa net.initial_guess().

    tol:
        Tolerância na norma do resíduo.

    update_network:
        Se True, grava a solução de volta em net.nodes[*].temperature.
    """
    if z0 is None:
        z0 = net.initial_guess()
    else:
        z0 = np.asarray(z0, dtype=float)

    if prefer_scipy:
        try:
            return _solve_with_scipy(
                net,
                z0=z0,
                tol=tol,
                max_iter=max_iter,
                update_network=update_network,
            )
        except ImportError:
            pass

    return _solve_with_newton_fd(
        net,
        z0=z0,
        tol=tol,
        max_iter=max_iter,
        update_network=update_network,
    )


def _solve_with_scipy(
    net,
    *,
    z0: np.ndarray,
    tol: float,
    max_iter: int,
    update_network: bool,
) -> SolverResult:
    from scipy.optimize import root

    sol = root(
        net.residual_steady,
        z0,
        method="hybr",
        options={"maxfev": max_iter * max(1, len(z0) + 1), "xtol": tol},
    )

    residual = net.residual_steady(sol.x)
    residual_norm = float(np.linalg.norm(residual))

    # Algumas vezes o scipy marca "not making good progress" mesmo com
    # resíduo numericamente pequeno. Para o método nodal, aceitamos a solução
    # se a norma do resíduo estiver abaixo do critério.
    success = residual_norm <= max(tol * 10.0, 1e-7)

    # Sempre grava a melhor solução encontrada (igual ao Newton FD).
    # O chamador pode verificar sol.success para saber se convergiu.
    if update_network:
        net.update_temperatures(sol.x)

    return SolverResult(
        success=success,
        message=str(sol.message),
        x=np.asarray(sol.x, dtype=float),
        residual=residual,
        residual_norm=residual_norm,
        iterations=getattr(sol, "nfev", None),
    )


def _solve_with_newton_fd(
    net,
    *,
    z0: np.ndarray,
    tol: float,
    max_iter: int,
    update_network: bool,
) -> SolverResult:
    z = np.asarray(z0, dtype=float).copy()

    damping_values = [1.0, 0.5, 0.25, 0.1, 0.05, 0.01]

    for iteration in range(max_iter):
        r = net.residual_steady(z)
        norm_r = float(np.linalg.norm(r))

        if norm_r <= tol:
            if update_network:
                net.update_temperatures(z)

            return SolverResult(
                success=True,
                message="Convergiu pelo Newton com diferenças finitas.",
                x=z,
                residual=r,
                residual_norm=norm_r,
                iterations=iteration,
            )

        J = finite_difference_jacobian(net.residual_steady, z)

        try:
            dz = np.linalg.solve(J, -r)
        except np.linalg.LinAlgError:
            dz, *_ = np.linalg.lstsq(J, -r, rcond=None)

        accepted = False

        for damping in damping_values:
            candidate = z + damping * dz
            r_candidate = net.residual_steady(candidate)

            if np.linalg.norm(r_candidate) < norm_r:
                z = candidate
                accepted = True
                break

        if not accepted:
            z = z + 0.01 * dz

    r = net.residual_steady(z)
    norm_r = float(np.linalg.norm(r))

    if update_network:
        net.update_temperatures(z)

    return SolverResult(
        success=False,
        message="Newton com diferenças finitas atingiu max_iter sem convergir.",
        x=z,
        residual=r,
        residual_norm=norm_r,
        iterations=max_iter,
    )


def finite_difference_jacobian(func, z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    f0 = np.asarray(func(z), dtype=float)

    n = len(z)
    m = len(f0)
    J = np.zeros((m, n), dtype=float)

    for k in range(n):
        step = 1e-6 * max(1.0, abs(z[k]))
        zp = z.copy()
        zm = z.copy()
        zp[k] += step
        zm[k] -= step

        fp = np.asarray(func(zp), dtype=float)
        fm = np.asarray(func(zm), dtype=float)

        J[:, k] = (fp - fm) / (2.0 * step)

    return J


def print_solution_summary(net, *, max_nodes: int = 30) -> None:
    """
    Imprime uma tabela simples de temperaturas por nó.
    """
    print("=== Temperaturas dos nós ===")

    for count, (node_id, node) in enumerate(net.nodes.items()):
        if count >= max_nodes:
            print(f"... {len(net.nodes) - max_nodes} nó(s) omitido(s)")
            break

        region = "-" if node.region is None else getattr(node.region, "name", str(node.region))
        material = "-" if node.material is None else getattr(node.material, "name", str(node.material))
        print(
            f"{node_id:4d} | {node.kind.value:10s} | "
            f"region={region:15s} | mat={material:10s} | "
            f"x={node.x:.6g} | y={node.y:.6g} | T={node.temperature:.6g}"
        )
