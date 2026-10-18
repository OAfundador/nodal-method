"""
solver.py

Solver para NodalNetwork — regime estacionário e transiente.

Regime permanente (solve_steady_state):
    R(T) = 0
    R_i = Q_i + soma(fluxos entrando no nó i)

Transiente (solve_transient):
    C_i · dT_i/dt = R_i(T)

    Integrado por Euler implícito:
        C_i/dt · (T_i^{n+1} - T_i^n) - R_i(T^{n+1}) = 0

    Nós ARITHMETIC têm C_i = 0 e são tratados como regime permanente
    a cada passo (quasi-estático). Nós BOUNDARY mantêm temperatura fixa.

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

    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol deve ser finita e positiva.")
    if not isinstance(max_iter, (int, np.integer)) or max_iter <= 0:
        raise ValueError("max_iter deve ser um inteiro positivo.")
    if z0.ndim != 1 or len(z0) != len(net.unknown_node_ids()):
        raise ValueError("z0 deve conter uma temperatura por nó desconhecido.")
    if not np.all(np.isfinite(z0)):
        raise ValueError("z0 deve conter apenas temperaturas finitas.")
    if len(z0) == 0:
        if update_network:
            net.update_temperatures(z0)
        return SolverResult(True, "Não há temperaturas desconhecidas.", z0.copy(),
                            np.array([], dtype=float), 0.0, 0)

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
        options={"maxfev": max_iter * max(1, len(z0) + 1), "xtol": min(tol, 1e-11)},
    )

    residual = net.residual_steady(sol.x)
    residual_norm = float(np.linalg.norm(residual))

    # Algumas vezes o scipy marca "not making good progress" mesmo com
    # resíduo numericamente pequeno. Para o método nodal, aceitamos a solução
    # se a norma do resíduo estiver abaixo do critério.
    success = bool(np.isfinite(residual_norm) and residual_norm <= tol)

    # xtol do SciPy mede o passo relativo em T, não o balanço em W.
    # Refina a solução quando ela ainda não satisfaz a tolerância física pedida.
    if not success and np.all(np.isfinite(sol.x)):
        refined = _solve_with_newton_fd(net, z0=sol.x, tol=tol, max_iter=max_iter,
                                        update_network=update_network)
        refined.message = "SciPy seguido de Newton: " + refined.message
        return refined

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
        success=bool(np.isfinite(norm_r) and norm_r <= tol),
        message=("Convergiu na última atualização do Newton."
                 if norm_r <= tol else
                 "Newton com diferenças finitas atingiu max_iter sem convergir."),
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


@dataclass
class TransientResult:
    """Resultado de solve_transient."""
    t_values:   list        # [float] — instantes de tempo salvos
    snapshots:  list        # [np.ndarray] — vetor T dos nós desconhecidos em cada instante
    node_ids:   list        # [int] — IDs dos nós desconhecidos (mesma ordem dos vetores)
    success:    bool
    message:    str

    def temperatures_at(self, step: int) -> dict:
        """Retorna {node_id: T} no passo step."""
        return dict(zip(self.node_ids, self.snapshots[step]))

    def t_max_history(self) -> list:
        """Temperatura máxima ao longo do tempo."""
        return [float(s.max()) for s in self.snapshots]


def _build_capacitance_vector(net, z=None) -> np.ndarray:
    """Vetor C com ρ·cp·V para cada nó desconhecido (0 para aritméticos)."""
    C = []
    temperatures = net.temperature_map_from_vector(z) if z is not None else None
    for node_id in net.unknown_node_ids():
        node = net.nodes[node_id]
        C.append(node.heat_capacity(temperature=temperatures[node_id] if temperatures is not None else None))
    return np.array(C, dtype=float)


def solve_transient(
    net,
    *,
    t_end: float,
    dt: float,
    z0: Optional[np.ndarray] = None,
    tol: float = 1e-6,
    max_newton_iter: int = 30,
    update_network: bool = True,
    save_every: int = 1,
) -> TransientResult:
    """
    Integração transiente por Euler implícito.

    Parâmetros
    ----------
    net        : NodalNetwork com materiais, volumes e condição inicial em node.temperature.
    t_end      : tempo final [s].
    dt         : passo de tempo [s].
    z0         : condição inicial (vetor de nós desconhecidos). Se None usa temperatures atuais.
    tol        : tolerância Newton por passo.
    max_newton_iter: iterações Newton máximas por passo.
    update_network : se True, atualiza node.temperature a cada passo salvo.
    save_every : salva snapshot a cada N passos.
    """
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("dt deve ser finito e positivo.")
    if not np.isfinite(t_end) or t_end < 0:
        raise ValueError("t_end deve ser finito e não negativo.")
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol deve ser finita e positiva.")
    if not isinstance(save_every, (int, np.integer)) or save_every <= 0:
        raise ValueError("save_every deve ser um inteiro positivo.")
    if not isinstance(max_newton_iter, (int, np.integer)) or max_newton_iter <= 0:
        raise ValueError("max_newton_iter deve ser um inteiro positivo.")
    if z0 is None:
        z0 = net.initial_guess()
    z = np.asarray(z0, dtype=float).copy()
    if z.ndim != 1 or len(z) != len(net.unknown_node_ids()) or not np.all(np.isfinite(z)):
        raise ValueError("z0 deve conter uma temperatura finita por nó desconhecido.")

    node_ids = list(net.unknown_node_ids())
    step_ratio = t_end / dt
    n_steps = int(np.ceil(step_ratio))
    # Um múltiplo exato pode aparecer como 45.00000000000001 em ponto
    # flutuante. Não criar um passo extra de ~1e-17 s nesse caso.
    nearest_steps = round(step_ratio)
    if nearest_steps > 0 and abs(step_ratio - nearest_steps) <= 8 * np.spacing(step_ratio):
        n_steps = nearest_steps

    t_values  = [0.0]
    snapshots = [z.copy()]
    n_diverge = 0

    for step in range(1, n_steps + 1):
        t = t_end if step == n_steps else step * dt
        step_dt = t - (step - 1) * dt
        z_old = z.copy()

        # Resíduo transiente: F(z_new) = C/dt*(z_new - z_old) - R(z_new) = 0
        def residual_tr(z_new):
            R = net.residual_steady(z_new)
            C = _build_capacitance_vector(net, z_new)
            return C / step_dt * (z_new - z_old) - R

        # Newton com diferenças finitas
        z_new = z_old.copy()
        converged = False
        for _ in range(max_newton_iter):
            F = residual_tr(z_new)
            norm_F = float(np.linalg.norm(F))
            if norm_F <= tol:
                converged = True
                break
            J = finite_difference_jacobian(residual_tr, z_new)
            try:
                dz = np.linalg.solve(J, -F)
            except np.linalg.LinAlgError:
                dz, *_ = np.linalg.lstsq(J, -F, rcond=None)
            # busca linear simples
            for alpha in (1.0, 0.5, 0.25, 0.1):
                cand = z_new + alpha * dz
                if np.linalg.norm(residual_tr(cand)) < norm_F:
                    z_new = cand
                    break
            else:
                z_new = z_new + 0.01 * dz

        # A última atualização de Newton também pode satisfazer o critério.
        converged = converged or float(np.linalg.norm(residual_tr(z_new))) <= tol
        if not converged:
            n_diverge += 1

        z = z_new

        if step % save_every == 0 or step == n_steps:
            t_values.append(t)
            snapshots.append(z.copy())
            if update_network:
                net.update_temperatures(z)

    msg = f"{len(t_values)} snapshots, dt={dt}s, t_end={t_end}s"
    if n_diverge:
        msg += f" | AVISO: {n_diverge} passo(s) nao convergiram"

    if update_network:
        net.update_temperatures(z)

    return TransientResult(
        t_values=t_values,
        snapshots=snapshots,
        node_ids=node_ids,
        success=(n_diverge == 0),
        message=msg,
    )


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
