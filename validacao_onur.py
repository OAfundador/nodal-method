"""Placa de Onur & Sivrioglu (1993), equações adimensionais 5–9.

Todos os valores desta rede são ADIMENSIONAIS: temperatura = theta = T/Ti,
tempo = Fo, comprimento = x/L. Os coeficientes unitários abaixo representam
a normalização da equação, não propriedades atribuídas a um material real.
O módulo não acessa a internet nem lê/envia o artigo durante a execução.
"""
from dataclasses import asdict, dataclass
from pathlib import Path
import csv
import json
import math

import numpy as np
from scipy.optimize import brentq

from materiais import Material
from nos import NodalNetwork, NodeKind, TransferKind
from solver import solve_steady_state, solve_transient

SOURCE = "https://doi.org/10.1007/BF01539532"


@dataclass(frozen=True)
class OnurCase:
    """Parâmetros da Fig. 2/5 por padrão; não requer escolher Ti ou k reais."""
    Bi: float = 1.
    Nr: float = 1.
    u: float = 1.
    fluid_ratio: float = .25
    environment_ratio: float = .25

    def __post_init__(self):
        values = asdict(self)
        if any(not math.isfinite(value) for value in values.values()):
            raise ValueError("Os parâmetros adimensionais devem ser finitos.")
        if self.Bi < 0 or self.Nr < 0 or self.Bi + self.Nr <= 0 or self.u < 0:
            raise ValueError("Bi e Nr devem ser não negativos, com uma troca ativa; u >= 0.")
        if self.fluid_ratio <= 0 or self.environment_ratio <= 0:
            raise ValueError("As razões de temperaturas absolutas devem ser positivas.")

    def outward_flux(self, surface_theta):
        return (self.Bi * (surface_theta - self.fluid_ratio) +
                self.Nr * (surface_theta ** 4 - self.environment_ratio ** 4))

    def analytical_steady(self, x):
        """Integração exata da EDP permanente + raiz escalar do balanço na face."""
        upper = max(1., self.fluid_ratio, self.environment_ratio)
        while self.outward_flux(upper) < self.u:
            upper *= 2
        surface = brentq(lambda theta: self.outward_flux(theta) - self.u,
                         0., upper, xtol=1e-14)
        return surface + .5 * self.u * (1 - np.asarray(x) ** 2)

    def paper_ordinate(self, theta):
        """Eixo das figuras: (T-Ti)/(Te-Ti), diferente de theta da Eq. 5."""
        if self.environment_ratio == 1:
            raise ValueError("O eixo das figuras é indefinido quando Te/Ti = 1.")
        return (np.asarray(theta) - 1) / (self.environment_ratio - 1)


class OnurSlab:
    """Volumes nodais com meia célula nas duas faces e radiação na mesma face da convecção."""
    def __init__(self, case=None, intervals=20):
        self.case = case or OnurCase()
        if not isinstance(intervals, int) or intervals < 2 or intervals > 200:
            raise ValueError("Use de 2 a 200 intervalos na placa.")
        self.x = np.linspace(0, 1, intervals + 1)
        dx = 1 / intervals
        self.weights = np.full(intervals + 1, dx)
        self.weights[[0, -1]] *= .5
        normalized = Material("Coeficientes adimensionais", k=1, rho=1, cp=1)
        self.net = NodalNetwork()
        self.ids = [self.net.add_node(f"theta_{i}", NodeKind.DIFFUSION, float(x), 0,
                    material=normalized, volume=float(weight),
                    source=self.case.u * float(weight), temperature=1.,
                    metadata={"units": "dimensionless"})
                    for i, (x, weight) in enumerate(zip(self.x, self.weights))]
        for a, b in zip(self.ids[:-1], self.ids[1:]):
            self.net.add_link(a, b, TransferKind.CONDUCTION, conductance=1 / dx)
        surface = self.ids[-1]
        fluid = self.net.add_node("Tf/Ti", NodeKind.BOUNDARY, 1, 0,
                                 fixed_temperature=self.case.fluid_ratio)
        environment = self.net.add_node("Te/Ti", NodeKind.BOUNDARY, 1, 0,
                                       fixed_temperature=self.case.environment_ratio)
        self.net.add_link(surface, fluid, TransferKind.CONVECTION, conductance=self.case.Bi)
        # Fatoração exata de theta_s^4 - theta_e^4. Não usar a rotina em Celsius.
        self.net.add_link(surface, environment, TransferKind.RADIATION,
                          conductance_func=lambda T: self.case.Nr *
                          (T[surface] + T[environment]) *
                          (T[surface] ** 2 + T[environment] ** 2))

    def steady(self):
        result = solve_steady_state(self.net, tol=1e-9, max_iter=100)
        if not result.success:
            raise RuntimeError(result.message)
        return result

    def transient(self, fo_end=2., dt=.02, *, z0=None, save_every=1):
        # Começa sempre da condição inicial publicada, salvo continuação explícita.
        initial = np.ones(len(self.ids)) if z0 is None else z0
        result = solve_transient(self.net, t_end=fo_end, dt=dt, z0=initial,
                                 tol=1e-9, save_every=save_every)
        if not result.success:
            raise RuntimeError(result.message)
        return result

    def transient_balance_error(self, history):
        """Erro de cada balanço global de Euler: dE/dFo = u - q_face."""
        times = np.asarray(history.t_values)
        values = np.asarray(history.snapshots)
        rate = np.diff(values @ self.weights) / np.diff(times)
        flux = self.case.outward_flux(values[1:, -1])
        return float(np.max(np.abs(rate - self.case.u + flux)))


def pure_convection_exact(x, fo, Bi=1., fluid_ratio=.25, terms=80):
    """Série obtida por separação de variáveis para Nr=u=0 (verificação independente)."""
    x = np.asarray(x)
    result = np.zeros_like(x, dtype=float)
    for n in range(terms):
        lower, upper = n * np.pi + 1e-12, n * np.pi + np.pi / 2 - 1e-12
        root = brentq(lambda mu: mu * np.tan(mu) - Bi, lower, upper)
        coefficient = 4 * np.sin(root) / (2 * root + np.sin(2 * root))
        result += coefficient * np.cos(root * x) * np.exp(-root ** 2 * fo)
    return fluid_ratio + (1 - fluid_ratio) * result


def combined_benchmark():
    """Caso com os três fenômenos, parâmetros publicados, sem ajustar à figura."""
    slab = OnurSlab(OnurCase(Bi=1, Nr=1, u=1), intervals=20)
    transient = slab.transient(fo_end=2, dt=.02)
    balance = slab.transient_balance_error(transient)
    numerical = slab.steady().x
    exact = slab.case.analytical_steady(slab.x)
    return {"name": "Onur & Sivrioglu (1993) · geração + convecção + radiação",
            "source": SOURCE, "parameters": asdict(slab.case),
            "theta_insulated": float(numerical[0]), "theta_surface": float(numerical[-1]),
            "steady_profile_error": float(np.max(np.abs(numerical - exact))),
            "transient_balance_error": balance, "Fo": 2., "dFo": .02,
            "transient_theta_insulated": float(transient.snapshots[-1][0]),
            "transient_theta_surface": float(transient.snapshots[-1][-1]),
            "passed": bool(np.max(np.abs(numerical - exact)) < 1e-7 and balance < 1e-7)}


def radiation_audit():
    """Fig. 8, u=5: leitura visual aproximada, sem apresentá-la como dados exatos."""
    rows = []
    for Nr, printed_approx in ((1, -3.5), (5, -4.3), (10, -4.7)):
        case = OnurCase(Bi=1, Nr=Nr, u=5)
        numerical = OnurSlab(case).steady().x
        exact = case.analytical_steady([0, 1])
        rows.append({"Bi": 1, "Nr": Nr, "u": 5, "Te_Ti": .25, "Tf_Ti": .25,
                     "theta_0": float(numerical[0]), "theta_1": float(numerical[-1]),
                     "Y_0_exact": float(case.paper_ordinate(exact[0])),
                     "Y_0_figure_8_visual_approx": printed_approx,
                     "steady_max_error": float(np.max(np.abs(numerical[[0, -1]] - exact)))})
    return rows


def reproduce_figures(output):
    """Parâmetros das Figs. 2 e 5; curvas próprias, sem digitalizar o artigo."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    colors = {1: "#087e8b", 5: "#c77136", 10: "#694c9f"}
    styles = {0: "-", 1: "--", 5: ":"}
    rows, steady_rows = [], []
    # Passos menores no início; após Fo=10 apenas a aproximação ao patamar.
    schedule = ((.001, .0001), (.01, .001), (.1, .002), (1., .01),
                (10., .05), (100., 5.), (1000., 50.))
    for Bi in (1, 5, 10):
        for u in (0, 1, 5):
            slab = OnurSlab(OnurCase(Bi=Bi, Nr=1, u=u), intervals=20)
            times, values, previous = [0.], [np.ones(21)], 0.
            for end, dt in schedule:
                history = slab.transient(end - previous, dt, z0=values[-1])
                times.extend(previous + t for t in history.t_values[1:])
                values.extend(history.snapshots[1:])
                previous = end
            values = np.asarray(values)
            ordinate = slab.case.paper_ordinate(values)
            for ax, index in zip(axes, (0, -1)):
                ax.semilogx(times[1:], ordinate[1:, index], color=colors[Bi],
                            linestyle=styles[u], label=f"Bi={Bi}, u={u}")
            numerical = slab.steady().x
            exact = slab.case.analytical_steady(slab.x)
            steady_rows.append({"Bi": Bi, "Nr": 1, "u": u, "Tf_Ti": .25, "Te_Ti": .25,
                                "theta_0": float(numerical[0]), "theta_1": float(numerical[-1]),
                                "Y_0": float(slab.case.paper_ordinate(numerical[0])),
                                "Y_1": float(slab.case.paper_ordinate(numerical[-1])),
                                "steady_max_error": float(np.max(np.abs(numerical - exact)))})
            for time, vector in zip(times, values):
                rows.append({"Bi": Bi, "Nr": 1, "u": u, "Fo": time,
                             "theta_0": float(vector[0]), "theta_1": float(vector[-1]),
                             "Y_0": float(slab.case.paper_ordinate(vector[0])),
                             "Y_1": float(slab.case.paper_ordinate(vector[-1]))})
    for ax, title in zip(axes, ("Face isolada · x/L = 0", "Face resfriada · x/L = 1")):
        ax.set(xlabel="Fo = αt/L²", ylabel="Y = (T − Ti)/(Te − Ti)", title=title, xlim=(.001, 1000))
        ax.grid(alpha=.22)
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, ncols=3, loc="lower center", bbox_to_anchor=(.5, -.13))
    fig.suptitle("Parâmetros de Onur & Sivrioglu (1993), Figs. 2 e 5 · Nr=1, Tf/Ti=Te/Ti=0,25")
    fig.savefig(output / "onur_curvas.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    for name, data in (("onur_transiente.csv", rows), ("onur_permanente.csv", steady_rows)):
        with (output / name).open("w", newline="", encoding="utf-8-sig") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    summary = {"source": SOURCE, "benchmark": combined_benchmark(), "steady": steady_rows,
               "radiation_audit": radiation_audit(),
               "note": "Curvas calculadas pelas equações publicadas. Não são pontos digitalizados das figuras."}
    (output / "onur_resumo.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    (output / "onur.html").write_text("""<!doctype html><html lang="pt-BR"><meta charset="utf-8">
<title>Placa com geração, convecção e radiação</title><style>body{font:16px/1.6 system-ui;max-width:1100px;margin:40px auto;padding:0 24px;color:#22323c}img{width:100%}code{background:#eef2f5;padding:2px 5px}</style>
<h1>Placa com geração, convecção e radiação</h1><p>Equações de Onur &amp; Sivrioglu (1993),
<a href="https://doi.org/10.1007/BF01539532">DOI 10.1007/BF01539532</a>. Parâmetros das figuras 2 e 5:
Nr=1, Bi=1/5/10, u=0/1/5 e Tf/Ti=Te/Ti=0,25. Não foi necessário escolher materiais dimensionais.</p>
<h2>O que muda entre os casos?</h2>
<p>São nove combinações do mesmo problema. <strong>Bi=hL/k</strong> compara a convecção na face com a
condução dentro da placa. Bi=1 indica escalas comparáveis; Bi=5 ou 10 intensifica a convecção, mantendo
os demais parâmetros. As cores separam esses três valores.</p>
<p><strong>u=q'''L²/(kTi)</strong> mede a geração interna de calor: u=0 significa <strong>sem geração</strong>
(linha contínua); u=1 é uma geração de referência (tracejada); u=5 representa cinco vezes essa geração
(pontilhada), para a mesma normalização. Esses números são adimensionais, não watts. Nr=1 mantém
a radiação ativa em todos os casos, inclusive u=0.</p>
<p>O gráfico da esquerda acompanha a face isolada; o da direita, a face que perde calor por convecção
e radiação. Fo é o tempo adimensional: avançar para a direita mostra a evolução até o regime permanente.
Y=0 é a temperatura inicial; Y=1 é a temperatura do ambiente; Y negativo indica temperatura acima da inicial.</p>
<img src="onur_curvas.png" alt="Curvas calculadas nas duas faces da placa">
<p>A integração usa a rede nodal do projeto, volumes com meia célula nas faces e Euler implícito.
O artigo usa diferenças finitas e Runge–Kutta de quarta ordem. As curvas acima são cálculos próprios;
não representam uma digitalização nem uma concordância quantitativa certificada com os gráficos impressos.</p>
<p>A variável das equações é θ=T/Ti; os gráficos do artigo usam Y=(T−Ti)/(Te−Ti).
Y negativo representa aquecimento acima de Ti quando Te&lt;Ti, e não temperatura absoluta negativa.</p>
<p>A verificação independente compara todo o perfil permanente com a parábola analítica e a raiz do
balanço da face: u=Bi(θs−Tf/Ti)+Nr(θs⁴−(Te/Ti)⁴). O transiente também é testado contra a série analítica
de convecção pura, refinamento temporal e fechamento da energia em cada passo.</p>
<h2>Divergência com os patamares impressos</h2>
<p>Na figura 8, com Bi=1, u=5 e ambas as razões de temperatura iguais a 0,25, uma leitura visual
aproximada dá Y(0)≈−3,5/−4,3/−4,7 para Nr=1/5/10. As equações publicadas resultam em
−3,8680/−3,2845/−3,0903. A leitura da imagem não substitui os dados numéricos originais, mas a
divergência de tendência impede declarar reprodução quantitativa dessas curvas. Não ajustamos
parâmetros nem atribuímos uma causa sem os dados/código dos autores.</p>
<p><a href="onur_transiente.csv">Históricos CSV</a> · <a href="onur_permanente.csv">Patamares permanentes CSV</a> ·
<a href="onur_resumo.json">Resultados e erros de verificação</a>. Consulte BENCHMARKS.md para critérios e limites.</p></html>""", encoding="utf-8")
    return summary


if __name__ == "__main__":
    summary = reproduce_figures(Path(__file__).resolve().parent / "saidas_comparacao")
    print(json.dumps(summary["benchmark"], ensure_ascii=False, indent=2))
