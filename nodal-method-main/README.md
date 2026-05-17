# Método Nodal em Python

Implementação em Python do método nodal para problemas estacionários de transferência de calor. Baseado no livro **Transferência de Calor Computacional — Método Nodal** (J. L. Ferraz Bastos & D. A. de Andrade, IPEN/CNEN-SP).

## Como funciona

O método nodal resolve o problema térmico em cinco passos: discretização → tipos de nó → condutâncias → balanços → solução de `R(T) = 0`. Cada passo está em um módulo:

| Módulo | Papel |
| --- | --- |
| `materiais.py` | `Material` com `k(T)`, `cp(T)`, `ρ(T)`, `μ(T)` por tabela, polinômio ou função |
| `geometria.py` | `Geometry2D` para descrever regiões, malha e BCs em 2D |
| `condutancias.py` | `conduction_G`, `convection_G`, `fluid_transport_G`, série/paralelo |
| `nos.py` | `NodalNetwork` — monta o resíduo `R(T)` da rede |
| `solver.py` | `solve_steady_state` — Newton-FD com fallback `scipy.optimize.root` |

```python
net = NodalNetwork()
T1 = net.add_node("T1", NodeKind.DIFFUSION, volume=V, source=Q)
T2 = net.add_node("T2", NodeKind.BOUNDARY, fixed_temperature=20.0)
net.add_link(T1, T2, TransferKind.CONDUCTION, conductance=G)

result = solve_steady_state(net, tol=1e-10)
print(net.nodes[T1].temperature)
```

Não-linearidades como `k(T)` e `h(T,P)` entram via `conductance_func` em vez de `conductance`.

## Quickstart

```bash
git clone https://github.com/<usuario>/<repo>.git
cd <repo>
pip install -r requirements.txt
python exemplos/exemplo_chip.py
```

## Exemplos

| Arquivo | O que faz |
| --- | --- |
| `exemplos/exemplo_chip.py` | Exemplo 1 do Cap. V do livro (componente eletrônico) — caminho lumped, 2-3 nós |
| `exemplos/exemplo_chip_geom.py` | Mesmo exemplo via `Geometry2D` com condutividade efetiva calibrada |
| `exemplos/exemplo_placa.py` | Placa 2D com fonte interna — caminho `Geometry2D` distribuído |
| `exemplos/exemplo_chip_interativo.py` | REPL para construir redes via comandos |

### Validação contra o livro

| Caso | Livro | Lumped (`exemplo_chip.py`) | Geometry2D (`exemplo_chip_geom.py`) |
| --- | --- | --- | --- |
| Item (a) sem convecção | 66.3 °C | 66.34 °C | 66.44 °C (T_méd) |
| Item (b) h=30 W/m²K | 60.266 °C | 60.27 °C | 60.32 °C (T_méd) |

## REPL interativo

```bash
python exemplos/exemplo_chip_interativo.py --interativo
```

```text
>>> g_cond k=380 A=7.854e-7 L=0.02
  G_cond = k·A/L = 0.0149226 W/K
>>> node Tc diffusion Q=4
>>> node Tp boundary T_fixed=44
>>> link Tc Tp cond G=0.17908
>>> solve
  Tc = 66.3364 °C
```

Comandos: `node`, `link`, `g_cond`, `g_conv` (caminho lumped); `domain`, `region`, `bc`, `source`, `mesh`, `build_from_geom` (caminho `Geometry2D`); `show`, `solve`, `reset`. Pode também rodar scripts de comandos:

```bash
python exemplos/exemplo_chip_interativo.py --script exemplos/comandos_chip.txt
```

## Estrutura

```
.
├── README.md
├── LICENSE
├── requirements.txt
├── materiais.py
├── geometria.py
├── condutancias.py
├── nos.py
├── solver.py
└── exemplos/
    ├── exemplo_chip.py
    ├── exemplo_chip_geom.py
    ├── exemplo_chip_interativo.py
    ├── exemplo_placa.py
    ├── comandos_chip.txt
    ├── comandos_chip_geom.txt
    └── comandos_placa.txt
```

## Stack

`numpy`, `scipy.optimize.root`. Apenas regime estacionário.

## Sobre o autor

Material desenvolvido enquanto aluno ouvinte da disciplina TNR5703 (PPGEN/USP, *Análise Termo-fluido-dinâmica de Reatores Nucleares*), com apoio de IA e revisão manual contra o livro de referência. Projeto de estudo, sem avaliação acadêmica.

Os exemplos, dados de entrada e a sequência didática derivam do livro de Bastos & Andrade. Se usar este código academicamente, cite o livro original.

## Licença

[MIT](LICENSE).
