# Método Nodal em Python

Implementação do método nodal para transferência de calor em regime estacionário. Resolve redes nodais genéricas com condução, convecção e transporte entálpico de fluido, com suporte a geometria 2D automática.

Baseado em **Transferência de Calor Computacional — Método Nodal** (J. L. Ferraz Bastos & D. A. de Andrade, IPEN/CNEN-SP).

---

## Filosofia

O REPL (`nodal_repl.py`) é **100% genérico** — não conhece nenhum domínio físico específico. O que resolver e como montar a rede é definido exclusivamente pelo arquivo `.txt` passado como argumento, usando as ferramentas genéricas do método nodal (materiais, nós, links, geometria 2D).

```
python exemplos/nodal_repl.py exemplos/meu_caso.txt
```

---

## Quickstart

```bash
git clone https://github.com/OAfundador/nodal-method.git
cd nodal-method
pip install -r requirements.txt
python exemplos/nodal_repl.py exemplos/interativo_chip.txt
```

---

## Módulos principais

| Módulo | Papel |
|---|---|
| `materiais.py` | `Material` com `k(T)`, `cp(T)`, `ρ(T)`, `μ(T)` — constante, polinômio ou expressão simbólica |
| `geometria.py` | `Geometry2D` — domínio retangular com regiões, fontes, BCs e malha |
| `condutancias.py` | `conduction_G`, `convection_G`, série/paralelo |
| `nos.py` | `NodalNetwork` — monta o resíduo `R(T)` e gerencia os links |
| `solver.py` | `solve_steady_state` — Newton-FD com fallback `scipy.optimize.root` |

### Uso direto via Python

```python
from nos import NodalNetwork, NodeKind, TransferKind
from solver import solve_steady_state

net = NodalNetwork()
Tc = net.add_node("Tc", NodeKind.DIFFUSION, volume=1e-6, source=4.0)
Tp = net.add_node("Tp", NodeKind.BOUNDARY,  fixed_temperature=44.0)
net.add_link(Tc, Tp, TransferKind.CONDUCTION, conductance=0.17908)

result = solve_steady_state(net, tol=1e-10)
print(net.nodes[Tc].temperature)  # → 66.34 °C
```

---

## Referência de comandos do REPL

### Materiais

```
material <nome> [phase=solid|fluid] [k=v] [rho=v] [cp=v] [mu=v]
                [k_poly=a0,a1,a2]       # k(T) = a0 + a1·T + a2·T²
                [k_expr="15+0.002*T"]   # k(T) expressão em T [°C]
materiais                               # lista todos os materiais
```

Quando `k_poly` ou `k_expr` são fornecidos, a condutância é reavaliada a cada iteração Newton com a temperatura local do link — sem nenhuma alteração no solver.

### Geometria 2D (montagem automática)

```
domain W=v H=v
region <nome> [mat=<material> | k=v] x0=v x1=v y0=v y1=v
source <regiao> q=v                      # fonte volumétrica [W/m³]
bc <face> <tipo> [T=v] [h=v] [T_inf=v]  # face: top/bottom/left/right
mesh nx=v ny=v [fix=auto]
show_geom
build_from_geom [tz=v]
```

### Rede nodal direta

```
node <nome> <tipo> [Q=v] [T_fixed=v] [T_init=v] [V=v]
     tipos: diffusion | arithmetic | fluid | boundary
link <a> <b> <tipo> G=v [dir=fwd|bwd|undirected]
     tipos: cond | conv | fluid | rad
fluid_chain n1 n2 [n3...] mdot=v cp=v [T_in=v]
show
solve
reset / reset_net
```

### Condutâncias prontas

```
g_cond k=v A=v L=v     →  G = k·A/L  [W/K]
g_conv h=v A=v          →  G = h·A    [W/K]
```

### Variáveis e loops

```
let <nome> = <expressao>       # define variável escalar ou lista
for <var>=<lo>..<hi>
  ...                          # linhas podem usar $var e $(expr)
end
```

### Expressões dinâmicas

Em qualquer campo `X_expr="..."` você pode usar:
- `math.*`, `T` (temperatura do nó em iteração)
- `dittus_h(material, T, P, mdot, A_flow, Dh)` — correlação de Dittus-Boelter
- `conduction_G(k,A,L)`, `convection_G(h,A)`
- `sum_over("var", lo, hi, "expr")`

### Visualização

```
viz  [png=arquivo.png] [titulo=Texto]
gif  [gif=arquivo.gif] [fps=10] [steps=40] [dpi=100]
```

---

## Exemplo validado — Chip eletrônico

Componente de Q = 4 W sobre substrato de cobre. Dois casos: (a) só condução, (b) com convecção h = 30 W/m²K.

**Arquivo:** `exemplos/interativo_chip.txt`

| Caso | Livro | Este código |
|---|---|---|
| (a) só condução | 66,3 °C | **66,34 °C** ✓ |
| (b) + convecção h=30 | 60,266 °C | **60,27 °C** ✓ |

---

## Estrutura do repositório

```
nodal-method/
├── condutancias.py
├── geometria.py
├── materiais.py
├── nos.py
├── solver.py
├── requirements.txt
└── exemplos/
    ├── nodal_repl.py               ← REPL genérico — ponto de entrada
    ├── interativo_chip.txt         ← chip eletrônico (validado)
    ├── interativo_reator.txt       ← EC tipo placa (referência)
    ├── interativo_projeto_final.txt← Projeto Final TNR5703 (referência)
    ├── exemplo_chip.py
    ├── exemplo_chip_geom.py
    └── exemplo_placa.py
```

> As pastas `saidas*/` são geradas localmente e não são rastreadas pelo git.

---

## Stack

`numpy`, `scipy`, `matplotlib`, `pillow`. Python 3.9+. Apenas regime **estacionário**.

---

## Sobre

Desenvolvido como material de estudo da disciplina **TNR5703 — Análise Termo-fluido-dinâmica de Reatores Nucleares** (PPGEN/USP), com apoio de IA e revisão manual contra o livro de referência e os enunciados da disciplina.

Se usar este código academicamente, cite o livro original: Bastos & Andrade, IPEN/CNEN-SP.

[MIT License](LICENSE)
