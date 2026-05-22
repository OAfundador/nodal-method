# Método Nodal em Python

Implementação do método nodal para transferência de calor em regime estacionário. Resolve redes nodais genéricas com condução, convecção e transporte entálpico de fluido, com suporte a geometria 2D automática.

Baseado em **Transferência de Calor Computacional — Método Nodal** (J. L. Ferraz Bastos & D. A. de Andrade, IPEN/CNEN-SP).

---

## Filosofia

O REPL (`nodal_repl.py`) é **100% genérico** — não conhece nenhum domínio físico específico. O que resolver e como montar a rede é definido exclusivamente pelo arquivo `.txt` passado como argumento, usando as ferramentas genéricas do método nodal: materiais, nós, links, geometria 2D.

---

## Como usar

### 1. Instalar

```bash
git clone https://github.com/OAfundador/nodal-method.git
cd nodal-method
pip install -r requirements.txt
```

### 2. Rodar um caso pronto

```bash
# Chip eletrônico — validado contra o livro
python exemplos/nodal_repl.py exemplos/interativo_chip.txt

# Canal de refrigeração com transporte entálpico
python exemplos/nodal_repl.py exemplos/comandos_canal.txt

# Dois canais acoplados por placa combustível
python exemplos/nodal_repl.py exemplos/comandos_canal_acoplado.txt

# Troca de calor em tubo cilíndrico
python exemplos/nodal_repl.py exemplos/comandos_tubo.txt

# Projeto Final TNR5703 — referência
python exemplos/nodal_repl.py exemplos/interativo_projeto_final.txt
```

### 3. Modo interativo (prompt `>>>`)

```bash
python exemplos/nodal_repl.py
```

Digite comandos diretamente e resolva qualquer problema passo a passo:

```
>>> material aluminio k=205 rho=2700 cp=900
>>> node T1 diffusion Q=100 V=1e-4
>>> node T2 boundary T_fixed=20
>>> link T1 T2 cond G=2.5
>>> solve
  convergiu=True, |R|=3.2e-11, iter=1
  T1 = 60.00 C
  T2 = 20.00 C
>>> viz png=minha_rede.png
```

---

## Arquivos disponíveis

| Arquivo | Problema | Comandos-chave |
|---|---|---|
| `interativo_chip.txt` | Chip eletrônico Q=4W — validado contra o livro | `domain`, `region`, `source`, `bc`, `mesh`, `build_from_geom` |
| `interativo_projeto_final.txt` | Projeto Final TNR5703 — núcleo completo | `material`, `node`, `link`, `for..end`, `solve` |
| `comandos_canal.txt` | Canal aquecido com transporte entálpico | `fluid_chain`, `node`, `link` |
| `comandos_canal_acoplado.txt` | Dois canais acoplados por placa combustível | `node`, `link`, `fluid_chain` |
| `comandos_tubo.txt` | Troca de calor em tubo cilíndrico | `node`, `link`, `g_cond`, `g_conv` |

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
                [k_expr="15+0.002*T"]   # expressão simbólica em T [°C]
materiais                               # lista todos os materiais
```

Quando `k_poly` ou `k_expr` são fornecidos, a condutância é reavaliada a cada iteração Newton com a temperatura local — sem nenhuma alteração no solver.

### Geometria 2D (montagem automática da rede)

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
let <nome> = <expressao>
for <var>=<lo>..<hi>
  ...                    # use $var e $(expr) nas linhas
end
```

### Expressões dinâmicas

Em qualquer campo `G_expr="..."` ou `Q_expr="..."`:

```
dittus_h(material, T, P, mdot, A_flow, Dh)   # Dittus-Boelter h [W/m²K]
conduction_G(k, A, L)
convection_G(h, A)
sum_over("i", lo, hi, "expr")
math.*                                         # sin, exp, log, etc.
T['nome_no']                                   # temperatura atual do nó
```

### Visualização

```
viz  [png=arquivo.png] [titulo=Texto]
gif  [gif=arquivo.gif] [fps=10] [steps=40] [dpi=100]
```

---

## Validação — Chip eletrônico

Componente de Q = 4 W sobre substrato de cobre. Dois casos do Capítulo V do livro:

| Caso | Livro | Este código |
|---|---|---|
| (a) só condução | 66,3 °C | **66,34 °C** ✓ |
| (b) + convecção h=30 W/m²K | 60,266 °C | **60,27 °C** ✓ |

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
    ├── nodal_repl.py                   ← REPL genérico — ponto de entrada
    ├── interativo_chip.txt             ← chip eletrônico (validado)
    ├── interativo_projeto_final.txt    ← Projeto Final TNR5703 (referência)
    ├── comandos_canal.txt              ← canal aquecido
    ├── comandos_canal_acoplado.txt     ← dois canais acoplados
    └── comandos_tubo.txt               ← troca de calor em tubo
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
