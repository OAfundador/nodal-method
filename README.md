# Método Nodal em Python

Implementação do método nodal para transferência de calor em regime estacionário. Resolve redes nodais genéricas com condução, convecção e transporte entálpico de fluido, com suporte a geometria 2D automática e modelo completo de elemento combustível tipo placa para reatores de pesquisa.

Baseado em **Transferência de Calor Computacional — Método Nodal** (J. L. Ferraz Bastos & D. A. de Andrade, IPEN/CNEN-SP).

---

## O que está pronto

**Três problemas completos**, cada um com arquivo de comandos interativo pronto para rodar:

| # | Problema | Arquivo interativo | Descrição |
|---|---|---|---|
| 1 | Chip eletrônico | `interativo_chip.txt` | Condução + convecção, validado contra o livro |
| 2 | Canal do reator | `interativo_reator.txt` | EC tipo placa, Dittus-Boelter, perfis axiais |
| 3 | Projeto Final TNR5703 | `interativo_projeto_final.txt` | Núcleo 5×5, canal mais quente, associação de canais |

---

## Quickstart

```bash
git clone https://github.com/OAfundador/nodal-method.git
cd nodal-method
pip install -r requirements.txt
```

### Rodar um dos três problemas prontos

```bash
# Problema 1 — chip eletrônico
python exemplos/nodal_repl.py exemplos/interativo_chip.txt

# Problema 2 — canal do reator (placa tipo EC)
python exemplos/nodal_repl.py exemplos/interativo_reator.txt

# Problema 3 — Projeto Final completo (núcleo 5×5)
python exemplos/nodal_repl.py exemplos/interativo_projeto_final.txt
```

### Abrir o REPL interativo

```bash
python exemplos/nodal_repl.py
```

No prompt `>>>`, você define materiais, geometria, rede nodal e resolve — sem escrever código Python.

---

## Como funciona

O método discretiza o domínio em nós, cada um representando um volume de controle com temperatura T. Links entre nós carregam condutâncias de condução, convecção ou transporte entálpico. O sistema não-linear `R(T) = 0` é resolvido por Newton-Raphson.

### Tipos de nó

| Tipo | Descrição |
|---|---|
| `DIFFUSION` | Volume sólido ou fluido com fonte de calor |
| `ARITHMETIC` | Temperatura média ponderada dos vizinhos (nó auxiliar) |
| `FLUID` | Nó de fluido com escoamento |
| `BOUNDARY` | Temperatura prescrita (condição de contorno) |

### Tipos de link

| Tipo | Equação |
|---|---|
| `CONDUCTION` | `q = G·(Ti − Tj)`,  G = k·A/L |
| `CONVECTION` | `q = G·(Ti − Tj)`,  G = h·A |
| `FLUID_TRANSPORT` | `q = ṁ·cp·(T_upstream − T_node)` — direcional |
| `RADIATION` | `q = G·(Ti⁴ − Tj⁴)` |

---

## Os três problemas

### Problema 1 — Chip eletrônico

Componente eletrônico com Q = 4 W sobre substrato de cobre. Dois casos: (a) só condução, (b) com convecção h = 30 W/m²K.

**Arquivo:** `exemplos/interativo_chip.txt`

```
# Trecho do arquivo — copie e cole no REPL ou rode direto:
material silicio k=150 rho=2330 cp=712
material pernas_cobre k=3.979
domain W=0.030 H=0.030
region camada_pernas mat=pernas_cobre x0=0.0 x1=0.030 y0=0.000 y1=0.020
region componente mat=silicio x0=0.0 x1=0.030 y0=0.020 y1=0.030
source componente q=4.444e5
bc bottom temperature T=44
mesh nx=12 ny=12 fix=auto
build_from_geom tz=0.030
solve
viz png=saidas/chip_rede.png
gif gif=saidas/chip.gif fps=8 steps=35
```

| Caso | Livro | Este código |
|---|---|---|
| (a) só condução | 66,3 °C | **66,34 °C** ✓ |
| (b) + convecção h=30 | 60,266 °C | **60,27 °C** ✓ |

---

### Problema 2 — Canal do reator (EC tipo placa)

Elemento combustível tipo placa (U₃Si₂-Al + revestimento de alumínio) refrigerado com água desmineralizada. Modelo nodal axial com 9 nós por camada: combustível, interface, revestimento e canal (topo e fundo), mais o nó de fluido central.

**Arquivo:** `exemplos/interativo_reator.txt`

Compara perfil axial **cosseno** vs **fluxo constante**.

| Parâmetro | Valor |
|---|---|
| Comprimento ativo (Lx) | 500 mm |
| Espessura combustível (df) | 1,0 mm |
| Espessura revestimento (dcl) | 0,5 mm |
| Espessura canal (dch) | 3,0 mm |
| Dh | ~5,5 mm |
| Camadas axiais | 40 |

Saídas geradas ao rodar: perfis axiais de temperatura (PNG), campo 2D (GIF), CSVs.

---

### Problema 3 — Projeto Final TNR5703

Núcleo 5×5 de um reator de pesquisa de 5 MW. O arquivo faz tudo em sequência:

**Arquivo:** `exemplos/interativo_projeto_final.txt`

**Etapa 1 — Mapa de potência radial**

```
nucleo_init rows=5 cols=5
nucleo_linha row=0 vals=1.321,1.563,0.981,1.628,1.030 tipos=EC,EC,EC,EC,EC
nucleo_linha row=1 vals=0.857,0.515,1.129,0.402,0.826 tipos=EC,CR,EC,CR,EC
nucleo_linha row=2 vals=1.050,1.877,0.000,1.914,0.979 tipos=EC,EC,CR,EC,EC
nucleo_linha row=3 vals=0.860,0.411,1.151,0.519,0.822 tipos=EC,CR,EC,CR,EC
nucleo_linha row=4 vals=0.906,1.028,0.878,1.044,0.867 tipos=EC,EC,EC,EC,EC
nucleo_mapa png=saidas_projeto_final/mapa_nucleo.png
```

**Etapa 2 — Canal mais quente** (fator 1,914, posição (2,3))

```
reator_geom Lx=0.6 Ly=0.0626 df=0.00076 dcl=0.00038 dch=0.00289 Lcanal=0.0671
reator_config n_axial=40 vazao=3.33e-4 potencia=21267 modo=cos T_in=30 P_in=160000 dP=10000
reator_solve
reator_show
```

**Etapa 3 — Mapa de temperaturas do núcleo completo**

```
nucleo_solve P_total=5e6 n_canais=19 png=saidas_projeto_final/mapa_temperaturas.png gif=saidas_projeto_final/mapa_temperaturas.gif fps=4
```

Resolve todos os 20 ECs ativos em sequência e gera mapa colorido (escala plasma) com Tf_max de cada posição.

**Etapa 4 — Associação de canais em paralelo**

```
ec_associacao n_int=17 n_ext=2 dch_int=0.00289 dch_ext=0.00452 potencia=21267 vazao=3.33e-4 png=saidas_projeto_final/ec_associacao.png gif=saidas_projeto_final/ec_associacao.gif
```

Distribui vazão entre 17 canais internos e 2 externos pelo mesmo ΔP (Blasius: `ṁ ∝ A·Dh^(5/7)`).

| Grandeza | Canal interno | Canal externo |
|---|---|---|
| Espessura (dch) | 2,89 mm | 4,52 mm |
| Faces aquecidas | 2 | 1 |
| Tf_max | ~63 °C | ~42 °C |

**Resultados — canal mais quente (modo cosseno, n_axial=40)**

| Grandeza | Cosseno | Constante |
|---|---|---|
| T_fluido entrada | 30,00 °C | 30,00 °C |
| T_fluido saída | 45,49 °C | 45,49 °C |
| T_superfície_max | 71,02 °C | 72,57 °C |
| T_combustível_max | **73,00 °C** | **74,41 °C** |
| Re (máx) | 16 094 | 16 094 |
| h (máx) | 10 455 W/m²K | 10 455 W/m²K |

---

## Referência de comandos do REPL

### Materiais

```
material <nome> [phase=solid|fluid] [k=v] [rho=v] [cp=v] [mu=v]
                [k_poly=a0,a1,a2]        # k(T) = a0 + a1·T + a2·T²  (T em °C)
                [k_expr="15+0.002*T"]    # k(T) expressão simbólica   (T em °C)
materiais                                # lista todos os materiais definidos
```

Quando `k_poly` ou `k_expr` são usados, a condutância de cada link de condução é reavaliada a cada iteração Newton com a temperatura local — sem modificação no solver.

**Correlações usadas nos exemplos:**

| Material | Expressão `k_expr` | Intervalo |
|---|---|---|
| UO2 (`interativo_reator.txt`) | `"1/(0.0452+2.46e-4*(T+273.15))"` | 8,4 W/mK @ 25°C → 3,5 W/mK @ 700°C |
| U₃Si₂-Al (`interativo_projeto_final.txt`) | `"1.73073*(3978.1/(724.61+1.8*T)+6.02366e-12*(1.8*T+492)**3)"` | 8,95 W/mK @ 25°C → 7,61 W/mK @ 100°C |

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
node <nome> <tipo> [Q=v] [T_fixed=v] [T_init=v]
link <a> <b> <tipo> G=v [dir=fwd|bwd|undirected]
fluid_chain n1 n2 [n3...] mdot=v cp=v [T_in=v]
show
solve
reset
```

### Condutâncias prontas

```
g_cond k=v A=v L=v     →  G = k·A/L  [W/K]
g_conv h=v A=v          →  G = h·A    [W/K]
```

### Visualização

```
viz  [png=arquivo.png] [titulo=Texto]
gif  [gif=arquivo.gif] [fps=10] [steps=40] [dpi=100]
```

### Mapa do núcleo

```
nucleo_init [rows=5] [cols=5]
nucleo_linha row=N vals=v1,v2,... [tipos=EC,CR,...]
nucleo_mapa  [png=mapa.png] [titulo=Texto]
nucleo_solve [P_total=5e6] [n_canais=19] [png=mapa_T.png]
             [gif=mapa_T.gif] [fps=4] [dir=saidas]
```

### Reator tipo placa

```
reator_geom   Lx=v Ly=v df=v dcl=v dch=v Lcanal=v
show_geom
reator_config n_axial=v vazao=v potencia=v modo=cos|constante
              T_in=v P_in=v dP=v
reator_solve
reator_show
reator_csv    [dir=.]
reator_plot   [dir=.]
reator_gif    [dir=.] [fps=12] [dpi=100]
```

### Associação de canais em paralelo

```
ec_associacao [n_int=17] [n_ext=2]
              [dch_int=0.00289] [dch_ext=0.00452]
              [potencia=21267] [vazao=3.33e-4]
              [png=ec_assoc.png] [gif=ec_assoc.gif] [fps=10] [steps=40]
```

---

## Módulos principais

| Módulo | Papel |
|---|---|
| `materiais.py` | `Material` com `k(T)`, `cp(T)`, `ρ(T)`, `μ(T)` — constante, polinômio ou expressão |
| `geometria.py` | `Geometry2D` — regiões retangulares, malha, BCs, montagem automática da rede |
| `condutancias.py` | `conduction_G`, `convection_G`, `fluid_transport_G`, série/paralelo |
| `nos.py` | `NodalNetwork` — monta o resíduo `R(T)` da rede nodal |
| `solver.py` | `solve_steady_state` — Newton-FD com fallback `scipy.optimize.root` |

### Exemplo de uso via API Python

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

## Estrutura do repositório

```
nodal-method/
├── README.md
├── LICENSE
├── requirements.txt
├── .gitignore
├── condutancias.py          # condutâncias de condução, convecção, transporte
├── geometria.py             # Geometry2D — domínio 2D com regiões e malha
├── materiais.py             # Material com k(T), cp(T), ρ(T), μ(T)
├── nos.py                   # NodalNetwork — monta R(T)
├── solver.py                # Newton-Raphson para regime estacionário
├── reator_placa/            # módulo legado do EC tipo placa (referência)
│   ├── geometria_reator.py
│   ├── modelo_nodal_reator.py
│   ├── exemplo_reator_placa.py
│   ├── propriedades_agua.py
│   ├── propriedades_combustivel.py
│   └── LEIA-ME.md
└── exemplos/
    ├── nodal_repl.py                  ← REPL genérico — ponto de entrada principal
    ├── interativo_chip.txt            ← Problema 1: chip eletrônico
    ├── interativo_reator.txt          ← Problema 2: canal do reator
    ├── interativo_projeto_final.txt   ← Problema 3: Projeto Final TNR5703
    ├── exemplo_chip.py
    ├── exemplo_chip_geom.py
    ├── exemplo_chip_interativo.py
    └── exemplo_placa.py
```

> As pastas `saidas*/` são geradas localmente ao rodar os scripts e não são rastreadas pelo git.

---

## Stack

| Biblioteca | Uso |
|---|---|
| `numpy` | álgebra linear, arrays nodais |
| `scipy` | fallback Newton (`scipy.optimize.root`) |
| `matplotlib` | figuras PNG e animações GIF |
| `pillow` | writer para GIF |

Python 3.9+. Apenas regime **estacionário**.

---

## Sobre o autor

Desenvolvido como material de estudo da disciplina **TNR5703 — Análise Termo-fluido-dinâmica de Reatores Nucleares** (PPGEN/USP). Os exemplos, dados de entrada e sequência didática derivam do livro de Bastos & Andrade e dos enunciados da disciplina.

Se usar este código academicamente, cite o livro original.

---

[MIT License](LICENSE)
