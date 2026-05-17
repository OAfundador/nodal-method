# Método Nodal em Python

Implementação em Python do método nodal para problemas estacionários de transferência de calor, com suporte a redes nodais genéricas, geometria 2D, escoamento interno e modelagem de elementos combustíveis tipo placa. Baseado no livro **Transferência de Calor Computacional — Método Nodal** (J. L. Ferraz Bastos & D. A. de Andrade, IPEN/CNEN-SP).

---

## Índice

- [Como funciona](#como-funciona)
- [Quickstart](#quickstart)
- [Módulos principais](#módulos-principais)
- [REPL interativo — `nodal_repl.py`](#repl-interativo--nodal_replpy)
- [Os três problemas resolvidos](#os-três-problemas-resolvidos)
  - [Problema 1 — Chip eletrônico](#problema-1--chip-eletrônico)
  - [Problema 2 — Canal do reator (EC tipo placa)](#problema-2--canal-do-reator-ec-tipo-placa)
  - [Problema 3 — Projeto Final (núcleo completo)](#problema-3--projeto-final-núcleo-completo)
- [Referência de comandos do REPL](#referência-de-comandos-do-repl)
- [Módulo `reator_placa`](#módulo-reator_placa)
- [Estrutura do repositório](#estrutura-do-repositório)
- [Validação](#validação)
- [Stack](#stack)
- [Sobre o autor](#sobre-o-autor)

---

## Como funciona

O método nodal resolve o problema térmico em cinco passos:

```
discretização → tipos de nó → condutâncias → balanços → R(T) = 0
```

Cada nó representa um volume de controle com temperatura T. Os links entre nós carregam condutâncias de condução, convecção, radiação ou transporte entálpico de fluido. O sistema não-linear `R(T) = 0` é resolvido por Newton-Raphson com diferenças finitas para o Jacobiano.

### Tipos de nó

| Tipo | Descrição |
|---|---|
| `DIFFUSION` | Volume sólido ou fluido com fonte de calor interna |
| `ARITHMETIC` | Nó auxiliar (temperatura é média ponderada dos vizinhos) |
| `FLUID` | Nó de fluido com escoamento |
| `BOUNDARY` | Temperatura prescrita — condição de contorno |

### Tipos de link

| Tipo | Equação |
|---|---|
| `CONDUCTION` | `q = G·(Ti − Tj)`,  `G = k·A/L` |
| `CONVECTION` | `q = G·(Ti − Tj)`,  `G = h·A` |
| `FLUID_TRANSPORT` | `q = ṁ·cp·(T_upstream − T_node)` — direcional |
| `RADIATION` | `q = G·(Ti⁴ − Tj⁴)` |

---

## Quickstart

```bash
git clone https://github.com/OAfundador/nodal-method.git
cd nodal-method
pip install -r requirements.txt

# Roda o exemplo do chip diretamente
python exemplos/exemplo_chip.py

# Abre o REPL interativo
python exemplos/nodal_repl.py

# Roda um dos três arquivos de comandos prontos
python exemplos/nodal_repl.py exemplos/interativo_chip.txt
python exemplos/nodal_repl.py exemplos/interativo_reator.txt
python exemplos/nodal_repl.py exemplos/interativo_projeto_final.txt
```

---

## Módulos principais

| Módulo | Papel |
|---|---|
| `materiais.py` | `Material` com `k(T)`, `cp(T)`, `ρ(T)`, `μ(T)` — constante, polinômio, tabela ou função |
| `geometria.py` | `Geometry2D` — regiões retangulares, malha, BCs, auto-montagem da rede nodal |
| `condutancias.py` | `conduction_G`, `convection_G`, `fluid_transport_G`, série/paralelo |
| `nos.py` | `NodalNetwork` — monta o resíduo `R(T)` da rede |
| `solver.py` | `solve_steady_state` — Newton-FD com fallback `scipy.optimize.root` |
| `reator_placa/` | Modelo completo do EC tipo placa: geometria, rede nodal axial, pós-processamento |

### Exemplo de uso via API

```python
from nos import NodalNetwork, NodeKind, TransferKind
from solver import solve_steady_state

net = NodalNetwork()
Tc = net.add_node("Tc", NodeKind.DIFFUSION,  volume=1e-6, source=4.0)
Tp = net.add_node("Tp", NodeKind.BOUNDARY,   fixed_temperature=44.0)
net.add_link(Tc, Tp, TransferKind.CONDUCTION, conductance=0.17908)

result = solve_steady_state(net, tol=1e-10)
print(net.nodes[Tc].temperature)  # -> 66.34 °C
```

Não-linearidades como `k(T)` entram via `conductance_func` em vez de `conductance`.

---

## REPL interativo — `nodal_repl.py`

O arquivo central do projeto é `exemplos/nodal_repl.py`. É um REPL genérico que permite construir e resolver qualquer rede nodal interativamente, sem escrever código Python.

```bash
# Modo interativo (prompt >>>)
python exemplos/nodal_repl.py

# Modo script (executa arquivo de comandos)
python exemplos/nodal_repl.py exemplos/interativo_chip.txt
```

Qualquer comando digitado no prompt pode também estar em um arquivo `.txt` — uma linha por comando, linhas com `#` são comentários.

---

## Os três problemas resolvidos

### Problema 1 — Chip eletrônico

Componente eletrônico com geração interna de calor `Q = 4 W`, montado sobre substrato de cobre. Dois itens: (a) apenas condução, (b) com convecção natural `h = 30 W/m²K`.

**Script pronto:** `exemplos/exemplo_chip.py`  
**Entrada interativa:** `exemplos/interativo_chip.txt`

```
python exemplos/nodal_repl.py exemplos/interativo_chip.txt
```

Saídas geradas em `saidas/`:
- `chip_rede.png` — visualização da rede nodal com temperaturas
- `chip_transiente.gif` — animação do campo de temperatura (frio → regime)

| Item | Livro | Este código |
|---|---|---|
| (a) só condução | 66.3 °C | **66.34 °C** |
| (b) + convecção | 60.266 °C | **60.27 °C** |

---

### Problema 2 — Canal do reator (EC tipo placa)

Elemento combustível (EC) tipo placa do reator de pesquisa. Placa de U₃Si₂-Al com revestimento de alumínio em canal de refrigeração com água desmineralizada. Modelo nodal axial com 9 nós por camada (Tf, Ti, Tcl, Ts, Tch — top e bottom + fluido central).

**Script pronto:** `exemplos/exemplo_placa.py`  
**Entrada interativa:** `exemplos/interativo_reator.txt`

```
python exemplos/nodal_repl.py exemplos/interativo_reator.txt
```

Geometria modelada:

| Parâmetro | Valor |
|---|---|
| Lx (comprimento ativo) | 500 mm |
| df (espessura combustível) | 0,76 mm |
| dcl (espessura revestimento) | 0,38 mm |
| dch (espessura canal) | 2,89 mm |
| Lcanal (largura do canal) | 67,1 mm |
| Dh | 5,54 mm |

Saídas geradas em `saidas_reator/`:
- `temperaturas_cos.png`, `temperaturas_constante.png` — perfis axiais
- `fluxo_cos.png`, `fluxo_constante.png` — distribuição de fluxo de calor
- `animacao_cos.gif`, `animacao_constante.gif` — evolução pseudo-transiente

Correlação de Dittus-Boelter para o coeficiente de convecção:  
`h = 0,023 · (k/Dh) · Re⁰·⁸ · Pr⁰·⁴`

---

### Problema 3 — Projeto Final (núcleo completo)

Modelagem do canal mais quente do núcleo 5×5 de um reator de pesquisa de 5 MW. O núcleo tem 25 posições com fatores de potência radial distintos; o canal mais quente pertence ao EC de fator **1,914**.

**Script pronto:** `exemplos/exemplo_placa.py` (com parâmetros do Projeto Final)  
**Entrada interativa:** `exemplos/interativo_projeto_final.txt`

```
python exemplos/nodal_repl.py exemplos/interativo_projeto_final.txt
```

O arquivo de comandos faz tudo em sequência:

1. **Mapa de potência radial** — insere a matriz 5×5 manualmente e gera o PNG
2. **Geometria e configuração** do EC mais quente
3. **Resolve o canal** nos modos cosseno e fluxo constante
4. **Mapa de temperaturas** — resolve todos os 20 ECs ativos e gera mapa colorido
5. **Associação de canais** — modela os 17 canais internos + 2 externos em paralelo

#### Inserção manual do mapa de potência

```
nucleo_init rows=5 cols=5
nucleo_linha row=0 vals=1.321,1.563,0.981,1.628,1.030 tipos=EC,EC,EC,EC,EC
nucleo_linha row=1 vals=0.857,0.515,1.129,0.402,0.826 tipos=EC,CR,EC,CR,EC
nucleo_linha row=2 vals=1.050,1.877,0.000,1.914,0.979 tipos=EC,EC,CR,EC,EC
nucleo_linha row=3 vals=0.860,0.411,1.151,0.519,0.822 tipos=EC,CR,EC,CR,EC
nucleo_linha row=4 vals=0.906,1.028,0.878,1.044,0.867 tipos=EC,EC,EC,EC,EC
nucleo_mapa png=saidas_projeto_final/mapa_nucleo.png
```

#### Geometria do EC do Projeto Final

| Parâmetro | Valor |
|---|---|
| Lx | 600 mm |
| df | 0,76 mm |
| dcl | 0,38 mm |
| dch | 2,89 mm |
| Lcanal | 67,1 mm |
| Dh | 5,54 mm |
| Vazão por canal | 3,33×10⁻⁴ m³/s |
| Potência (canal mais quente) | 21 267 W |

#### Resultados

| Modo | Tf_max | Ts_max | Tch_saída |
|---|---|---|---|
| Cosseno | **73,00 °C** | 71,02 °C | 45,49 °C |
| Constante | **74,41 °C** | 72,57 °C | 45,49 °C |

Saídas geradas em `saidas_projeto_final/`:
- `mapa_nucleo.png` — mapa radial de potência (salmon = EC, azul = CR)
- `mapa_temperaturas.png` — mapa 5×5 com Tf_max de cada EC (escala plasma)
- `mapa_temperaturas.gif` — animação revelando o núcleo EC por EC
- `temperaturas_cos.png`, `temperaturas_constante.png`
- `animacao_cos.gif`, `animacao_constante.gif`
- `ec_associacao.png` — perfis axiais do canal interno vs externo lado a lado
- `ec_associacao.gif` — pseudo-transiente da associação de canais
- `resultado_cos.csv`, `resultado_constante.csv`

---

## Referência de comandos do REPL

### Materiais

```
material <nome> [phase=solid|fluid] [k=v] [rho=v] [cp=v] [mu=v]
                [k_poly=a0,a1,a2] [k_expr="15.0+0.002*T"]
materiais
```

`k` pode ser constante, polinômio em T (`k_poly`) ou expressão simbólica (`k_expr`).

### Geometria 2D (caminho automático)

```
domain W=v H=v
region <nome> [mat=<material> | k=v] x0=v x1=v y0=v y1=v
source <regiao> q=v
bc <face> <tipo> [T=v] [h=v] [T_inf=v]
mesh nx=v ny=v
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

### Condutâncias

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
             [plots=no] [gif=mapa_T.gif] [fps=4] [dir=saidas_nucleo]
```

### Reator tipo placa

```
reator_geom   Lx=v Ly=v df=v dcl=v dch=v Lcanal=v
reator_config n_axial=v vazao=v potencia=v modo=cos|constante
              T_in=v P_in=v dP=v
reator_solve
reator_show
reator_csv    [dir=.]
reator_plot   [dir=.]
reator_gif    [dir=.] [fps=12] [dpi=100]
```

### Associação de canais

```
ec_associacao [n_int=17] [n_ext=2]
              [dch_int=0.00289] [dch_ext=0.00452]
              [potencia=21267] [vazao=3.33e-4]
              [png=ec_assoc.png] [gif=ec_assoc.gif] [fps=10] [steps=40]
```

Distribui a vazão total entre canais internos e externos impondo o mesmo ΔP (Blasius turbulento: `ṁ ∝ A · Dh^(5/7)`). Canal interno aquecido por 2 faces, externo por 1 face.

---

## Módulo `reator_placa`

```
reator_placa/
├── geometria_reator.py      # GeometriaReator: dataclass com Lx, Ly, df, dcl, dch, Lcanal → Dh, A_flow
├── modelo_nodal_reator.py   # construir_rede, gerar_chute_inicial, extrair_resultados
└── exemplo_reator_placa.py  # texto_resumo, salvar_csv, plot_temperaturas, plot_fluxo
```

A rede nodal do EC tem `9 × n_axial + 1` nós:

```
Tf_top → Ti_top → Tcl_top → Ts_top
                                ↕  (convecção Dittus-Boelter)
                              Tch[j] → Tch[j+1]   (transporte entálpico)
                                ↕
Tf_bot → Ti_bot → Tcl_bot → Ts_bot
```

---

## Estrutura do repositório

```
.
├── README.md
├── LICENSE
├── requirements.txt
├── condutancias.py
├── geometria.py
├── materiais.py
├── nos.py
├── solver.py
├── reator_placa/
│   ├── geometria_reator.py
│   ├── modelo_nodal_reator.py
│   └── exemplo_reator_placa.py
└── exemplos/
    ├── nodal_repl.py              ← REPL genérico principal
    ├── interativo_chip.txt        ← comandos: problema do chip
    ├── interativo_reator.txt      ← comandos: canal do reator
    ├── interativo_projeto_final.txt ← comandos: Projeto Final completo
    ├── exemplo_chip.py
    ├── exemplo_chip_geom.py
    └── exemplo_placa.py
```

---

## Validação

### Chip eletrônico (Cap. V, Bastos & Andrade)

| Caso | Livro | `nodal_repl.py` |
|---|---|---|
| Item (a) — só condução | 66,3 °C | **66,34 °C** |
| Item (b) — + convecção h=30 | 60,266 °C | **60,27 °C** |

### EC tipo placa (dados do Projeto Final TNR5703)

| Grandeza | Cosseno | Constante |
|---|---|---|
| T_fluido_entrada | 30,00 °C | 30,00 °C |
| T_fluido_saída | 45,49 °C | 45,49 °C |
| T_superfície_max | 71,02 °C | 72,57 °C |
| T_combustível_max | **73,00 °C** | **74,41 °C** |
| Re (máx) | 16 094 | 16 094 |
| h (máx) | 10 455 W/m²K | 10 455 W/m²K |

---

## Stack

| Biblioteca | Uso |
|---|---|
| `numpy` | álgebra linear, arrays nodais |
| `scipy` | fallback para o solver Newton (`scipy.optimize.root`) |
| `matplotlib` | visualizações PNG e GIF |
| `pillow` | writer para animações GIF |

Apenas regime **estacionário**. Python 3.9+.

---

## Sobre o autor

Material desenvolvido enquanto aluno ouvinte da disciplina **TNR5703 — Análise Termo-fluido-dinâmica de Reatores Nucleares** (PPGEN/USP), com apoio de IA e revisão manual contra o livro de referência. Projeto de estudo, sem avaliação acadêmica.

Os exemplos, dados de entrada e a sequência didática derivam do livro de Bastos & Andrade e dos enunciados da disciplina. Se usar este código academicamente, cite o livro original.

---

## Licença

[MIT](LICENSE).
