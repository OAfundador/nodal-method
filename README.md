# Método Nodal em Python

Ferramenta em Python para montar e resolver redes térmicas com condução, convecção, radiação e transporte entálpico de fluido. Inclui editor visual de regiões retangulares, montagem automática da malha 2D e um núcleo para regime permanente e transiente.

A versão atual é a branch padrão [`codex/revisao-reproducao`](https://github.com/OAfundador/nodal-method/tree/codex/revisao-reproducao). A implementação anterior está preservada em [`legado`](https://github.com/OAfundador/nodal-method/tree/legado).

Baseado em **Transferência de Calor Computacional — Método Nodal** (J. L. Ferraz Bastos & D. A. de Andrade, IPEN/CNEN-SP).

---

## Editor visual local

Instale as dependências uma vez:

```bash
python -m pip install -r requirements.txt
```

No Windows, abra **`iniciar_visual.bat`**, ou execute:

```bash
python app_visual.py
```

O navegador abre em `http://127.0.0.1:8765`. Mantenha o terminal aberto enquanto usa a ferramenta. Os cálculos usam o solver Python deste repositório; a interface não depende de serviços externos, CDNs ou publicação.

1. Em **Materiais**, cadastre nome, fase e propriedades. Para condução estacionária, `k` é obrigatório; `rho`, `cp` e `mu` são opcionais. Também é possível informar os coeficientes de `k(T)`.
2. Selecione um material e use **Retângulo**: clique e arraste para desenhar. Também pode arrastar um cartão de material para a área de desenho. Use **Selecionar** para mover; os quatro cantos permitem redimensionar. O painel à direita aceita coordenadas e dimensões exatas.
3. Encoste as regiões, sem sobreposição. O encaixe na grade ajuda a alinhar as bordas. Defina a geração de calor de cada região, se houver.
4. Em **Contornos**, configure temperatura, convecção, radiação, convecção + radiação, fluxo imposto ou isolamento nas faces externas do domínio. Para radiação, informe emissividade, fator de vista e temperatura radiativa; no caso combinado, informe também h e a temperatura do fluido. Em **Malha**, escolha `nx` e `ny` e clique em **Gerar malha**. As faces verdes indicam o contato entre regiões.
5. Clique em **Resolver modelo** para obter o mapa e o perfil de temperaturas. **Exportar CSV** salva as temperaturas dos centros das células. **Salvar modelo** e **Abrir modelo** preservam/restauram os dados em JSON.

A interface mostra comprimentos em **mm**, temperaturas em **°C** e propriedades em SI. O JSON e o núcleo numérico usam comprimentos em **m**. Áreas vazias são omitidas e suas bordas internas são isoladas; os contornos configuráveis pertencem às faces externas do domínio. O editor calcula regime permanente, com até 900 células e 1600 incógnitas. Uma região desconectada de qualquer referência térmica é rejeitada.

Para regiões fluidas, informe `cp`, vazão **mássica**, direção e temperatura de entrada no painel da região. Na interface sólido-fluido, informe também `h` interno. O transporte é distribuído entre as cadeias paralelas da região. Não há cálculo automático de vazão, queda de pressão ou coeficiente de convecção a partir de correlações no editor.

O menu oferece **parede uniforme**, **parede composta**, a **geometria 2D do chip** e o **NAFEMS T4** com ponto de avaliação na superfície. O botão **Benchmarks** executa os quatro casos do livro, os casos públicos NAFEMS/TMAP8 e a placa com geração, convecção e radiação de Onur & Sivrioglu, descritos em [BENCHMARKS.md](BENCHMARKS.md). As equações desse artigo foram verificadas, mas há divergências com alguns patamares de seus gráficos. Atalhos: `V` seleciona, `R` desenha, `Del` exclui a região selecionada e `Ctrl+Z` desfaz uma alteração. A ferramenta também permite apagar ou redimensionar uma região pelo painel de propriedades.

Se a porta estiver ocupada, use `python app_visual.py --porta 8766`. Para iniciar sem abrir o navegador, acrescente `--sem-abrir`.

## Executar os benchmarks e testes

```bash
python -m pip install -r requirements.txt
python validacao_livro.py
python validacao_online.py
python validacao_onur.py
python -m unittest discover -s tests -v
```

Os benchmarks verificam os balanços e a redução do erro ao refinar a malha ou o passo de tempo. A suíte também cobre os exemplos antigos, entradas inválidas, exportação e a API local.

`validacao_onur.py` gera `saidas_comparacao/onur.html`, curvas e históricos adimensionais das duas faces, com os parâmetros das figuras 2 e 5 de Onur & Sivrioglu (1993). Também é possível gerar o relatório pelo botão **Benchmarks** do editor. A análise explicita diferenças entre as equações e a figura 8, sem ajustar parâmetros para forçar uma reprodução.

## Exemplo de saída: placa de reator

Resultados calculados pela rede nodal para um caso local com 51 pontos axiais, no canal 14. A saída corresponde ao plano z=−0,300 m.

| Perfil de geração | Máxima no combustível [K] | Fluido na saída [°C] | Combustível na saída [°C] |
| --- | ---: | ---: | ---: |
| Cossenoidal | 359,3492 | 45,0149 | 45,0149 |
| Constante | 347,5215 | 45,0200 | 74,3715 |

O [CSV dos resultados](resultados/placa_reator.csv) inclui as máximas e as temperaturas de saída dos quatro componentes. Os arquivos de entrada desse caso permanecem locais. A tabela registra a saída calculada sob as hipóteses desse modelo; as verificações reproduzíveis com os dados publicados estão em [BENCHMARKS.md](BENCHMARKS.md).

## Melhorias desta versão

- Editor com cadastro de materiais, desenho, movimento e redimensionamento de regiões, malha conectada e mapa de temperaturas.
- Contornos de temperatura, fluxo, isolamento, convecção, radiação e convecção + radiação na mesma superfície.
- Execução de arquivos organizada em `CaseRunner`, com mensagens de erro por linha e exportação CSV.
- Opções de tolerância e iteração do REPL aplicadas ao solver; convergência conferida pelo resíduo do balanço.
- Transiente com último passo ajustado ao tempo final, tratamento do arredondamento decimal e capacitâncias reavaliadas durante a solução.
- Testes de regressão, balanço de energia e refinamento contra soluções analíticas e referências públicas.

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

Use `--quiet` para reduzir a saída, `--csv arquivo.csv` para exportar os nós e `--listar` para listar os exemplos. Erros em arquivos interrompem a execução, com nome do arquivo e número da linha, e retornam código de saída diferente de zero.

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
| `visual/` | Modelo validado, servidor local e interface do editor |
| `app_visual.py` | Inicia o editor no navegador |
| `validacao_*.py` | Casos de referência e geração dos relatórios |

### Uso direto via Python

```python
from nos import NodalNetwork, NodeKind, TransferKind
from solver import solve_steady_state

net = NodalNetwork()
Tc = net.add_node("Tc", NodeKind.DIFFUSION, x=0, y=0, volume=1e-6, source=4.0)
Tp = net.add_node("Tp", NodeKind.BOUNDARY, x=1, y=0, fixed_temperature=44.0)
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
solve tol=1e-8 max_iter=500 linhas=20
export csv=temperaturas.csv
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

Componente de Q = 4 W sobre substrato de cobre. Três casos do Capítulo V do livro:

| Caso | Livro | Este código |
|---|---|---|
| (a) só condução, modelo concentrado | 66,3 °C | **66,3375 °C** |
| (b) + convecção h=30 W/m²K | 60,3 °C | **60,2663 °C** |
| (e) + radiação | 59,1 °C | **59,0782 °C** |

Esses valores usam a rede concentrada do exercício 1 do livro. O arquivo `interativo_chip.txt` e o exemplo visual distribuem a geração em uma geometria 2D e resultam em **66,4854 °C**: não são a mesma discretização do modelo concentrado. O quarto benchmark do livro é a esfera em resfriamento (exercício 2), comparada à solução exponencial exata e com refinamento temporal.

---

## Estrutura do repositório

```
nodal-method/
├── app_visual.py                      ← abre o editor visual
├── iniciar_visual.bat                 ← atalho para Windows
├── condutancias.py
├── geometria.py
├── materiais.py
├── nos.py
├── solver.py
├── requirements.txt
├── BENCHMARKS.md                      ← fontes, resultados e critérios
├── validacao_livro.py
├── validacao_online.py
├── validacao_onur.py
├── visual/                            ← interface e API local
├── tests/                             ← suíte de testes
├── resultados/                        ← exemplos de saídas calculadas
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

`numpy`, `scipy`, `matplotlib`, `pillow`. Python 3.9+.

O editor e os exemplos geométricos usam regime estacionário. A API `solve_transient` e o comando `solve_transient` do REPL também oferecem Euler implícito, com capacitâncias reavaliadas durante o cálculo.

---

## Sobre

Desenvolvido como material de estudo da disciplina **TNR5703 — Análise Termo-fluido-dinâmica de Reatores Nucleares** (PPGEN/USP), com apoio de IA e revisão manual contra o livro de referência e os enunciados da disciplina.

Se usar este código academicamente, cite o livro original: Bastos & Andrade, IPEN/CNEN-SP.

[MIT License](LICENSE)
