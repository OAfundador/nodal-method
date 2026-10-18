# Referências públicas e verificação independente

Casos de referência para verificar o solver por soluções publicadas, balanços de energia e refinamento espacial/temporal. Os programas executam os casos localmente, sem fazer consultas externas.

## Casos já executáveis

| Caso | Fenômenos | Referência e resultado desta implementação |
| --- | --- | --- |
| NAFEMS T4 | Condução 2D + convecção | Placa de 0,6 × 1 m, k=52 W/(m·K), base a 100 °C, esquerda isolada, topo/direita com h=750 W/(m²·K) para 0 °C. No ponto E=(0,6;0,2) m: referência 18,25 °C; malha nodal 18×30: **18,33846 °C**. |
| NAFEMS T2 | Condução + radiação não linear | Barra de 0,1 m, k=55,6 W/(m·K), ponta A a 1000 K, ponta B radiando para 300 K, emissividade 0,98. Referência em B: 927 K; calculado: **927,00761 K**. |
| INL/TMAP8 ver-1fa | Condução com geração volumétrica | Placa 1D de 1,6 m, k=10 W/(m·K), geração de 10⁴ W/m³, uma face isolada e outra a 300 K. O perfil permanente é comparado célula a célula com a solução parabólica exata. Com 40 células, erro RMS de **0,2 °C**. |

Fontes primárias que fornecem os dados e as referências:

- [COMSOL: benchmark NAFEMS T4, geometria e ponto E](https://www.comsol.com/model/download/1329071/models.mph.heat_convection_2d.pdf). A documentação do [Abaqus T4](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-t4.htm) também apresenta o caso, com a referência arredondada para 18,3 °C.
- [Abaqus: NAFEMS T2](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-t2.htm).
- [Idaho National Laboratory: TMAP8 ver-1fa](https://tmap8.inl.gov/verification_and_validation/ver-1fa.html), incluindo o arquivo de entrada público.

Os resultados vêm de `validacao_online.py` e entram no botão **Benchmarks** do editor. O NAFEMS T4 também está no menu de exemplos e identifica o ponto E no desenho. As tolerâncias declaradas são 0,15 °C no T4, 0,1 °C no T2 e 0,21 °C no erro RMS do caso com geração. São critérios para estas malhas, e não limites de engenharia.

### Refinamento do NAFEMS T4

| Malha | T em E [°C] | Diferença para 18,25 °C |
| --- | ---: | ---: |
| 6 × 10 | 18,69861 | +0,44861 |
| 12 × 20 | 18,43414 | +0,18414 |
| 18 × 30 | 18,33846 | +0,08846 |

E fica na superfície convectiva, não no centro de uma célula. A avaliação usa a resistência de meia célula até a superfície e interpola ao longo das faces. Os testes verificam a redução do erro ao refinar a malha.

No caso TMAP8, a seção transversal é unitária e cancela nos balanços do perfil 1D. Verifica-se o limite permanente da solução; o programa de origem também apresenta a integração transiente. Refinar 10 → 20 → 40 células reduz o erro RMS por um fator 4 a cada etapa. A comparação usa todo o perfil, pois conferir apenas uma máxima pode ocultar erro espacial.

## Artigo com geração, convecção e radiação simultâneas

**Onur, N.; Sivrioglu, M. (1993).** *Transient heat conduction with uniform heat generation in a slab subjected to convection and radiation cooling*. Wärme- und Stoffübertragung, 28, 345–349. [Publicação e DOI](https://link.springer.com/article/10.1007/BF01539532).

O problema reúne os mecanismos solicitados: geração uniforme numa placa, isolamento em uma face e resfriamento simultâneo por convecção e radiação na outra, com evolução transiente. É uma boa continuação para avaliar não linearidade radiativa e refinamento temporal.

`validacao_onur.py` implementa as equações adimensionais 5–9 do artigo usando a rede e o integrador do projeto. O cálculo usa diretamente θ=T/Ti, x/L e Fo, sem escolher um material dimensional. A geração aparece nos volumes nodais, incluindo os meios volumes nas faces; convecção e radiação atuam na mesma temperatura da face resfriada.

O comando `python validacao_onur.py` gera `saidas_comparacao/onur.html`, gráfico e históricos CSV para os parâmetros das figuras 2 e 5: Nr=1; Bi=1,5,10; u=0,1,5; Tf/Ti=Te/Ti=0,25. As curvas usam 20 intervalos e Euler implícito com passo variável por faixas de Fo; o artigo usa Runge–Kutta de quarta ordem. Não se afirma igualdade dos dois métodos nem erro transiente certificado para todo esse gráfico. As figuras do artigo usam **Y=(T−Ti)/(Te−Ti)**, apesar de a variável das equações ser θ=T/Ti. Para Te<Ti, Y negativo representa aquecimento acima de Ti, não temperatura absoluta negativa.

As verificações independentes incluem:

- Perfil permanente exato: θ(x)=θs+u(1−x²)/2, com θs obtido pela raiz escalar de `u = Bi(θs−Tf/Ti) + Nr(θs⁴−(Te/Ti)⁴)`. Os 45 casos que combinam razões de 0,25/0,5/0,75, u=0/1/5 e os valores de Bi/Nr das figuras são conferidos em todo o perfil, com erro máximo menor que 10⁻⁸.
- Balanço global em cada passo transiente: armazenamento = geração − perdas. No caso Bi=Nr=u=1, erro máximo de aproximadamente 4,6×10⁻¹⁰, adimensional.
- Convecção pura (Nr=u=0): comparação com a série analítica obtida por separação de variáveis, cujas raízes satisfazem μ tan μ=Bi. Em Fo=0,4, refinar ΔFo=0,04 → 0,02 → 0,01 reduz o erro máximo de θ de 0,003364 → 0,001642 → 0,000812.

Para Bi=Nr=u=1 e Tf/Ti=Te/Ti=0,25, o limite permanente calculado é **θ(0)=1,31427633** e **θ(1)=0,81427633**. Os erros são adimensionais; não devem ser lidos como graus Celsius.

Para ler as nove curvas: a cor distingue Bi=1/5/10 (intensidade relativa da convecção), e o traço distingue u=0/1/5 (geração interna). Linha contínua: sem geração; tracejada: u=1; pontilhada: u=5. Bi=hL/k e u=q'''L²/(kTi), com L sendo a distância entre a face isolada e a resfriada no domínio calculado. A radiação permanece ativa, Nr=1. O texto do relatório explica essa leitura junto ao gráfico.

### Divergência que impede declarar reprodução quantitativa das figuras

A solução das equações e a comparação com os gráficos impressos são verificações distintas. Na figura 8, os patamares da face isolada com Bi=1, u=5 e razões iguais a 0,25 apresentam a seguinte diferença:

| Nr | Y(0) pelas equações, confirmado pela rede | Y(0) aproximado por leitura visual da figura 8 |
| --- | ---: | ---: |
| 1 | −3,8680 | −3,5 |
| 5 | −3,2845 | −4,3 |
| 10 | −3,0903 | −4,7 |

A terceira coluna é uma leitura aproximada da imagem, sem dados tabulados ou digitalização de precisão. A tendência também difere: pelas equações, aumentar Nr reduz θ quando a placa está mais quente que o ambiente. Isso decorre diretamente do balanço da face. Não há base para ajustar o código às curvas ou atribuir uma causa específica à diferença. **As equações foram implementadas e verificadas; a concordância quantitativa com essas figuras permanece não demonstrada.** Os dados/código originais dos autores seriam necessários para esclarecer o ponto.

O botão **Benchmarks** apresenta esse caso combinado. No editor de geometria, os contornos também permitem radiação isolada ou convecção + radiação, exigindo emissividade, fator de vista e temperaturas radiativa/do fluido explícitos. O editor geométrico continua calculando regime permanente; os transientes do artigo são executados pelo módulo de validação.
