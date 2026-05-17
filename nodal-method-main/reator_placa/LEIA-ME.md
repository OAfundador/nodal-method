# Aplicação: EC tipo placa — LEIA ANTES DE USAR

## Aviso importante

Esta pasta contém um exemplo **didático/demonstrativo** de aplicação do
método nodal a um elemento combustível (EC) tipo placa.

**Todos os dados numéricos aqui são FICTÍCIOS.** Dimensões, potência,
vazão, propriedades do combustível e correlações foram escolhidas com
valores redondos e plausíveis apenas para exercitar a estrutura do
método nodal e produzir resultados fisicamente coerentes.

## O que isto NÃO é

- Não é a solução de nenhum reator real.
- Não é validado contra dados experimentais ou de projeto.
- Não substitui um cálculo de termo-hidráulica para licenciamento ou
  análise de segurança.
- Não foi avaliado academicamente.

## O que isto É

Um exemplo da **metodologia** completa do método nodal aplicada a um
problema multi-físico realista (condução acoplada com k(T), convecção
forçada com h(T,P) por Dittus-Boelter, transporte entálpico
direcional), usando inteiramente a infraestrutura do framework:

- `materiais.py` para `k(T)` tabelado e propriedades da água;
- `geometria.py` (`Geometry2D`) para documentar a seção transversal;
- `nos.py` (`NodalNetwork`) para a topologia 9 nós/camada × N camadas;
- `solver.py` para Newton sobre `R(T) = 0` em ~370 incógnitas.

## Para usar com dados reais

Para aplicar a um problema concreto:

1. Edite `geometria_reator.py` com as suas dimensões.
2. Edite `propriedades_combustivel.py` com a tabela `k(T)` do seu
   combustível.
3. Edite `propriedades_agua.py` se precisar de outras correlações de
   propriedades (ou troque o fluido inteiro).
4. Edite os parâmetros em `exemplo_reator_placa.py`
   (`P_PLACA_W`, `VAZAO_CANAL_M3_S`) com os do seu caso.
5. Rode e valide contra solução analítica ou referência conhecida.

A topologia da rede nodal (em `modelo_nodal_reator.py`) deve continuar
adequada para qualquer EC tipo placa com refrigeração entre placas.

## Como rodar

A partir da pasta raiz do framework (uma acima desta):

```bash
python reator_placa/exemplo_reator_placa.py
```

Saídas em `saidas_reator_placa/`.
