# 3. Método de solução e implementação

Há três camadas distintas:

1. **Modelo matemático** — programação linear contínua (Seção 2): variáveis reais, objetivo e restrições lineares.
2. **Algoritmo** — método simplex (implementação revisada do solver HiGHS), que encontra o ótimo global do LP.
3. **Software** — Python + Pyomo constroem o modelo; o HiGHS o resolve. Nenhum algoritmo próprio (heurística, gradiente etc.) foi implementado.

Não há solução inicial a informar: o simplex parte de uma base construída pelo próprio solver.

## 3.1 Método de otimização

O problema é um **programa linear (PL)** de horizonte diário (24 variáveis de cada tipo, restrições de balanço e de SOC). Como a função objetivo e as restrições são lineares e as variáveis são contínuas, o conjunto viável é um poliedro e o ótimo, se existir, está em um vértice. Isso justifica o **método simplex**, que percorre bases viáveis até a otimalidade primal-dual.

O tamanho da instância é pequeno (dezenas de variáveis e restrições), de modo que um solver de PL encontra a solução exata em fração de segundo. Não se usam metaheurísticas nem métodos de ponto interior de forma explícita: o notebook chama `SolverFactory("highs")` com a configuração padrão. O HiGHS, para PL, emprega em geral o **simplex dual revisado** [1, 2].

**Passos essenciais (o que o código faz):**

1. Ler os dados: demanda (caso 0 ou 1), perfis FV e eólico, parâmetros da bateria e tarifas.
2. Declarar o `ConcreteModel` do Pyomo e o conjunto de horas \(t = 0,\ldots,23\).
3. Criar as variáveis \(E_t^{\mathrm{carga,bat}}\), \(E_t^{\mathrm{descarga,bat}}\), \(E_t^{\mathrm{compra}}\), \(E_t^{\mathrm{venda}}\) com domínio \(\mathbb{R}_+\) e limites superiores.
4. Definir a FOB: minimizar \(\sum_t c_t E_t^{\mathrm{compra}}\).
5. Incluir a dinâmica e os limites do SOC e o balanço de potência em cada hora.
6. Resolver com HiGHS até o status de otimalidade (critérios internos do solver).
7. Recuperar o valor da FOB e os perfis horários (compra, venda, carga/descarga, SOC) para os gráficos.

**Algoritmo 1 — Resolução do PL da microrrede (simplex / HiGHS)**

Entrada: \(n\), \(P_t^{\mathrm{sol}}\), \(P_t^{\mathrm{eol}}\), \(P_t^{\mathrm{dem}}\), \(c_t\), \(\mathrm{SOC}^{\mathrm{ini}}\), \(\mathrm{SOC}^{\min}_t\), \(\mathrm{SOC}^{\max}_t\), \(\mu_t^{\mathrm{carga}}\), \(\mu_t^{\mathrm{descarga}}\), limites de carga, descarga, compra e venda.

1. Construir o PL: variáveis, FOB (1) e restrições (2a)–(8).
2. Entregar o modelo ao HiGHS (simplex dual revisado, opções padrão).
3. Iterar: escolher coluna/linha da base, pivotear, atualizar os custos reduzidos e as viabilidades primal e dual, até que não haja custo reduzido melhorador (ótimo) ou se detecte ilimitado/inviável.
4. Se o status for *optimal*, ler \(E_t^{\mathrm{compra}}\), \(E_t^{\mathrm{venda}}\), \(E_t^{\mathrm{carga,bat}}\), \(E_t^{\mathrm{descarga,bat}}\) e \(\mathrm{SOC}_t\).

Saída: solução ótima das variáveis de decisão, trajetória do SOC e valor da FOB.

O passo 3 ocorre **dentro do solver**; o código do trabalho não implementa o pivoteamento. A adaptação em relação à literatura de EMS em microrredes [3] é só de modelagem (PL contínuo, sem diesel nesta instância, venda sem receita na FOB), não do algoritmo de solução.

Texto corrido para colar em 3.1:

> O modelo da Seção 2 é um programa linear contínuo. Por isso adota-se o método simplex, que garante ótimo global em um número finito de pivôs quando o PL é viável e limitado. A instância (24 horas, poucas dezenas de variáveis) é resolvida de forma exata pelo solver HiGHS, acionado via Pyomo, sem heurísticas e sem ajuste manual de parâmetros. O procedimento consiste em montar o modelo, invocá-lo o solver até a otimalidade e extrair as decisões horárias e o custo de compra.

## 3.2 Implementação computacional

Tabela 2 — Informações mínimas para reprodutibilidade

| Item | Informação a registrar |
|------|------------------------|
| Linguagem/ambiente | Python 3, notebook Jupyter (`codigo_linear_microrede.ipynb`). Desenvolvimento no Google Colab (metadado do arquivo). |
| Solver/bibliotecas | **Pyomo** (`import pyomo.environ as pyo`) para modelagem. **HiGHS** via `pyo.SolverFactory("highs")`, sem opções extras (`solver.solve(model)`). Apoio: NumPy (dados) e Matplotlib (gráficos). |
| Hardware/sistema | Não é crítico: o PL é pequeno. Registrar a máquina em que o notebook foi executado (Colab ou PC local). Em execução local típica: Windows 10. Preencher processador e RAM se o template exigir. |
| Código e dados | Repositório: https://github.com/Jakson-Almeida/Optimization — arquivo `T1/codigo_linear_microrede.ipynb`. Os dados (demanda, FV, eólica, bateria, tarifas) estão no próprio notebook, célula “Dados do problema”. Execução: abrir o `.ipynb`, instalar Pyomo e HiGHS, correr as células em ordem. Instância reportada: `caso = 1` (demanda alta + 345 W fixos). |

Trechos que sustentam a tabela:

- Modelagem: `model = pyo.ConcreteModel()` e `pyo.Var(..., domain=pyo.NonNegativeReals, bounds=(0, lim))`.
- Solver: `solver = pyo.SolverFactory("highs")` e `resultado = solver.solve(model)`.
- Indicador retornado: `print("FOB =", pyo.value(model.fob))`.

Antes de fechar o resumo, vale colar no notebook (e no texto) as versões:

```python
import sys, pyomo
print(sys.version)
print("pyomo", pyomo.__version__)
```

e, no terminal, `highs --version` (ou a versão impressa no log do `solve`).

## 3.3 Parâmetros e critérios de parada

Não há metaheurística: não se aplicam tamanho de população, número de execuções independentes nem sementes.

**Parâmetros da instância (não do algoritmo)** — os que de fato alteram a solução:

| Parâmetro | Valor usado |
|-----------|-------------|
| Horizonte | \(n = 24\) h |
| Caso de demanda | `caso = 1` (alta), mais \(345\,\mathrm{W}\) fixos |
| \(C^{\mathrm{bat}}\) | \(16 \times 672 = 10752\,\mathrm{Wh}\) |
| \(\mathrm{SOC}^{\mathrm{ini}}\) | \(0{,}6\,C^{\mathrm{bat}}\) |
| \(\mathrm{SOC}^{\min}\), \(\mathrm{SOC}^{\max}\) | \(0{,}5\,C^{\mathrm{bat}}\) e \(0{,}9\,C^{\mathrm{bat}}\) |
| \(\mu^{\mathrm{carga}}\), \(\mu^{\mathrm{descarga}}\) | 1 |
| Limites carga/descarga | \(1000\,\mathrm{W}\) |
| Limites compra/venda | \(1350\,\mathrm{W}\) |
| Tarifas \(c_t\) | 1 fora do pico; 2 nas horas \(t \in \{6,\ldots,18\}\) (índice 1–24) |

**Critérios de parada do solver:** nenhum `tee`, `timelimit`, gap ou tolerância foi passado no código. Vale o padrão do HiGHS (otimalidade quando as folgas primal e dual estão abaixo das tolerâncias internas; em PL isso corresponde a uma solução ótima, não a um gap MIP). Não há limite de iterações nem de tempo definido pelo usuário.

**Método de ajuste:** não houve calibração. Os dados operacionais vêm da instância; as opções do HiGHS permaneceram as de fábrica.

Texto corrido para colar em 3.3:

> Como o método é um solver de PL, os únicos parâmetros numéricos do estudo são os da instância (horizonte, perfis, limites da bateria e tarifas). A parada ocorre quando o HiGHS declara otimalidade, com tolerâncias padrão. Não se impôs limite de tempo nem de iterações, e não há parâmetros estocásticos a registrar.

## Referências do método

[1] G. B. Dantzig, *Linear Programming and Extensions*. Princeton University Press, 1963. (método simplex)

[2] Q. Huangfu and J. A. J. Hall, “Parallelizing the dual revised simplex method,” *Mathematical Programming Computation*, 2018. Solver HiGHS: https://highs.dev/

[3] A. C. Luna et al., “Mixed-Integer-Linear-Programming-Based Energy Management System for Hybrid PV-Wind-Battery Microgrids…,” *IEEE TIE*, 2017. (formulação EMS; neste trabalho o modelo é PL contínuo)
