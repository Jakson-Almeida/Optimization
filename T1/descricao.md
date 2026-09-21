# 2. Descrição e formulação do problema

A modelagem matemática adota **programação linear (PL)**, visando à **minimização do custo de compra de energia da concessionária** ao longo de um horizonte diário de 24 horas. A microrrede é representada como um sistema de **barra única**, com geração fotovoltaica e eólica (não despacháveis), armazenamento em bateria e intercâmbio com a rede principal (compra e venda).

Nesta instância, o gerador a diesel **não entra no modelo** (a variável correspondente permanece comentada no código). O custo operacional reduz-se, portanto, à energia adquirida da concessionária, valorizada por uma tarifa horária.

A formulação segue a linha de sistemas de gerenciamento de energia (EMS) para microrredes híbridas FV–eólica–bateria, em especial Luna et al. e revisões de operação ótima de microrredes [1, 2].

## 2.1 Descrição do problema

### Sistema considerado

A microrrede opera em intervalo horário $\Delta t = 1\,\mathrm{h}$, no horizonte $T = \{1,2,\ldots,24\}$. Em cada hora $t \in T$ estão disponíveis:

| Recurso | Papel no modelo | Dados de entrada |
|---------|-----------------|------------------|
| Geração fotovoltaica | Parâmetro (perfil conhecido) | $P_t^{\mathrm{sol}}$ (W) |
| Geração eólica | Parâmetro (perfil conhecido) | $P_t^{\mathrm{eol}}$ (W) |
| Bateria | Recurso despachável | capacidade, SOC inicial, limites de SOC, carga e descarga, eficiências |
| Concessionária | Compra e venda | limites de intercâmbio e tarifa horária $c_t$ |
| Carga | Demanda a atender | $P_t^{\mathrm{dem}}$ (W), perfil de demanda baixa ou alta, acrescido de uma parcela fixa |

A geração renovável é **não despachável**: o modelo deve absorver $P_t^{\mathrm{sol}}$ e $P_t^{\mathrm{eol}}$ em cada hora, seja na carga, na bateria ou na venda à rede.

**Hipóteses:** (i) operação determinística (perfis de geração e demanda conhecidos); (ii) barra única, sem perdas de rede nem restrições de tensão/fluxo; (iii) passo horário, de modo que, numericamente, potência em W coincide com energia em Wh no intervalo; (iv) eficiências de carga e descarga unitárias nesta instância; (v) a venda à rede **não gera receita** na função objetivo — serve como folga para o excedente renovável que a bateria não consegue armazenar.

### Dados de entrada da instância

- Horizonte: $n_{\mathrm{horas}} = 24$.
- Demanda: perfil variável (caso 0 = demanda baixa; caso 1 = demanda alta, usado na implementação) somado a uma carga fixa de $345\,\mathrm{W}$.
- Capacidade da bateria: $C^{\mathrm{bat}} = 16 \times 672 = 10752\,\mathrm{Wh}$.
- Estado de carga inicial (e final): $\mathrm{SOC}^{\mathrm{ini}} = 0{,}6\,C^{\mathrm{bat}}$.
- Faixa de SOC: $\mathrm{SOC}^{\min} = 0{,}5\,C^{\mathrm{bat}}$, $\mathrm{SOC}^{\max} = 0{,}9\,C^{\mathrm{bat}}$.
- Limites de carga e descarga: $E^{\mathrm{carga,max}} = E^{\mathrm{descarga,max}} = 1000\,\mathrm{W}$.
- Limites de compra e venda: $E^{\mathrm{compra,max}} = E^{\mathrm{venda,max}} = 1350\,\mathrm{W}$.
- Tarifas: $c_t = 1$ nas horas $t \in \{1,\ldots,5\} \cup \{19,\ldots,24\}$ e $c_t = 2$ nas horas $t \in \{6,\ldots,18\}$ (unidades monetárias por Wh).
- Eficiências: $\mu^{\mathrm{carga}}_t = \mu^{\mathrm{descarga}}_t = 1$.

### Decisões, objetivo prático e restrições

**Decisões (por hora):** quanto comprar e vender da/para a concessionária; quanto carregar e descarregar a bateria.

**Objetivo prático:** reduzir o gasto com energia da rede, deslocando carga da bateria para horas de tarifa baixa e usando o excedente renovável em vez de comprar nas horas caras.

**Restrições relevantes:** balanço de potência em cada hora; dinâmica e limites do SOC; limites de carga/descarga da bateria; limites de compra/venda; condição cíclica $\mathrm{SOC}_{\mathrm{final}} = \mathrm{SOC}_{\mathrm{inicial}}$.

## 2.2 Formulação matemática

### Conjuntos, índices, parâmetros e variáveis

**Índice:** $t \in T = \{1,2,\ldots,n_{\mathrm{horas}}\}$, $n_{\mathrm{horas}} = 24$.

**Parâmetros**

| Símbolo | Significado | Unidade |
|---------|-------------|---------|
| $P_t^{\mathrm{sol}}$ | Geração fotovoltaica disponível | W |
| $P_t^{\mathrm{eol}}$ | Geração eólica disponível | W |
| $P_t^{\mathrm{dem}}$ | Demanda da microrrede | W |
| $c_t$ | Tarifa de compra da concessionária | u.m./Wh |
| $\mathrm{SOC}^{\mathrm{ini}}$ | Estado de carga no início do horizonte | Wh |
| $\mathrm{SOC}^{\min}_t$, $\mathrm{SOC}^{\max}_t$ | Limites mínimo e máximo do SOC | Wh |
| $\mu^{\mathrm{carga}}_t$, $\mu^{\mathrm{descarga}}_t$ | Eficiências de carga e descarga | — |
| $E^{\mathrm{carga,min}}_t$, $E^{\mathrm{carga,max}}_t$ | Limites de carga da bateria | W |
| $E^{\mathrm{descarga,min}}_t$, $E^{\mathrm{descarga,max}}_t$ | Limites de descarga da bateria | W |
| $E^{\mathrm{compra,min}}_t$, $E^{\mathrm{compra,max}}_t$ | Limites de compra da rede | W |
| $E^{\mathrm{venda,min}}_t$, $E^{\mathrm{venda,max}}_t$ | Limites de venda à rede | W |

Nesta instância, os limites inferiores de carga, descarga, compra e venda são nulos.

**Variáveis de decisão** (domínio: reais não negativos)

| Símbolo | Significado | Unidade | Domínio |
|---------|-------------|---------|---------|
| $E_t^{\mathrm{compra}}$ | Energia comprada da concessionária | W | $[0,\,E^{\mathrm{compra,max}}_t]$ |
| $E_t^{\mathrm{venda}}$ | Energia vendida à concessionária | W | $[0,\,E^{\mathrm{venda,max}}_t]$ |
| $E_t^{\mathrm{carga,bat}}$ | Energia destinada à carga da bateria | W | $[0,\,E^{\mathrm{carga,max}}_t]$ |
| $E_t^{\mathrm{descarga,bat}}$ | Energia retirada da bateria | W | $[0,\,E^{\mathrm{descarga,max}}_t]$ |

**Variável de estado** (determinada pelas decisões)

| Símbolo | Significado | Unidade |
|---------|-------------|---------|
| $\mathrm{SOC}_t$ | Estado de carga da bateria no início da hora $t$ | Wh |

### Modelo

Nas equações abaixo, $n = n_{\mathrm{horas}} = 24$. As restrições (2b) e (3)–(8) valem no domínio indicado.

**(1)**

$$
\min \sum_{t=1}^{n} c_t E_t^{\mathrm{compra}}
$$

sujeito a:

**(2a)**

$$
\mathrm{SOC}_{1} = \mathrm{SOC}^{\mathrm{ini}}
$$

**(2b)**  $t = 1, 2, \dots, n-1$

$$
\mathrm{SOC}_{t+1}
= \mathrm{SOC}_{t}
+ \mu_{t}^{\mathrm{carga}} E_{t}^{\mathrm{carga,bat}}
- \mu_{t}^{\mathrm{descarga}} E_{t}^{\mathrm{descarga,bat}}
$$

**(2c)**

$$
\mathrm{SOC}_{n} = \mathrm{SOC}^{\mathrm{ini}}
$$

**(3)**  $t \in T$

$$
E_{t}^{\mathrm{compra}} + P_{t}^{\mathrm{eol}} + P_{t}^{\mathrm{sol}} + E_{t}^{\mathrm{descarga,bat}}
= E_{t}^{\mathrm{venda}} + E_{t}^{\mathrm{carga,bat}} + P_{t}^{\mathrm{dem}}
$$

**(4)**  $t \in T$

$$
\mathrm{SOC}_{t}^{\min} \le \mathrm{SOC}_{t} \le \mathrm{SOC}_{t}^{\max}
$$

**(5)**  $t \in T$

$$
E_{t}^{\mathrm{carga,bat,min}} \le E_{t}^{\mathrm{carga,bat}} \le E_{t}^{\mathrm{carga,bat,max}}
$$

**(6)**  $t \in T$

$$
E_{t}^{\mathrm{descarga,bat,min}} \le E_{t}^{\mathrm{descarga,bat}} \le E_{t}^{\mathrm{descarga,bat,max}}
$$

**(7)**  $t \in T$

$$
E_{t}^{\mathrm{compra,min}} \le E_{t}^{\mathrm{compra}} \le E_{t}^{\mathrm{compra,max}}
$$

**(8)**  $t \in T$

$$
E_{t}^{\mathrm{venda,min}} \le E_{t}^{\mathrm{venda}} \le E_{t}^{\mathrm{venda,max}}
$$

LaTeX para colar no Word (Alt+= → colar, sem cifrões):

- (2b) `SOC_{t+1} = SOC_t + \mu_t^{carga} E_t^{carga,bat} - \mu_t^{descarga} E_t^{descarga,bat}`
- (2c) `SOC_n = SOC^{ini}`
- (4) `SOC_t^{min} \le SOC_t \le SOC_t^{max}`

- **(1)** minimiza o custo de compra. Diferentemente de um modelo que também precifique o diesel ou a venda, só $E_t^{\mathrm{compra}}$ entra na função objetivo, ponderada por $c_t$.
- **(2a)–(2c)** descrevem a dinâmica do SOC, com condição inicial e fechamento cíclico do dia ($\mathrm{SOC}_{\mathrm{final}} = \mathrm{SOC}_{\mathrm{inicial}}$).
- **(3)** é o balanço de potência: geração renovável, descarga e compra cobrem demanda, carga da bateria e venda.
- **(4)–(8)** são os limites operacionais da bateria e do intercâmbio com a rede.

O modelo é linear, contínuo e resolvido com o solver HiGHS via Pyomo.

### Referências da formulação

[1] A. C. Luna et al., “Mixed-Integer-Linear-Programming-Based Energy Management System for Hybrid PV-Wind-Battery Microgrids: Modeling, Design, and Experimental Verification,” *IEEE Transactions on Industrial Electronics*, 2017. https://ieeexplore.ieee.org/document/7492611/

[2] Revisões de operação e gerenciamento de energia em microrredes. https://www.sciencedirect.com/science/article/pii/S1364032115000696
