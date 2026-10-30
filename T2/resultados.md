# Resultados — PNL cinemática do manipulador

Instância: dois elos, tempo mínimo, um disco, três pontos no elo 2, \(K=10\) intervalos, aceleração constante por peça. A colisão é exigida nos nós e nos pontos médios. Cinco partidas por método (semente 0).

Jacobiana analítica conferida por diferenças finitas no chute de desvio: erro relativo \(1{,}0\times 10^{-7}\) (igualdades) e \(2{,}2\times 10^{-7}\) (desigualdades).

Sucesso = ponto **viável** e parada de ótimo local (`gtol` / KKT / raio mínimo da região de confiança), não o teto de iterações.

## Bateria (5 chutes × 3 métodos)

Estatísticas de \(t_f\) só entre as corridas **bem-sucedidas** (viáveis e convergidas). CPU em todas as cinco.

| Método | Sucesso | \(t_f\) min (s) | \(t_f\) média (s) | \(t_f\) máx (s) | CPU média (s) |
|--------|---------|-----------------|-------------------|-----------------|---------------|
| Pontos interiores (`trust-constr` + SQP) | 4/5 | 1,192 | 1,193 | 1,194 | 6,16 |
| Pattern search (COBYLA) | 4/5 | 2,367 | 2,637 | 2,762 | 20,42 |
| Newton-SQP (`SLSQP`) | 4/5 | 1,194 | 1,194 | 1,194 | 0,75 |

**Melhor:** pontos interiores a partir do chute `desvio_ruido_1`, \(t_f=1{,}192\) s (violações de igualdade \(10^{-14}\), desigualdades nulas; mensagem `Optimization terminated successfully`).

A única partida que falha nos três métodos é `reta_juntas`: a reta no espaço das juntas atravessa o disco, e não há direção viável de descida a partir desse chute.

Detalhe dos pontos interiores:

| Chute | \(t_f\) (s) | conv | viável | iter | CPU (s) | parada |
|-------|-------------|------|--------|------|---------|--------|
| reta_juntas | 0,667 | não | não | 366 | 7,05 | maxiter |
| desvio_esquerda | 1,194 | sim | sim | 202 | 11,35 | KKT (SQP) |
| desvio_ruido_1 | 1,192 | sim | sim | 174 | 6,37 | KKT (SQP) |
| desvio_ruido_2 | 1,194 | sim | sim | 163 | 3,17 | KKT (SQP) |
| desvio_ruido_3 | 1,194 | sim | sim | 169 | 2,86 | KKT (SQP) |

## Comparação com o paper

Massaro et al. (cinemática, um disco, \(N=3\)): \(t_f=1{,}180\) s, colocação LGR e malha \(h\)–\(p\).

O melhor valor da PNL simplificada, \(1{,}192\) s, está na mesma ordem de grandeza. Não se espera coincidência: aqui há 10 intervalos com aceleração constante, o paper usa dezenas de intervalos de Gauss–Radau. O caminho ótimo local contorna o disco pela esquerda no C-space, como no caso \(N=3\) do artigo.

## Leitura dos métodos

**Pontos interiores.** Barreira Byrd–Omojokun (`trust-constr`) com Jacobiana analítica das restrições, seguida de polimento SQP. A barreira sozinha neste problema bang-bang tende a esbarrar no teto de iterações mesmo já viável; o polimento fecha o KKT. As quatro partidas viáveis convergem para o mesmo poço, \(t_f\approx 1{,}19\) s.

**Pattern search.** COBYLA (Powell, 1994): busca direta sem derivadas, restrições linearizadas numa região de confiança. Nas quatro partidas viáveis a parada é o raio mínimo da região — ótimo local do modelo linear —, não o máximo de avaliações. Sem gradiente o \(t_f\) local fica maior (\(2{,}37\) a \(2{,}76\) s): a busca direta não atravessa o feixe fino de restrições ativas que o SQP segue.

**Newton-SQP.** `SLSQP` é Newton no sistema KKT com Hessiana quasi-Newton do Lagrangiano. Fecha as quatro partidas viáveis em \(t_f=1{,}194\) s, com CPU médio abaixo de 1 s. A penalidade quadrática que se tentou antes não tem ponto estacionário quando as desigualdades estão folgadas (\(\nabla\Phi=e_{t_f}\)); o SQP trata as restrições de forma explícita e resolve isso.

## Figuras

Geradas em `T2/figuras/` a partir da melhor trajetória viável:

- `melhor.gif` — animação
- `melhor_estroboscopio.png` — configurações igualmente espaçadas no tempo
- `melhor_perfis.png` — \(q(t)\), \(\dot q(t)\), \(a(t)\)
- `melhor_cspace.png` — caminho no espaço das juntas sobre a região de colisão
- `cspace_regiao.png` — a mesma região, sem trajetória (notebook 02)
