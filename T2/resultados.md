# Resultados — PNL cinemática do manipulador

Instância: dois elos, tempo mínimo, um disco, três pontos no elo 2, \(K=50\) intervalos, aceleração constante por peça. A colisão é exigida nos nós e nos pontos médios. Cinco partidas por método (semente 0). Decisão \(\mathbf{z}\in\mathbb{R}^{101}\), 4 igualdades e 507 desigualdades.

Jacobiana analítica conferida por diferenças finitas no chute de desvio: erro relativo \(1{,}0\times 10^{-8}\) (igualdades) e \(2{,}2\times 10^{-7}\) (desigualdades).

Sucesso = ponto **viável** e parada de ótimo local (`gtol` / KKT / raio mínimo da região de confiança), não o teto de iterações.

## Bateria (5 chutes × 3 métodos)

Estatísticas de \(t_f\) só entre as corridas **bem-sucedidas** (viáveis e convergidas). CPU em todas as cinco.

| Método | Sucesso | \(t_f\) min (s) | \(t_f\) média (s) | \(t_f\) máx (s) | CPU média (s) |
|--------|---------|-----------------|-------------------|-----------------|---------------|
| Pontos interiores (`trust-constr` + SQP) | 4/5 | 1,181 | 1,181 | 1,181 | 25,99 |
| Pattern search (COBYLA) | 0/5 | — | — | — | 35,00 |
| Newton-SQP (`SLSQP`) | 4/5 | 1,181 | 1,181 | 1,181 | 11,06 |

**Melhor:** pontos interiores a partir do chute `desvio_ruido_2`, \(t_f=1{,}181\) s (violações de igualdade \(10^{-15}\), desigualdades nulas; mensagem `Optimization terminated successfully`).

A única partida que falha nos três métodos é `reta_juntas`: a reta no espaço das juntas atravessa o disco, e não há direção viável de descida a partir desse chute.

Detalhe dos pontos interiores:

| Chute | \(t_f\) (s) | conv | viável | iter | CPU (s) | parada |
|-------|-------------|------|--------|------|---------|--------|
| reta_juntas | 0,566 | não | não | 453 | 35,94 | maxiter |
| desvio_esquerda | 1,181 | sim | sim | 364 | 18,91 | KKT (SQP) |
| desvio_ruido_1 | 1,181 | sim | sim | 393 | 18,68 | KKT (SQP) |
| desvio_ruido_2 | 1,181 | sim | sim | 881 | 37,61 | KKT (SQP) |
| desvio_ruido_3 | 1,181 | sim | sim | 396 | 18,81 | KKT (SQP) |

## Comparação com o paper

Massaro et al. (cinemática, um disco, \(N=3\)): \(t_f=1{,}180\) s, colocação LGR e malha \(h\)–\(p\).

Com \(K=50\) o melhor valor desta transcrição, \(1{,}181\) s, coincide com o paper na milésima. A malha mais fina (antes \(K=10\), \(t_f\approx 1{,}192\) s) aproxima o bang-bang contínuo o bastante para apagar a folga da discretização. O caminho ótimo local contorna o disco pela esquerda no C-space, como no caso \(N=3\) do artigo.

## Leitura dos métodos

**Pontos interiores.** Barreira Byrd–Omojokun (`trust-constr`) com Jacobiana analítica, seguida de polimento SQP. Nas quatro partidas viáveis fecha KKT em \(t_f\approx 1{,}181\) s.

**Pattern search.** COBYLA em \(\mathbb{R}^{101}\) com 507 restrições não chega ao raio mínimo no orçamento da bateria: permanece viável perto do chute (\(\approx 5{,}4\) s) e para em `MAXFUN`. Sem derivadas, a malha fina torna o feixe de restrições ativas estreito demais para a busca direta.

**Newton-SQP.** `SLSQP` fecha as quatro partidas viáveis no mesmo poço, \(t_f=1{,}181\) s, com CPU médio de 11 s. A reta inviável de novo falha (derivada direcional positiva na busca linear).

## Figuras

Geradas em `T2/figuras/` a partir da melhor trajetória viável:

- `melhor.gif` — animação
- `melhor_estroboscopio.png` — configurações igualmente espaçadas no tempo
- `melhor_perfis.png` — \(q(t)\), \(\dot q(t)\), \(a(t)\)
- `melhor_cspace.png` — caminho no espaço das juntas sobre a região de colisão
- `cspace_regiao.png` — a mesma região, sem trajetória (notebook 02)
