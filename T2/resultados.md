# Resultados — PNL cinemática do manipulador

Instância: dois elos, tempo mínimo, um disco, três pontos no elo 2, \(K=10\) intervalos, aceleração constante por peça. A colisão é exigida nos nós e nos pontos médios. Cinco partidas por método (semente 0).

Jacobiana analítica conferida por diferenças finitas no chute de desvio: erro relativo \(1{,}0\times 10^{-7}\) (igualdades) e \(2{,}2\times 10^{-7}\) (desigualdades).

## Bateria (5 chutes × 3 métodos)

Estatísticas de \(t_f\) só entre as corridas **viáveis**. CPU em todas as cinco.

| Método | Viáveis | \(t_f\) min (s) | \(t_f\) média (s) | \(t_f\) máx (s) | CPU média (s) |
|--------|---------|-----------------|-------------------|-----------------|---------------|
| Pontos interiores (`trust-constr`) | 5/5 | 1,400 | 2,039 | 2,484 | 9,91 |
| Pattern search (Hooke–Jeeves) | 0/5 | — | — | — | 1,41 |
| Newton–Levenberg na penalidade | 0/5 | — | — | — | 1,32 |

**Melhor viável:** pontos interiores a partir do chute `desvio_esquerda`, \(t_f=1{,}400\) s (violações de igualdade \(1{,}6\times 10^{-7}\), desigualdades nulas).

Detalhe dos pontos interiores (todos viáveis):

| Chute | \(t_f\) (s) | iter | CPU (s) | parada |
|-------|-------------|------|---------|--------|
| reta_juntas | 2,484 | 500 | 9,73 | maxiter |
| desvio_esquerda | 1,400 | 500 | 12,67 | maxiter |
| desvio_ruido_1 | 1,976 | 184 | 3,44 | gtol |
| desvio_ruido_2 | 2,098 | 500 | 11,60 | maxiter |
| desvio_ruido_3 | 2,240 | 500 | 12,13 | maxiter |

## Comparação com o paper

Massaro et al. (cinemática, um disco, \(N=3\)): \(t_f=1{,}180\) s, colocação LGR e malha \(h\)–\(p\).

O melhor valor da PNL simplificada, \(1{,}40\) s, está na mesma ordem de grandeza. Não se espera coincidência: aqui há 10 intervalos com aceleração constante, o paper usa dezenas de intervalos de Gauss–Radau. O caminho ótimo local contorna o disco pela esquerda no C-space, como no caso \(N=3\) do artigo; a reta no espaço das juntas permanece inviável.

A dispersão dos pontos interiores (\(1{,}40\) a \(2{,}48\) s) é o sintoma da **não convexidade** medida no notebook 02: o complementar do disco, puxado pela cinemática direta, parte o plano \((q_1,q_2)\) em uma região não convexa. Chutes diferentes caem em mínimos locais diferentes, todos viáveis.

## Leitura dos métodos

**Pontos interiores.** Único método que devolveu as cinco soluções viáveis. A Jacobiana analítica das restrições entra no `trust-constr`. O objetivo é linear, então a Hessiana de \(J\) é nula; o trabalho numérico está nas restrições não lineares (integração bilinear em \(t_f\) e folga trigonométrica do disco).

**Pattern search.** Hooke–Jeeves na penalidade quadrática, sem derivada, em \(\mathbb{R}^{21}\). Com o orçamento da bateria não restaurou viabilidade a partir da reta que atravessa o disco, nem a partir dos ruídos, e não melhorou o desvio viável o bastante sem voltar a violar o obstáculo. É o comportamento esperado de busca direta em dimensão moderada com restrições ativas estreitas.

**Newton.** A Hessiana analítica de \(\Phi=t_f+\rho(\|h\|^2+\|\max(g,0)\|^2)\) é correta (conferida por diferenças centrais na peça de integração). Mesmo assim o método não fechou as cinco partidas em pontos viáveis: quando as igualdades estão quase satisfeitas, \(\nabla^2\Phi\) tem posto baixo no espaço tangente; \(\max(g,0)\) não é \(C^2\) na fronteira do disco. O Levenberg regulariza o passo, mas a penalidade quadrática não substitui um tratamento explícito das desigualdades (barreira ou SQP).

## Figuras

Geradas em `T2/figuras/` a partir da melhor trajetória viável:

- `melhor.gif` — animação
- `melhor_estroboscopio.png` — configurações igualmente espaçadas no tempo
- `melhor_perfis.png` — \(q(t)\), \(\dot q(t)\), \(a(t)\)
- `melhor_cspace.png` — caminho no espaço das juntas sobre a região de colisão
- `cspace_regiao.png` — a mesma região, sem trajetória (notebook 02)
