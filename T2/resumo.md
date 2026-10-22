# Resumo — planejamento de trajetória de tempo mínimo por controle ótimo

**Fonte.** Massaro, M.; Lovato, S.; Bottin, M.; Rosati, G. *An Optimal Control Approach to the Minimum-Time Trajectory Planning of Robotic Manipulators*. Robotics 2023, 12, 64. https://doi.org/10.3390/robotics12030064

O problema é levar um manipulador da configuração inicial à final no menor tempo, com o caminho, a lei temporal e os controles obtidos numa única otimização. A formulação é um problema de controle ótimo (OCP) não linear na forma de Bolza. Um método direto transcreve esse OCP para um problema de programação não linear (PNL), resolvido por colocação pseudoespectral e um solver de ponto interior.

O texto abaixo concentra a modelagem, as restrições e a PNL. Os tempos de manobra aparecem só para ligar cada conjunto de restrições ao efeito sobre a solução.

## 1. Controle ótimo na forma de Bolza

O sistema é descrito por equações diferenciais ordinárias

$$
\dot{\mathbf{x}}(t) = \mathbf{f}\big(\mathbf{x}(t), \mathbf{u}(t)\big), \qquad t_i \le t \le t_f,
$$

com estado $\mathbf{x}$, controle $\mathbf{u}$ e horizonte $[t_i, t_f]$. O índice de desempenho é

$$
J = M\big(\mathbf{x}(t_i), \mathbf{x}(t_f)\big) + \int_{t_i}^{t_f} \ell\big(\mathbf{x}(t), \mathbf{u}(t)\big)\, dt,
$$

em que $M$ é o custo de Mayer e $\ell$ é o custo de Lagrange. Minimiza-se $J$ sujeito à dinâmica, às restrições de caminho

$$
\mathbf{c}\big(\mathbf{x}(t), \mathbf{u}(t)\big) \le \mathbf{0},
$$

às condições de contorno

$$
\mathbf{b}\big(\mathbf{x}(t_i), \mathbf{x}(t_f)\big) = \mathbf{0}
$$

e aos limites de controle $|\mathbf{u}(t)| \le \mathbf{u}_{\max}$.

Neste trabalho, $J$ é o tempo de manobra ou a distância percorrida no espaço das juntas. As restrições de caminho limitam velocidade, aceleração ou torque e afastam o robô dos obstáculos. A solução é o caminho, o cronograma do movimento e os controles.

Para tempo mínimo com $t$ como variável independente usa-se $M = 0$ e $\ell = 1$, ou, de forma equivalente, $M = t_f$ e $\ell = 0$. O instante final $t_f$ é livre.

## 2. Transcrição para programação não linear

O OCP contínuo é transcrito para uma PNL por um método direto. A trajetória é aproximada por polinômios e as equações do problema são impostas em pontos de colocação. Os pontos são raízes de polinômios de Legendre e a aproximação usa polinômios de Lagrange, no esquema de Legendre–Gauss–Radau (LGR). Estados e controles de todos os intervalos são obtidos de uma só vez, como solução de um sistema algébrico.

A transcrição é feita no GPOPS-II. A PNL resultante é resolvida pelo IPOPT, um algoritmo de ponto interior com busca linear e filtro, adequado a PNL esparsa de grande porte.

### 2.1 Por que o problema é não linear

A não linearidade entra por quatro vias, conforme a variante:

1. **Horizonte livre.** O tempo final é variável de decisão. Na colocação, o intervalo físico é ligado a um intervalo normalizado, e a dinâmica escalada pelo comprimento do horizonte produz produtos entre $t_f$ e $\mathbf{f}(\mathbf{x}, \mathbf{u})$.
2. **Cinemática e obstáculos.** A posição cartesiana de um ponto do elo depende de senos e cossenos dos ângulos. A folga em relação a um disco é uma desigualdade quadrática nessas coordenadas, não convexa (o complementar de um disco).
3. **Dinâmica do manipulador.** A matriz de inércia $\mathbf{M}(\mathbf{q})$ depende de $\cos q_2$ e $\sin q_2$. O vetor de Coriolis $\mathbf{c}(\mathbf{q}, \dot{\mathbf{q}})$ é quadrático nas velocidades. A aceleração $\mathbf{M}(\mathbf{q})^{-1}(\mathbf{Q} - \mathbf{c})$ é uma função racional e trigonométrica do estado.
4. **Regularização de caminho.** O termo de distância usa $\sqrt{\dot{q}_1^2 + \dot{q}_2^2}$. No planejamento em distância percorrida, a tangente unitária $(q_1')^2 + (q_2')^2 = 1$ é quadrática.

O IPOPT é um solver local. Com obstáculos, a PNL tem vários mínimos locais, e a solução numérica depende da estimativa inicial. O artigo usa várias estimativas para conferir as soluções apresentadas e não afirma otimalidade global.

### 2.2 Malha $h$–$p$ e derivadas

Em cada intervalo de malha de comprimento $h$ há $p$ pontos de colocação. Métodos $h$ movem as fronteiras dos intervalos com $p$ fixo. Métodos $p$ variam a ordem do polinômio numa malha fixa. O esquema usado varia os dois: se o erro num intervalo é grande e distribuído, aumenta-se $p$; se o erro se concentra num ponto, o intervalo é subdividido. Todas as otimizações partem de 50 intervalos iguais, com 4 pontos de colocação em cada um. A solução típica converge em 2 a 5 refinos de malha, com tolerância de malha $10^{-5}$ no GPOPS-II e tolerância $10^{-7}$ no IPOPT.

As derivadas da função objetivo e das restrições em relação a estados e controles são obtidas por diferenciação automática (ADiGator), com precisão de máquina, sem derivação simbólica manual nem diferenças finitas.

### 2.3 Escala numérica

A convergência depende da escala das variáveis. Todas são adimensionalizadas com uma massa de referência $m_0 = m_1 + m_2 + M$, um comprimento $l_0 = l_1 + l_2$ e uma grandeza de referência da variante:

- restrições cinemáticas: aceleração $\ddot{q}_0 = \max(\ddot{q}_{1,\max}, \ddot{q}_{2,\max})$, e o tempo escala com $1/\sqrt{\ddot{q}_0}$;
- restrições dinâmicas: torque $Q_0 = \max(Q_{1,\max}, Q_{2,\max})$, e o tempo escala com $\sqrt{m_0 l_0^2 / Q_0}$.

### 2.4 Porte da PNL

| Variante | Equações (ordem de grandeza) | Tempo de CPU |
|----------|------------------------------|--------------|
| Só cinemática | 1000–1600 | 1–2 s |
| Com dinâmica | 1000–5000 | 2–30 s |

Máquina de referência: PC com Intel i7 e 16 GB de RAM. O tempo de CPU fica entre 1 e 30 vezes o tempo da manobra. Obstáculos múltiplos, via-points e o robô de quatro elos ocupam o extremo superior. Monitorar $N = 10$ pontos de colisão custa cerca de 50% a mais do que $N = 3$. Com obstáculos, a estimativa inicial precisa ser quase viável, isto é, quase sem colisão.

## 3. Modelo do manipulador de dois elos

O exemplo base é um manipulador planar de dois graus de liberdade. $q_1$ é a rotação do primeiro elo em relação à base e $q_2$ é a rotação do segundo elo em relação ao primeiro. Os dados coincidem com os de Ghasemi, Kashiri e Dardel (2012).

| Símbolo | Descrição | Valor |
|---------|-----------|-------|
| $l_1$, $l_2$ | Comprimentos | 0,4 m e 0,25 m |
| $m_1$, $m_2$ | Massas dos elos | 29,58 kg e 15 kg |
| $I_1$, $I_2$ | Inércias em relação ao centro de massa | 0,417 kg·m² e 0,206 kg·m² |
| $(b_1, h_1)$ | Centro de massa do elo 1, no elo | $(l_1/2,\ 0)$ |
| $(b_2, h_2)$ | Centro de massa do elo 2, no elo | $(l_2/2,\ 0)$ |
| $M$ | Massa da carga paga | 6 kg |

Nesta instância $h_1 = h_2 = 0$, então os termos em $h$ da dinâmica se anulam. O modelo geral abaixo os mantém.

### 3.1 Cinemática

Um ponto $P$ do segundo elo, com parâmetro $\zeta \in [0, 1]$, tem coordenadas

$$
\begin{aligned}
x_P(\mathbf{q}, \zeta) &= l_1 \cos q_1 + \zeta\, l_2 \cos(q_1 + q_2), \\
y_P(\mathbf{q}, \zeta) &= l_1 \sin q_1 + \zeta\, l_2 \sin(q_1 + q_2).
\end{aligned}
$$

$\zeta = 0$ é a junta entre os elos e $\zeta = 1$ é o efetuador. A mesma expressão, com $\zeta = i/N$, fornece os pontos em que a restrição de colisão é escrita.

As condições de contorno são impostas nas juntas, $\mathbf{q}(t_i) = \mathbf{q}_i$ e $\mathbf{q}(t_f) = \mathbf{q}_f$. Se a tarefa viesse em coordenadas cartesianas, as juntas correspondentes sairiam da cinemática inversa. Fixar $\mathbf{q}_i$ e $\mathbf{q}_f$ escolhe ainda a configuração do manipulador entre as que colocam o efetuador no mesmo ponto. No exemplo,

$$
q_{1i} = q_{2i} = 0, \qquad q_{1f} = 1\ \mathrm{rad}, \qquad q_{2f} = -0{,}5\ \mathrm{rad},
$$

com manobras repouso–repouso: $\dot{\mathbf{q}}(t_i) = \dot{\mathbf{q}}(t_f) = \mathbf{0}$.

### 3.2 Dinâmica

As equações vêm do Lagrangeano $L = T - V$. O movimento é planar e horizontal, então $V = 0$. A energia cinética é

$$
T = \frac{1}{2} m_1 |\mathbf{v}_{G_1}|^2 + \frac{1}{2} I_1 \dot{q}_1^2 + \frac{1}{2} m_2 |\mathbf{v}_{G_2}|^2 + \frac{1}{2} I_2 (\dot{q}_1 + \dot{q}_2)^2 + \frac{1}{2} M |\mathbf{v}_{P_2}|^2.
$$

Os torques nas juntas são $\mathbf{Q} = \{Q_1, Q_2\}^\mathsf{T}$. As equações de Lagrange ficam na forma compacta

$$
\mathbf{M}(\mathbf{q})\, \ddot{\mathbf{q}} + \mathbf{c}(\mathbf{q}, \dot{\mathbf{q}}) = \mathbf{Q},
$$

com

$$
\mathbf{M}(\mathbf{q}) =
\begin{bmatrix}
m_{11} & m_{12} \\
m_{12} & m_{22}
\end{bmatrix}.
$$

Os elementos, com o centro de massa do elo 2 em $(b_2, h_2)$, são

$$
\begin{aligned}
m_{11} &= I_1 + I_2 + m_1(b_1^2 + h_1^2) \\
&\quad + m_2\big(l_1^2 + b_2^2 + h_2^2 + 2 l_1 b_2 \cos q_2 - 2 l_1 h_2 \sin q_2\big) \\
&\quad + M\big(l_1^2 + l_2^2 + 2 l_1 l_2 \cos q_2\big), \\[4pt]
m_{12} &= I_2 + m_2\big(b_2^2 + h_2^2 + l_1 b_2 \cos q_2 - l_1 h_2 \sin q_2\big) \\
&\quad + M\big(l_2^2 + l_1 l_2 \cos q_2\big), \\[4pt]
m_{22} &= I_2 + m_2(b_2^2 + h_2^2) + M l_2^2.
\end{aligned}
$$

O vetor $\mathbf{c}$ reúne os termos quadráticos de velocidade. Com $h_2 = 0$ na instância numérica, ele se reduz a

$$
\mathbf{c}(\mathbf{q}, \dot{\mathbf{q}}) =
l_1 (m_2 b_2 + M l_2) \sin q_2
\begin{Bmatrix}
-\dot{q}_2 (2\dot{q}_1 + \dot{q}_2) \\
\dot{q}_1^2
\end{Bmatrix}.
$$

Na forma usada pelo OCP dinâmico, as acelerações são o controle inverso da dinâmica,

$$
\begin{Bmatrix} \dot{v}_{q_1} \\ \dot{v}_{q_2} \end{Bmatrix}
= \mathbf{M}(\mathbf{q})^{-1} \big(\mathbf{Q} - \mathbf{c}(\mathbf{q}, \dot{\mathbf{q}})\big),
$$

acopladas às cadeias de integradores $\dot{q}_1 = v_{q_1}$ e $\dot{q}_2 = v_{q_2}$. É essa expressão, e não um limite de aceleração constante, que torna o OCP dinâmico não linear no estado.

## 4. Restrições e variantes do problema

Todas as variantes compartilham o repouso–repouso e os limites de velocidade

$$
|\dot{q}_1| \le \dot{q}_{1,\max} = 3\ \mathrm{rad/s}, \qquad
|\dot{q}_2| \le \dot{q}_{2,\max} = 8\ \mathrm{rad/s}.
$$

O que muda é o vetor de estado, o controle, a dinâmica $\mathbf{f}$, o custo e o conjunto de desigualdades.

### 4.1 Planejamento cinemático de tempo mínimo

| Item | Definição |
|------|-----------|
| Estado | $\mathbf{x} = \{q_1, q_2, v_{q_1}, v_{q_2}\}^\mathsf{T}$ |
| Controle | $\mathbf{u} = \{a_{q_1}, a_{q_2}\}^\mathsf{T}$ (acelerações) |
| Dinâmica | $\dot{q}_i = v_{q_i}$, $\dot{v}_{q_i} = a_{q_i}$ |
| Caminho | limites de velocidade |
| Controle | $\|a_{q_i}\| \le \ddot{q}_{i,\max}$ |
| Contorno | posições inicial e final, velocidades nulas |
| Custo | $\displaystyle \int_{t_i}^{t_f} 1\, dt$ |

Os limites de aceleração vêm do coeficiente de velocidade $c_v = \dot{q}_{\max} / (\Delta q / T)$. Para um perfil trapezoidal simétrico, $c_v = T / (T - T_a)$ e

$$
\ddot{q}_{\max} = \frac{\dot{q}_{\max}^2}{(c_v - 1)\, \Delta q}.
$$

Com o valor típico $c_v = 1{,}5$ obtêm-se 18 rad/s² na junta 1 e 128 rad/s² na junta 2. Adota-se o menor dos dois para as duas juntas: $\ddot{q}_{1,\max} = \ddot{q}_{2,\max} = 18$ rad/s². O trapézio serve apenas para calibrar esse número. A otimização recupera o trapézio na junta gargalo; ele não é imposto como lei de movimento.

Há infinitas soluções de tempo mínimo. A junta mais exigida, aqui $q_1$, satura velocidade e aceleração e fixa $t_f = 0{,}500$ s. A outra junta pode seguir qualquer perfil cuja integral de velocidade feche $\Delta q_2$ exatamente nesse instante. A solução clássica é o segmento de reta no plano $(q_1, q_2)$, com velocidades trapezoidais. Outras estimativas iniciais produzem caminhos curvos no espaço das juntas, com o mesmo $t_f$ e com o trapézio preservado só em $q_1$.

### 4.2 Distância mínima no espaço das juntas

Troca-se a variável independente para a distância $s$ percorrida no plano $(q_1, q_2)$, com custo de Mayer $M = s_f$. As derivadas passam a ser em relação a $s$ (indicadas por $'$) e a dinâmica continua uma cadeia de integradores. Impõe-se a tangente unitária

$$
(q_1')^2 + (q_2')^2 = 1,
$$

além dos limites de velocidade. A solução numérica é a reta entre os extremos. Esse problema planeja só o caminho: o tempo de execução fica para uma etapa posterior de timing.

Para escolher, entre os tempos mínimos, o de menor distância, soma-se um termo penalizado ao custo temporal:

$$
J = \int_0^{t_f} \Big(1 + w \sqrt{\dot{q}_1^2 + \dot{q}_2^2}\Big)\, dt.
$$

Com $w = 0$ recupera-se o tempo mínimo puro. Com $w > 0$ o solver fica com a reta no espaço das juntas. Essa forma é a usada sempre que há obstáculo e se quer evitar a multiplicidade de soluções de mesmo $t_f$.

### 4.3 Obstáculos

O obstáculo de referência é um disco no plano cartesiano, centrado em $(x_0, y_0) = (0{,}45,\ 0{,}25)$ m, com raio $R = 0{,}1$ m. A folga é exigida em $N$ pontos do segundo elo:

$$
\big(x_P(\mathbf{q}, \zeta) - x_0\big)^2 + \big(y_P(\mathbf{q}, \zeta) - y_0\big)^2 > R^2,
\qquad \zeta = \frac{i}{N},\quad i \in \{1, \ldots, N\}.
$$

São desigualdades de caminho, não lineares e não convexas, porque $x_P$ e $y_P$ são trigonométricas em $\mathbf{q}$. Formas mais complexas são uniões de discos ou métodos de colisão mais gerais. O primeiro elo pode receber as mesmas desigualdades; no exemplo elas ficam no segundo.

O inteiro $N$ muda a solução:

- $N = 1$ vigia só o efetuador. O caminho no espaço das juntas permanece a reta, $t_f = 0{,}500$ s, e o corpo do elo 2 atravessa o disco. Isso só seria aceitável num SCARA se o obstáculo fosse mais baixo que o elo.
- $N = 3$ desvia o elo pela esquerda do disco, o caminho em $(q_1, q_2)$ deixa de ser reto e o tempo sobe para $t_f = 1{,}180$ s.

Para $C$ discos a mesma desigualdade vale para cada centro $(x_{0j}, y_{0j})$ e cada raio $R_j$, com $j \in \{1, \ldots, C\}$. No exemplo de três discos os centros e raios são $(0{,}52,\ 0{,}12)$ m com $0{,}1$ m, $(0{,}45,\ 0{,}2)$ m com $0{,}075$ m e $(0{,}45,\ 0{,}3)$ m com $0{,}125$ m.

### 4.4 Controle pelo jerk

Os tempos mínimos com aceleração como controle têm saltos de aceleração no início, no fim e nas quebras do trapézio. Para eliminá-los, o controle passa a ser o jerk $\mathbf{u} = \{j_{q_1}, j_{q_2}\}^\mathsf{T}$ e o estado ganha as acelerações:

$$
\mathbf{x} = \{q_1, q_2, v_{q_1}, v_{q_2}, a_{q_1}, a_{q_2}\}^\mathsf{T},
\qquad \dot{a}_{q_i} = j_{q_i},
$$

junto com a cadeia de integradores já usada. Mantêm-se os limites de velocidade e de aceleração e acrescentam-se

$$
|j_{q_1}| \le 500\ \mathrm{rad/s^3}, \qquad |j_{q_2}| \le 200\ \mathrm{rad/s^3},
$$

com acelerações nulas nos dois extremos, $\ddot{\mathbf{q}}(t_i) = \ddot{\mathbf{q}}(t_f) = \mathbf{0}$. O custo continua sendo o de tempo com peso na distância. Com o obstáculo e $N = 3$, o perfil fica próximo do caso anterior, arredondado onde o jerk satura, e o tempo passa a $t_f = 1{,}286$ s.

### 4.5 Planejamento com a dinâmica

Os limites de aceleração são substituídos por limites de torque, e a velocidade continua limitada:

$$
|Q_1| \le Q_{1,\max} = 25\ \mathrm{N \cdot m}, \qquad
|Q_2| \le Q_{2,\max} = 9\ \mathrm{N \cdot m}.
$$

| Item | Definição |
|------|-----------|
| Estado | $\mathbf{x} = \{q_1, q_2, v_{q_1}, v_{q_2}\}^\mathsf{T}$ |
| Controle | $\mathbf{u} = \{Q_1, Q_2\}^\mathsf{T}$ |
| Dinâmica | integradores de $\mathbf{q}$ e $\dot{\mathbf{v}}_q = \mathbf{M}(\mathbf{q})^{-1}(\mathbf{Q} - \mathbf{c})$ |
| Caminho | limites de velocidade e de torque |
| Contorno | posições e velocidades de repouso–repouso |
| Custo | tempo de manobra |

A estrutura bang–bang dos torques é resultado da otimização, não uma hipótese. Sem obstáculo, nenhuma junta atinge a velocidade máxima. Há dois arcos em $Q_1$ e três em $Q_2$, o caminho em $(q_1, q_2)$ é curvo e $t_f = 1{,}002$ s. As acelerações locais chegam a cerca de 18 rad/s² em $q_1$ e 59 rad/s² em $q_2$, acima do teto cinemático da seção 4.1, porque aqui quem satura é o torque.

Forçar a reta de distância mínima no espaço das juntas aumenta o tempo para $1{,}081$ s. Só a junta rápida permanece bang–bang, e as acelerações máximas caem para cerca de 3,5 rad/s² e 1,8 rad/s². A reta, ótima no problema cinemático, deixa de sê-lo quando $\mathbf{M}(\mathbf{q})$ e $\mathbf{c}(\mathbf{q}, \dot{\mathbf{q}})$ entram no modelo.

Sem a carga ($M = 0$), os mesmos torques máximos movem menos inércia: $t_f = 0{,}843$ s, com acelerações de pico próximas de 18 rad/s² e 69 rad/s². Sobre a reta no espaço das juntas o tempo sobe para $0{,}921$ s.

Com o disco e $N = 1$, ao contrário do caso cinemático, a solução já muda em relação ao problema sem obstáculo: $t_f = 1{,}046$ s, e o elo 2 ainda colide. Com $N = 3$ o efetuador passa pela esquerda, $t_f = 1{,}098$ s, e a estrutura passa a três arcos em cada torque.

### 4.6 Limite na taxa de variação do torque

O bang–bang tem saltos de torque. Para removê-los, o controle passa a ser a taxa $\mathbf{u} = \{J_1, J_2\}^\mathsf{T}$ e o torque entra no estado:

$$
\mathbf{x} = \{q_1, q_2, v_{q_1}, v_{q_2}, Q_1, Q_2\}^\mathsf{T},
\qquad \dot{Q}_i = J_i,
$$

com a dinâmica de $\mathbf{M}$ e $\mathbf{c}$ inalterada. Mantêm-se os limites de velocidade e de torque e acrescentam-se

$$
|J_1| \le 250\ \mathrm{N \cdot m/s}, \qquad |J_2| \le 100\ \mathrm{N \cdot m/s},
$$

além de $\mathbf{Q}(t_i) = \mathbf{Q}(t_f) = \mathbf{0}$. A taxa satura no início, no fim e nas antigas comutações, o que suaviza a velocidade e impõe aceleração nula nos extremos. Sem obstáculo, $t_f = 1{,}106$ s. Com o disco e $N = 3$, $t_f = 1{,}216$ s.

Com os três discos e torque como controle, $t_f = 1{,}362$ s, com quatro arcos em $Q_1$ e cinco em $Q_2$. Com limite de taxa, $t_f = 1{,}491$ s. Nesse exemplo, subir de $N = 3$ para $N = 10$ não altera a trajetória.

### 4.7 Via-points em fases

O caminho entre pontos de passagem continua livre. Exigem-se $V$ configurações intermediárias. O OCP é dividido em $V + 1$ fases, resolvidas em conjunto: da configuração inicial ao primeiro via-point, entre via-points consecutivos, e do último via-point ao destino.

Em cada fase $v \in \{1, \ldots, V\}$,

$$
\mathbf{q}(t_{f,v}) = \mathbf{q}_v.
$$

A continuidade de estado e de tempo entre fases é

$$
\mathbf{x}(t_{f,v}) = \mathbf{x}(t_{i,v+1}), \qquad t_{f,v} = t_{i,v+1}.
$$

Os instantes de comutação e as velocidades nos via-points são variáveis da PNL. O via-point pode ser dado em coordenadas cartesianas e convertido por cinemática inversa.

No exemplo, $\mathbf{q}_{v,1} = \{1{,}0,\ -2{,}0\}^\mathsf{T}$ e $\mathbf{q}_{v,2} = \{1{,}4,\ -1{,}1\}^\mathsf{T}$, com limite de taxa de torque e os três discos. O tempo sobe para $1{,}771$ s, com comutações em $0{,}682$ s e $1{,}177$ s. O torque de $q_1$ conserva o padrão min–max–min–max; o de $q_2$ ganha um par min–max adicional, e os trechos de torque constante mudam de duração.

### 4.8 Manipulador redundante de quatro elos

O mesmo OCP é aplicado a um planar de quatro juntas, obtido dividindo cada elo do modelo base de modo que, com $q_2 = 0$ e $q_4 = 0$, massa, centro de massa e inércia dos pares coincidam com os elos originais. Dois graus de liberdade bastariam para a posição cartesiana do efetuador; os outros dois são redundância.

A energia cinética soma as contribuições dos quatro elos e da carga no efetuador $P_4$, com a velocidade angular do elo $k$ igual à soma $\dot{q}_1 + \cdots + \dot{q}_k$. O torque é $\mathbf{Q} = \{Q_1, Q_2, Q_3, Q_4\}^\mathsf{T}$ e a dinâmica volta à forma $\mathbf{M}(\mathbf{q})\ddot{\mathbf{q}} + \mathbf{c}(\mathbf{q}, \dot{\mathbf{q}}) = \mathbf{Q}$, agora com $\mathbf{M}$ de ordem $4 \times 4$. O artigo não explicita essas matrizes.

O controle é a taxa de torque $\mathbf{J} \in \mathbb{R}^4$, para evitar o bang–bang. Há 12 estados, $\{\mathbf{q}, \dot{\mathbf{q}}, \mathbf{Q}\}$, e 4 controles. As desigualdades são

$$
\begin{aligned}
|\dot{q}_1|, |\dot{q}_2| &\le 3\ \mathrm{rad/s}, &
|\dot{q}_3|, |\dot{q}_4| &\le 8\ \mathrm{rad/s}, \\[2pt]
|Q_1| &\le 25,\ 
|Q_2| \le 15,\ 
|Q_3| \le 9,\ 
|Q_4| \le 5 \quad (\mathrm{N \cdot m}), \\[2pt]
|J_1| &\le 250,\ 
|J_2| \le 150,\ 
|J_3| \le 100,\ 
|J_4| \le 50 \quad (\mathrm{N \cdot m/s}).
\end{aligned}
$$

Contorno: $\mathbf{q}(t_i) = \mathbf{0}$ e $\mathbf{q}(t_f) = \{1,\ 0,\ -0{,}5,\ 0\}^\mathsf{T}$ rad, de modo que os pares de elos comecem e terminem alinhados com os elos do modelo de dois corpos. A colisão é vigiada nos elos 2, 3 e 4, com $N = 2$ pontos por elo, nos três discos da seção 4.6. O custo é o tempo de manobra. A redundância reduz esse tempo para $t_f = 1{,}029$ s, abaixo dos $1{,}491$ s do manipulador de dois elos com os mesmos discos e com limite de taxa.

## 5. Efeito das restrições sobre o tempo ótimo

Tempos obtidos com os dados da seção 3. Servem como leitura da PNL: cada linha muda o conjunto viável ou o custo, e o $t_f$ é o valor ótimo local reportado.

| Formulação | Restrições que definem o caso | $t_f$ (s) |
|------------|-------------------------------|-----------|
| Cinemática, tempo mínimo | velocidade, aceleração 18 rad/s² | 0,500 |
| Cinemática, obstáculo, $N = 1$ | idem, disco só no efetuador | 0,500 |
| Cinemática, obstáculo, $N = 3$ | disco em três pontos do elo 2 | 1,180 |
| Cinemática, jerk, obstáculo, $N = 3$ | jerk limitado, aceleração nula nos extremos | 1,286 |
| Dinâmica, sem obstáculo | torque 25 N·m e 9 N·m | 1,002 |
| Dinâmica, caminho reto imposto | dinâmica mais a reta em $(q_1, q_2)$ | 1,081 |
| Dinâmica, sem carga | $M = 0$ | 0,843 |
| Dinâmica, sem carga, caminho reto | $M = 0$ e reta imposta | 0,921 |
| Dinâmica, obstáculo, $N = 1$ | torque e disco no efetuador | 1,046 |
| Dinâmica, obstáculo, $N = 3$ | torque e três pontos | 1,098 |
| Dinâmica, taxa de torque | $\|J\|$ limitado, torque nulo nos extremos | 1,106 |
| Dinâmica, taxa de torque, $N = 3$ | idem, com o disco | 1,216 |
| Dinâmica, três discos | $C = 3$, $N = 3$, torque como controle | 1,362 |
| Dinâmica, três discos, taxa de torque | $C = 3$ e $\|J\|$ limitado | 1,491 |
| Via-points | duas fases internas, três discos, taxa de torque | 1,771 |
| Quatro elos, redundante | 12 estados, três discos, taxa de torque | 1,029 |

Duas comparações do artigo situam essas soluções. O método PRM, com caminho amostrado e timing por splines 4-3-4 separado da busca do caminho, produziu tempos entre cerca de 1,4 s e 4,0 s (média 2,0 s, mediana 1,9 s, em 1000 resoluções). O melhor caso encontrado, 1,418 s, fica acima dos 1,286 s do OCP cinemático com jerk e o mesmo obstáculo. E o caminho de distância mínima no espaço das juntas, ótimo sob limites cinemáticos, aumenta o tempo assim que a dinâmica de $\mathbf{M}$ e $\mathbf{c}$ é incluída.

## 6. Síntese da modelagem

O OCP de tempo mínimo fica definido por cinco escolhas, todas explícitas nas seções anteriores:

1. **Custo.** Integral do tempo, com ou sem peso $w$ na distância no espaço das juntas. Com variável independente $s$, o custo passa a ser o comprimento do caminho e entra a restrição de tangente unitária.
2. **Dinâmica $\mathbf{f}$.** Cadeia de integradores até a ordem do controle escolhido (aceleração, jerk, torque ou taxa de torque). Quando há dinâmica do robô, a última equação de aceleração é $\mathbf{M}(\mathbf{q})^{-1}(\mathbf{Q} - \mathbf{c})$.
3. **Desigualdades de caminho.** Limites de velocidade e, conforme o caso, de aceleração ou de torque, mais uma desigualdade quadrática-trigonométrica por ponto monitorado e por disco.
4. **Contorno e fases.** Repouso–repouso. Com jerk, aceleração nula nos extremos; com taxa de torque, torque nulo nos extremos. Cada via-point abre uma fase, com posição prescrita e continuidade de estado e de tempo.
5. **PNL.** Colocação LGR em malha $h$–$p$, derivadas por diferenciação automática, variáveis adimensionalizadas, IPOPT sobre a PNL esparsa. A não convexidade vem do horizonte livre, dos obstáculos e, no caso dinâmico, de $\mathbf{M}(\mathbf{q})$ e $\mathbf{c}(\mathbf{q}, \dot{\mathbf{q}})$.
