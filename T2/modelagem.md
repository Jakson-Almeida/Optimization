# Modelagem — tempo mínimo cinemático de um manipulador planar

Instância cinemática de Massaro et al. (*Robotics* 2023, seção 4.2): dois elos, repouso–repouso, um obstáculo circular, três pontos monitorados no segundo elo. A dinâmica de Lagrange fica de fora. A transcrição usa aceleração constante por intervalo (não a colocação LGR do paper).

---

## 1. Problema contínuo (OCP)

Estado \(\mathbf{x}(t)=\{q_1(t),q_2(t),v_1(t),v_2(t)\}^\mathsf{T}\) e controle \(\mathbf{u}(t)=\{a_1(t),a_2(t)\}^\mathsf{T}\), \(t\in[0,t_f]\):

\[
\dot q_i = v_i,\qquad \dot v_i = a_i,\qquad i=1,2.
\]

Função objetivo (Bolza degenerada):

\[
\min_{t_f,\,\mathbf{u}(\cdot)} \quad J = t_f = \int_0^{t_f} 1\,\mathrm{d}t.
\]

Condições de ponta, limites de caminho e obstáculo:

\[
\begin{aligned}
\mathbf{q}(0)&=\mathbf{q}_0, &
\mathbf{v}(0)&=\mathbf{0}, &
\mathbf{q}(t_f)&=\mathbf{q}_f, &
\mathbf{v}(t_f)&=\mathbf{0}, \\
|v_i(t)|&\le v_i^{\max}, &
|a_i(t)|&\le a^{\max}, \\
\bigl(x_P(\mathbf{q}(t),\zeta)-x_0\bigr)^2
+\bigl(y_P(\mathbf{q}(t),\zeta)-y_0\bigr)^2
&\ge R^2,
\qquad \zeta\in\bigl\{\tfrac13,\tfrac23,1\bigr\}.
\end{aligned}
\]

A cinemática direta de um ponto do elo 2 (eq. 6 do paper) é

\[
x_P(\mathbf{q},\zeta)=l_1\cos q_1+\zeta\,l_2\cos(q_1+q_2),\qquad
y_P(\mathbf{q},\zeta)=l_1\sin q_1+\zeta\,l_2\sin(q_1+q_2).
\]

\(\zeta=0\) é a junta entre os elos; \(\zeta=1\) é o efetuador.

---

## 2. PNL discreta (o que se otimiza)

O horizonte \([0,t_f]\) é partido em \(K\) intervalos iguais de comprimento \(h=t_f/K\). Nós \(k=0,\ldots,K\):

\[
t_k = k\,h.
\]

A aceleração de cada junta é **constante** no intervalo \(k\):

\[
a_i(t)=a_{i,k},\qquad t\in[t_k,t_{k+1}),\quad k=0,\ldots,K-1.
\]

Posição e velocidade **não** são variáveis de decisão: integram-se em forma fechada a partir de \(\mathbf{q}_0\), \(\mathbf{v}_0=\mathbf{0}\) e das acelerações,

\[
\begin{aligned}
v_{i,k+1} &= v_{i,k} + h\,a_{i,k}, \\
q_{i,k+1} &= q_{i,k} + h\,v_{i,k} + \tfrac12 h^2 a_{i,k}.
\end{aligned}
\]

No meio do intervalo \(k\), \(t_{k+1/2}=t_k+h/2\),

\[
q_{i,k+1/2}
= q_{i,k} + \tfrac{h}{2} v_{i,k} + \tfrac12\bigl(\tfrac{h}{2}\bigr)^2 a_{i,k}.
\]

### Variável de decisão

\[
\mathbf{z}
=
\bigl(
t_f,\;
a_{1,0},\ldots,a_{1,K-1},\;
a_{2,0},\ldots,a_{2,K-1}
\bigr)
\in\mathbb{R}^{1+2K}.
\]

Com \(K=50\), \(\mathbf{z}\in\mathbb{R}^{101}\).

### Função objetivo (FOB)

\[
\min_{\mathbf{z}} \quad f(\mathbf{z}) = t_f.
\]

### Sujeito a

**Igualdades** (repouso no alvo), \(h(\mathbf{z})=\mathbf{0}\in\mathbb{R}^4\):

\[
\begin{aligned}
q_{1,K}(\mathbf{z}) - q_{1,f} &= 0, \\
q_{2,K}(\mathbf{z}) - q_{2,f} &= 0, \\
v_{1,K}(\mathbf{z}) &= 0, \\
v_{2,K}(\mathbf{z}) &= 0.
\end{aligned}
\]

**Caixas:**

\[
t_f^{\min} \le t_f \le t_f^{\max},
\qquad
\bigl|a_{i,k}\bigr| \le a^{\max},
\quad i=1,2,\quad k=0,\ldots,K-1.
\]

**Velocidade** nos nós \(k=0,\ldots,K\):

\[
\bigl|v_{i,k}(\mathbf{z})\bigr| \le v_i^{\max},
\quad i=1,2.
\]

**Obstáculo** nos nós \(k=0,\ldots,K\) **e** nos pontos médios \(k+1/2\), para cada \(\zeta\in\{1/3,2/3,1\}\). Equivalente a \(g\le 0\):

\[
g_{\mathrm{obs}}(\mathbf{q},\zeta)
=
R^2 - \bigl(x_P(\mathbf{q},\zeta)-x_0\bigr)^2 - \bigl(y_P(\mathbf{q},\zeta)-y_0\bigr)^2
\le 0.
\]

Contagem: \(4(K+1)\) desigualdades de velocidade e \(3(2K+1)\) de obstáculo, no total \(N_{\mathrm{ineq}}=507\) quando \(K=50\). As condições iniciais \(q_{i,0}=q_{i,0}^{\mathrm{dado}}\) e \(v_{i,0}=0\) estão **fixas na integração**, não entram em \(\mathbf{z}\).

---

## 3. Forma compacta

\[
\begin{aligned}
\min_{\mathbf{z}}\quad
& t_f \\
\mathrm{s.a.}\quad
& q_{i,K}(\mathbf{z}) = q_{i,f},
&& i=1,2, \\
& v_{i,K}(\mathbf{z}) = 0,
&& i=1,2, \\
& \bigl|v_{i,k}(\mathbf{z})\bigr| \le v_i^{\max},
&& i=1,2,\; k=0,\ldots,K, \\
& g_{\mathrm{obs}}\bigl(\mathbf{q}_k(\mathbf{z}),\zeta\bigr) \le 0,
&& k\in\{0,\ldots,K\}\cup\{1/2,\ldots,K-1/2\},\; \zeta\in\{\tfrac13,\tfrac23,1\}, \\
& t_f^{\min} \le t_f \le t_f^{\max}, \\
& \bigl|a_{i,k}\bigr| \le a^{\max},
&& i=1,2,\; k=0,\ldots,K-1.
\end{aligned}
\]

A PNL **não é convexa**: \(g_{\mathrm{obs}}\) composta com a cinemática direta tem Hessiana indefinida no espaço das juntas.

---

## 4. Variáveis e parâmetros

### Decisão (\(\mathbf{z}\))

| Símbolo | Unidade | Papel |
|---------|---------|--------|
| \(t_f\) | s | Tempo total da manobra; é a FOB. |
| \(a_{i,k}\) | rad/s² | Aceleração (constante) da junta \(i\) no intervalo \(k\). |

### Estado discreto (dependente de \(\mathbf{z}\))

| Símbolo | Unidade | Papel |
|---------|---------|--------|
| \(h=t_f/K\) | s | Passo de tempo; muda com \(t_f\). |
| \(t_k=kh\) | s | Instante do nó \(k\). |
| \(q_{i,k}\) | rad | Ângulo da junta \(i\) em \(t_k\). |
| \(v_{i,k}\) | rad/s | Velocidade da junta \(i\) em \(t_k\). |
| \(q_{i,k+1/2}\) | rad | Ângulo no ponto médio do intervalo \(k\) (só para colisão). |

### Dados do robô e da manobra

| Símbolo | Valor | Papel |
|---------|-------|--------|
| \(l_1\) | \(0{,}4\) m | Comprimento do primeiro elo. |
| \(l_2\) | \(0{,}25\) m | Comprimento do segundo elo. |
| \(\mathbf{q}_0\) | \((0,\,0)\) rad | Configuração inicial. |
| \(\mathbf{q}_f\) | \((1{,}0,\,-0{,}5)\) rad | Configuração final. |
| \(\mathbf{v}_0\) | \(\mathbf{0}\) | Partida do repouso (fixo). |
| \(v_1^{\max}\) | \(3\) rad/s | Limite de velocidade da junta 1. |
| \(v_2^{\max}\) | \(8\) rad/s | Limite de velocidade da junta 2. |
| \(a^{\max}\) | \(18\) rad/s² | Limite de aceleração (as duas juntas). |
| \(t_f^{\min},\,t_f^{\max}\) | \(0{,}2\) s, \(8\) s | Caixa de \(t_f\). |
| \(K\) | \(50\) | Número de intervalos iguais. |

### Obstáculo e amostragem do elo 2

| Símbolo | Valor | Papel |
|---------|-------|--------|
| \((x_0,y_0)\) | \((0{,}45,\,0{,}25)\) m | Centro do disco no plano de trabalho. |
| \(R\) | \(0{,}1\) m | Raio do disco. |
| \(\zeta\) | \(1/3,\;2/3,\;1\) | Fração do segundo elo (junta 2 → efetuador). |
| \(x_P,\,y_P\) | — | Coordenadas cartesianas do ponto \(\zeta\) no elo 2. |
| \(g_{\mathrm{obs}}\) | — | Folga: \(\le 0\) fora ou na borda do disco; \(>0\) é colisão. |

### Índices

| Símbolo | Conjunto | Papel |
|---------|----------|--------|
| \(i\) | \(\{1,2\}\) | Junta. |
| \(k\) | \(\{0,\ldots,K-1\}\) para \(a\); \(\{0,\ldots,K\}\) para \(q,v\) | Intervalo / nó. |

---

## 5. O que esta transcrição **não** inclui

- Torque, inércia ou a equação de Lagrange (só cinemática).
- Colocação de Legendre–Gauss–Radau e malha \(h\)–\(p\) do paper: aqui \(a(t)\) é constante em cada \(h\), daí o perfil em degraus.
- Colisão em tempo contínuo: o disco é exigido só nos \(K+1\) nós e nos \(K\) meios. Entre amostras a trajetória pode ainda cortar o obstáculo.
- Elos 1 e a base: só três pontos do elo 2 são monitorados, como \(N=3\) no artigo.
