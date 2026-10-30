"""Transcrição da PNL: decisão z, integração, restrições, Jacobiana e Hessiana.

z = (tf, a1_0..a1_{K-1}, a2_0..a2_{K-1}) in R^{1+2K}

A aceleração é constante em cada intervalo. A cadeia q̈ = a integra-se em
forma fechada. q e v não são variáveis de decisão.
"""

from dataclasses import dataclass

import numpy as np

from .cinematica import min_folga, ponto_elo2
from .parametros import N_A, N_EQ, N_INEQ, N_OBS, N_VEL, N_Z, P, ZETA


def n_decisao():
    return N_Z


def unpack(z):
    z = np.asarray(z, dtype=float).reshape(-1)
    if z.size != N_Z:
        raise ValueError(
            f"z tem {z.size} componentes, esperado {N_Z} (K={P.K}). "
            "A bateria foi gravada com outro K; apague saida/bateria.pkl e rode de novo."
        )
    tf = float(z[0])
    a = np.empty((2, P.K), dtype=float)
    a[0] = z[1 : 1 + P.K]
    a[1] = z[1 + P.K : 1 + 2 * P.K]
    return tf, a


def pack(tf, a):
    a = np.asarray(a, dtype=float)
    return np.concatenate(([float(tf)], a[0], a[1]))


def limites():
    lb = np.empty(N_Z)
    ub = np.empty(N_Z)
    lb[0], ub[0] = P.tf_min, P.tf_max
    lb[1:], ub[1:] = -P.a_max, P.a_max
    return lb, ub


def projetar(z):
    lb, ub = limites()
    return np.clip(np.asarray(z, dtype=float), lb, ub)


@dataclass
class Trajetoria:
    tf: float
    t: np.ndarray
    q: np.ndarray  # (2, K+1)
    v: np.ndarray  # (2, K+1)
    a: np.ndarray  # (2, K)


def integrar(z) -> Trajetoria:
    tf, a = unpack(z)
    h = tf / P.K
    q = np.zeros((2, P.K + 1))
    v = np.zeros((2, P.K + 1))
    q[:, 0] = P.q0
    for k in range(P.K):
        v[:, k + 1] = v[:, k] + h * a[:, k]
        q[:, k + 1] = q[:, k] + h * v[:, k] + 0.5 * h * h * a[:, k]
    t = np.linspace(0.0, tf, P.K + 1)
    return Trajetoria(tf=tf, t=t, q=q, v=v, a=a)


def objetivo(z):
    return float(np.asarray(z, dtype=float).reshape(-1)[0])


def grad_objetivo(z):
    g = np.zeros(N_Z)
    g[0] = 1.0
    return g


def igualdades(z):
    tr = integrar(z)
    return np.array(
        [
            tr.q[0, -1] - P.qf[0],
            tr.q[1, -1] - P.qf[1],
            tr.v[0, -1],
            tr.v[1, -1],
        ]
    )


def _amostras_q(tr: Trajetoria):
    """q nos nós e nos pontos médios dos intervalos (colisão)."""
    qs = [tr.q[:, k].copy() for k in range(P.K + 1)]
    dt = 0.5 * (tr.tf / P.K)
    for k in range(P.K):
        qs.append(tr.q[:, k] + dt * tr.v[:, k] + 0.5 * dt * dt * tr.a[:, k])
    return qs


def desigualdades(z):
    """g(z) <= 0. Velocidade em ± e obstáculo (R² - dist²) em nós, meios e ζ."""
    tr = integrar(z)
    g = np.empty(N_INEQ)
    i = 0
    v_max = np.array(P.v_max)
    for k in range(P.K + 1):
        for j in range(2):
            g[i] = tr.v[j, k] - v_max[j]
            g[i + 1] = -tr.v[j, k] - v_max[j]
            i += 2
    assert i == N_VEL
    for qk in _amostras_q(tr):
        for zeta in ZETA:
            x, y = ponto_elo2(qk[0], qk[1], zeta)
            g[i] = P.R_obs ** 2 - (x - P.x0_obs) ** 2 - (y - P.y0_obs) ** 2
            i += 1
    assert i == N_INEQ
    return g


def violacao(z):
    h = igualdades(z)
    g = desigualdades(z)
    return float(np.max(np.abs(h))), float(np.max(np.maximum(g, 0.0)))


def viavel(z, eq_tol=1e-4, ineq_tol=1e-4):
    eq, ineq = violacao(z)
    return eq <= eq_tol and ineq <= ineq_tol


def avaliar(z):
    tr = integrar(z)
    eq = igualdades(z)
    g = desigualdades(z)
    eq_max, ineq_max = violacao(z)
    return {
        "tf": tr.tf,
        "trajetoria": tr,
        "igualdades": eq,
        "desigualdades": g,
        "violacao_eq": eq_max,
        "violacao_ineq": ineq_max,
        "viavel": eq_max <= 1e-4 and ineq_max <= 1e-4,
        "colide": bool(np.any(g[N_VEL:] > 1e-12)),
    }


# ---------------------------------------------------------------------------
# Sensibilidades da cadeia de integradores (analíticas)
# ---------------------------------------------------------------------------

@dataclass
class Sensibilidades:
    dq_dtf: np.ndarray  # (2, K+1)
    dv_dtf: np.ndarray
    dq_da: np.ndarray  # (2, K+1, K)  — só a junta correspondente
    dv_da: np.ndarray
    d2q_dtf2: np.ndarray
    d2q_dtf_da: np.ndarray
    d2v_dtf_da: np.ndarray


def sensibilidades(z) -> tuple[Trajetoria, Sensibilidades]:
    """Derivadas primeiras e segundas de q, v em relação a tf e a_{i,j}."""
    tf, a = unpack(z)
    h = tf / P.K
    invK = 1.0 / P.K
    tr = integrar(z)

    dq_dtf = np.zeros((2, P.K + 1))
    dv_dtf = np.zeros((2, P.K + 1))
    dq_da = np.zeros((2, P.K + 1, P.K))
    dv_da = np.zeros((2, P.K + 1, P.K))
    d2q_dtf2 = np.zeros((2, P.K + 1))
    d2q_dtf_da = np.zeros((2, P.K + 1, P.K))
    d2v_dtf_da = np.zeros((2, P.K + 1, P.K))

    for k in range(P.K):
        # primeiras, tf
        dv_dtf[:, k + 1] = dv_dtf[:, k] + invK * a[:, k]
        dq_dtf[:, k + 1] = (
            dq_dtf[:, k]
            + invK * tr.v[:, k]
            + h * dv_dtf[:, k]
            + (h * invK) * a[:, k]
        )
        # primeiras, a_j da mesma junta
        dv_da[:, k + 1, :] = dv_da[:, k, :]
        dv_da[:, k + 1, k] += h
        dq_da[:, k + 1, :] = dq_da[:, k, :] + h * dv_da[:, k, :]
        dq_da[:, k + 1, k] += 0.5 * h * h
        # segundas
        d2q_dtf2[:, k + 1] = (
            d2q_dtf2[:, k] + (2.0 * invK) * dv_dtf[:, k] + (invK * invK) * a[:, k]
        )
        d2v_dtf_da[:, k + 1, :] = d2v_dtf_da[:, k, :]
        d2v_dtf_da[:, k + 1, k] += invK
        d2q_dtf_da[:, k + 1, :] = (
            d2q_dtf_da[:, k, :]
            + invK * dv_da[:, k, :]
            + h * d2v_dtf_da[:, k, :]
        )
        d2q_dtf_da[:, k + 1, k] += h * invK

    sens = Sensibilidades(
        dq_dtf=dq_dtf,
        dv_dtf=dv_dtf,
        dq_da=dq_da,
        dv_da=dv_da,
        d2q_dtf2=d2q_dtf2,
        d2q_dtf_da=d2q_dtf_da,
        d2v_dtf_da=d2v_dtf_da,
    )
    return tr, sens


def _coluna_q_meio(sens: Sensibilidades, tr: Trajetoria, junta: int, k: int):
    """Gradiente de q_junta no ponto médio do intervalo k."""
    h = tr.tf / P.K
    dt = 0.5 * h
    g = np.zeros(N_Z)
    g[0] = (
        sens.dq_dtf[junta, k]
        + (0.5 / P.K) * tr.v[junta, k]
        + dt * sens.dv_dtf[junta, k]
        + dt * (0.5 / P.K) * tr.a[junta, k]
    )
    sl = slice(1 + junta * P.K, 1 + (junta + 1) * P.K)
    g[sl] = sens.dq_da[junta, k, :] + dt * sens.dv_da[junta, k, :]
    g[1 + junta * P.K + k] += 0.5 * dt * dt
    return g


def _coluna_q_amostra(sens, tr, junta, idx):
    if idx <= P.K:
        return _coluna_junta(sens, junta, idx, "q")
    return _coluna_q_meio(sens, tr, junta, idx - (P.K + 1))


def _coluna_junta(sens: Sensibilidades, junta: int, knot: int, qual: str):
    """Gradiente  (N_Z,)  de q_junta[knot] ou v_junta[knot]."""
    g = np.zeros(N_Z)
    if qual == "q":
        g[0] = sens.dq_dtf[junta, knot]
        sl = slice(1 + junta * P.K, 1 + (junta + 1) * P.K)
        g[sl] = sens.dq_da[junta, knot, :]
    else:
        g[0] = sens.dv_dtf[junta, knot]
        sl = slice(1 + junta * P.K, 1 + (junta + 1) * P.K)
        g[sl] = sens.dv_da[junta, knot, :]
    return g


def _hessiana_q(sens: Sensibilidades, junta: int, knot: int):
    """Hessiana (N_Z, N_Z) de q_junta[knot]. Só blocos tf-tf e tf-a_junta."""
    H = np.zeros((N_Z, N_Z))
    H[0, 0] = sens.d2q_dtf2[junta, knot]
    sl = slice(1 + junta * P.K, 1 + (junta + 1) * P.K)
    mixed = sens.d2q_dtf_da[junta, knot, :]
    H[0, sl] = mixed
    H[sl, 0] = mixed
    return H


def _hessiana_v(sens: Sensibilidades, junta: int, knot: int):
    H = np.zeros((N_Z, N_Z))
    sl = slice(1 + junta * P.K, 1 + (junta + 1) * P.K)
    mixed = sens.d2v_dtf_da[junta, knot, :]
    H[0, sl] = mixed
    H[sl, 0] = mixed
    return H


def jac_igualdades(z):
    _, sens = sensibilidades(z)
    J = np.zeros((N_EQ, N_Z))
    J[0] = _coluna_junta(sens, 0, P.K, "q")
    J[1] = _coluna_junta(sens, 1, P.K, "q")
    J[2] = _coluna_junta(sens, 0, P.K, "v")
    J[3] = _coluna_junta(sens, 1, P.K, "v")
    return J


def _fk_derivadas(q1, q2, zeta):
    """g = R² - (x-x0)² - (y-y0)² e derivadas em (q1, q2)."""
    theta = q1 + q2
    c1, s1 = np.cos(q1), np.sin(q1)
    ct, st = np.cos(theta), np.sin(theta)
    x = P.l1 * c1 + zeta * P.l2 * ct
    y = P.l1 * s1 + zeta * P.l2 * st
    dx = x - P.x0_obs
    dy = y - P.y0_obs
    x_q1 = -P.l1 * s1 - zeta * P.l2 * st
    x_q2 = -zeta * P.l2 * st
    y_q1 = P.l1 * c1 + zeta * P.l2 * ct
    y_q2 = zeta * P.l2 * ct
    g = P.R_obs ** 2 - dx * dx - dy * dy
    g_q1 = -2.0 * dx * x_q1 - 2.0 * dy * y_q1
    g_q2 = -2.0 * dx * x_q2 - 2.0 * dy * y_q2
    x_11 = -P.l1 * c1 - zeta * P.l2 * ct
    x_12 = -zeta * P.l2 * ct
    x_22 = -zeta * P.l2 * ct
    y_11 = -P.l1 * s1 - zeta * P.l2 * st
    y_12 = -zeta * P.l2 * st
    y_22 = -zeta * P.l2 * st
    g_11 = -2.0 * (x_q1 * x_q1 + dx * x_11 + y_q1 * y_q1 + dy * y_11)
    g_12 = -2.0 * (x_q1 * x_q2 + dx * x_12 + y_q1 * y_q2 + dy * y_12)
    g_22 = -2.0 * (x_q2 * x_q2 + dx * x_22 + y_q2 * y_q2 + dy * y_22)
    Hg = np.array([[g_11, g_12], [g_12, g_22]])
    return g, np.array([g_q1, g_q2]), Hg


def jac_desigualdades(z):
    tr, sens = sensibilidades(z)
    J = np.zeros((N_INEQ, N_Z))
    i = 0
    for k in range(P.K + 1):
        for j in range(2):
            col = _coluna_junta(sens, j, k, "v")
            J[i] = col
            J[i + 1] = -col
            i += 2
    qs = _amostras_q(tr)
    for idx, qk in enumerate(qs):
        Jq1 = _coluna_q_amostra(sens, tr, 0, idx)
        Jq2 = _coluna_q_amostra(sens, tr, 1, idx)
        for zeta in ZETA:
            _, dg_dq, _ = _fk_derivadas(qk[0], qk[1], zeta)
            J[i] = dg_dq[0] * Jq1 + dg_dq[1] * Jq2
            i += 1
    return J


def _hessiana_obstaculo(tr, sens, idx, zeta, qs=None):
    if qs is None:
        qs = _amostras_q(tr)
    qk = qs[idx]
    _, dg_dq, Hg_q = _fk_derivadas(qk[0], qk[1], zeta)
    Jq1 = _coluna_q_amostra(sens, tr, 0, idx)
    Jq2 = _coluna_q_amostra(sens, tr, 1, idx)
    Jq = np.vstack([Jq1, Jq2])
    H = Jq.T @ Hg_q @ Jq
    if idx <= P.K:
        H += dg_dq[0] * _hessiana_q(sens, 0, idx)
        H += dg_dq[1] * _hessiana_q(sens, 1, idx)
    return H


def penalidade(z, rho=20.0):
    z = np.asarray(z, dtype=float).reshape(-1)
    h = igualdades(z)
    g = desigualdades(z)
    s = np.maximum(g, 0.0)
    return float(z[0] + rho * (h @ h + s @ s))


def grad_hess_penalidade(z, rho=20.0):
    """Gradiente e Hessiana analíticos de Φ = tf + ρ (||h||² + ||max(g,0)||²)."""
    z = np.asarray(z, dtype=float).reshape(-1)
    tr, sens = sensibilidades(z)
    h = igualdades(z)
    g = desigualdades(z)
    Jh = jac_igualdades(z)
    Jg = jac_desigualdades(z)

    grad = grad_objetivo(z) + 2.0 * rho * (Jh.T @ h)
    Hess = 2.0 * rho * (Jh.T @ Jh)

    H_eq = [
        _hessiana_q(sens, 0, P.K),
        _hessiana_q(sens, 1, P.K),
        _hessiana_v(sens, 0, P.K),
        _hessiana_v(sens, 1, P.K),
    ]
    for hi, Hi in zip(h, H_eq):
        Hess += 2.0 * rho * hi * Hi

    ativo = g > 0.0
    if np.any(ativo):
        ga = g[ativo]
        Ja = Jg[ativo]
        grad += 2.0 * rho * (Ja.T @ ga)
        Hess += 2.0 * rho * (Ja.T @ Ja)
        idx_ativos = np.flatnonzero(ativo)
        for ia, gi in zip(idx_ativos, ga):
            if ia < N_VEL:
                knot = ia // 4
                resto = ia % 4
                junta = resto // 2
                sinal = 1.0 if (resto % 2 == 0) else -1.0
                Hess += 2.0 * rho * gi * sinal * _hessiana_v(sens, junta, knot)
            else:
                j = ia - N_VEL
                knot = j // P.n_zeta
                zeta = ZETA[j % P.n_zeta]
                Hess += 2.0 * rho * gi * _hessiana_obstaculo(tr, sens, knot, zeta)
    return grad, Hess


def hess_igualdades(z, v):
    """Σ v_i ∇² h_i, para o Lagrangian do trust-constr."""
    _, sens = sensibilidades(z)
    H = np.zeros((N_Z, N_Z))
    blocos = [
        _hessiana_q(sens, 0, P.K),
        _hessiana_q(sens, 1, P.K),
        _hessiana_v(sens, 0, P.K),
        _hessiana_v(sens, 1, P.K),
    ]
    for vi, Hi in zip(np.asarray(v, dtype=float).reshape(-1), blocos):
        if abs(vi) > 1e-16:
            H += float(vi) * Hi
    return H


def hess_desigualdades(z, v):
    """Σ v_i ∇² g_i, só nas componentes com multiplicador não nulo."""
    tr, sens = sensibilidades(z)
    H = np.zeros((N_Z, N_Z))
    qs = None
    for ia, vi in enumerate(np.asarray(v, dtype=float).reshape(-1)):
        if abs(vi) < 1e-12:
            continue
        if ia < N_VEL:
            knot = ia // 4
            resto = ia % 4
            junta = resto // 2
            sinal = 1.0 if (resto % 2 == 0) else -1.0
            H += float(vi) * sinal * _hessiana_v(sens, junta, knot)
        else:
            if qs is None:
                qs = _amostras_q(tr)
            j = ia - N_VEL
            idx = j // P.n_zeta
            zeta = ZETA[j % P.n_zeta]
            H += float(vi) * _hessiana_obstaculo(tr, sens, idx, zeta, qs=qs)
    return H


def merita_l1(z, rho=20.0):
    h = igualdades(z)
    g = desigualdades(z)
    return float(z[0] + rho * (np.sum(np.abs(h)) + np.sum(np.maximum(g, 0.0))))


def checar_jacobiana(z, eps=1e-6):
    """Diferenças finitas vs. Jacobiana analítica. Devolve erros relativos."""
    z = np.asarray(z, dtype=float).reshape(-1)
    Jh = jac_igualdades(z)
    Jg = jac_desigualdades(z)
    Jh_fd = np.zeros_like(Jh)
    Jg_fd = np.zeros_like(Jg)
    h0 = igualdades(z)
    g0 = desigualdades(z)
    for j in range(N_Z):
        zp = z.copy()
        zp[j] += eps
        Jh_fd[:, j] = (igualdades(zp) - h0) / eps
        Jg_fd[:, j] = (desigualdades(zp) - g0) / eps

    def rel(a, b):
        den = max(np.linalg.norm(b), 1e-12)
        return float(np.linalg.norm(a - b) / den)

    return rel(Jh, Jh_fd), rel(Jg, Jg_fd)


# ---------------------------------------------------------------------------
# Chutes iniciais
# ---------------------------------------------------------------------------

def _acel_rest_discreto(dq, n_int, h, amax):
    """Bang-coast-bang discreto com v_final=0 e deslocamento dq (se |α|≤amax)."""
    a = np.zeros(n_int)
    if abs(dq) < 1e-14 or n_int <= 0 or h <= 0.0:
        return a
    w = n_int - np.arange(n_int) - 0.5
    melhor = None
    n_acc_0 = max(1, n_int // 4)
    n_acc_1 = max(n_acc_0, n_int // 2)
    for n_acc in range(n_acc_0, n_acc_1 + 1):
        n_dec = n_acc
        n_coast = n_int - n_acc - n_dec
        if n_coast < 0:
            continue
        coef = w[:n_acc].sum() - w[n_int - n_dec :].sum()
        if abs(coef) < 1e-14:
            continue
        alpha = dq / (h * h * coef)
        cand = np.zeros(n_int)
        cand[:n_acc] = alpha
        cand[n_int - n_dec :] = -alpha
        melhor = cand
        if np.max(np.abs(cand)) <= amax + 1e-12:
            return cand
    if melhor is None:
        return a
    return np.clip(melhor, -amax, amax)


def fechar_terminal(z):
    """Ajusta a_{K-2}, a_{K-1} de cada junta para q(tf)=qf e v(tf)=0 (forma fechada)."""
    tf, a = unpack(z)
    h = tf / P.K
    tr = integrar(z)
    for j in range(2):
        qkm2 = tr.q[j, P.K - 2]
        vkm2 = tr.v[j, P.K - 2]
        s = -vkm2 / h
        a[j, -2] = (P.qf[j] - qkm2 - 2.0 * h * vkm2) / (h * h) - 0.5 * s
        a[j, -1] = s - a[j, -2]
    return projetar(pack(tf, a))


def _acel_from_q(q, h):
    """Reproduz q[0..K-2] com v_0=0; os dois últimos a fecham q_K e v_K=0."""
    a = np.zeros(P.K)
    v = 0.0
    q_now = q[0]
    for k in range(P.K - 2):
        a[k] = 2.0 * (q[k + 1] - q_now - h * v) / (h * h)
        q_now = q_now + h * v + 0.5 * h * h * a[k]
        v = v + h * a[k]
    # q_now = q_{K-2}, v = v_{K-2}
    S = -v / h
    a[-2] = (q[-1] - q_now - 2.0 * h * v) / (h * h) - 0.5 * S
    a[-1] = S - a[-2]
    return a


def caminho_cspace_livre(n_grid=90):
    """Polilinha em (q1, q2) de q0 a qf que contorna o disco no C-space."""
    from collections import deque

    q1s = np.linspace(P.q1_lim[0], min(P.q1_lim[1], 1.45), n_grid)
    q2s = np.linspace(P.q2_lim[0], P.q2_lim[1], n_grid)
    livre = np.empty((n_grid, n_grid), dtype=bool)
    for i, qq2 in enumerate(q2s):
        for j, qq1 in enumerate(q1s):
            livre[i, j] = min_folga(qq1, qq2) <= -8e-3

    def idx(q1, q2):
        j = int(np.clip(np.argmin(np.abs(q1s - q1)), 0, n_grid - 1))
        i = int(np.clip(np.argmin(np.abs(q2s - q2)), 0, n_grid - 1))
        return i, j

    s = idx(*P.q0)
    g = idx(*P.qf)
    vis = np.zeros_like(livre)
    pai = {}
    dq = deque([s])
    vis[s] = True
    achou = False
    while dq:
        cur = dq.popleft()
        if cur == g:
            achou = True
            break
        i, j = cur
        for di, dj in (
            (1, 0), (-1, 0), (0, 1), (0, -1),
            (1, 1), (1, -1), (-1, 1), (-1, -1),
        ):
            ni, nj = i + di, j + dj
            if 0 <= ni < n_grid and 0 <= nj < n_grid and livre[ni, nj] and not vis[ni, nj]:
                vis[ni, nj] = True
                pai[(ni, nj)] = cur
                dq.append((ni, nj))
    if not achou:
        raise RuntimeError("Não há caminho livre no C-space com a malha usada.")
    nos = [g]
    cur = g
    while cur != s:
        cur = pai[cur]
        nos.append(cur)
    nos.reverse()
    pts = np.array([[q1s[j], q2s[i]] for i, j in nos], dtype=float)
    pts[0] = P.q0
    pts[-1] = P.qf
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s_cum = np.concatenate([[0.0], np.cumsum(seg)])
    s_cum /= max(s_cum[-1], 1e-12)
    s_nos = np.linspace(0.0, 1.0, P.K + 1)
    q = np.column_stack([
        np.interp(s_nos, s_cum, pts[:, 0]),
        np.interp(s_nos, s_cum, pts[:, 1]),
    ]).T
    return q


def chute_reta(tf=0.55):
    """Perfil rest–repouso direto no espaço das juntas (atravessa o disco)."""
    h = tf / P.K
    a1 = _acel_rest_discreto(P.qf[0] - P.q0[0], P.K, h, P.a_max)
    a2 = _acel_rest_discreto(P.qf[1] - P.q0[1], P.K, h, P.a_max)
    return pack(tf, np.vstack([a1, a2]))


def _suavizar_q(q, npassos=8):
    """Média móvel nas juntas, extremos fixos em q0 e qf."""
    q = np.asarray(q, dtype=float).copy()
    for _ in range(npassos):
        q[:, 1:-1] = 0.25 * q[:, :-2] + 0.5 * q[:, 1:-1] + 0.25 * q[:, 2:]
    q[:, 0] = P.q0
    q[:, -1] = P.qf
    return q


def chute_desvio(tf=None, via=None):
    """Segue um caminho livre no C-space (contorna o disco pela esquerda)."""
    del via
    q = _suavizar_q(caminho_cspace_livre())
    tf0 = 2.5 if tf is None else float(tf)

    def _a(tf_):
        h = float(tf_) / P.K
        return np.vstack([_acel_from_q(q[0], h), _acel_from_q(q[1], h)])

    a_ref = _a(tf0)
    if tf is None:
        need = float(np.max(np.abs(a_ref)))
        if need > P.a_max:
            tf0 = tf0 * np.sqrt(need / P.a_max) * 1.06
        tf0 = float(np.clip(tf0, P.tf_min, P.tf_max))
        a_ref = _a(tf0)
    z = pack(tf0, a_ref)
    if float(np.max(np.abs(a_ref))) > P.a_max + 1e-9:
        z = fechar_terminal(z)
    return z


def construir_chutes(semente=0):
    """Cinco partidas compartilhadas pelos três métodos."""
    rng = np.random.default_rng(semente)
    z_reta = chute_reta()
    z_desvio = chute_desvio()
    chutes = [z_reta, z_desvio]
    nomes = ["reta_juntas", "desvio_esquerda"]
    for i in range(3):
        z = z_desvio.copy()
        for _ in range(40):
            ruido = rng.normal(0.0, 0.12, size=N_A)
            cand = z_desvio.copy()
            cand[1:] = np.clip(cand[1:] + ruido, -P.a_max, P.a_max)
            cand[0] = float(np.clip(cand[0] + rng.normal(0.0, 0.05), P.tf_min, P.tf_max))
            cand = fechar_terminal(cand)
            if viavel(cand):
                z = cand
                break
        chutes.append(z)
        nomes.append(f"desvio_ruido_{i + 1}")
    return list(zip(nomes, chutes))
