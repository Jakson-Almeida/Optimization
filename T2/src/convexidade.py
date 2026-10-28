"""Jacobiana e Hessiana analíticas (SymPy) das peças da PNL."""

import numpy as np
import sympy as sp

from .parametros import P, ZETA


q1, q2, zeta = sp.symbols("q1 q2 zeta", real=True)
tf, h, a_k, v_k, v_kp1 = sp.symbols("t_f h a_k v_k v_{k+1}", real=True)
l1, l2, x0, y0, R = sp.symbols("l1 l2 x0 y0 R", positive=True)

_subs_num = {
    l1: P.l1,
    l2: P.l2,
    x0: P.x0_obs,
    y0: P.y0_obs,
    R: P.R_obs,
}


def mapa_fk():
    theta = q1 + q2
    x = l1 * sp.cos(q1) + zeta * l2 * sp.cos(theta)
    y = l1 * sp.sin(q1) + zeta * l2 * sp.sin(theta)
    return x, y


def folga_simbolica():
    """g(q, ζ) = R² - ||r - r0||². Viável quando g ≤ 0."""
    x, y = mapa_fk()
    g = R ** 2 - (x - x0) ** 2 - (y - y0) ** 2
    return sp.simplify(g)


def jacobiana_folga():
    g = folga_simbolica()
    J = sp.Matrix([[sp.diff(g, q1), sp.diff(g, q2)]])
    return sp.simplify(J)


def hessiana_folga():
    g = folga_simbolica()
    H = sp.hessian(g, [q1, q2])
    return sp.simplify(H)


def residual_bilinear():
    """v_{k+1} - v_k - (tf/K) a_k  (h = tf/K)."""
    K = sp.Integer(P.K)
    return v_kp1 - v_k - (tf / K) * a_k


def hessiana_residual():
    r = residual_bilinear()
    vars_ = [tf, a_k, v_k, v_kp1]
    return sp.hessian(r, vars_), vars_


def objetivo_simbolico():
    return tf


def hessiana_objetivo():
    return sp.hessian(tf, [tf])


def classificar_autovalores(H_num, atol=1e-10):
    w = np.linalg.eigvalsh(np.asarray(H_num, dtype=float))
    pos = np.all(w > atol)
    neg = np.all(w < -atol)
    nulo = np.all(np.abs(w) <= atol)
    indef = (np.any(w > atol) and np.any(w < -atol))
    if pos:
        tipo = "definida positiva (convexidade local estrita)"
    elif neg:
        tipo = "definida negativa (concavidade local estrita)"
    elif nulo:
        tipo = "nula (linear / afim)"
    elif indef:
        tipo = "indefinida (sela, não convexa nem côncava)"
    else:
        tipo = "semidefinida"
    return w, tipo


def numerica_folga(q1_val, q2_val, zeta_val=1.0):
    """g, Jacobiana 1x2 e Hessiana 2x2 em um ponto, a partir do SymPy."""
    g = folga_simbolica()
    J = jacobiana_folga()
    H = hessiana_folga()
    subs = dict(_subs_num)
    subs.update({q1: q1_val, q2: q2_val, zeta: zeta_val})
    g_n = float(g.subs(subs))
    J_n = np.array(J.subs(subs), dtype=float)
    H_n = np.array(H.subs(subs), dtype=float)
    return g_n, J_n, H_n


def grade_cspace(n1=181, n2=181):
    """Maior g (pior violação) no retângulo (q1, q2) do paper. Positivo = colisão."""
    q1g = np.linspace(P.q1_lim[0], P.q1_lim[1], n1)
    q2g = np.linspace(P.q2_lim[0], P.q2_lim[1], n2)
    Q1, Q2 = np.meshgrid(q1g, q2g)
    G = np.full(Q1.shape, -np.inf)
    for zeta in ZETA:
        theta = Q1 + Q2
        x = P.l1 * np.cos(Q1) + zeta * P.l2 * np.cos(theta)
        y = P.l1 * np.sin(Q1) + zeta * P.l2 * np.sin(theta)
        g = P.R_obs ** 2 - (x - P.x0_obs) ** 2 - (y - P.y0_obs) ** 2
        G = np.maximum(G, g)
    return q1g, q2g, G


def latex_pecas():
    """Dicionário de expressões em LaTeX para o notebook."""
    x, y = mapa_fk()
    g = folga_simbolica()
    J = jacobiana_folga()
    H = hessiana_folga()
    r = residual_bilinear()
    Hr, vars_ = hessiana_residual()
    return {
        "fk_x": sp.latex(x),
        "fk_y": sp.latex(y),
        "g": sp.latex(g),
        "J_g": sp.latex(J),
        "H_g": sp.latex(H),
        "residual": sp.latex(r),
        "H_residual": sp.latex(Hr),
        "vars_residual": [sp.latex(v) for v in vars_],
        "H_obj": sp.latex(hessiana_objetivo()),
    }
