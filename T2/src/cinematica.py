"""Cinemática direta do elo 2 e folga em relação ao disco."""

import numpy as np

from .parametros import P, ZETA


def ponto_elo2(q1, q2, zeta, l1=P.l1, l2=P.l2):
    """Coordenadas de um ponto do segundo elo (eq. 6 do paper).

    zeta=0 é a junta entre os elos; zeta=1 é o efetuador.
    """
    theta = q1 + q2
    x = l1 * np.cos(q1) + zeta * l2 * np.cos(theta)
    y = l1 * np.sin(q1) + zeta * l2 * np.sin(theta)
    return x, y


def juntas_elo2(q1, q2, l1=P.l1, l2=P.l2):
    """Origem, junta 2 e efetuador, para desenhar o manipulador."""
    x1 = l1 * np.cos(q1)
    y1 = l1 * np.sin(q1)
    x2 = x1 + l2 * np.cos(q1 + q2)
    y2 = y1 + l2 * np.sin(q1 + q2)
    return np.array([0.0, x1, x2]), np.array([0.0, y1, y2])


def distancia2_obstaculo(q1, q2, zeta, x0=P.x0_obs, y0=P.y0_obs):
    x, y = ponto_elo2(q1, q2, zeta)
    return (x - x0) ** 2 + (y - y0) ** 2


def folga_obstaculo(q1, q2, zeta, R=P.R_obs):
    """R² - dist². Negativo ou nulo = fora do disco (viável)."""
    return R ** 2 - distancia2_obstaculo(q1, q2, zeta)


def colide(q1, q2, zetas=ZETA, R=P.R_obs):
    """True se algum dos pontos monitorados está dentro do disco."""
    for zeta in zetas:
        if folga_obstaculo(q1, q2, zeta, R=R) > 0.0:
            return True
    return False


def min_folga(q1, q2, zetas=ZETA, R=P.R_obs):
    """Maior violação (mais positiva) entre os pontos do elo 2."""
    return max(folga_obstaculo(q1, q2, zeta, R=R) for zeta in zetas)
