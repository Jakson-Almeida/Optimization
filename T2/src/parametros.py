"""Dados numéricos do paper (Massaro et al., Robotics 2023, seção 4.2)."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Parametros:
    l1: float = 0.4
    l2: float = 0.25
    q0: tuple = (0.0, 0.0)
    qf: tuple = (1.0, -0.5)
    v_max: tuple = (3.0, 8.0)
    a_max: float = 18.0
    x0_obs: float = 0.45
    y0_obs: float = 0.25
    R_obs: float = 0.1
    n_zeta: int = 3
    K: int = 10
    tf_min: float = 0.2
    tf_max: float = 5.0
    # C-space do paper (seção 4.4), para a figura de convexidade.
    q1_lim: tuple = (0.0, 0.5 * np.pi)
    q2_lim: tuple = (-2.0 * np.pi / 3.0, 0.0)


P = Parametros()

ZETA = tuple((i / P.n_zeta) for i in range(1, P.n_zeta + 1))  # 1/3, 2/3, 1
N_A = 2 * P.K
N_Z = 1 + N_A  # 21
N_EQ = 4  # q1(tf), q2(tf), v1(tf), v2(tf)
N_VEL = 2 * 2 * (P.K + 1)  # duas juntas, ±, em cada nó
N_NOS_OBS = P.K + 1 + P.K  # nós e pontos médios dos intervalos
N_OBS = N_NOS_OBS * P.n_zeta
N_INEQ = N_VEL + N_OBS
