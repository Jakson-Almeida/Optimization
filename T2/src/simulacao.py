"""Desenho do manipulador, obstáculo, perfis e animação."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation, patches
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from .cinematica import juntas_elo2, ponto_elo2
from .convexidade import grade_cspace
from .parametros import P
from .transcricao import Trajetoria


def _salvar(fig: Figure, caminho):
    FigureCanvasAgg(fig)
    fig.tight_layout()
    fig.savefig(caminho, dpi=140)


def _eixos_workspace(ax):
    raio = P.l1 + P.l2 + 0.08
    ax.set_aspect("equal")
    ax.set_xlim(-0.15, raio)
    ax.set_ylim(-raio * 0.35, raio)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    circ = patches.Circle(
        (P.x0_obs, P.y0_obs),
        P.R_obs,
        facecolor="#c44e52",
        edgecolor="#841c1c",
        alpha=0.85,
        zorder=2,
    )
    ax.add_patch(circ)
    ax.plot(0, 0, "ks", ms=7, zorder=4)
    ax.grid(True, alpha=0.3)


def desenhar_configuracao(ax, q1, q2, cor="#4c72b0", alpha=1.0, lw=3.0):
    xs, ys = juntas_elo2(q1, q2)
    ax.plot(xs, ys, "-", color=cor, lw=lw, alpha=alpha, solid_capstyle="round", zorder=3)
    ax.plot(xs, ys, "o", color=cor, ms=6, alpha=alpha, zorder=4)
    return xs, ys


def figura_estroboscopio(tr: Trajetoria, caminho, n_frames=10):
    fig = Figure(figsize=(6.2, 5.4))
    ax = fig.add_subplot(111)
    _eixos_workspace(ax)
    idx = np.linspace(0, P.K, n_frames).round().astype(int)
    cmap = plt.cm.Blues
    for i, k in enumerate(idx):
        cor = cmap(0.35 + 0.6 * i / max(n_frames - 1, 1))
        desenhar_configuracao(ax, tr.q[0, k], tr.q[1, k], cor=cor, alpha=0.95, lw=2.4)
    xe, ye = ponto_elo2(tr.q[0], tr.q[1], 1.0)
    ax.plot(xe, ye, color="#2a2a2a", lw=1.2, alpha=0.7, label="efetuador")
    ax.set_title(f"Estroboscópio ($t_f$ = {tr.tf:.3f} s)")
    _salvar(fig, caminho)
    return caminho


def figura_perfis(tr: Trajetoria, caminho):
    fig = Figure(figsize=(7.2, 7.8))
    axes = fig.subplots(3, 1, sharex=False)
    t = tr.t
    axes[0].plot(t, tr.q[0], label=r"$q_1$")
    axes[0].plot(t, tr.q[1], label=r"$q_2$")
    axes[0].axhline(P.qf[0], ls=":", color="C0", alpha=0.5)
    axes[0].axhline(P.qf[1], ls=":", color="C1", alpha=0.5)
    axes[0].set_ylabel("posição (rad)")
    axes[0].legend(frameon=False)
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, tr.v[0], label=r"$\dot q_1$")
    axes[1].plot(t, tr.v[1], label=r"$\dot q_2$")
    axes[1].axhline(P.v_max[0], ls="--", color="C0", alpha=0.4)
    axes[1].axhline(-P.v_max[0], ls="--", color="C0", alpha=0.4)
    axes[1].axhline(P.v_max[1], ls="--", color="C1", alpha=0.4)
    axes[1].axhline(-P.v_max[1], ls="--", color="C1", alpha=0.4)
    axes[1].set_ylabel("velocidade (rad/s)")
    axes[1].legend(frameon=False)
    axes[1].grid(True, alpha=0.3)

    axes[2].step(np.concatenate([[0.0], t[1:]]), np.concatenate([tr.a[0], tr.a[0, -1:]]),
                 where="post", label=r"$a_1$")
    axes[2].step(np.concatenate([[0.0], t[1:]]), np.concatenate([tr.a[1], tr.a[1, -1:]]),
                 where="post", label=r"$a_2$")
    axes[2].axhline(P.a_max, ls="--", color="k", alpha=0.3)
    axes[2].axhline(-P.a_max, ls="--", color="k", alpha=0.3)
    axes[2].set_xlabel("t (s)")
    axes[2].set_ylabel("aceleração (rad/s²)")
    axes[2].legend(frameon=False)
    axes[2].grid(True, alpha=0.3)
    fig.suptitle(f"Perfis da trajetória ($t_f$ = {tr.tf:.3f} s)")
    _salvar(fig, caminho)
    return caminho


def figura_cspace(tr: Trajetoria | None, caminho, q1g=None, q2g=None, G=None):
    if G is None:
        q1g, q2g, G = grade_cspace(n1=161, n2=161)
    fig = Figure(figsize=(6.4, 5.6))
    ax = fig.add_subplot(111)
    ax.contourf(q1g, q2g, G, levels=[-1, 0, 1], colors=["#d9f0d3", "#f4c2c2"], alpha=0.9)
    ax.contour(q1g, q2g, G, levels=[0.0], colors=["#841c1c"], linewidths=1.4)
    ax.plot(P.q0[0], P.q0[1], "o", color="#2a2a2a", ms=8, label="inicial")
    ax.plot(P.qf[0], P.qf[1], "s", color="#2a2a2a", ms=8, label="final")
    ax.plot([P.q0[0], P.qf[0]], [P.q0[1], P.qf[1]], "k--", lw=1.0, alpha=0.6, label="reta (inviável)")
    if tr is not None:
        _, qd, _ = interpolar_trajetoria(tr, n=200)
        ax.plot(qd[0], qd[1], color="#4c72b0", lw=2.0, label="trajetória")
    ax.set_xlabel(r"$q_1$ (rad)")
    ax.set_ylabel(r"$q_2$ (rad)")
    ax.set_title("Espaço das juntas: verde viável, vermelho colisão")
    ax.legend(frameon=False, loc="lower left")
    ax.set_xlim(P.q1_lim)
    ax.set_ylim(P.q2_lim)
    ax.grid(True, alpha=0.25)
    _salvar(fig, caminho)
    return caminho


def interpolar_trajetoria(tr: Trajetoria, n=80):
    """Amostras densas no tempo, aceleração constante por intervalo."""
    t_denso = np.linspace(0.0, tr.tf, n)
    q = np.zeros((2, n))
    v = np.zeros((2, n))
    h = tr.tf / P.K
    for i, t in enumerate(t_denso):
        if t >= tr.tf:
            q[:, i] = tr.q[:, -1]
            v[:, i] = tr.v[:, -1]
            continue
        k = min(int(np.floor(t / h + 1e-12)), P.K - 1)
        tau = t - k * h
        v[:, i] = tr.v[:, k] + tau * tr.a[:, k]
        q[:, i] = tr.q[:, k] + tau * tr.v[:, k] + 0.5 * tau * tau * tr.a[:, k]
    return t_denso, q, v


def animar(tr: Trajetoria, caminho, n_frames=80, fps=20):
    t, q, _ = interpolar_trajetoria(tr, n=n_frames)
    xe, ye = ponto_elo2(q[0], q[1], 1.0)
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    _eixos_workspace(ax)
    (linha,) = ax.plot([], [], "-", color="#4c72b0", lw=3.2, zorder=3)
    (juntas,) = ax.plot([], [], "o", color="#4c72b0", ms=7, zorder=4)
    (rastro,) = ax.plot([], [], color="#2a2a2a", lw=1.1, alpha=0.65)
    titulo = ax.set_title("")

    def init():
        linha.set_data([], [])
        juntas.set_data([], [])
        rastro.set_data([], [])
        return linha, juntas, rastro, titulo

    def update(i):
        xs, ys = juntas_elo2(q[0, i], q[1, i])
        linha.set_data(xs, ys)
        juntas.set_data(xs, ys)
        rastro.set_data(xe[: i + 1], ye[: i + 1])
        titulo.set_text(f"t = {t[i]:.3f} s / {tr.tf:.3f} s")
        return linha, juntas, rastro, titulo

    ani = animation.FuncAnimation(
        fig, update, frames=n_frames, init_func=init, blit=True, interval=1000 / fps
    )
    caminho = Path(caminho)
    caminho.parent.mkdir(parents=True, exist_ok=True)
    ani.save(str(caminho), writer=animation.PillowWriter(fps=fps))
    plt.close(fig)
    return str(caminho)


def salvar_tudo(tr: Trajetoria, pasta: Path, prefixo="melhor"):
    pasta = Path(pasta)
    pasta.mkdir(parents=True, exist_ok=True)
    out = {
        "estroboscopio": str(figura_estroboscopio(tr, pasta / f"{prefixo}_estroboscopio.png")),
        "perfis": str(figura_perfis(tr, pasta / f"{prefixo}_perfis.png")),
        "cspace": str(figura_cspace(tr, pasta / f"{prefixo}_cspace.png")),
        "gif": animar(tr, pasta / f"{prefixo}.gif"),
    }
    return out
