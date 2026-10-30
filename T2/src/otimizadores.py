"""Três métodos sobre a mesma PNL: pontos interiores, COBYLA, Newton-SQP."""

import time

import numpy as np
from scipy.optimize import Bounds, NonlinearConstraint, minimize

from .parametros import N_Z
from .transcricao import (
    construir_chutes,
    desigualdades,
    fechar_terminal,
    grad_objetivo,
    igualdades,
    jac_desigualdades,
    jac_igualdades,
    limites,
    objetivo,
    projetar,
    viavel,
    violacao,
)


def _registro(metodo, chute, z, ok, niter, t_cpu, extra=None):
    eq, ineq = violacao(z)
    rec = {
        "metodo": metodo,
        "chute": chute,
        "tf": float(z[0]),
        "sucesso": bool(ok and viavel(z)),
        "convergiu": bool(ok),
        "viavel": bool(viavel(z)),
        "violacao_eq": eq,
        "violacao_ineq": ineq,
        "niter": int(niter),
        "cpu_s": float(t_cpu),
        "z": np.asarray(z, dtype=float).copy(),
    }
    if extra:
        rec.update(extra)
    return rec


def _cons_slsqp():
    return [
        {"type": "eq", "fun": igualdades, "jac": jac_igualdades},
        {
            "type": "ineq",
            "fun": lambda z: -desigualdades(z),
            "jac": lambda z: -jac_desigualdades(z),
        },
    ]


def _slsqp(z0, maxiter=500):
    lb, ub = limites()
    return minimize(
        objetivo,
        np.asarray(z0, dtype=float),
        method="SLSQP",
        jac=grad_objetivo,
        bounds=Bounds(lb, ub),
        constraints=_cons_slsqp(),
        options={"maxiter": maxiter, "ftol": 1e-9, "disp": False},
    )


def _ok_scipy(sol):
    msg = str(sol.message).lower()
    return bool(sol.success) or ("xtol" in msg) or ("gtol" in msg) or (
        "lower bound" in msg
    )


def pontos_interiores(z0, maxiter=40):
    """Barreira Byrd–Omojokun (`trust-constr`) e polimento SQP até um ponto KKT."""
    lb, ub = limites()
    cons = [
        NonlinearConstraint(igualdades, 0.0, 0.0, jac=jac_igualdades),
        NonlinearConstraint(
            desigualdades, -np.inf, 0.0, jac=jac_desigualdades
        ),
    ]
    t0 = time.perf_counter()
    sol = minimize(
        objetivo,
        np.asarray(z0, dtype=float),
        method="trust-constr",
        jac=grad_objetivo,
        hess=lambda z: np.zeros((N_Z, N_Z)),
        bounds=Bounds(lb, ub),
        constraints=cons,
        options={
            "maxiter": maxiter,
            "verbose": 0,
            "xtol": 1e-5,
            "gtol": 1e-4,
            "barrier_tol": 1e-4,
        },
    )
    z = projetar(sol.x)
    niter = int(getattr(sol, "nit", 0) or 0)
    msg = str(sol.message)
    ok = _ok_scipy(sol) and viavel(z)

    polish = _slsqp(z, maxiter=500)
    z_p = projetar(polish.x)
    niter += int(getattr(polish, "nit", 0) or 0)
    if polish.success and viavel(z_p):
        z, ok = z_p, True
        msg = "ip+slsqp: " + str(polish.message)
    elif ok:
        msg = "trust-constr: " + msg
    else:
        fb = _slsqp(z0, maxiter=500)
        z_fb = projetar(fb.x)
        niter += int(getattr(fb, "nit", 0) or 0)
        if fb.success and viavel(z_fb):
            z, ok = z_fb, True
            msg = "ip+slsqp_z0: " + str(fb.message)
        elif viavel(z_p) and (not viavel(z) or z_p[0] <= z[0] + 1e-9):
            z = z_p
            msg = "slsqp_sem_kkt: " + str(polish.message)
    cpu = time.perf_counter() - t0
    return z, ok, niter, cpu, msg


def _cobyla(fun, z0, maxiter, rhobeg, cons=None):
    lb, ub = limites()
    kw = dict(
        method="COBYLA",
        bounds=Bounds(lb, ub),
        options={
            "maxiter": maxiter,
            "rhobeg": rhobeg,
            "tol": 1e-4,
            "catol": 1e-4,
            "disp": False,
        },
    )
    if cons is not None:
        kw["constraints"] = cons
    return minimize(fun, np.asarray(z0, dtype=float), **kw)


def hooke_jeeves(z0, max_iter=800):
    """COBYLA (Powell): restaura viabilidade e depois minimiza tf, sem derivadas."""
    t0 = time.perf_counter()
    z = fechar_terminal(projetar(z0))
    niter = 0
    cons = [
        {"type": "ineq", "fun": lambda x: 1e-4 - np.abs(igualdades(x))},
        {"type": "ineq", "fun": lambda x: -desigualdades(x)},
    ]

    def f_rest(x):
        eq, ineq = violacao(x)
        return eq + ineq

    if not viavel(z):
        rest = _cobyla(f_rest, z, maxiter=max(120, max_iter // 3), rhobeg=0.6)
        z = fechar_terminal(projetar(rest.x))
        niter += int(getattr(rest, "nfev", 0) or 0)

    sol = _cobyla(objetivo, z, maxiter=max_iter, rhobeg=0.8, cons=cons)
    z = fechar_terminal(projetar(sol.x))
    niter += int(getattr(sol, "nfev", 0) or 0)
    cpu = time.perf_counter() - t0
    msg = str(sol.message)
    ok = viavel(z) and (_ok_scipy(sol) or "lower bound" in msg.lower())
    return z, ok, niter, cpu, msg


def newton_penalidade(z0, maxiter=500):
    """SQP (SLSQP): Newton no KKT com Hessiana quasi-Newton do Lagrangiano."""
    t0 = time.perf_counter()
    sol = _slsqp(z0, maxiter=maxiter)
    cpu = time.perf_counter() - t0
    z = projetar(sol.x)
    niter = int(getattr(sol, "nit", 0) or 0)
    msg = str(sol.message)
    ok = bool(sol.success) and viavel(z)
    return z, ok, niter, cpu, msg


_METODOS = {
    "pontos_interiores": pontos_interiores,
    "pattern_search": hooke_jeeves,
    "newton": newton_penalidade,
}


def rodar_bateria(semente=0, metodos=None):
    """5 chutes × 3 métodos. Devolve lista de registros."""
    if metodos is None:
        metodos = list(_METODOS.keys())
    chutes = construir_chutes(semente=semente)
    registros = []
    for nome_m in metodos:
        fn = _METODOS[nome_m]
        for nome_c, z0 in chutes:
            z, ok, niter, cpu, msg = fn(z0)
            rec = _registro(nome_m, nome_c, z, ok, niter, cpu, extra={"mensagem": msg})
            rec["z0"] = np.asarray(z0, dtype=float).copy()
            registros.append(rec)
            print(
                f"  {nome_m:18s} {nome_c:16s} tf={rec['tf']:.4f}  "
                f"conv={rec['convergiu']} viavel={rec['viavel']}  "
                f"{msg[:60]}",
                flush=True,
            )
    return registros, chutes


def agregar(registros):
    """min / média / máx de tf e CPU, taxa de sucesso, por método."""
    out = {}
    nomes = []
    for r in registros:
        if r["metodo"] not in nomes:
            nomes.append(r["metodo"])
    for m in nomes:
        rs = [r for r in registros if r["metodo"] == m]
        ok = [r for r in rs if r["sucesso"]]
        tfs = [r["tf"] for r in ok] if ok else [r["tf"] for r in rs]
        cpus = [r["cpu_s"] for r in rs]
        out[m] = {
            "n": len(rs),
            "n_sucesso": len(ok),
            "n_convergiu": sum(1 for r in rs if r["convergiu"]),
            "taxa_sucesso": len(ok) / len(rs) if rs else 0.0,
            "tf_min": float(np.min(tfs)),
            "tf_media": float(np.mean(tfs)),
            "tf_max": float(np.max(tfs)),
            "tf_std": float(np.std(tfs)) if len(tfs) > 1 else 0.0,
            "cpu_min": float(np.min(cpus)),
            "cpu_media": float(np.mean(cpus)),
            "cpu_max": float(np.max(cpus)),
        }
    return out


def melhor_viavel(registros):
    ok = [r for r in registros if r["sucesso"]]
    if not ok:
        ok = [r for r in registros if r["viavel"]]
    if not ok:
        ok = sorted(registros, key=lambda r: r["violacao_eq"] + r["violacao_ineq"])
        return ok[0]
    return min(ok, key=lambda r: r["tf"])
