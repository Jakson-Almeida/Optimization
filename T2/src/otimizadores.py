"""Três métodos sobre a mesma PNL: pontos interiores, Hooke–Jeeves, Newton."""

import time

import numpy as np
from scipy.optimize import Bounds, NonlinearConstraint, minimize

from .parametros import N_INEQ, N_Z, P
from .transcricao import (
    construir_chutes,
    desigualdades,
    grad_hess_penalidade,
    igualdades,
    jac_desigualdades,
    jac_igualdades,
    limites,
    objetivo,
    penalidade,
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
        "sucesso": bool(viavel(z)),
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


def pontos_interiores(z0, maxiter=500):
    lb, ub = limites()
    bounds = Bounds(lb, ub)
    cons = [
        NonlinearConstraint(igualdades, 0.0, 0.0, jac=jac_igualdades),
        NonlinearConstraint(
            desigualdades,
            np.full(N_INEQ, -np.inf),
            np.zeros(N_INEQ),
            jac=jac_desigualdades,
        ),
    ]
    t0 = time.perf_counter()
    sol = minimize(
        objetivo,
        np.asarray(z0, dtype=float),
        method="trust-constr",
        jac=lambda z: np.array([1.0] + [0.0] * (N_Z - 1)),
        hess=lambda z: np.zeros((N_Z, N_Z)),
        bounds=bounds,
        constraints=cons,
        options={
            "maxiter": maxiter,
            "verbose": 0,
            "xtol": 1e-8,
            "gtol": 1e-6,
            "barrier_tol": 1e-6,
        },
    )
    cpu = time.perf_counter() - t0
    niter = sol.nit if hasattr(sol, "nit") else 0
    return sol.x, bool(sol.success), niter, cpu, sol.message


def hooke_jeeves(z0, rho=150.0, delta0=0.25, tol=5e-3, max_iter=80):
    """Hooke–Jeeves em R^{21} sobre a penalidade, com projeção nas caixas."""
    z = projetar(z0)
    delta = float(delta0)
    niter = 0
    t0 = time.perf_counter()

    def phi(x):
        return penalidade(x, rho=rho)

    def explorar(x, d):
        melhor = x
        f_mel = phi(x)
        for j in range(N_Z):
            for s in (1.0, -1.0):
                cand = x.copy()
                cand[j] += s * d
                cand = projetar(cand)
                f_c = phi(cand)
                if f_c < f_mel:
                    melhor, f_mel = cand, f_c
        return melhor

    while delta > tol and niter < max_iter:
        niter += 1
        x_exp = explorar(z, delta)
        if phi(x_exp) >= phi(z) - 1e-16:
            delta *= 0.5
            continue
        direcao = x_exp - z
        z = x_exp
        while True:
            padrao = projetar(z + direcao)
            x_exp = explorar(padrao, delta)
            if phi(x_exp) < phi(z) - 1e-16:
                direcao = x_exp - z
                z = x_exp
            else:
                break
    cpu = time.perf_counter() - t0
    return z, True, niter, cpu, f"delta={delta:.3e}"


def newton_penalidade(z0, rho=150.0, maxiter=80, gtol=1e-6):
    """Newton–Levenberg em Φ, com Hessiana analítica da penalidade."""
    z = projetar(z0)
    niter = 0
    msg = "maxiter"
    lam = 1e-2
    t0 = time.perf_counter()
    for niter in range(1, maxiter + 1):
        grad, Hess = grad_hess_penalidade(z, rho=rho)
        gn = float(np.linalg.norm(grad, ord=np.inf))
        if gn < gtol:
            msg = "gtol"
            break
        H_lm = Hess + lam * np.eye(N_Z)
        wmin = float(np.linalg.eigvalsh(H_lm)[0])
        if wmin < 1e-8:
            H_lm = H_lm + (1e-8 - wmin) * np.eye(N_Z)
        try:
            p = -np.linalg.solve(H_lm, grad)
        except np.linalg.LinAlgError:
            p = -grad
        pn = float(np.linalg.norm(p))
        if pn > 4.0:
            p = p * (4.0 / pn)
        phi0 = penalidade(z, rho=rho)
        alpha = 1.0
        aceitou = False
        for _ in range(16):
            z_new = projetar(z + alpha * p)
            passo = z_new - z
            if penalidade(z_new, rho=rho) <= phi0 + 1e-4 * alpha * float(grad @ passo):
                z = z_new
                aceitou = True
                lam = max(lam * 0.4, 1e-8)
                break
            alpha *= 0.5
        if not aceitou:
            lam = min(lam * 6.0, 1e5)
            z = projetar(z - (1e-3 / max(gn, 1e-12)) * grad)
            msg = "backtracking falhou"
    cpu = time.perf_counter() - t0
    ok = msg in {"gtol", "maxiter"}
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
        ok = sorted(registros, key=lambda r: r["violacao_eq"] + r["violacao_ineq"])
        return ok[0]
    return min(ok, key=lambda r: r["tf"])
