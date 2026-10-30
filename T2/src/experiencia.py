"""Roda a bateria, grava pickle/JSON e gera as figuras da melhor trajetória."""

import json
import pickle

import matplotlib

matplotlib.use("Agg")
import numpy as np

from .otimizadores import agregar, melhor_viavel, rodar_bateria
from .parametros import N_Z, P
from .paths import figuras_dir, saida_dir  # Path via pathlib inside paths.py
from .simulacao import salvar_tudo
from .transcricao import checar_jacobiana, construir_chutes, integrar


def _json_seguro(obj):
    if isinstance(obj, dict):
        return {k: _json_seguro(v) for k, v in obj.items() if k not in {"z", "z0"}}
    if isinstance(obj, (list, tuple)):
        return [_json_seguro(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        return float(obj)
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def executar(semente=0):
    chutes = construir_chutes(semente=semente)
    err_h, err_g = checar_jacobiana(chutes[1][1])
    registros, _ = rodar_bateria(semente=semente)
    resumo = agregar(registros)
    melhor = melhor_viavel(registros)
    tr = integrar(melhor["z"])
    figs = salvar_tudo(tr, figuras_dir(), prefixo="melhor")

    saida = saida_dir()
    payload = {
        "K": int(P.K),
        "n_z": int(N_Z),
        "registros": registros,
        "resumo": resumo,
        "melhor": melhor,
        "jac_err_eq": err_h,
        "jac_err_ineq": err_g,
        "figuras": figs,
    }
    with open(saida / "bateria.pkl", "wb") as f:
        pickle.dump(payload, f)
    with open(saida / "resumo.json", "w", encoding="utf-8") as f:
        json.dump(
            {
                "K": int(P.K),
                "n_z": int(N_Z),
                "resumo": _json_seguro(resumo),
                "melhor": _json_seguro(melhor),
                "jac_err_eq": err_h,
                "jac_err_ineq": err_g,
                "figuras": figs,
                "corridas": _json_seguro(registros),
            },
            f,
            indent=2,
            ensure_ascii=False,
        )
    return payload


def bateria_compativel(payload):
    """True se o pickle foi gravado com o K atual."""
    if not isinstance(payload, dict) or "melhor" not in payload:
        return False
    if int(payload.get("K", -1)) != int(P.K):
        return False
    z = np.asarray(payload["melhor"].get("z", []), dtype=float).reshape(-1)
    return z.size == N_Z


def carregar_ou_executar(semente=0, forcar=False):
    """Lê saida/bateria.pkl se for do K atual; senão rerroda a bateria."""
    pkl = saida_dir() / "bateria.pkl"
    if not forcar and pkl.exists():
        with open(pkl, "rb") as f:
            payload = pickle.load(f)
        if bateria_compativel(payload):
            return payload
        print(
            f"bateria.pkl é de K={payload.get('K', '?')} "
            f"(z com {np.asarray(payload.get('melhor', {}).get('z', [])).reshape(-1).size} "
            f"componentes); o modelo atual tem K={P.K}, n_z={N_Z}. Recalculando."
        )
    return executar(semente=semente)


if __name__ == "__main__":
    dados = executar()
    print("jac eq/ineq", dados["jac_err_eq"], dados["jac_err_ineq"])
    print("resumo", dados["resumo"])
    print("melhor", dados["melhor"]["metodo"], dados["melhor"]["tf"], dados["melhor"]["sucesso"])
    print("figuras", dados["figuras"])
