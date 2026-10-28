"""Localiza a pasta T2 para os notebooks, independente do cwd."""

import sys
from pathlib import Path


def t2_root() -> Path:
    here = Path.cwd().resolve()
    candidatos = [here, here.parent, here / "T2", here.parent / "T2"]
    for cand in candidatos:
        if (cand / "src" / "parametros.py").exists():
            return cand
    return Path(__file__).resolve().parent.parent


def preparar():
    root = t2_root()
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    from src.paths import figuras_dir, saida_dir

    figuras_dir()
    saida_dir()
    return root
