"""Pastas do trabalho T2, independentes do cwd."""

from pathlib import Path


def t2_root() -> Path:
    return Path(__file__).resolve().parent.parent


def figuras_dir() -> Path:
    d = t2_root() / "figuras"
    d.mkdir(parents=True, exist_ok=True)
    return d


def saida_dir() -> Path:
    d = t2_root() / "saida"
    d.mkdir(parents=True, exist_ok=True)
    return d
