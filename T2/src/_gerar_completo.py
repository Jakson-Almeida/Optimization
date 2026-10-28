"""Gera T2_completo.ipynb: um arquivo autossuficiente, sem pasta src/."""

from __future__ import annotations

import re
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"


def strip_relativos(texto: str) -> str:
    saida = []
    pulando = False
    for linha in texto.splitlines():
        if pulando:
            if ")" in linha:
                pulando = False
            continue
        if re.match(r"^\s*from \.", linha):
            if "(" in linha and ")" not in linha:
                pulando = True
            continue
        saida.append(linha)
    return "\n".join(saida).strip() + "\n"


def ler(nome: str) -> str:
    return strip_relativos((SRC / nome).read_text(encoding="utf-8"))


def md(nb, texto):
    nb.cells.append(nbf.v4.new_markdown_cell(texto))


def py(nb, texto):
    nb.cells.append(nbf.v4.new_code_cell(texto))


def construir():
    nb = nbf.v4.new_notebook()
    nb.metadata["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }

    md(
        nb,
        r"""# Tempo mínimo cinemático de um manipulador planar — PNL

Notebook **autossuficiente** (Massaro et al., *Robotics* 2023, seção 4.2).
Não depende da pasta `src/` nem do restante do repositório: basta este `.ipynb`.

Dois elos, repouso–repouso, um disco, três pontos no segundo elo. A dinâmica de
Lagrange e a colocação LGR ficam de fora; a transcrição é aceleração constante
por intervalo.

**No Colab:** Runtime → Run all. A célula de instalação é inofensiva se os
pacotes já existirem. A bateria de otimização leva da ordem de um minuto.
""",
    )
    py(
        nb,
        """\
# Colab / máquina nova: instala o que faltar. No ambiente local, pode ignorar.
import sys
if "google.colab" in sys.modules:
    %pip install -q numpy scipy sympy matplotlib pillow
""",
    )
    py(
        nb,
        """\
from collections import deque
from dataclasses import dataclass
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from IPython.display import Image, display
from matplotlib import animation, patches
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from scipy.optimize import Bounds, NonlinearConstraint, minimize

SAIDA = Path.cwd() / "t2_saida"
SAIDA.mkdir(parents=True, exist_ok=True)
print("figuras em", SAIDA.resolve())
""",
    )
    md(nb, "## Código do modelo\n\nParâmetros, cinemática, transcrição da PNL, convexidade, solvers e desenho.")
    py(nb, ler("parametros.py"))
    py(nb, ler("cinematica.py"))
    py(nb, ler("transcricao.py"))
    py(nb, ler("convexidade.py"))
    py(nb, ler("otimizadores.py"))
    py(nb, ler("simulacao.py"))

    md(
        nb,
        r"""## 1. Modelagem

OCP contínuo: estado \(\mathbf{x}=\{q_1,q_2,v_{q_1},v_{q_2}\}^\mathsf{T}\),
controle \(\mathbf{u}=\{a_{q_1},a_{q_2}\}^\mathsf{T}\),
\(\dot q_i=v_{q_i}\), \(\dot v_{q_i}=a_{q_i}\), \(J=t_f\).

PNL: \(\mathbf{z}=(t_f,a_{1,0},\ldots,a_{1,K-1},a_{2,0},\ldots,a_{2,K-1})\in\mathbb{R}^{21}\),
\(K=10\), \(h=t_f/K\),

\[
v_{i,k+1}=v_{i,k}+h a_{i,k},\qquad
q_{i,k+1}=q_{i,k}+h v_{i,k}+\tfrac12 h^2 a_{i,k}.
\]

Igualdades \(q_{i,K}=q_{i,f}\), \(v_{i,K}=0\). Obstáculo nos nós e nos pontos médios.
""",
    )
    py(
        nb,
        """\
print("l1, l2 =", P.l1, P.l2)
print("q0, qf =", P.q0, P.qf)
print("v_max, a_max =", P.v_max, P.a_max)
print("obstáculo =", (P.x0_obs, P.y0_obs, P.R_obs), "zeta =", ZETA)
print("K =", P.K, "dimensão de z =", n_decisao(), "igualdades =", N_EQ, "desigualdades =", N_INEQ)

z_reta = chute_reta()
z_desvio = chute_desvio()


def relatorio(nome, a):
    tr = a["trajetoria"]
    print(
        f"{nome:18s} tf={tr.tf:.3f}s  eq={a['violacao_eq']:.2e}  "
        f"ineq={a['violacao_ineq']:.2e}  viavel={a['viavel']}  colide={a['colide']}"
    )


relatorio("reta_juntas", avaliar(z_reta))
relatorio("desvio_esquerda", avaliar(z_desvio))
print("cinco chutes:")
for nome, z in construir_chutes():
    relatorio(nome, avaliar(z))

err_eq, err_ineq = checar_jacobiana(z_desvio)
print(f"erro relativo da Jacobiana  eq={err_eq:.2e}  ineq={err_ineq:.2e}")
""",
    )

    md(
        nb,
        r"""## 2. Convexidade

Derivadas pelo SymPy, sem diferenças finitas. A figura do C-space é a cara da
região viável em \((q_1,q_2)\).

| Peça | Forma | Convexidade |
|------|--------|-------------|
| \(J=t_f\) | linear | convexa (Hessiana nula) |
| caixas de \(t_f\) e \(a\) | intervalos | convexas |
| residual \(v_{k+1}-v_k-(t_f/K)a_k\) | bilinear | nem convexo nem côncavo |
| \(g(\mathbf{q},\zeta)\) | trigonométrica | Hessiana indefinida |
| \(\{g\le 0\}\) no plano das juntas | complementar de um disco via FK | **não convexo** |
""",
    )
    py(
        nb,
        """\
print("Hessiana de J = tf:")
print(hessiana_objetivo())
Hr, vars_r = hessiana_residual()
print("Hessiana do residual bilinear em", vars_r)
print(Hr)

g, Jg, Hg = numerica_folga(0.7, -0.4, zeta_val=1.0)
w, tipo = classificar_autovalores(Hg)
print("Em (q1,q2)=(0.7,-0.4), zeta=1: g =", g)
print("autovalores", w, "→", tipo)

q1g, q2g, G = grade_cspace(n1=161, n2=161)
cspace_png = SAIDA / "cspace_regiao.png"
figura_cspace(None, cspace_png, q1g=q1g, q2g=q2g, G=G)
print("fração do retângulo em colisão:", float((G > 0).mean()))
display(Image(filename=str(cspace_png)))
print("A PNL não é convexa; as cinco partidas abaixo não são cosmética.")
""",
    )

    md(
        nb,
        r"""## 3. Otimização

Mesma PNL, três métodos, **cinco partidas** (semente 0): reta inviável, desvio
pela esquerda, e três ruídos.

| Método | Implementação |
|--------|----------------|
| Pontos interiores | `trust-constr` com Jacobiana analítica |
| Pattern search | Hooke–Jeeves na penalidade quadrática |
| Newton | Newton–Levenberg na mesma penalidade, Hessiana analítica |
""",
    )
    py(
        nb,
        """\
print("Rodando 5 chutes × 3 métodos...")
registros, _ = rodar_bateria(semente=0)
resumo = agregar(registros)
melhor = melhor_viavel(registros)

print(f"{'método':20s} {'chute':18s} {'tf':>8s} {'viável':>7s} {'eq':>10s} {'ineq':>10s} {'iter':>6s} {'CPU':>8s}")
for r in registros:
    print(
        f"{r['metodo']:20s} {r['chute']:18s} {r['tf']:8.4f} {str(r['viavel']):>7s} "
        f"{r['violacao_eq']:10.2e} {r['violacao_ineq']:10.2e} {r['niter']:6d} {r['cpu_s']:8.3f}"
    )
print()
print(f"{'método':20s} {'ok':>6s} {'tf min':>8s} {'tf méd':>8s} {'tf máx':>8s} {'std':>8s} {'CPU méd':>8s}")
for m, s in resumo.items():
    print(
        f"{m:20s} {s['n_sucesso']}/{s['n']:<3d} {s['tf_min']:8.4f} {s['tf_media']:8.4f} "
        f"{s['tf_max']:8.4f} {s['tf_std']:8.4f} {s['cpu_media']:8.3f}"
    )
print(f"\\nmelhor viável: {melhor['metodo']} / {melhor['chute']}  tf={melhor['tf']:.4f} s")
""",
    )

    md(
        nb,
        r"""## 4. Simulação

Estroboscópio, C-space, perfis e GIF da melhor trajetória viável.
""",
    )
    py(
        nb,
        """\
tr = integrar(melhor["z"])
print(f"método={melhor['metodo']}  chute={melhor['chute']}  tf={tr.tf:.4f} s  viável={melhor['viavel']}")
figs = salvar_tudo(tr, SAIDA, prefixo="melhor")
for k, p in figs.items():
    print(k, p)
display(Image(filename=str(SAIDA / "melhor_estroboscopio.png")))
display(Image(filename=str(SAIDA / "melhor_cspace.png")))
display(Image(filename=str(SAIDA / "melhor_perfis.png")))
display(Image(filename=str(SAIDA / "melhor.gif")))
""",
    )

    md(
        nb,
        r"""## 5. Resultados

Comparação com o paper (\(t_f\approx 1{,}180\) s, obstáculo com \(N=3\)).
""",
    )
    py(
        nb,
        """\
print(f"{'método':20s} {'ok':>6s} {'tf min':>8s} {'tf méd':>8s} {'tf máx':>8s} {'CPU méd':>8s}")
for m, s in resumo.items():
    print(
        f"{m:20s} {s['n_sucesso']}/{s['n']:<3d} {s['tf_min']:8.4f} {s['tf_media']:8.4f} "
        f"{s['tf_max']:8.4f} {s['cpu_media']:8.3f}"
    )
print(f"\\nmelhor: {melhor['metodo']} / {melhor['chute']}  tf={melhor['tf']:.4f} s")
print("paper, cinemática com obstáculo N=3: tf ≈ 1.180 s")
""",
    )
    md(
        nb,
        r"""**Leitura.** Pontos interiores (`trust-constr`) são o método que devolve
soluções viáveis com \(t_f\) competitivo; chutes diferentes caem em mínimos
locais distintos porque o C-space não é convexo. Pattern search e Newton, na
penalidade em \(\mathbb{R}^{21}\), não fecham viabilidade com o orçamento
usado. O \(t_f\) do paper (1,18 s) usa colocação LGR bem mais fina; o valor
deste notebook deve ser lido na mesma ordem de grandeza. O caminho contorna
o disco pela esquerda, como no caso \(N=3\).
""",
    )
    return nb


def main():
    nb = construir()
    destino = ROOT / "T2_completo.ipynb"
    with open(destino, "w", encoding="utf-8") as f:
        nbf.write(nb, f)
    print("escreveu", destino, "células", len(nb.cells))


if __name__ == "__main__":
    main()
