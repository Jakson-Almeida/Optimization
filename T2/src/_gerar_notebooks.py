"""Gera os cinco notebooks do T2. Executar da raiz do repositório ou de T2."""

from pathlib import Path

import nbformat as nbf

from .ambiente import t2_root

BOOT = """\
from pathlib import Path
import sys
here = Path.cwd().resolve()
for cand in [here, here.parent, here / "T2", here.parent / "T2"]:
    if (cand / "src" / "parametros.py").exists():
        if str(cand) not in sys.path:
            sys.path.insert(0, str(cand))
        break
from src.ambiente import preparar
T2 = preparar()
print("T2 =", T2)
"""


def md(nb, text):
    nb.cells.append(nbf.v4.new_markdown_cell(text))


def py(nb, text):
    nb.cells.append(nbf.v4.new_code_cell(text))


def nb01():
    nb = nbf.v4.new_notebook()
    md(
        nb,
        r"""# 1. Modelagem — tempo mínimo cinemático

Problema alinhado à seção 4.2 de Massaro et al., *Robotics* 2023: levar um
manipulador planar de dois elos de \(\mathbf{q}(0)=\mathbf{0}\) a
\(\mathbf{q}(t_f)=(1,-0{,}5)\) rad, em **repouso–repouso**, no menor \(t_f\),
sem colidir com um disco no plano de trabalho.

**O que entra.** Cinemática direta, limites de velocidade e de aceleração,
um obstáculo circular, três pontos monitorados no segundo elo.

**O que fica de fora (simplificação em relação ao paper).** Dinâmica de
Lagrange, jerk, colocação LGR/GPOPS-II. A transcrição é aceleração constante
por intervalo com integração exata da cadeia \(\ddot q=a\).
""",
    )
    py(nb, BOOT)
    md(
        nb,
        r"""## OCP contínuo (Bolza)

Estado \(\mathbf{x}=\{q_1,q_2,v_{q_1},v_{q_2}\}^\mathsf{T}\), controle
\(\mathbf{u}=\{a_{q_1},a_{q_2}\}^\mathsf{T}\),

\[
\dot q_i=v_{q_i},\qquad \dot v_{q_i}=a_{q_i},\qquad
J=\int_0^{t_f} 1\,dt=t_f.
\]

Caminho: \(|\dot q_1|\le 3\), \(|\dot q_2|\le 8\) rad/s, \(|a_i|\le 18\) rad/s²,
e para \(\zeta\in\{1/3,2/3,1\}\)

\[
(x_P(\mathbf{q},\zeta)-x_0)^2+(y_P-y_0)^2\ge R^2,
\]

com \(x_P=l_1\cos q_1+\zeta l_2\cos(q_1+q_2)\) (eq. 6 do paper).
""",
    )
    py(
        nb,
        """\
from src.parametros import P, ZETA, N_Z, N_EQ, N_INEQ
from src.cinematica import ponto_elo2, juntas_elo2
from src.transcricao import (
    chute_desvio, chute_reta, construir_chutes, avaliar, checar_jacobiana, n_decisao
)

print("l1, l2 =", P.l1, P.l2)
print("q0, qf =", P.q0, P.qf)
print("v_max, a_max =", P.v_max, P.a_max)
print("obstáculo =", (P.x0_obs, P.y0_obs, P.R_obs), "zeta =", ZETA)
print("K =", P.K, "dimensão de z =", n_decisao(), "igualdades =", N_EQ, "desigualdades =", N_INEQ)
""",
    )
    md(
        nb,
        r"""## PNL discreta

Decisão \(\mathbf{z}=(t_f,\,a_{1,0},\ldots,a_{1,K-1},\,a_{2,0},\ldots,a_{2,K-1})\in\mathbb{R}^{21}\),
\(K=10\), \(h=t_f/K\):

\[
v_{i,k+1}=v_{i,k}+h a_{i,k},\qquad
q_{i,k+1}=q_{i,k}+h v_{i,k}+\tfrac12 h^2 a_{i,k},
\]

com \(q_{i,0}\) e \(v_{i,0}=0\) fixos. Igualdades: \(q_{i,K}=q_{i,f}\), \(v_{i,K}=0\).
Obstáculo nos nós **e** nos pontos médios dos intervalos.
""",
    )
    py(
        nb,
        """\
z_reta = chute_reta()
z_desvio = chute_desvio()
a_reta = avaliar(z_reta)
a_desvio = avaliar(z_desvio)

def relatorio(nome, a):
    tr = a["trajetoria"]
    print(f"{nome:18s} tf={tr.tf:.3f}s  eq={a['violacao_eq']:.2e}  ineq={a['violacao_ineq']:.2e}  "
          f"viavel={a['viavel']}  colide={a['colide']}  qf={tr.q[:,-1]}  vf={tr.v[:,-1]}")

relatorio("reta_juntas", a_reta)
relatorio("desvio_esquerda", a_desvio)
print()
print("cinco chutes da bateria:")
for nome, z in construir_chutes():
    relatorio(nome, avaliar(z))
""",
    )
    py(
        nb,
        """\
err_eq, err_ineq = checar_jacobiana(z_desvio)
print(f"erro relativo da Jacobiana  igualdades={err_eq:.2e}  desigualdades={err_ineq:.2e}")
print("A Jacobiana analítica é a que os pontos interiores recebem.")
""",
    )
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    return nb


def nb02():
    nb = nbf.v4.new_notebook()
    md(
        nb,
        r"""# 2. Convexidade — Jacobiana e Hessiana analíticas

Antes de otimizar, classifica-se cada peça da PNL. As derivadas saem do SymPy,
sem diferenças finitas. A figura do C-space é a “cara” da região viável em
\((q_1,q_2)\); a PNL em \(\mathbb{R}^{21}\) é um levantamento dessa geometria
com \(t_f\) e as acelerações.
""",
    )
    py(nb, BOOT)
    py(
        nb,
        """\
import numpy as np
from IPython.display import display, Math
from src.convexidade import (
    classificar_autovalores, grade_cspace, hessiana_folga, hessiana_objetivo,
    hessiana_residual, jacobiana_folga, latex_pecas, numerica_folga,
)
from src.parametros import P
from src.paths import figuras_dir
from src.simulacao import figura_cspace

lx = latex_pecas()
display(Math(r"\\text{FK }x_P=" + lx["fk_x"]))
display(Math(r"\\text{FK }y_P=" + lx["fk_y"]))
display(Math(r"g(\\mathbf{q},\\zeta)=" + lx["g"]))
display(Math(r"\\nabla_q g=" + lx["J_g"]))
display(Math(r"\\nabla^2_q g=" + lx["H_g"]))
""",
    )
    md(
        nb,
        r"""## Classificação das peças

| Peça | Forma | Convexidade |
|------|--------|-------------|
| \(J=t_f\) | linear | convexa (Hessiana nula) |
| caixas de \(t_f\) e \(a\) | intervalos | convexas |
| residual \(v_{k+1}-v_k-(t_f/K)a_k\) | bilinear | nem convexo nem côncavo |
| \(g(\mathbf{q},\zeta)=R^2-\|\mathbf{r}-\mathbf{r}_0\|^2\) | trigonométrica | Hessiana indefinida |
| \(\{g\le 0\}\) no plano das juntas | complementar de um disco via FK | **não convexo** |
""",
    )
    py(
        nb,
        """\
print("Hessiana de J = tf:")
print(hessiana_objetivo())

Hr, vars_r = hessiana_residual()
print("\\nHessiana do residual bilinear nas variáveis", vars_r)
print(Hr)

g, J, H = numerica_folga(0.7, -0.4, zeta_val=1.0)
w, tipo = classificar_autovalores(H)
print("\\nEm (q1,q2)=(0.7,-0.4), zeta=1:")
print("g =", g)
print("J =", J)
print("H =\\n", H)
print("autovalores", w, "→", tipo)
""",
    )
    py(
        nb,
        """\
from src.transcricao import integrar
q1g, q2g, G = grade_cspace(n1=201, n2=201)
caminho = figuras_dir() / "cspace_regiao.png"
figura_cspace(None, caminho, q1g=q1g, q2g=q2g, G=G)
print("figura:", caminho)
print("fração do retângulo em colisão:", float((G > 0).mean()))
print("Conclusão: a PNL não é convexa. Mínimos locais distintos são esperados;")
print("as cinco partidas do notebook 03 não são cosmética.")
""",
    )
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    return nb


def nb03():
    nb = nbf.v4.new_notebook()
    md(
        nb,
        r"""# 3. Otimização — bateria de métodos

Mesma PNL, três métodos, **cinco partidas** compartilhadas (semente 0):

1. reta no espaço das juntas (inviável: atravessa o disco);
2. desvio pela esquerda, viável;
3–5. ruído gaussiano sobre o desvio, reamostrado até ficar viável.

| Método | Implementação |
|--------|----------------|
| Pontos interiores | `trust-constr` (barreira) com Jacobiana analítica, polido por SQP até um ponto KKT |
| Pattern search | COBYLA (Powell), sem derivadas; para quando o raio da região de confiança atinge o mínimo |
| Newton | SQP (`SLSQP`): Newton no sistema KKT com Hessiana quasi-Newton do Lagrangiano |

Sucesso = ponto **viável** e parada de ótimo local (`gtol` / `xtol` / KKT / raio mínimo), não o teto de iterações.

Se `T2/saida/bateria.pkl` já existe, o notebook o recarrega. Apague o arquivo
para rerrodar os 15 problemas.
""",
    )
    py(nb, BOOT)
    py(
        nb,
        """\
import pickle
import numpy as np
from src.experiencia import executar
from src.otimizadores import agregar
from src.paths import saida_dir

pkl = saida_dir() / "bateria.pkl"
if pkl.exists():
    with open(pkl, "rb") as f:
        payload = pickle.load(f)
    print("carregado", pkl)
else:
    payload = executar()
    print("bateria executada e gravada")

registros = payload["registros"]
resumo = payload["resumo"]
print("erro da Jacobiana eq/ineq:", payload["jac_err_eq"], payload["jac_err_ineq"])
""",
    )
    py(
        nb,
        """\
print(f"{'método':20s} {'chute':18s} {'tf':>8s} {'conv':>5s} {'viável':>7s} {'eq':>10s} {'ineq':>10s} {'iter':>6s} {'CPU':>8s}")
for r in registros:
    print(f"{r['metodo']:20s} {r['chute']:18s} {r['tf']:8.4f} {str(r['convergiu']):>5s} {str(r['viavel']):>7s} "
          f"{r['violacao_eq']:10.2e} {r['violacao_ineq']:10.2e} {r['niter']:6d} {r['cpu_s']:8.3f}")
""",
    )
    py(
        nb,
        """\
print(f"{'método':20s} {'ok':>6s} {'tf min':>8s} {'tf méd':>8s} {'tf máx':>8s} {'std':>8s} {'CPU méd':>8s}")
for m, s in resumo.items():
    print(f"{m:20s} {s['n_sucesso']}/{s['n']:<3d} {s['tf_min']:8.4f} {s['tf_media']:8.4f} "
          f"{s['tf_max']:8.4f} {s['tf_std']:8.4f} {s['cpu_media']:8.3f}")
melhor = payload["melhor"]
print(f"\\nmelhor (viável e convergido): {melhor['metodo']} / {melhor['chute']}  tf={melhor['tf']:.4f} s")
""",
    )
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    return nb


def nb04():
    nb = nbf.v4.new_notebook()
    md(
        nb,
        r"""# 4. Simulação do manipulador

Animação, estroboscópio (tipo Fig. 6 do paper) e perfis \(q(t)\), \(\dot q(t)\)
da melhor trajetória viável da bateria. Figuras em `T2/figuras/`.
""",
    )
    py(nb, BOOT)
    py(
        nb,
        """\
import pickle
from IPython.display import Image, display
from src.paths import figuras_dir, saida_dir
from src.simulacao import salvar_tudo
from src.transcricao import integrar

with open(saida_dir() / "bateria.pkl", "rb") as f:
    payload = pickle.load(f)
melhor = payload["melhor"]
tr = integrar(melhor["z"])
print(f"método={melhor['metodo']}  chute={melhor['chute']}  tf={tr.tf:.4f} s  viável={melhor['viavel']}")
figs = salvar_tudo(tr, figuras_dir(), prefixo="melhor")
for k, p in figs.items():
    print(k, p)
""",
    )
    py(
        nb,
        """\
display(Image(filename=str(figuras_dir() / "melhor_estroboscopio.png")))
display(Image(filename=str(figuras_dir() / "melhor_cspace.png")))
display(Image(filename=str(figuras_dir() / "melhor_perfis.png")))
""",
    )
    py(
        nb,
        """\
from IPython.display import HTML
HTML(f'<img src=\"../figuras/melhor.gif\" width=\"480\">')
""",
    )
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    return nb


def nb05():
    nb = nbf.v4.new_notebook()
    md(
        nb,
        r"""# 5. Resultados da programação não linear

Leitura da bateria: o que cada método consegue na PNL cinemática, e como isso
se compara ao paper (\(t_f\approx 1{,}180\) s, obstáculo com \(N=3\)).
""",
    )
    py(nb, BOOT)
    py(
        nb,
        """\
import pickle
from src.paths import saida_dir

with open(saida_dir() / "bateria.pkl", "rb") as f:
    payload = pickle.load(f)
resumo = payload["resumo"]
melhor = payload["melhor"]

print(f"{'método':20s} {'ok':>6s} {'tf min':>8s} {'tf méd':>8s} {'tf máx':>8s} {'CPU méd':>8s}")
for m, s in resumo.items():
    print(f"{m:20s} {s['n_sucesso']}/{s['n']:<3d} {s['tf_min']:8.4f} {s['tf_media']:8.4f} "
          f"{s['tf_max']:8.4f} {s['cpu_media']:8.3f}")
print(f"\\nmelhor: {melhor['metodo']} / {melhor['chute']}  tf={melhor['tf']:.4f} s")
print("paper, cinemática com obstáculo N=3: tf ≈ 1.180 s")
""",
    )
    md(
        nb,
        r"""## Leitura

- **Pontos interiores** (`trust-constr` + polimento SQP) param em um ponto KKT
  nas quatro partidas viáveis. A Jacobiana analítica das restrições entra no
  solver; o polimento fecha o critério de ótimo local que a barreira sozinha
  não atinge neste problema bang-bang.
- **Pattern search** (COBYLA) também converge nas quatro partidas viáveis: a
  parada é o raio mínimo da região de confiança, não o teto de avaliações. Sem
  derivadas o \(t_f\) local fica maior (\(\approx 2{,}4\)–\(2{,}8\) s).
- **Newton-SQP** (`SLSQP`) fecha KKT nas mesmas quatro partidas, com
  \(t_f\approx 1{,}194\) s. A reta no espaço das juntas atravessa o disco e os
  três métodos falham ali (não há direção viável de descida a partir desse chute).
- **Paper.** \(t_f\approx 1{,}18\) s com colocação LGR e malha bem mais fina.
  O melhor desta transcrição, \(t_f\approx 1{,}192\) s, está na mesma ordem.
  A passagem é pela esquerda do disco no C-space, como no caso \(N=3\).
""",
    )
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    return nb


def main():
    pasta = t2_root() / "notebooks"
    pasta.mkdir(parents=True, exist_ok=True)
    arquivos = {
        "01_modelagem.ipynb": nb01(),
        "02_convexidade.ipynb": nb02(),
        "03_otimizacao.ipynb": nb03(),
        "04_simulacao.ipynb": nb04(),
        "05_resultados.ipynb": nb05(),
    }
    for nome, nb in arquivos.items():
        path = pasta / nome
        with open(path, "w", encoding="utf-8") as f:
            nbf.write(nb, f)
        print("escreveu", path)


if __name__ == "__main__":
    main()
