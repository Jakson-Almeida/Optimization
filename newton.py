# Max f(x) = sen(x) + cos(8x)
# Ponto estacionário por Newton: x <- x - f'(x) / f''(x)

import numpy as np


def f(x):
    return np.sin(x) + np.cos(8 * x)


def df(x):
    return np.cos(x) - 8 * np.sin(8 * x)


def d2f(x):
    return -np.sin(x) - 64 * np.cos(8 * x)


def newton(x0, tol=1e-10, max_iter=100):
    x = float(x0)
    for cont in range(1, max_iter + 1):
        derivada = df(x)
        segunda = d2f(x)
        if abs(derivada) <= tol:
            return x, cont - 1
        if abs(segunda) < 1e-14:
            raise ZeroDivisionError(f"f''({x}) é nula; Newton não pode continuar.")
        x = x - derivada / segunda
        print(cont, x)
    raise RuntimeError(f"Newton não convergiu a partir de x0 = {x0}.")


def classificar(x):
    segunda = d2f(x)
    if segunda < 0:
        return "máximo local"
    if segunda > 0:
        return "mínimo local"
    return "ponto inconclusivo"


if __name__ == "__main__":
    for x0 in (3, 45):
        x, cont = newton(x0)
        print(
            f"x0 = {x0}: x* = {x:.10f}, f(x*) = {f(x):.10f}, "
            f"f'(x*) = {df(x):.3e}, f''(x*) = {d2f(x):.6f}, "
            f"{classificar(x)}, iterações = {cont}"
        )
