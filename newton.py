# Max f(x)=sen(x)+cos(8x)

import numpy as np

cont = 0
x = 10
xH = x + 1e7

def deri2(x):
    return -np.sin(x)-64*np.cos(8*x)

def deltaY(x, xH):
    return np.cos(x)-8*np.sin(8*x)-(np.cos(xH)-8*np.sin(8*xH))

dY = deltaY(x, xH)

while np.abs(deltaY(x, xH)) > 10e-6:
    cont += 1
    dX = np.power(np.abs(deri2(x)), -1) * deltaY(x, xH)
    x_novo = x + dX
    dY = deltaY(x, x_novo)
    print(cont)

print(x)