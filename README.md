# Modelo de la formación de radios en la aleta caudal del pez cebra

Para el estudio, y por tanto la modelización, de las formaciones de radios en la aleta caudal del pez cebra, es sabido que sigue un sistema de ecuaciones de derivadas parciales como el que sigue:
<p align="center"><img src="https://rawgit.com/rafafrdz/Pez-Cebra/master/svgs/3fd1cf893dd111654a8d539bf87cc0d0.svg?invert_in_darkmode" align=middle width=350.29005pt height=35.777445pt/></p>

<p align="center"><img src="https://rawgit.com/rafafrdz/Pez-Cebra/master/svgs/23b7c9de8e03fd79f519af67b349900e.svg?invert_in_darkmode" align=middle width=239.6757pt height=35.777445pt/></p>

donde:

* La unidad de la variable <img src="https://rawgit.com/rafafrdz/Pez-Cebra/master/svgs/178e553be50b61dae0734edfd972a567.svg?invert_in_darkmode" align=middle width=15.33114pt height=20.22207pt/> es la longitud característica de una célula
* La unidad de la variable <img src="https://rawgit.com/rafafrdz/Pez-Cebra/master/svgs/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode" align=middle width=5.9361555pt height=20.22207pt/> es el tiempo característico en el que se produce una división celular
* <img src="https://rawgit.com/rafafrdz/Pez-Cebra/master/svgs/d97f9c2920507770c422858be1cbe2ef.svg?invert_in_darkmode" align=middle width=43.29534pt height=24.6576pt/> es la concentración de una sustancia activadora de la formación de células de radio en <img src="https://rawgit.com/rafafrdz/Pez-Cebra/master/svgs/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.3951pt height=14.15535pt/> en el instante <img src="https://rawgit.com/rafafrdz/Pez-Cebra/master/svgs/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode" align=middle width=5.9361555pt height=20.22207pt/>, de forma que, una célula pasa a ser del tipo radio cuando el valor de r es mayor que un cierto valor crítico <img src="https://rawgit.com/rafafrdz/Pez-Cebra/master/svgs/45afd16d417643a14faef303967f027b.svg?invert_in_darkmode" align=middle width=7.8730245pt height=19.41522pt/>
* <img src="https://rawgit.com/rafafrdz/Pez-Cebra/master/svgs/883d91f1b04662f98e1ea7399dd6fe8e.svg?invert_in_darkmode" align=middle width=43.132815pt height=24.6576pt/> es la concentración de una sustancia inhibidora de la formación de células de radio en <img src="https://rawgit.com/rafafrdz/Pez-Cebra/master/svgs/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.3951pt height=14.15535pt/> en el instante <img src="https://rawgit.com/rafafrdz/Pez-Cebra/master/svgs/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode" align=middle width=5.9361555pt height=20.22207pt/>
* <img src="https://rawgit.com/rafafrdz/Pez-Cebra/master/svgs/3a6198dde5b4af22c79a43f92c59e0b4.svg?invert_in_darkmode" align=middle width=136.61076pt height=26.94318pt/> son constantes positivas

Se trata por tanto de un modelo de reacción-difusión: Las sustancias <img src="https://rawgit.com/rafafrdz/Pez-Cebra/master/svgs/89f2e0d2d24bcf44db73aab8fc03252c.svg?invert_in_darkmode" align=middle width=7.8730245pt height=14.15535pt/> y <img src="https://rawgit.com/rafafrdz/Pez-Cebra/master/svgs/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710483pt height=21.68331pt/> son creadas y destruidas en las células y pasan de unas a otras a través de sus membranas por difusión.

```python
# -*- coding: utf-8 -*-
from numpy import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D

# Apartado c

def u0(x):
    return 1.0 * (10 <= x) * (x <= 20)

def v0(x):
    return -1.0 * (10 <= x) * (x <= 20)

def programa(L, T, Nx, u0, v0, delta, gamma, l):
    L = float(L)
    T = float(T)
    Nx = int(Nx)
    delta = float(delta)
    gamma = float(gamma)
    l = float(l)

    x = linspace(0, L, Nx + 1)  # partición del intervalo [0,L]
    dx = float(L / Nx)  # diferencia entre los distintos nodos
    dx2 = dx * dx

    # diferencia entre los distintos niveles de tiempo de manera que esté bajo la condicion CFL y así sea consistente y estable
    dt = min(0.39*dx2/(2*l), dx2/2)
    Mt = int(T / dt)

    # Los distintos valores para u^n (comenzando con los que proporciona u(x,0), es decir u0(x) "u^0")
    solu0 = u0(x)
    solv0 = v0(x)

    # Los distintos valores para u^(n+1)
    solu1 = 0 * solu0
    solv1 = 0 * solv0

    for n in range(1, Mt):
        tn = n * dt
        for i in range(1, Nx):

            # G(xi,tn) = solu0[i]*(1 - pow(solu0[i],2))-solv0[i]
            # H(xi,tn) = gamma * solu0[i] - delta * solv0[i]
            funG = solu0[i] * (1 - pow(solu0[i], 2)) - solv0[i]
            funH = gamma * solu0[i] - delta * solv0[i]
            solu1[i] = solu0[i] + dt * funG + dt / dx2 * (solu0[i - 1] - 2 * solu0[i] + solu0[i + 1])
            solv1[i] = solv0[i] + dt*funH + dt*l/dx2 * (solv0[i-1]-2*solv0[i] + solv0[i+1])

            # Condiciones frontera para u
            funG0 = solu0[0] * (1 - pow(solu0[0], 2)) - solv0[0]
            funGNx = solu0[Nx] * (1 - pow(solu0[Nx], 2)) - solv0[Nx]
            solu1[0] = solu0[0] + dt * funG0 + 2*dt/dx2*(solu0[1]-solu0[0])
            solu1[Nx] = solu0[Nx] + dt * funGNx + 2 * dt / dx2 * (solu0[Nx-1] - solu0[Nx])
            # Condiciones frontera para v
            funH0 = gamma * solu0[0] - delta * solv0[0]
            funHNx = gamma * solu0[Nx] - delta * solv0[Nx]
            solv1[0] = solv0[0] + dt*funH0 + 2*dt*l/dx2 * (solv0[1]-solv0[0])
            solv1[Nx] = solv0[Nx] + dt * funHNx + 2 * dt * l / dx2 * (solv0[Nx-1] - solv0[Nx])

        solu0 = solu1
        solv0 = solv1
    print('Para u: ' + str(solu0))
    print('Para v: ' + str(solv0))
    plot(x,solu0)
    plot(x,solv0)
    show()

# Apartado d
# Hemos tomado T = 10 y una partición del intervalo [0,L] en 70 nodos
programa(30,10,70,u0,v0,1.78,1.8,12)

```

