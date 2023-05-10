import numpy as npy
import matplotlib.pyplot as mat
import pandas as pds


def derivada(xi, yi, M, E):
    dy_dx = 1.2 * yi * (1 -(yi/M)) - E * yi
    return dy_dx

def algoritmo_runge_kutta(x0, y0, h, sol, M, E):
    x = npy.zeros(sol)
    y = npy.zeros(sol)
    x[0] = x0
    y[0] = y0

    for i in range(sol-1):
        k1 = derivada(x[i], y[i], M, E)
        k2 = derivada(x[i] + (1/2) * h, y[i] + (1/2) * h * k1, M, E)
        k3 = derivada(x[i] + (1/2) * h, y[i] + (1/2) * h * k2, M, E)
        k4 = derivada(x[i] + h, y[i] + h * k3, M, E)

        x[i+1] = x[i] + h
        y[i+1] = y[i] + (1/6) * h * (k1+ 2 * k2 + 2 * k3 + k4)

    return x,y

def resultado(M, E, x0, y0, h, sol):
    for i in range(len(M)):
        solx, soly = algoritmo_runge_kutta(x0, y0, h, sol, M[i], E[i])
        print("Set", i+1, ":")
        print("xn")
        print(solx)
        print("yn")
        print(soly)
        mat.scatter(solx,soly)
        mat.show()

x0 = 0
y0 = 10
h = 0.05
xfin = 12
sol = int((xfin - x0)/h)

print("Escenario A")
M_A = [10.5, 10, 10.5]
E_A = [0.55, 0.60, 0.50]
resultado(M_A, E_A, x0, y0, h, sol)

print("Escenario B")
M_B = [6, 5.5, 6.5]
E_B = [0.10, 0.10, 0.10]
#resultado(M_B, E_B, x0, y0, h, sol)

print("Escenario C")
M_C = [9, 8.5, 9]
E_C = [0.4, 0.45, 0.4]
#resultado(M_A, E_A, x0, y0, h, sol)
