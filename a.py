import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# El 1.2 se puede cambiar por el K, el K es a criterio propio
def funtion(xi, yi, M, E):
    dydx = 1.2*yi*(1 -(yi/M)) - E*yi
    
    return dydx

def runge_kutta(x0, y0, h, num_sol, M, E):
    x = np.zeros(num_sol)
    y = np.zeros(num_sol)
    x[0] = x0
    y[0] = y0
    
    for i in range(num_sol-1):
        k1 = funtion(x[i], y[i], M, E)
        k2 = funtion(x[i] + (1/2)*h, y[i] + (1/2)*h*k1, M, E)
        k3 = funtion(x[i] + (1/2)*h, y[i] + (1/2)*h*k2, M, E)
        k4 = funtion(x[i] + h, y[i] + h*k3, M, E)
        
        x[i+1] = x[i] + h
        y[i+1] = y[i] + (1/6)*h*(k1+ 2*k2 + 2*k3 + k4)
        
    return x,y
    
def stage(M, E, x0, y0, h, num_sol):
    for i in range(len(M)):
        solx, soly = runge_kutta(x0, y0, h, num_sol, M[i], E[i])
        print("Set", i+1, ":")
       #print("xn")
        print(solx)
        #print("yn")
        print(soly)
        plt.scatter(solx,soly)
        plt.show()
        
x0=0
y0=10
h=0.05
xfin = 12
num_sol = int((xfin- x0)/h)
M_A = [10.5, 10, 10.5]
E_A = [0.55, 0.60, 0.50]




#stage(M_A, E_A, x0, y0, h, num_sol)
#La tabla muestra como han disminuido los mosquitos dado un insecticida

M_B = [6, 5.5, 6.5]
E_B = [0.10, 0.10, 0.10]
stage(M_B, E_B, x0, y0, h, num_sol)

M_C = [9, 8.5, 9]
E_C = [0.4, 0.45, 0.4]
#stage(M_A, E_A, x0, y0, h, num_sol)