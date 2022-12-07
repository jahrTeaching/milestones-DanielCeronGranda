from numpy import array, zeros, reshape #Importación de funciones numéricas
from numpy.linalg import norm
import matplotlib.pyplot as plt

#Ecuaciones de la órbita de Kepler
def Kepler(U,t):
    d = ((U[0])**2+(U[1])**2)**1.5 
    x = U[2]
    y = U[3]
    vx = -U[0]/d
    vy = -U[1]/d
    return array([x,y,vx,vy])

def Linear_oscillator(U,t):
    x = U[1]
    y = -U[0]
    return array([x,y])

 #Función de representación del oscilador lineal
def Wrapped_Repr_Linear_Oscillator(t, x, Scheme):
    plt.plot(t,x)
    plt.xlabel("t")
    plt.ylabel("x")
    i = "Oscilador lineal con esquema: "+Scheme.__name__
    plt.title(i)
    plt.show()
    return i

