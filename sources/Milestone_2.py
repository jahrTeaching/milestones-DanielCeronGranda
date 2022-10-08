"""
CÁLCULO NUMÉRICO: MILESTONE 2
Prototypes to integrate orbits with functions

1) Euler Method
2) Crank-Nicholson Method
3) RK4 Method
4) Inverse Euler Method
5) Increase and decrease time step
"""
# IMPORTACIÓN
from numpy import array, zeros #Importación de funciones numéricas
from scipy.optimize import newton #Importación del método de Newton
import matplotlib.pyplot as plt #Importación de gráficas

#Funciones utilizadas en todos los apartados

#Función para integrar un paso de Euler
def Euler(U0, deltat, t, F0):
    U = U0 + deltat*F0
    return U

#Función para integrar un paso de Crank-Nicholson
def CN(U0, deltat, t, F0):
    def Res(X):
        return X-U0-deltat*0.5*(F0+ Kepler(X,t))
    return newton(Res, U0)

#Función para integrar un paso de Euler inverso
def Euler_inverso(U0, deltat, t, F0):
    def Res2(Y):
        return Y-U0-deltat*Kepler(Y,t)
    return newton(Res2, U0,tol = 1e-1)

#Función para integrar un paso de RK4
def RK4(U0, deltat, t, F0):
    k1 = F0 #Primer coeficiente
    U2 = U0 +  k1*deltat/2
    k2 = Kepler(U2,t) #Segundo coeficiente
    U3 = U0 +  k2*deltat/2
    k3 = Kepler(U3,t) #Tercer coeficiente
    U4 = U0 +  k3*deltat
    k4 = Kepler(U4,t) #Cuarto coeficiente
    U = U0+ deltat*(k1+2*k2+2*k3+k4)/6 #Método RK4 
    return U

#Ecuaciones de la órbita de Kepler
def Kepler(U,t):
    d = ((U[0])**2+(U[1])**2)**1.5 
    x = U[2]
    y = U[3]
    vx = -U[0]/d
    vy = -U[1]/d
    F = array([x,y,vx,vy])
    return F

#Función para integrar un problema de Cauchy con un esquema arbitrario
def Cauchy_problem(Temporal_scheme, U0, n, deltat, F0):
    t=0
    for i in range (n):
        t = t+deltat
        U = Temporal_scheme(U0, deltat, t, F0)
        matrizU[i,:] = U 
        U0 = U
        F0 = Kepler(U,t) #Obtengo mi nuevo valor de F(U,t)
        continue
    return (matrizU)

#Definición de condiciones iniciales
n = 200 # Número de particiones
deltat = 0.1 #Intervalo de tiempo entre iteraciones
U0 = array([1,0,0,1]) #Condiciones iniciales
matrizU = array(zeros([n,4])) #Matriz que contendrá en sus columnas los valores de las variables, cada fila es un paso de tiempo
F0 = array([0,1,-1,0])#Valor de F en el instante inicial
U = zeros(4)
x = zeros(n)
y = zeros(n)

# APARTADO 1: Euler
matrizU = Cauchy_problem(Euler, U0, n, deltat, F0)

#Las dos primeras columnas de la matriz son x e y de la órbita
x = matrizU[:,0]
y = matrizU[:,1]

#Representación gráfica
plt.plot(x,y)
plt.show()

#APARTADO 2: Crank-Nicholson
matrizU = Cauchy_problem(CN, U0, n, deltat, F0)

#Las dos primeras columnas de la matriz son x e y de la órbita
x = matrizU[:,0]
y = matrizU[:,1]

#Representación gráfica
plt.plot(x,y)
plt.show()

#APARTADO 3: RK4
matrizU = Cauchy_problem(RK4, U0, n, deltat, F0)

#Las dos primeras columnas de la matriz son x e y de la órbita
x = matrizU[:,0]
y = matrizU[:,1]

#Representación gráfica
plt.plot(x,y)
plt.show()

#APARTADO 4: Euler inverso
deltat=0.08
matrizU = Cauchy_problem(Euler_inverso, U0, n, deltat, F0)
#Las dos primeras columnas de la matriz son x e y de la órbita
x = matrizU[:,0]
y = matrizU[:,1]

#Representación gráfica
plt.plot(x,y)
plt.show()

#APARTADO 5: Aumentar y disminuir delta t
#Definición de condiciones iniciales
n = 100 # Número de particiones
deltat = 0.2 #Intervalo de tiempo entre iteraciones
U0 = array([1,0,0,1]) #Condiciones iniciales
matrizU = array(zeros([n,4])) #Matriz que contendrá en sus columnas los valores de las variables, cada fila es un paso de tiempo
F0 = array([0,1,-1,0])#Valor de F en el instante inicial
U = zeros(4)
x = zeros(n)
y = zeros(n)

Esquemas = [Euler, CN, RK4, Euler_inverso]

for i in Esquemas:
    matrizU = Cauchy_problem(i, U0, n, deltat, F0)
    #Las dos primeras columnas de la matriz son x e y de la órbita
    x = matrizU[:,0]
    y = matrizU[:,1]

    #Representación gráfica
    plt.plot(x,y)
    plt.show()
    continue

n = 400 # Número de particiones
deltat = 0.05 #Intervalo de tiempo entre iteraciones
U0 = array([1,0,0,1]) #Condiciones iniciales
matrizU = array(zeros([n,4])) #Matriz que contendrá en sus columnas los valores de las variables, cada fila es un paso de tiempo
F0 = array([0,1,-1,0])#Valor de F en el instante inicial
U = zeros(4)
x = zeros(n)
y = zeros(n)

for i in Esquemas:
    matrizU = Cauchy_problem(i, U0, n, deltat, F0)
    #Las dos primeras columnas de la matriz son x e y de la órbita
    x = matrizU[:,0]
    y = matrizU[:,1]

    #Representación gráfica
    plt.plot(x,y)
    plt.show()
    continue