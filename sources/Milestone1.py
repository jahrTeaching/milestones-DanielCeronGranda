""" 
CÁLCULO NUMÉRICO: MILESTONE 1
Prototypes to integrate orbits without functions

1)Euler method
2)Range-Kutta fourth order
3)Change time steps of the orbits
"""
# IMPORTACIÓN
from numpy import array, zeros #Importación de funciones numéricas
import matplotlib.pyplot as plt #Importación de gráficas

#APARTADO 1: EULER

#Definición de condiciones iniciales
n = 200 # Número de particiones
deltat = 0.1 #Intervalo de tiempo entre iteraciones
U0 = array([1,0,0,1]) #Condiciones iniciales
matrizU = array(zeros([n,4])) #Matriz que contendrá en sus columnas los valores de las variables, cada fila es un paso de tiempo
F = array([0,1,-1,0])#Valor de F en el instante inicial
U = zeros(4)
x = zeros(n)
y = zeros(n)

#Aplicación del método de Euler
for i in range(n):
    for j in range(4):
         U[j] = U0[j]+ deltat*F[j] #Método de Euler
         continue
    matrizU[i,:] = U #Almacena en una columna de la matriz el valor de las variables para el paso de tiempo analizado
    U0 = U #Cambio las condiciones iniciales del paso
    d = ((U[0])**2+(U[1])**2)**1.5
    F = array([U[2], U[3], -U[0]/d, -U[1]/d]) #Obtengo mi nuevo valor de F(U,t)
    continue

#Las dos primeras columnas de la matriz son x e y de la órbita
x = matrizU[:,0]
y = matrizU[:,1]

#Representación gráfica
plt.plot(x,y)
plt.show()

#APARTADO 2: RK4

#Definición de condiciones iniciales
n = 200 # Número de particiones
deltat = 0.1 #Intervalo de tiempo entre iteraciones
U0 = array([1,0,0,1]) #Condiciones iniciales
matrizU = array(zeros([n,4])) #Matriz que contendrá en sus columnas los valores de las variables, cada fila es un paso de tiempo
F = array([0,1,-1,0])#Valor de F en el instante inicial
U = zeros(4)
x = zeros(n)
y = zeros(n)
k1 = zeros(4)
k2 = zeros(4)
k3 = zeros(4)
k4 = zeros(4)

#Aplicación del método RK4
for i in range(n):
    k1 = F #Primer coeficiente
    U2 = U0 +  k1*deltat/2
    d2 = ((U2[0])**2+(U2[1])**2)**1.5
    k2 = array([U2[2], U2[3], -U2[0]/d2, -U2[1]/d2]) #Segundo coeficiente
    U3 = U0 +  k2*deltat/2
    d3 = ((U3[0])**2+(U3[1])**2)**1.5
    k3 = array([U3[2], U3[3], -U3[0]/d3, -U2[1]/d3]) #Tercer coeficiente
    U4 = U0 +  k3*deltat
    d4 = ((U4[0])**2+(U4[1])**2)**1.5
    k4 = array([U4[2], U4[3], -U4[0]/d4, -U4[1]/d4]) #Cuarto coeficiente
    U = U0+ deltat*(k1+2*k2+2*k3+k4)/6 #Método RK4 
    matrizU[i,:] = U #Almacena en una columna de la matriz el valor de las variables para el paso de tiempo analizado
    U0 = U #Cambio las condiciones iniciales del paso
    d = ((U[0])**2+(U[1])**2)**1.5
    F = array([U[2], U[3], -U[0]/d, -U[1]/d]) #Obtengo mi nuevo valor de F(U,t)
    continue

#Las dos primeras columnas de la matriz son x e y de la órbita
x = matrizU[:,0]
y = matrizU[:,1]

#Representación gráfica
plt.plot(x,y)
plt.show()

#APARTADO 3: Efecto del cambio del intervalo de delta de t
#3.1 Para el método de Euler

#Definición de condiciones iniciales
n = 400 # Número de particiones
deltat = 0.05 #Intervalo de tiempo entre iteraciones
U0 = array([1,0,0,1]) #Condiciones iniciales
matrizU = array(zeros([n,4])) #Matriz que contendrá en sus columnas los valores de las variables, cada fila es un paso de tiempo
F = array([0,1,-1,0])#Valor de F en el instante inicial
U = zeros(4)
x = zeros(n)
y = zeros(n)

#Aplicación del método de Euler
for i in range(n):
    for j in range(4):
         U[j] = U0[j]+ deltat*F[j] #Método de Euler
         continue
    matrizU[i,:] = U #Almacena en una columna de la matriz el valor de las variables para el paso de tiempo analizado
    U0 = U #Cambio las condiciones iniciales del paso
    d = ((U[0])**2+(U[1])**2)**1.5
    F = array([U[2], U[3], -U[0]/d, -U[1]/d]) #Obtengo mi nuevo valor de F(U,t)
    continue

#Las dos primeras columnas de la matriz son x e y de la órbita
x = matrizU[:,0]
y = matrizU[:,1]

#Representación gráfica
plt.plot(x,y)
plt.show()

#3.2 Para el método RK4
#Definición de condiciones iniciales
n = 50 # Número de particiones
deltat = 0.4 #Intervalo de tiempo entre iteraciones
U0 = array([1,0,0,1]) #Condiciones iniciales
matrizU = array(zeros([n,4])) #Matriz que contendrá en sus columnas los valores de las variables, cada fila es un paso de tiempo
F = array([0,1,-1,0])#Valor de F en el instante inicial
U = zeros(4)
x = zeros(n)
y = zeros(n)
k1 = zeros(4)
k2 = zeros(4)
k3 = zeros(4)
k4 = zeros(4)

#Aplicación del método RK4
for i in range(n):
    k1 = F #Primer coeficiente
    U2 = U0 +  k1*deltat/2
    d2 = ((U2[0])**2+(U2[1])**2)**1.5
    k2 = array([U2[2], U2[3], -U2[0]/d2, -U2[1]/d2]) #Segundo coeficiente
    U3 = U0 +  k2*deltat/2
    d3 = ((U3[0])**2+(U3[1])**2)**1.5
    k3 = array([U3[2], U3[3], -U3[0]/d3, -U2[1]/d3]) #Tercer coeficiente
    U4 = U0 +  k3*deltat
    d4 = ((U4[0])**2+(U4[1])**2)**1.5
    k4 = array([U4[2], U4[3], -U4[0]/d4, -U4[1]/d4]) #Cuarto coeficiente
    U = U0+ deltat*(k1+2*k2+2*k3+k4)/6 #Método RK4 
    matrizU[i,:] = U #Almacena en una columna de la matriz el valor de las variables para el paso de tiempo analizado
    U0 = U #Cambio las condiciones iniciales del paso
    d = ((U[0])**2+(U[1])**2)**1.5
    F = array([U[2], U[3], -U[0]/d, -U[1]/d]) #Obtengo mi nuevo valor de F(U,t)
    continue

#Las dos primeras columnas de la matriz son x e y de la órbita
x = matrizU[:,0]
y = matrizU[:,1]

#Representación gráfica
plt.plot(x,y)
plt.show()
