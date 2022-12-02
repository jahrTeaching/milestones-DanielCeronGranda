# Importacion de funciones de otros programas 
from Funciones.Cauchy_Problem import Cauchy_problem, Cauchy_problem_Leap_Frog
from Funciones.Esquemas_temporales import Euler, CN, RK4, Euler_inverso, Leap_Frog
from Funciones.Funciones_a_integrar import Linear_oscillator, Wrapped_Repr_Linear_Oscillator
from Funciones.Stability import Stability_Region, Wrapped_Stability_Zones
# IMPORTACIÓN
from numpy import array, zeros, linspace, transpose #Importación de funciones numéricas
from scipy.optimize import newton #Importación del método de Newton
import matplotlib.pyplot as plt #Importación de gráficas
import math

#APARTAD0 1: Integrar el oscilador lineal con varios esquemas
#Definición de condiciones iniciales
n = 2000 # Número de particiones
deltat = 0.1 #Intervalo de tiempo entre iteraciones
U0 = array([1,0]) #Condiciones iniciales
matrizU = array(zeros([n,4])) #Matriz que contendrá en sus columnas los valores de las variables, cada fila es un paso de tiempo
F0 = array([0,-1])#Valor de F en el instante inicial
U = zeros(2)
x = zeros(n)
Esquemas = [Euler, Euler_inverso, CN, RK4]

#Solución analítica
t = linspace(0, n*deltat, n)
x = array(zeros(len(t)))
for i in range (len(t)):
    x[i] = math.cos(t[i])
    continue
plt.plot(t,x)
plt.xlabel("t(s)")
plt.ylabel("x")
plt.title("Solución analítica del problema del oscilador lineal")
plt.grid()
plt.show()

for j in Esquemas:
    matrizU = Cauchy_problem(j, U0, n, deltat, Linear_oscillator, F0, 2)
    x = matrizU[:,0]
    Wrapped_Repr_Linear_Oscillator(t, x, j)
    continue

matriz_U1 = Cauchy_problem(Euler, U0, 10, deltat, Linear_oscillator, F0, 2)
U1 = matriz_U1[0,:]
matrizU = Cauchy_problem_Leap_Frog(Leap_Frog, U0, n, deltat, Linear_oscillator, U1, 2)
x = matrizU[:,0]
Wrapped_Repr_Linear_Oscillator(t, x, Leap_Frog)

#APARTADO 2: Regiones de estabilidad absoluta
Esquemas = [Euler, Euler_inverso, CN, RK4]
for i in Esquemas:
    rho, Re, Im  = Stability_Region(i, n/10, -4, 6, -5, 5)
    Wrapped_Stability_Zones(rho, Re, Im, i)
    continue

rho, Re, Im  = Stability_Region(Leap_Frog, n/10, -1, 1, -1, 1)
Wrapped_Stability_Zones(rho, Re, Im, Leap_Frog)