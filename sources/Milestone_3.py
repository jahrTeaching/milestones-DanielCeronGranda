# Importación de funciones de otros programas 
from Funciones.Cauchy_Problem import Cauchy_problem
from Funciones.Esquemas_temporales import Euler, CN, RK4, Euler_inverso
from Funciones.Funciones_a_integrar import Kepler
# IMPORTACIÓN
from numpy import array, zeros, sqrt, log, abs, linspace, log10 #Importación de funciones numéricas
from scipy.optimize import newton #Importación del método de Newton
import matplotlib.pyplot as plt #Importación de gráficas

#Definición de condiciones iniciales
deltat = 0.1 #Intervalo de tiempo entre iteraciones
U0 = array([1,0,0,1]) #Condiciones iniciales
F0 = array([0,1,-1,0])#Valor de F en el instante inicial
T = 100 #Tiempo para el que se va a integrar el error 
Error =array(zeros([int(T/deltat)]))

# Definición de funciones

#Función para calcular el error de un esquema temporal para un tiempo y delta de t dados
def Error_esquema (Cauchy_problem, Temporal_scheme, U0, deltat, F, F0, t, q):
   N = len(t)
   C = len(U0)
   E = array(zeros([N,C]))
   matriz_Un = array(Cauchy_problem(Temporal_scheme, U0, N, deltat, F, F0, len(U0)))
   matriz_U2n = array(Cauchy_problem(Temporal_scheme, U0, 2*N, deltat/2, F, F0, len(U0)))
   for i in range (N):
       E[i-1,:] = (matriz_U2n[2*i-1,:]-matriz_Un[i-1,:])/(1-0.5**q)
       continue
   return E

#Función para calcular el orden de un esquema numérico para un tiempo T y dando un número de puntos 
#de la gráfica N-E llamado K. deltatmin y deltatmax son los deltat máximo y mínimo para los que se 
#calculará la gráfica N-E. Hay que darlos con cuidado, según el esquema que se esté usando. Darlos 
#incorrectamente puede provocar que falle la función.
def Orden_esquema(Cauchy_problem, Temporal_scheme, U0, F, F0, T, K, deltatmin, deltatmax):
   h = linspace(deltatmax, deltatmin, K)
   n = array(zeros([K]))
   C = len(U0)
   U = array(zeros([K,C]))
   norma_U = array(zeros([K]))
   norma_cuad_U = array(zeros([K]))
   for k in range(K):
       deltat = h[k]
       N = int(T/deltat)
       matriz_Un = array(Cauchy_problem(Temporal_scheme, U0, N, deltat, F, F0, len(U0)))
       matriz_U2n = array(Cauchy_problem(Temporal_scheme, U0, 2*N, deltat/2, F, F0, len(U0)))
       U[k,:] = matriz_U2n[2*N-1,:]-matriz_Un[N-1,:]
       n[k] = N
       for j in range(len(U0)):
               norma_cuad_U[k] = norma_cuad_U[k]+U[k,j-1]**2
               continue
       norma_U[k]=sqrt(norma_cuad_U[k])
       continue
   logU = log10(norma_U)
   logN = log10(n)
   vector_q = array(zeros([int(K/2)]))
   for i in range (int(K/2)):
       vector_q[i] = (logU[i+int(K/4)]-logU[i-10+int(K/4)])/(logN[i+int(K/4)]-logN[i-10+int(K/4)])
       continue
   q = 0
   for i in range(int(K/2)-1):
       q= q + vector_q[i]
       continue
   media_q = q/(K/2)
   Orden = (-1)*round(media_q)  
   print("El orden del esquema es", Orden)
   logE = logU - log10(1-0.5**Orden)
   plt.plot(logN, logE)
   plt.xlabel("log(N)")
   plt.ylabel("log(E)")
   plt.show()
   return Orden

Esquemas = [Euler, Euler_inverso, CN, RK4]
t = linspace(0, T, num=int(T/deltat))
for j in Esquemas:
    Orden = Orden_esquema(Cauchy_problem, j, U0, Kepler, F0, 1, 100, 0.01, 0.001)
    matriz_Error = Error_esquema (Cauchy_problem, j, U0, deltat, Kepler, F0, t, Orden)
    norma_Error = array(zeros([len(t)]))
    for i in range(int(T/deltat)):
        norma_Error[i] = sqrt(matriz_Error[i,0]**2+matriz_Error[i,1]**2)
        continue

    plt.plot(t, norma_Error)
    plt.xlabel("Tiempo de integración(s)")
    plt.ylabel("Error")
    plt.show()
    continue

