from numpy import array, zeros, sqrt, log, abs, linspace, log10 #Importación de funciones numéricas
from scipy.optimize import newton #Importación del método de Newton
import matplotlib.pyplot as plt #Importación de gráficas

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
   print(media_q)
   Orden = (-1)*round(media_q)  
   print("El orden del esquema es", Orden)
   logE = logU - log10(1-0.5**Orden)
   plt.plot(logN, logE)
   plt.xlabel("log(N)")
   plt.ylabel("log(E)")
   plt.show()
   return Orden
