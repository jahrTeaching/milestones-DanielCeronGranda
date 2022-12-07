# Importación de funciones de otros programas 
from Funciones.Cauchy_Problem import Cauchy_problem
from Funciones.Esquemas_temporales import  RK4
# IMPORTACIÓN
from numpy import array, zeros, reshape #Importación de funciones numéricas
from numpy.linalg import norm
import matplotlib.pyplot as plt #Importación de gráficas

#Función del problema de los N cuerpos
def F_NBody(U,t): 
  Nb = 4
  Nc = 2  
    
  Us = reshape( U, (Nb, Nc, 2) )      
  F =  zeros(len(U))   
  Fs = reshape(F, (Nb, Nc, 2))  

  r = reshape(Us[:, :, 0], (Nb, Nc))     
  v = reshape(Us[:, :, 1], (Nb, Nc)) 
 
  drdt = reshape(Fs[:, :, 0], (Nb, Nc)) 
  dvdt = reshape(Fs[:, :, 1], (Nb, Nc))
  dvdt[:,:] = 0   

  for i in range(Nb):   
     drdt[i,:] = v[i,:]
     for j in range(Nb): 
         if j != i:  
           d = r[j,:] - r[i,:]
           dvdt[i,:] = dvdt[i,:] + d[:]/norm(d)**3 
   
  return F

#Función de integración del problema de los N-cuerpos
def Int_N_Body(Nb, Nc, N, t0, tf, U0, Scheme):
  deltat = (tf-t0)/(N+1)
  F0 = array(zeros(Nc*2*Nb))
  U = Cauchy_problem(Scheme, U0, N, deltat, F_NBody, F0, 2*Nb*Nc)
  Us  = reshape(U, (N, Nb, Nc, 2)) 
  r   = reshape(Us[:, :, :, 0], (N, Nb, Nc)) 
  #Meter gráficos, títulos... 
  for i in range(Nb):
    plt.plot(r[:, i, 0], r[:, i, 1])
    
  plt.axis('equal')
  plt.grid()
  plt.show()
  return U

#Simulación de ejemplo del problema de los N cuerpos
Nb = 4 #Define el número de cuerpos
Nc = 2 #Define el número de coordenadas
t0 = 0
tf = 3 #Tiempo final de simulación

#Condiciones iniciales (se pueden variar)
U0 = array(zeros([Nc*2*Nb]))
U0 = array([2,0,2,0, 1,-1,3,1, 1,-1,0,-2, 0,-1,0,1]) #Para cada cuerpo el orden es (x0, vx0, y0, vy0)

matrizU = Int_N_Body(Nb, Nc, 100*tf, t0, tf, U0, RK4)