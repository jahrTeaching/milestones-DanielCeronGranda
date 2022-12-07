#Funciones relacionadas con la estabilidad
# Importacion de funciones de otros programas 
# IMPORTACIÓN
from numpy import array, zeros, linspace, transpose #Importación de funciones numéricas
import matplotlib.pyplot as plt #Importación de gráficas

#Región de estabilidad
def Stability_Region(Esquema, N, x0, xf, y0, yf): 

    M = int(N)
    Re = linspace(x0, xf, M)
    Im = linspace(y0, yf, M)
    rho = array(zeros([M, M]))

    for i in range(M): 
      for j in range(M):

          w = complex(Re[i], Im[j])
          r = Esquema( 1., 1., 0., w,lambda u, t:w*u)
          rho[i, j] = abs(r) 

    return rho, Re, Im 

#Función de representación de regiones de estabilidad
def Wrapped_Stability_Zones(rho, Re, Im, Scheme):
    plt.contour( Re, Im, transpose(rho), linspace(0, 1, 10))
    plt.xlabel("Re(|r|)")
    plt.ylabel("Im(|r|)")
    j = "Estabilidad "+Scheme.__name__
    plt.title(j)
    plt.axis('equal')
    plt.grid()
    plt.show()
    return j