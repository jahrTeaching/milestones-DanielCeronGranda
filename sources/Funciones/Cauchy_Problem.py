from numpy import array, zeros #Importación de funciones numéricas
#Función para integrar un problema de Cauchy con un esquema arbitrario durante todo un tiempo n*deltat
def Cauchy_problem(Temporal_scheme, U0, n, deltat, F, F0, d): #d = dimensión de U
    t=0
    c = int(d)
    matrizU = array(zeros([n,c]))
    for i in range (n):
        t = t+deltat
        U = Temporal_scheme(U0, deltat, t, F0, F)
        matrizU[i,:] = U 
        U0 = U
        F0 = F(U,t) #Obtengo mi nuevo valor de F(U,t)
        continue
    return (matrizU)
