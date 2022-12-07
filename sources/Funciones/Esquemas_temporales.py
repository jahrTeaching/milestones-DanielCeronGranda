from scipy.optimize import newton #Importación del método de Newton

#Función para integrar un paso de Euler
def Euler(U0, deltat, t, F0, F):
    U = U0 + deltat*F0
    return U

#Función para integrar un paso de Crank-Nicholson
def CN(U0, deltat, t, F0, F):
    def Res(X):
        return X-U0-deltat*0.5*(F0+ F(X,t))
    return newton(Res, U0)

#Función para integrar un paso de Euler inverso
def Euler_inverso(U0, deltat, t, F0, F):
    def Res2(Y):
        return Y-U0-deltat*F(Y,t)
    return newton(Res2, U0,tol = 1e-1)

#Función para integrar un paso de RK4
def RK4(U0, deltat, t, F0, F):
    k1 = F0 #Primer coeficiente
    U2 = U0 +  k1*deltat/2
    k2 = F(U2,t) #Segundo coeficiente
    U3 = U0 +  k2*deltat/2
    k3 = F(U3,t) #Tercer coeficiente
    U4 = U0 +  k3*deltat
    k4 = F(U4,t) #Cuarto coeficiente
    U = U0+ deltat*(k1+2*k2+2*k3+k4)/6 #Método RK4 
    return U

#Función para integrar un paso de Leap Frog
#Cuidado al usarla, aquí hay que meter dos condiciones de U en vez de
#una condición de U y otra de F
def Leap_Frog(U0, deltat, t, U1, F):
    F1 = F(U1, t)
    U2 = U0 + 2*deltat*F1
    return U2
