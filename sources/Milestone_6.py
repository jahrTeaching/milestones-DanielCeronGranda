# Importacion de funciones de otros programas 
from Funciones.Cauchy_Problem import Cauchy_problem, Cauchy_problem_Leap_Frog
from Funciones.Esquemas_temporales import Euler, CN, RK4, Euler_inverso, Leap_Frog
from Funciones.Funciones_a_integrar import Linear_oscillator, Wrapped_Repr_Linear_Oscillator
from Funciones.Stability import Stability_Region, Wrapped_Stability_Zones
from Funciones.Lagrange_Points import Emb_RK_DOPRI, RK_DOPRI, R3BP, deltatmax, Butcher_matriz_DOPRI, Lagrange_Points, Jacobian, stb_Lp
# IMPORTACIÓN
from numpy import array, zeros, linspace, transpose, matmul, sqrt, shape, around #Importación de funciones numéricas
from numpy.linalg import norm, eig
from scipy.optimize import newton, fsolve #Importación del método de Newton
import matplotlib.pyplot as plt #Importación de gráficas

#Condiciones iniciales
mu_earth_Moon = 1.22741e-2
N = 10000
t = linspace(0, 100, N)

# Condiciones iniciales para la obtención de los puntos de Lagrange
U0 = zeros( [5,4] ) 
U0[0,:] = [0.8, 0.6, 0 , 0]
U0[1,:] = [0.8, -0.6, 0, 0]
U0[2,:] = [-0.1, 0, 0, 0]
U0[3,:] = [0.1, 0, 0, 0]
U0[4,:] = [1.01, 0, 0, 0]

#Obtención de los puntos de Lagrange
L_p = Lagrange_Points(U0, 5, mu_earth_Moon)
print("\n" + str(L_p) + "\n")

#Estabilidad de los puntos de Lagrange
L_p_stb = zeros( 4 )

for i in range(5):
    L_p_stb[:2] = L_p[i, :]
    estab = around( stb_Lp(L_p_stb, mu_earth_Moon), 5)
    print(str(estab) + "\n")

#Órbitas alrededor de los puntos de Lagrange
U_0LPO = zeros( [shape(U0)[0], 4] )
eps = 1e-5 
for i in range( shape(U0)[0] ):
    U_0LPO[i, :2] = L_p[i, :] + eps
    U_0LPO[i, 2:] = eps

def F(U,t):
   return R3BP(U, t, mu_earth_Moon)

F0 = zeros(4)
Scheme = [ Euler, Euler_inverso, CN, RK4, Emb_RK_DOPRI]
for j in range (len(Scheme)):
    for i in range(shape(U0)[0]):    
        U_LP = Cauchy_problem(Scheme[j], U_0LPO[i,:], N, 100/N, F, F0, len(U_0LPO)) 
        LP = ["L4", "L5", "L3", "L1", "L2"]
    
        fig, ax = plt.subplots( figsize = (10, 10) )

        for k in range( len(L_p) ):
            ax.plot( L_p[k, 0],L_p[k, 1], 'o', label = LP[i] )
        
        ax.plot( U_LP[0,:], U_LP[1,:], color = 'b', label = "Orbit" )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("Puntos de Lagrange del sistema Tierra-Luna y órbita alrededor de" + LP[j])
        plt.legend()
        ax.grid()
        if LP[j] == "L4" or LP[j] == "L5": 
            fig, ay = plt.subplots( figsize = (7,7) )
            ay.plot( L_p[j, 0],L_p[j, 1], 'o', label = LP[j] )
            ay.plot( U_LP[0,:], U_LP[1,:], color = 'b', label = "Orbit"  )
            ay.set_xlabel("x")
            ay.set_ylabel("y")
            ay.set_title("Orbit arround" + LP[j])
            plt.legend()
            ay.grid()