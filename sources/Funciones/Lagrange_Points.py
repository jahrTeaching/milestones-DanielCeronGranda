from numpy import array, zeros, matmul #Importación de funciones numéricas
from numpy.linalg import norm, eig
from scipy.optimize import fsolve #Importación del método de Newton

#Función de un paso de Runge-Kutta embebido: Dormand-Prince
def Emb_RK_DOPRI(U0, deltat, t, F0, F):
    Ord4 = RK_DOPRI("4", U0, t, deltat, F)
    Ord5 = RK_DOPRI("5", U0, t, deltat, F)
    
    (a, b, bs, c, q, Ne) = Butcher_matriz_DOPRI()
    
    tolerance = 1e-9
    h = min( deltat, deltatmax(Ord4 - Ord5, tolerance, min(q), deltat) )
    
    N = int( deltat / h ) + 1
    h = deltat / N
    
    for i in range(N):
        time = t + i * deltat / N
        Ord4 = Ord5
        Ord5 = RK_DOPRI("5", Ord4, time, h, F)
    
    U = Ord5
    
    return U

#Función de un paso del esquema Dormand Prince
def RK_DOPRI(Order, U0, t, deltat, F):
    
    (a, b, bs, c, q, Ne) = Butcher_matriz_DOPRI()
    N = len(U0)
    k = zeros( [Ne, N] )
    k[0,:] = F( U0, t + c[0]*deltat )
    
    if Order == "4":
        
        for i in range(1,Ne):
            Up = U0
            
            for j in range(i):
                Up = Up + deltat * a[i,j]*k[j,:]
                
            k[i,:] = F( Up, t + c[i]*deltat )
        
        U1 = U0 + deltat * matmul(b,k)
        
    elif Order == "5":
        
        for i in range(1,Ne):
            Up = U0
            
            for j in range(i):
                Up = Up + deltat * a[i,j]*k[j,:]
                
            k[i,:] = F(Up, t + c[i] * deltat)
            
        U1 = U0 + deltat * matmul(bs,k)
    
    return U1

#Definición de paso de tiempo utilizado
def deltatmax(dU, tolerance, q, dt):
    
    if( norm(dU) > tolerance ):
        deltat =  dt *(tolerance/norm(dU))**(1/(q+1)) 
    else:
        deltat = dt 
        
    return deltat

#Definición de la matriz de Butcher para el esquema Dormand-Prince
def Butcher_matriz_DOPRI():
    q = [5,4]
    Ne = 7 

    a = zeros( [Ne, Ne-1] )
    b = zeros( [Ne] )
    bs = zeros( [Ne] )
    c = zeros( [Ne] )
    
    c[:] = [ 0., 1./5, 3./10, 4./5, 8./9, 1., 1. ]

    a[0,:] = [          0.,           0.,           0.,         0.,           0.,     0. ]
    a[1,:] = [      1./5  ,           0.,           0.,         0.,           0.,     0. ]
    a[2,:]= [      3./40 ,        9./40,           0.,         0.,           0.,     0. ]
    a[3,:] = [     44./45 ,      -56./15,        32./9,         0.,           0.,     0. ]
    a[4,:] = [ 19372./6561, -25360./2187,  64448./6561,  -212./729,           0.,     0. ]
    a[5,:] = [  9017./3168,    -355./33 ,  46732./5247,    49./176, -5103./18656,     0. ]
    a[6,:]= [    35./384 ,           0.,    500./1113,   125./192, -2187./6784 , 11./84 ]

    b[:]  = [ 35./384   , 0.,   500./1113,  125./192,  -2187./6784  ,  11./84  ,     0.]
    bs[:] = [5179./57600, 0., 7571./16695,  393./640, -92097./339200, 187./2100, 1./40 ]
    
    return (a, b, bs, c, q, Ne)

def R3BP(U0, t, mu): 
    
    r = U0[:2] 
    v = U0[2:] 

    D1 = ( ( r[0] + mu )**2 + r[1]**2 )**(3/2)
    D2 = ( (r[0] - 1 + mu)**2 + r[1]**2 )**(3/2)

    dvdt_1 =  r[0]+2*v[1]-(1-mu)*(r[0]+mu)/D1-mu*(r[0]-1+mu)/D2
    dvdt_2 = r[1]-2*v[0]-(1-mu)*r[1]/D1-mu*r[1]/D2

    return array( [ v[0], v[1], dvdt_1, dvdt_2] )

def Lagrange_Points(U0, N_LP, mu):

    LP = zeros( [5,2] )

    def F(Y):
        X = zeros(4)
        X[:2] = Y
        X[2:] = 0
        Z = R3BP(X,0,mu)
        return Z[2:4]
   
    for i in range(N_LP):
        LP[i,:] = fsolve(F,U0[i,:2])
        
    return LP

def Jacobian(F, x0):
    N = len(x0)
    dx = 1e-9
    Jacobiano = zeros([N,N])
    for j in range(N):
        x = zeros(N)
        x[j] = dx
        Jacobiano[:,j] = ( F(x0 + x ) - F(x0 - x) ) /(2*dx)
    return Jacobiano

def stb_Lp(U0, mu):
    def F(X):
        return R3BP(X, 0 , mu)

    J = Jacobian(F, U0)
    Autovalor, autovector = eig(J)
    return Autovalor