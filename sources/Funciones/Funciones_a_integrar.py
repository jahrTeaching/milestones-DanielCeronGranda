from numpy import array, zeros #Importación de funciones numéricas
#Ecuaciones de la órbita de Kepler
def Kepler(U,t):
    d = ((U[0])**2+(U[1])**2)**1.5 
    x = U[2]
    y = U[3]
    vx = -U[0]/d
    vy = -U[1]/d
    return array([x,y,vx,vy])