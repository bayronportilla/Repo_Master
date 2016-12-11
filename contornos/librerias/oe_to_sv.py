#!/usr/bin/env python

############################################################
# De elementos orbitales a vector de estado. Dados los elementos
# orbitales de un objeto para uns instante de tiempo determinado
# calcula las componenetes cartesianas en el plano de referencia
# de la ecliptica en forma de array [x,y,z]. Para llamar a la funcion
# usar: elements(a=xxx , e=xxx, i=xxx, Omega=xxx, omega=xxx, M=xxx)
# donde todos los angulos deben introducirse en grados y las
# distancias en unidades canonicas de longitud. El orden en la
# entrada de los argumentos es irrelevante PERO DEBEN SER
# ESPECIFICAMENTE LOS QUE ESTAN AHI. 

# @bayronportilla 2016


import numpy as np
from numpy import cos, sin
from scipy.optimize import newton



InDeg = 180.0 / np.pi
InRad = np.pi / 180.0
G = 1.0


############################################################
# Calculo de la posicion en el plano fundamental

def pos_ref(**kwargs):
    for key in kwargs:
        if( key == 'm_1' ):
            m_1 = kwargs[key]
        elif( key == 'm_2' ):
            m_2 = kwargs[key]
        elif( key == 'a' ):
            a = kwargs[key]
        elif( key == 'e' ):
            e = kwargs[key]
        elif( key == 'inc'):
            inc = kwargs[key]*InRad
        elif( key == 'Omega'):
            Omega = kwargs[key]*InRad
        elif( key == 'omega'):
            omega = kwargs[key]*InRad
        elif( key == 'M'):
            M = kwargs[key]*InRad
            
    if ( inc == 0.0 ):
        Omega = 0.0

    if ( e == 0.0 ):
        omega = 0.0
    
            
    ############################################################
    # Inicia calculo de cantidades derivadas
    
    def kepler(E):
        return E-e*np.sin(E) - M 
    E = newton(kepler,M)

    # True anomaly
    f = 2*np.arctan( np.sqrt((1.0 + e) / (1.0 - e)) * np.tan(E/2.0) )

    # Radial distance
    r = a*(1.0-e**2) / (1.0 + e*np.cos(f))


    # Posicion en plano orbital
    
    r_orb = [ r*np.cos(f) , r*np.sin(f) , 0.0 ]    
    
        
            
    matrix=np.matrix([ [cos(Omega)*cos(omega) - sin(Omega)*sin(omega)*cos(inc), -cos(Omega)*sin(omega) - sin(Omega)*cos(inc)*cos(omega), sin(inc)*sin(Omega)], \
                       [sin(Omega)*cos(omega) + cos(Omega)*cos(inc)*sin(omega), -sin(Omega)*sin(omega) + cos(Omega)*cos(inc)*cos(omega), -cos(Omega)*sin(inc)], \
                       [sin(inc)*sin(omega), sin(inc)*cos(omega), cos(inc)] ])

    
    r_ref = np.array( np.dot(matrix , r_orb) )
    r_ref = [r_ref[0][0] , r_ref[0][1] , r_ref[0][2]]

    
    return r_ref





############################################################
# Calculo de la velocidad en el plano fundamental

def vel_ref(**kwargs):
    for key in kwargs:
        if( key == 'm_1' ):
            m_1 = kwargs[key]
        elif( key == 'm_2' ):
            m_2 = kwargs[key]
        elif( key == 'a' ):
            a = kwargs[key]
        elif( key == 'e' ):
            e = kwargs[key]
        elif( key == 'inc'):
            inc = kwargs[key]*InRad
        elif( key == 'Omega'):
            Omega = kwargs[key]*InRad
        elif( key == 'omega'):
            omega = kwargs[key]*InRad
        elif( key == 'M'):
            M = kwargs[key]*InRad

    if ( inc == 0.0 ):
        Omega = 0.0

    if ( e == 0.0 ):
        omega = 0.0

    ############################################################
    # Inicia calculo de cantidades derivadas
    
    def kepler(E):
        return E-e*np.sin(E) - M 
    E = newton(kepler,M)

    # True anomaly
    f = 2*np.arctan( np.sqrt((1.0 + e) / (1.0 - e)) * np.tan(E/2.0) )

    # Radial distance
    r = a*(1.0-e**2) / (1.0 + e*np.cos(f))

    
    coef = ( G*(m_1 + m_2) / (a*(1.0 - e**2)) )**0.5
    v_orb = [-coef * np.sin(f) , coef*(e + np.cos(f)) , 0.0  ]

            
    
    matrix=np.matrix([ [cos(Omega)*cos(omega) - sin(Omega)*sin(omega)*cos(inc), -cos(Omega)*sin(omega) - sin(Omega)*cos(inc)*cos(omega), sin(inc)*sin(Omega)], \
                       [sin(Omega)*cos(omega) + cos(Omega)*cos(inc)*sin(omega), -sin(Omega)*sin(omega) + cos(Omega)*cos(inc)*cos(omega), -cos(Omega)*sin(inc)], \
                       [sin(inc)*sin(omega), sin(inc)*cos(omega), cos(inc)] ])

    
    v_ref = np.array( np.dot(matrix , v_orb) )
    v_ref = [v_ref[0][0] , v_ref[0][1] , v_ref[0][2]]
    
    
    return v_ref










    
    

    

    




