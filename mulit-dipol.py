# -*- coding: utf-8 -*-
"""
Spyder Editor

Dies ist eine temporäre Skriptdatei.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.special import jve
from scipy.special import spherical_in

def D(r,u,s,w):
    x = (1+1j)*(r/np.sqrt((2/w*u*s)))
    x = np.sqrt((1j*w*u*s))*r
    D = ((2*u+1)*x-((2*u+1)+x**2)*np.tanh(x))/(((u-1)*x+((1-u)+x**2)*np.tanh(x)))
    return D

def H(rl,I,z):
    H = (I*rl**2)/(2*(rl**2+z**2)**(3/2))
    return H


def Dipol2(s,r,u,u2,w,R,Q,H,D):
    E = 1j*w*u2/2*H*r**3*R**(-2)*np.sin(Q)*D
    return E


def Multipol4(s,u,w,u2,I,Q,r,R):
    k = np.sqrt(1j*s*u*w)
    E_m = -1j*w*u2*I/2*(((np.sin(Q)/2*(r/R)**3*(-np.sin(Q))**2*(jve(1+3/2,1j*k*r)/jve(1-1/2,1j*k*r))))
            +(np.sin(Q)/6*(r/R)**5*(-3*np.sin(Q)*np.cos(Q))**2*(jve(2+3/2,1j*k*r)/jve(2-1/2,1j*k*r)))
            +np.sin(Q)/12*(r/R)**7*(3/2*(1-(np.cos(Q))**2)**(1/2)*(5*(np.cos(Q))**2-1))**2*(jve(3+3/2,1j*k*r)/jve(3-1/2,1j*k*r))
            +np.sin(Q)/20*(r/R)**9*(5/2*(1-(np.cos(Q))**2)**(1/2)*(7*(np.cos(Q))**3-3*np.cos(Q)))**2*(jve(4+3/2,1j*k*r)/jve(4-1/2,1j*k*r)))
    return E_m




def Comp_Phy_Ex01():
    r = 0.01
    rl = 0.10
    z = 0.15
    z = np.linspace(0.01,1,100)
    R = np.sqrt(rl**2+z**2)
    Q = np.arctan(rl/z)
    u = 1
    u2 = 1
    I = 1 #Verhältnis unabhängig von I
    s = 37e6
    w = 2*np.pi*20000
    U_rel = abs(Multipol4(s,u,w,u2,I,Q,r,R))/abs(Dipol2(s,r,u,u2,w,R,Q,H(rl,I,z),D(r,u,s,w)))
    plt.rcParams["figure.figsize"] = (9,10)
    plt.xlabel('Tiefe z [m]')
    plt.ylabel('$|U_{Multipol}| / |U_{Dipol}| $')
    plt.subplot(2, 1, 1)
    plt.plot(z,(U_rel), 'b',label='Aluminium: $r_K$ = 1 cm, r = 10 cm, $\mu_r$ = 1, $\sigma = 37~MS/m$')
    r = 0.03
    rl = 0.15
    z = 0.15
    z = np.linspace(0.01,1,100)
    R = np.sqrt(rl**2+z**2)
    Q = np.arctan(rl/z)
    u = 1
    u2 = 1
    I = 1 #Verhältnis unabhängig von I
    s = 37e6
    w = 2*np.pi*20000
    U_rel = abs(Multipol4(s,u,w,u2,I,Q,r,R))/abs(Dipol2(s,r,u,u2,w,R,Q,H(rl,I,z),D(r,u,s,w)))
    plt.plot(z,(U_rel), 'c',label='Aluminium: $r_K$ = 3 cm, $r$ = 20 cm, $\mu_r$ = 1, $\sigma = 37~MS/m$')
    r = 0.05
    rl = 0.20
    z = 0.15
    z = np.linspace(0.01,1,100)
    R = np.sqrt(rl**2+z**2)
    Q = np.arctan(rl/z)
    u = 1
    u2 = 1
    I = 1 #Verhältnis unabhängig von I
    s = 37e6
    w = 2*np.pi*20000
    U_rel = abs(Multipol4(s,u,w,u2,I,Q,r,R))/abs(Dipol2(s,r,u,u2,w,R,Q,H(rl,I,z),D(r,u,s,w)))
    #plt.plot(z,(U_rel), 'c',label='Aluminium: $r_K$ = 5 cm, $r$ = 20 cm, $\mu_r$ = 1, $\sigma = 37~MS/m$')
    r = 0.01
    rl = 0.10
    z = 0.15
    z = np.linspace(0.01,1,100)
    R = np.sqrt(rl**2+z**2)
    Q = np.arctan(rl/z)
    u = 4000
    u2 = 1
    I = 1 #Verhältnis unabhängig von I
    s = 8.6e6
    w = 2*np.pi*20000
    U_rel = abs(Multipol4(s,u,w,u2,I,Q,r,R))/abs(Dipol2(s,r,u,u2,w,R,Q,H(rl,I,z),D(r,u,s,w)))
    #plt.plot(z,(U_rel), 'k',label='Eisen: $r_K$ = 1 cm, r = 10 cm, $\mu_r$ = 4000, $\sigma = 8.6~MS/m$')
    r = 0.03
    rl = 0.20
    z = 0.15
    z = np.linspace(0.01,1,100)
    R = np.sqrt(rl**2+z**2)
    Q = np.arctan(rl/z)
    u = 4000
    u2 = 1
    I = 1 #Verhältnis unabhängig von I
    s = 8.6e6
    w = 2*np.pi*20000
    U_rel = abs(Multipol4(s,u,w,u2,I,Q,r,R))/abs(Dipol2(s,r,u,u2,w,R,Q,H(rl,I,z),D(r,u,s,w)))
    plt.plot(z,(U_rel), 'r',label='Eisen: $r_K$ = 3 cm, r = 20 cm, $\mu_r$ = 4000, $\sigma = 8.6~MS/m$')
    r = 0.01
    rl = 0.10
    z = 0.15
    z = np.linspace(0.01,1,100)
    R = np.sqrt(rl**2+z**2)
    Q = np.arctan(rl/z)
    u = 1
    u2 = 1
    I = 1 #Verhältnis unabhängig von I
    s = 10e6
    w = 2*np.pi*20000
    U_rel = abs(Multipol4(s,u,w,u2,I,Q,r,R))/abs(Dipol2(s,r,u,u2,w,R,Q,H(rl,I,z),D(r,u,s,w)))
    plt.subplot(2, 1, 2)
    plt.plot(z,(U_rel), 'b',label='$r_K$ = 1 cm, r = 10 cm, $\mu_r$ = 1, $\sigma = 10~MS/m$')
    r = 0.01
    rl = 0.10
    z = 0.15
    z = np.linspace(0.01,1,100)
    R = np.sqrt(rl**2+z**2)
    Q = np.arctan(rl/z)
    u = 10
    u2 = 1
    I = 1 #Verhältnis unabhängig von I
    s = 10e6
    w = 2*np.pi*20000
    U_rel = abs(Multipol4(s,u,w,u2,I,Q,r,R))/abs(Dipol2(s,r,u,u2,w,R,Q,H(rl,I,z),D(r,u,s,w)))
    plt.plot(z,(U_rel), 'r',label='$r_K$ = 1 cm, r = 10 cm, $\mu_r$ = 10, $\sigma = 10~MS/m$')
    r = 0.01
    rl = 0.1
    z = 0.15
    z = np.linspace(0.01,1,100)
    R = np.sqrt(rl**2+z**2)
    Q = np.arctan(rl/z)
    u = 100
    u2 = 1
    I = 1 #Verhältnis unabhängig von I
    s = 10e6
    w = 2*np.pi*20000
    U_rel = abs(Multipol4(s,u,w,u2,I,Q,r,R))/abs(Dipol2(s,r,u,u2,w,R,Q,H(rl,I,z),D(r,u,s,w)))
    plt.plot(z,(U_rel), 'g',label='$r_K$ = 1 cm, r = 10 cm, $\mu_r$ = 100, $\sigma = 10~MS/m$')
    r = 0.01
    rl = 0.10
    z = 0.15
    z = np.linspace(0.01,1,100)
    R = np.sqrt(rl**2+z**2)
    Q = np.arctan(rl/z)
    u = 100
    u2 = 1
    I = 1 #Verhältnis unabhängig von I
    s = 20e6
    w = 2*np.pi*20000
    U_rel = abs(Multipol4(s,u,w,u2,I,Q,r,R))/abs(Dipol2(s,r,u,u2,w,R,Q,H(rl,I,z),D(r,u,s,w)))
    plt.plot(z,(U_rel), 'k',label='$r_K$ = 1 cm, r = 10 cm, $\mu_r$ = 100, $\sigma = 20~MS/m$')
    plt.legend()
    plt.show
    #plt.savefig('1.png')
Comp_Phy_Ex01()