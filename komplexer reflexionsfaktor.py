# -*- coding: utf-8 -*-
"""
Spyder Editor

Dies ist eine temporäre Skriptdatei.
"""

import numpy as np
import matplotlib.pyplot as plt


def D_Stahlreal(r,u,s,w):
    x = (1+1j)*(r/np.sqrt((2/w*u*s)))
    D = ((2*u+1)*x-((2*u+1)+x**2)*np.tanh(x))/(((u-1)*x+((1-u)+x**2)*np.tanh(x)))
    return D

def D_Stahlrealverbesserung(r,u,s,w):
    x = np.sqrt((1j*w*u*s))*r
    D = ((2*u+1)*x-((2*u+1)+x**2)*np.tanh(x))/(((u-1)*x+((1-u)+x**2)*np.tanh(x)))
    return D


def komplRef():
    r = np.logspace(-7,-1,10000)
    u1 = 10
    s1 = 10e6
    w = 20000
    ######
    #verbesserung
    D0_real = D_Stahlrealverbesserung(r,u1,s1,w).real
    D0_imag = D_Stahlrealverbesserung(r,u1,s1,w).imag
    delta0 = r/np.sqrt(2/(w*u1*s1))
    ######
    D1 = []
    D1i = []
    for i in range(len(r)):
        D1.append(D_Stahlreal(r[i],u1,s1,w).real)
    for i in range(len(r)):
        D1i.append(D_Stahlreal(r[i],u1,s1,w).imag)
    u2 = 100
    s1 = 10e6
    w = 20000    
    ######
    #verbesserung
    D1_real = D_Stahlrealverbesserung(r,u2,s1,w).real
    D1_imag = D_Stahlrealverbesserung(r,u2,s1,w).imag
    delta1 = r/np.sqrt(2/(w*u2*s1))
    ######
    D2 = []
    D2i = []
    for i in range(len(r)):
        D2.append(D_Stahlreal(r[i],u2,s1,w).real)
    for i in range(len(r)):
        D2i.append(D_Stahlreal(r[i],u2,s1,w).imag)
    u3 = 1000
    s1 = 10e6
    w = 20000
    ######
    #verbesserung
    D2_real = D_Stahlrealverbesserung(r,u3,s1,w).real
    D2_imag = D_Stahlrealverbesserung(r,u3,s1,w).imag
    delta2 = r/np.sqrt(2/(w*u3*s1))
    ######
    D3 = []
    D3i = []
    for i in range(len(r)):
        D3.append(D_Stahlreal(r[i],u3,s1,w).real)
    for i in range(len(r)):
        D3i.append(D_Stahlreal(r[i],u3,s1,w).imag)
    u4 = 1
    s2 = 37e6
    w = 20000
    ######
    #verbesserung
    D3_real = D_Stahlrealverbesserung(r,u4,s2,w).real
    D3_imag = D_Stahlrealverbesserung(r,u4,s2,w).imag
    delta3 = r/np.sqrt(2/(w*u4*s2))
    ######
    D4 = []
    D4i = []
    for i in range(len(r)):
        D4.append(D_Stahlreal(r[i],u4,s2,w).real)
    for i in range(len(r)):
        D4i.append(D_Stahlreal(r[i],u4,s2,w).imag)
    u5 = 4000
    s3 = 8.6e6
    w = 20000
    ######
    #verbesserung
    D4_real = D_Stahlrealverbesserung(r,u5,s3,w).real
    D4_imag = D_Stahlrealverbesserung(r,u5,s3,w).imag
    delta4 = r/np.sqrt(2/(w*u5*s3))
    ######
    u6 = 1
    s4 = 58e6
    w = 20000
    ######
    #verbesserung
    D5_real = D_Stahlrealverbesserung(r,u6,s4,w).real
    D5_imag = D_Stahlrealverbesserung(r,u6,s4,w).imag
    delta5 = r/np.sqrt(2/(w*u6*s4))
    ######
    plt.rcParams["figure.figsize"] = (9,5)
    plt.xlabel('Induktionszahl $r_K / \delta$')
    plt.ylabel('Reflexionsfaktor D')
    ########
    plt.plot(delta0,D0_real,'b-',label='$\mu_r$ = 10, $\sigma = 10~MS/m$')
    plt.plot(delta0,D0_imag,'b--')
    plt.plot(delta1,D1_real,'g-',label='$\mu_r$ = 100, $\sigma = 10~MS/m$')
    plt.plot(delta1,D1_imag,'g--')
    plt.plot(delta2,D2_real,'r-',label='$\mu_r$ = 1000, $\sigma = 10~MS/m$')
    plt.plot(delta2,D2_imag,'r--')
    plt.plot(delta3,D3_real,'c-',label='Aluminium: $\mu_r$ = 1, $\sigma = 37~MS/m$')
    plt.plot(delta3,D3_imag,'c--')
    plt.plot(delta4,D4_real,'k-',label='Eisen: $\mu_r$ = 4000, $\sigma = 8.6~MS/m$')
    plt.plot(delta4,D4_imag,'k--')
    #plt.plot(delta5,D5_real,'r-',label='Kupfer: $\mu$ = 1, $\sigma = 58~MS/m$')
    #plt.plot(delta5,D5_imag,'r--')
    #######
    #plt.plot(r,D1,'b-',label='$\mu$ = 10, $\sigma = 8,6~MS/m$')
    #plt.plot(r,D1i,'b--')
    #plt.plot(r,D2,'g-',label='$\mu$ = 100, $\sigma = 8,6~MS/m$')
    #plt.plot(r,D2i,'g--')
    #plt.plot(r,D3,'y-',label='$\mu$ = 1000, $\sigma = 8,6~MS/m$')
    #plt.plot(r,D3i,'y--')
    #plt.plot(r,D4,'r-',label='$\mu$ = 1, $\sigma = 37~MS/m$')
    #plt.plot(r,D4i,'r--')
    plt.xscale('log')
    plt.legend()
    plt.show
    #plt.savefig('D.png')
komplRef()