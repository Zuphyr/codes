# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 10:49:09 2023

@author: zuphyr
"""


from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import math
import time
from scipy import *
from multiprocessing import Pool
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


import matplotlib as mpl
ti=time.time() 


"QFI for non-hermitian"


# t=np.linspace(0,3000, 1000)

# # F=(np.arcsin(t))**2/(1-t**2)

# # G=(np.arcsin(t))**2

# Delta_0=0.005
# Delta_1=0.007
# Delta_2=0.01
# g=10
# kappa=10

# I_0=(np.sin(Delta_0*t)-Delta_0*t)**2/(Delta_0*Delta_0*Delta_0*Delta_0*Delta_0*Delta_0*100000000000000)
# I_1=(np.sin(Delta_1*t)-Delta_1*t)**2/(Delta_1*Delta_1*Delta_1*Delta_1*Delta_1*Delta_1*100000000000000)
# I_2=(np.sin(Delta_2*t)-Delta_2*t)**2/(Delta_2*Delta_2*Delta_2*Delta_2*Delta_2*Delta_2*100000000000000)

# plt.figure()

# # plt.plot(r, np.log(G))

# # plt.plot(r, np.log(F))

# plt.plot(t, (I_0),'r')
# plt.plot(t, (I_1),'b')
# plt.plot(t, (I_2),'g')


# plt.rcParams['savefig.dpi'] = 600 #图片像素
# plt.rcParams['figure.dpi'] = 600 #分辨率
# plt.rcParams['figure.figsize']=(3.26,2.36)
# # plt.rcParams['figure.figsize']=(1.63,1.18)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)

# plt.show()







"QFI for COMag"

"parameters"

omega_m=2*np.pi*10
Delta_BF=omega_m*1
Delta_h=0
Delta_v=Delta_BF
Delta_m=0

kappa_h_in=omega_m/40
kappa_h_ex=1*kappa_h_in #critical coupling
kappa_h=kappa_h_in+kappa_h_ex

kappa_v_in=1*kappa_h_in  #inequal H-V decay
kappa_v_ex=kappa_h_ex
kappa_v=kappa_v_in+kappa_v_ex

# kappa_v=1*kappa_h  #equal H-V decay
# kappa_v_ex=kappa_h_ex

kappa_m=2.74*omega_m/1000

m_c=Delta_BF/(2*np.sqrt(6*2.6))
ratio_0=0.999
ratio_1=0.9
ratio_2=0.8
ratio_3=1.1

m=ratio_3*m_c

g_h=2*np.pi*6*m
g_v=-2*np.pi*2.6*m

# theta_p=np.pi*0
# eta=kappa_h*0.1   # eta=kappa_ex 是否需要修改？
# eta_h=kappa_h*0.5
# eta_h=kappa_h*0.01*cos(theta_p)
# eta_v=kappa_h*0.01*sin(theta_p)
# pm_in=kappa_m*0

eta=0*kappa_h_ex


"operators and states"
N_o=10

a=destroy(N_o)
I_o=qeye(N_o)

a_h=tensor(a,I_o)
a_v=tensor(I_o,a)


psi_i=(tensor(basis(N_o,2),basis(N_o,0))+tensor(basis(N_o,0),basis(N_o,2)))/np.sqrt(2)
rho_i=psi_i*psi_i.trans()


i=0
t=np.linspace(0, 0.01, 1000)
I_3=np.zeros(1000)

def QFI(dt):
    U=-1j*dt*(Delta_h*a_h.dag()*a_h+Delta_v*a_v.dag()*a_v+g_h*a_h.dag()*a_v+g_v*a_v.dag()*a_h+eta*(a_h.dag()+a_h))
    U_t=U.expm()
    f=U_t*(a_h.dag()*a_v+a_h*a_v.dag())*U_t.dag()
    
    I_3=(f*f*rho_i).tr()-((f*rho_i).tr())**2

    return I_3

while i<1000:
    dt=t[i]*kappa_h/np.pi
    I_3[i]=real(QFI(dt))
    i=i+1

plt.figure()

plt.plot(t, (I_0), 'r')
plt.plot(t, (I_1), 'b')
plt.plot(t, (I_2), 'g')
plt.plot(t, (I_3), 'm')


plt.rcParams['savefig.dpi'] = 600 #图片像素
plt.rcParams['figure.dpi'] = 600 #分辨率
plt.rcParams['figure.figsize']=(3.26,2.36)
# plt.rcParams['figure.figsize']=(1.63,1.18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.show()



tf=time.time() 
print('cost',tf-ti,'s')