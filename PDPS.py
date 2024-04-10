# -*- coding: utf-8 -*-
"""
Created on Tue May 17 16:07:16 2022

@author: zuphyr
"""


from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import math
import time
from scipy import *
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

ti=time.time() 


"parameters"

omega_m=2*np.pi*1
Delta_h=0
Delta_v=0
Delta_m=0

kappa_h_in=omega_m/40
kappa_h_ex=1*kappa_h_in #critical coupling
kappa_h=kappa_h_in+kappa_h_ex

kappa_v_in=2*kappa_h_in  #inequal H-V decay
kappa_v_ex=kappa_h_ex
kappa_v=kappa_v_in+kappa_v_ex

# kappa_v=1*kappa_h  #equal H-V decay
# kappa_v_ex=kappa_h_ex

kappa_m=2.74*omega_m/1000
g0=0.5*kappa_h




# theta_p=np.pi*0
# eta=kappa_h*0.1   
# eta_h=kappa_h*0.5
# eta_h=kappa_h*0.01*cos(theta_p)
# eta_v=kappa_h*0.01*sin(theta_p)
# pm_in=kappa_m*0
pm_in=kappa_m*1

"operators and states"
N_o=4
N_m=8

a=destroy(N_o)
b=destroy(N_m)
I_o=qeye(N_o)
I_m=qeye(N_m)

a_h=tensor(a,I_o,I_m)
a_v=tensor(I_o,a,I_m)
bb=tensor(I_o,I_o,b)


# psi_i=tensor(basis(N_o,0),basis(N_o,0),basis(N_m,0))
# rho_i=psi_i*psi_i.trans()

n_m=1

# b_h_in=eta_h/sqrt(kappa_h_ex)
# b_v_in=eta_v/sqrt(kappa_v_ex)

# bh=b_h_in*I_o-np.sqrt(kappa_h_ex)*a
# bv=b_v_in*I_o-np.sqrt(kappa_v_ex)*a

# b_h=b_h_in-1j*sqrt(kappa_h_ex)*a_h
# b_v=b_v_in-1j*sqrt(kappa_v_ex)*a_v

# br=(bh+sqrt(2)*(0.9817+0.1739*1j)*bv)/sqrt(2)
# bl=(bh-1j*bv)/sqrt(2)

# b_h=tensor(bh,I_o,I_m)
# b_v=tensor(I_o,bv,I_m)

# b_r=tensor(br,I_o,I_m)
# b_l=tensor(I_o,bl,I_m)


"output photons blockade with fixed different incident angle"

"HH"
# def output_g20_for_HH(delta):
#     Delta_m=0
#     Delta_h=delta
#     Delta_v=delta
    
#     angle=np.pi/2
#     th=np.cos(angle)
#     tv=np.sin(angle)
    
#     H_free=-Delta_h*a_h.dag()*a_h-Delta_v*a_v.dag()*a_v+Delta_m*bb.dag()*bb
#     H_i=-g0*a_h.dag()*a_v*(bb.dag())-g0*a_v.dag()*a_h*(bb)

#     # H_p=eta_h*a_h.dag()+conjugate(eta_h)*a_h+eta_v*a_v.dag()+conjugate(eta_v)*a_v #HV probe
#     H_p=eta*tv*a_v.dag()+conjugate(eta)*tv*a_v+eta*th*a_h.dag()+conjugate(eta)*th*a_h+pm_in*(bb.dag()+bb) #V probe
#     H=H_free+H_i+H_p

#     c_list=[]
#     c_list.append(sqrt(kappa_h)*a_h)
#     c_list.append(sqrt(kappa_v)*a_v)
#     c_list.append(sqrt(kappa_m*(n_m+1))*bb)
#     c_list.append(sqrt(kappa_m*(n_m))*bb.dag())
    
#     rho_ss=steadystate(H, c_list)
    
#     p_h_n=expect(a_h.dag()*a_h,rho_ss)
#     # p_v_n=expect(a_v.dag()*a_v,rho_ss)
    
#     # p_h_nn=expect(a_h.dag()*a_h.dag()*a_h*a_h,rho_ss)
#     # p_v_nn=expect(a_v.dag()*a_v.dag()*a_v*a_v,rho_ss)

#     # p_h_plus=expect((a_h+a_h.dag()), rho_ss)
#     # p_v_plus=expect((a_v+a_v.dag()), rho_ss)
#     p_h_minus=expect((a_h.dag()-a_h),rho_ss)
#     # p_v_minus=expect((a_v.dag()-a_v),rho_ss)
#     # # p_cross=expect((a_h.dag()*a_v+a_v.dag()*a_h),rho_ss)
#     # p_cross_minus=expect((a_h.dag()*a_v-a_h*a_v.dag()), rho_ss)

#     I_h=th*th+1j*th*p_h_minus+p_h_n
#     # I_v=tv*tv+1j*tv*p_v_minus+p_v_n
    
#     p_1op=expect((a_h.dag()-a_h), rho_ss)
#     p_2op=expect((-a_h.dag()*a_h.dag()+4*a_h.dag()*a_h-a_h*a_h),rho_ss)
#     p_3op=expect((a_h.dag()*(a_h.dag()-a_h)*a_h), rho_ss)
#     p_4op=expect(a_h.dag()*a_h.dag()*a_h*a_h,rho_ss)
    
    
#     # g20_hh=p_h_nn/(p_h_n*p_h_n)
#     # g20_v=p_v_nn/(p_v_n*p_v_n)
#     g20_hh_output=(th*th*th*th+2j*th*th*th*p_1op+th*th*p_2op+2j*th*p_3op+p_4op)/(I_h*I_h)

#     return g20_hh_output

# output_g20_for_HH(10*g0)

# i=0
# detune=5*g0
# detuning=np.zeros(79)
# g20_HH_out=np.zeros(79)

# while i<=39:
#     detuning[i]=-detune
#     detuning[78-i]=detune
#     g20_HH_out[i]=output_g20_for_HH(detuning[i])
#     g20_HH_out[78-i]=output_g20_for_HH(detuning[78-i])
#     print(i)
#     i=i+1
#     # detune=detune-0.02*kappa_c
#     detune=detune/1.2
#     while i==39:
#         detuning[i]=0
#         g20_HH_out[i]=output_g20_for_HH(detuning[i])
#         print('finish')
#         break
    
# fig=plt.figure()
# # plt.plot(detuning)

# plt.plot(detuning/g0,g20_HH_out,'b')
# # plt.plot(detuning/g0,np.log(g20_HH_out),'b')
# # plt.plot(-1*Delta/omega_m,S1,'b')
# # plt.xlabel('$\Delta /\omega _m$')
# # plt.xlim(-2, 2)
# # plt.ylim(-0.5,2)
# plt.title('$g^2_{HH,out}(0)$ with $\pi /2$ incident')
# plt.xlabel('$\Delta /g0 $')
# plt.rcParams['savefig.dpi'] = 600 
# plt.rcParams['figure.dpi'] = 600 
# plt.rcParams['figure.figsize']=(3.26,2.36)
# plt.show()


"---------------------------versus Delta---------------------------------"


# "VV"

# def output_g20_for_HH(delta):
#        Delta_m=0
#        Delta_h=delta
#        Delta_v=delta
      
#        "i-o angles"
#        theta_in=np.pi/2
#        theta_out=0
#        # theta_o=np.pi*105/180
      

#        "decay rates"
#        kappa_h_ex=1*kappa_h_in
#        # kappa_h_ex=kappa_ex
#        kappa_h=kappa_h_in+kappa_h_ex

#        kappa_v_ex=kappa_h_ex
#        kappa_v=kappa_v_in+kappa_v_ex
      
#        "input amplitudes"
#        eta_h=kappa_h_in*0.6*np.cos(theta_in)
#        eta_v=kappa_h_in*0.6*np.sin(theta_in)
      
#        "output operators"
#        b_h_in=eta_h/sqrt(kappa_h_ex)
#        b_v_in=eta_v/sqrt(kappa_v_ex)

#        b_h=b_h_in-1j*sqrt(kappa_h_ex)*a_h
#        b_v=b_v_in-1j*sqrt(kappa_v_ex)*a_v

#        H_free=-Delta_h*a_h.dag()*a_h-Delta_v*a_v.dag()*a_v-Delta_m*bb.dag()*bb
#        H_i=+g0*a_h.dag()*a_v*(bb.dag())+g0*a_v.dag()*a_h*(bb)

#        # H_p=eta_h*a_h.dag()+conjugate(eta_h)*a_h+eta_v*a_v.dag()+conjugate(eta_v)*a_v #HV probe
#        H_p=eta_v*a_v.dag()+conjugate(eta_v)*a_v+eta_h*a_h.dag()+conjugate(eta_h)*a_h+pm_in*(bb.dag()+bb) #V probe
#        H=H_free+H_i+H_p

#        c_list=[]
#        c_list.append(sqrt(kappa_h)*a_h)
#        c_list.append(sqrt(kappa_v)*a_v)
#        c_list.append(sqrt(kappa_m*(n_m+1))*bb)
#        c_list.append(sqrt(kappa_m*(n_m))*bb.dag())
      
#        rho_ss=steadystate(H, c_list)
      
#        # rho_out=tensor(np.cos(theta_o)*I_o*bh,np.sin(theta_o)*I_o*bv,I_m)*rho_ss*tensor(np.cos(theta_o)*I_o*bh,np.sin(theta_o)*I_o*bv,I_m).dag()
      
#        # b_out_h=(np.cos(theta_o)**2*I_o)*a+(0.5*np.sin(2*theta_o)*I_o)*a
#        # b_out_v=(0.5*np.sin(2*theta_o)*I_o)*a+(np.sin(theta_o)**2*I_o)*a
#        # b_h=tensor(b_out_h,I_o,I_m)
#        # b_v=tensor(I_o,b_out_v,I_m)
#        b_out=tensor(np.cos(theta_out)*I_o,I_o,I_m)*b_h+tensor(I_o,np.sin(theta_out)*I_o,I_m)*b_v
#        # a_out=tensor(np.cos(theta_out)*I_o,I_o,I_m)*a_h+tensor(I_o,np.sin(theta_out)*I_o,I_m)*a_v

#        # p_h_n=expect(b_h.dag()*b_h,rho_ss)
#        # p_v_n=expect(b_v.dag()*b_v,rho_ss)
      
#        # p_int_n=expect((b_h.dag()*b_v+b_h*b_v.dag()),rho_ss)
      
#        # p_v_nn=expect(b_v.dag()*b_v.dag()*b_v*b_v,rho_ss)
#        # p_h_nn=expect(b_h.dag()*b_h.dag()*b_h*b_h,rho_ss)
      
      
#        # p_hv_nn=expect(b_h.dag()*b_h.dag()*b_v*b_v,rho_ss)
#        # p_vh_nn=expect(b_v.dag()*b_v.dag()*b_h*b_h,rho_ss)
#        # p_3hv_nn=expect((b_h.dag()*b_h.dag()*b_h*b_v+b_h.dag()*b_h*b_h*b_v.dag()),rho_ss)
#        # p_h3v_nn=expect((b_h*b_v.dag()*b_v.dag()*b_v+b_h.dag()*b_v.dag()*b_v*b_v),rho_ss)


      
#        G_out=expect(b_out.dag()*b_out.dag()*b_out*b_out,rho_ss)
#        I_out=expect(b_out.dag()*b_out,rho_ss)
#        # I=expect(a_out.dag()*a_out,rho_ss)

#        g20_out=(G_out)/(I_out*I_out)
      
#        return g20_out



# # output_g20_for_VV(0*g0)

# i=0
# detune=1*omega_m
# g20_HH_out_4=np.zeros(201)
# detuning=np.zeros(201)

# while i<=99:
#     detuning[i]=-detune
#     detuning[200-i]=detune
#     g20_HH_out_4[i]=np.real(output_g20_for_HH(detuning[i]))
#     g20_HH_out_4[200-i]=np.real(output_g20_for_HH(detuning[200-i]))
#     print(i)
#     i=i+1
#     # detune=detune-0.2*g0
#     detune=detune-0.01*omega_m
#     while i==100:
#         detuning[i]=0
#         g20_HH_out_4[i]=np.real(output_g20_for_HH(detuning[i]))
#         print('finish')
#         break
    
# fig=plt.figure()
# # plt.plot(detuning)

# # plt.plot(deg/np.pi,g2_0_h,'b')
# plt.plot(detuning/kappa_h,((g20_HH_out_4)),'b')
# # plt.plot(-1*Delta/omega_m,S1,'b')
# # plt.xlabel('$\Delta /\omega _m$')
# # plt.xlim(-20, 20)
# # plt.title('$g^2_{LL,out}(0)$ with $90$ incident')
# plt.xlabel('$\Delta /\kappa_h $')
# # plt.xlim(-20,20)
# plt.rcParams['savefig.dpi'] = 600
# plt.rcParams['figure.dpi'] = 600
# plt.rcParams['figure.figsize']=(3.26,2.36)

# plt.show()


"---------------------------versus input-output polarizer---------------------------------"

def output_g20_tuning_OPZR(theta_i,theta_o):
    Delta_m=0
    Delta_h=0
    Delta_v=0

    
    "i-o angles"
    theta_in=theta_i
    theta_out=theta_o
    # theta_o=np.pi*105/180
    
    
    "decay rates"
    kappa_h_ex=1*kappa_h_in
    # kappa_h_ex=kappa_ex
    kappa_h=kappa_h_in+kappa_h_ex

    kappa_v_ex=kappa_h_ex
    kappa_v=kappa_v_in+kappa_v_ex
    
    "input amplitudes"
    eta_h=kappa_h_in*0.6*np.cos(theta_in)
    eta_v=kappa_h_in*0.6*np.sin(theta_in)
    
    "output operators"
    b_h_in=eta_h/sqrt(kappa_h_ex)
    b_v_in=eta_v/sqrt(kappa_v_ex)

    b_h=b_h_in-1j*sqrt(kappa_h_ex)*a_h
    b_v=b_v_in-1j*sqrt(kappa_v_ex)*a_v
    
    

    H_free=-Delta_h*a_h.dag()*a_h-Delta_v*a_v.dag()*a_v-Delta_m*bb.dag()*bb
    H_i=+g0*a_h.dag()*a_v*(bb.dag())+g0*a_v.dag()*a_h*(bb)

    # H_p=eta_h*a_h.dag()+conjugate(eta_h)*a_h+eta_v*a_v.dag()+conjugate(eta_v)*a_v #HV probe
    H_p=eta_v*a_v.dag()+conjugate(eta_v)*a_v+eta_h*a_h.dag()+conjugate(eta_h)*a_h+pm_in*(bb.dag()+bb) #V probe
    H=H_free+H_i+H_p

    c_list=[]
    c_list.append(sqrt(kappa_h)*a_h)
    c_list.append(sqrt(kappa_v)*a_v)
    c_list.append(sqrt(kappa_m*(n_m+1))*bb)
    c_list.append(sqrt(kappa_m*(n_m))*bb.dag())
    
    rho_ss=steadystate(H, c_list)
    
    # rho_out=tensor(np.cos(theta_o)*I_o*bh,np.sin(theta_o)*I_o*bv,I_m)*rho_ss*tensor(np.cos(theta_o)*I_o*bh,np.sin(theta_o)*I_o*bv,I_m).dag()
    
    # b_out_h=(np.cos(theta_o)**2*I_o)*a+(0.5*np.sin(2*theta_o)*I_o)*a
    # b_out_v=(0.5*np.sin(2*theta_o)*I_o)*a+(np.sin(theta_o)**2*I_o)*a
    # b_h=tensor(b_out_h,I_o,I_m)
    # b_v=tensor(I_o,b_out_v,I_m)
    b_out=tensor(np.cos(theta_out)*I_o,I_o,I_m)*b_h+tensor(I_o,np.sin(theta_out)*I_o,I_m)*b_v
    b_out_ort=tensor(np.cos(theta_out)*I_o,I_o,I_m)*b_h-tensor(I_o,np.sin(theta_out)*I_o,I_m)*b_v

    # a_out=tensor(np.cos(theta_o)*I_o,I_o,I_m)*a_h+tensor(I_o,np.sin(theta_o)*I_o,I_m)*a_v

    # p_h_n=expect(b_h.dag()*b_h,rho_ss)
    # p_v_n=expect(b_v.dag()*b_v,rho_ss)
    
    # p_int_n=expect((b_h.dag()*b_v+b_h*b_v.dag()),rho_ss)
    
    # p_v_nn=expect(b_v.dag()*b_v.dag()*b_v*b_v,rho_ss)
    # p_h_nn=expect(b_h.dag()*b_h.dag()*b_h*b_h,rho_ss)
    
    
    # p_hv_nn=expect(b_h.dag()*b_h.dag()*b_v*b_v,rho_ss)
    # p_vh_nn=expect(b_v.dag()*b_v.dag()*b_h*b_h,rho_ss)
    # p_3hv_nn=expect((b_h.dag()*b_h.dag()*b_h*b_v+b_h.dag()*b_h*b_h*b_v.dag()),rho_ss)
    # p_h3v_nn=expect((b_h*b_v.dag()*b_v.dag()*b_v+b_h.dag()*b_v.dag()*b_v*b_v),rho_ss)


    
    G_out=expect(b_out.dag()*b_out.dag()*b_out*b_out,rho_ss)
    I_out=expect(b_out.dag()*b_out,rho_ss)
    
    # G_out_ort=expect(b_out_ort.dag()*b_out_ort.dag()*b_out_ort*b_out_ort,rho_ss)
    # I_out_ort=expect(b_out_ort.dag()*b_out_ort,rho_ss)

    # I=expect(a_out.dag()*a_out,rho_ss)

    # G_out_cross=expect(b_out.dag()*I_out_ort.dag()*I_out_ort*I_out,rho_ss)
    # I_out_cross=expect(b_out.dag()*b_out,rho_ss)

    g20_out=(G_out)/(I_out*I_out)
    # g20_out_ort=(G_out_ort)/(I_out_ort*I_out_ort)

    # g20_out_cross=(G_out_cross)/(I_out*I_out_ort)
    return g20_out

# F=output_g20_tuning_OPZR(45/180*np.pi,174/180*np.pi)
# print(F)



"__________________g20 versus polarizations_____________________"
# angle_in=np.zeros(90)


# angle_out=np.zeros(101)


# angle_o=80*np.pi/180
# i=0
# g20_out=np.zeros(101)

# while i<=99:
#     angle_out[i]=angle_o
#     g20_out[i]=np.real(output_g20_tuning_OPZR(angle_out[i]))
#     angle_o=angle_o+np.pi*(0.2/180)
#     print(angle_out[i]*180/np.pi)
#     i=i+1
#     if i==100:
#         angle_out[i]=100*np.pi/180
#         g20_out[i]=np.real(output_g20_tuning_OPZR(angle_out[i]))
#         print(angle_out[i]*180/np.pi)
#         break


# fig=plt.figure()

# plt.plot(angle_out*180/np.pi,np.log10((g20_out)),'b')
# plt.xlim(80,100)
# plt.xlabel('$\\theta_{out}$')

# plt.rcParams['savefig.dpi'] = 600
# plt.rcParams['figure.dpi'] = 600
# plt.rcParams['figure.figsize']=(3.26,2.36)

# plt.show()


"_______________________3D_______________________"

x=np.linspace(0, 210, 36)   #input angle
y=np.linspace(0, 210, 36)  #output angle
z=np.zeros((36,36))
g20=np.zeros((36,36))

# g20_xx=np.zeros((36,36))
# g20_xx=g20

# g20_yy=np.zeros((36,36))
# g20_yy=g20

# g20_xy=np.zeros((36,36))
# g20_xy=g20

for i, xx in enumerate(x):
    theta_in=np.pi*(xx+0)/180
    for j, yy in enumerate(y):
        theta_out=np.pi*(yy+0)/180
        tx=np.cos(theta_in)
        ty=np.sin(theta_out)
        z[i,j]=tx**2+ty**2
        g20[i,j]=output_g20_tuning_OPZR(theta_in,theta_out)

        print('input angle=',theta_in/np.pi,'pi,','output angle=',theta_out/np.pi,'pi')
        

fig=plt.figure()
p=plt.pcolor(y,x,np.log10(g20),cmap=cm.jet)
# p=plt.pcolor(y,x,np.log10(g20),cmap=cm.gist_heat) # photon intensity

plt.ylabel('$\\theta_{in}$')
plt.xlabel('$\\theta_{out}$')
plt.title('$log_{10}[g^{(2)}_{out}(0)]$')
# plt.title('$log_{10}[n_{out}]$')


# plt.xticks([0,30,60,90,120,150,180,210])
# plt.yticks([0,30,60,90,120,150,180,210])

# plt.xlim(0,180)
# plt.ylim(0,180)

# plt.xlim(60,180)
# plt.ylim(90,210)

# plt.xlim(90,210)
# plt.ylim(0,120)

fig.colorbar(p)
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
plt.rcParams['figure.figsize']=(3.26,2.36)

fig.show()

"_________________cross-correlation criterion_______________"

# G_cross=np.zeros((36,36))
# xx=0
# yy=0
# xy=0

# for i,xx in enumerate(x):
#     for j,yy in enumerate(y):
#         xx=g20_xx[i,j]
#         yy=g20_yy[i,j]
#         xy=g20_xy[i,j]
#         G_cross[i,j]=(xy**2)/(xx*yy)



# # g20_in_172=np.zeros(200)
# # g20_in_172=g20

# x=np.linspace(0, 210, 200)   #input angle
# g20=np.zeros(200)


# for j, xx in enumerate(x):
#     theta_in=np.pi*xx/180
#     theta_out=np.pi*172/180
#     tx=np.cos(theta_in)
#     ty=np.sin(theta_out)
#     # z[i,j]=tx**2+ty**2
#     g20[j]=output_g20_tuning_OPZR(theta_in,theta_out)

#     print('input angle=',theta_in/np.pi,'pi,','output angle=',theta_out/np.pi,'pi')












# for i in range (0,100):
#     # Delta[i]=float(nos[i])
#     g20_out[i]=float(nos[i])
    
# # print(Delta/omega_m)
# # print(angle_out)
# print(g20_out)




tf=time.time()
print(tf-ti,'s')