#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show

eiz1x = np.loadtxt('data/eiz_x_51.txt') # LDS in perpendicular direction
eiz1z = np.loadtxt('data/eiz_z_51.txt') # LDS in parallel direction

eiz2x = np.loadtxt('data/eiz_x_290.txt') # LDS in perpendicular direction
eiz2z = np.loadtxt('data/eiz_z_290.txt') # LDS in parallel direction

eiz_w = np.loadtxt('data/eiz_w.txt') # LDS of water, intervening medium

# Constants
r_1 = 0.5e-9
r_2 = 0.5e-9
theta = np.pi/2
c = 2.99e8              # [m/s]
coeff = 2.411e14        # [rad/s]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23 # [J]

# Matsubara frequencies
#ns = z/(0.159)#

ns = np.arange(0.,500.) 
zs = ns * coeff         
Ls = np.arange(1e-9,100e-9,1e-9)  # separation distance between 2 cyclinders
thetas = np.arange(np.pi/12, np.pi/2,0.01)

#Integration vars
r = np.linspace(0.,2.**17, 1.+2.**17)
phi =  np.linspace(0.,2.*np.pi, 1.+2.**17)

#T  = np.linspace(0.,2.**17, 1.+2.**17)
#U  = np.linspace(0.,2.**17, 1.+2.**17)

def rho1(Z,Eiz1z):
    return np.sqrt(r*r + Z*Z*Eiz1z/(c*c))

def rho2(Z,Eiz1z):
    return np.sqrt(r*r + Z*Z*Eiz2z/(c*c))

def rho3(Z,Eiz1z):
    return np.sqrt(r*r + Z*Z*Eiz3/(c*c))

def rhotilda1(Z,Eiz1x,Eiz1z):
    return np.sqrt(r*r+(Eiz1z/Eiz1x-1.)*r*np.cos(phi)*r*np.cos(phi)+Z*Z*Eiz1z/(c*c))

def rhotilda2(Z,Eiz2x,Eiz2z,theta): 
    return np.sqrt(r*r+(Eiz2z/Eiz2x-1.)r*np.cos(phi-theta)*r*np.cos(phi-theta)+Z*Z*Eiz2z/(c*c))

def exponential(Rho3,L): 
    return np.expi(-2.*Rho3*L)

def rsinphi(phi): 
    return r*np.sin(phi)

def gamma(Rho1,Rho2,Rho3,Rhotil1,Rhotil2,Eiz1x,Eiz2x,Eiz3,theta):
    return (Rho1 + Rho3)*(Rho2 + Rho3)*(((Eiz3*Rho1+Eiz1x*Rho3)\
            -(Eiz1x*(Rhotil1-Rho1)*(rsinphi*rsinphi-Rho1*Rho3))/(Rho1*Rho1-\
            rsinphi*rsinphi))*((Eiz3*Rho2 + Eiz2x*Rho3)\
            -(Eiz2x*(Rhotil2-Rho2)*r*np.sin(phi+theta)*r*np.sin(phi+theta)-Rho2*Rho3)\
            /(Rho2*Rho2-r*np.sin(phi+theta)*r*np.sin(phi+theta))))

def A(Rho1,Rho2,Rho3,Rhotil2,Eiz1x,Eiz2x,Eiz3,exponential):
    return ((Rho1 + Rho3)*(Rho2+Rho3)-(Rho1-Rho3)*(Rho2-Rho3)*exponential)\
            *((Eiz3*Rho1+Eiz1x*Rho3)*(Eiz3*Rho2 +\
            Eiz2x*Rho3)-(Eiz3*Rho1-Eiz1x*Rho3)*(Eiz3*Rho2 -\
            Eiz2x*Rho3)*exponential)\
            -((Eiz2x*(Rhotil2-Rho2))/(Rho1*Rho1 - rsinphi*rsinphi))\
            *((rsinphi*rsinphi-Rho1*Rho3)*((Eiz3*Rho2 +\
            Eiz2x*Rho3)*(Rho2+Rho3)*(Rho1+Rho3)+2.*(Eiz2x-Eiz3)\
            *(rsinphi*rsinphi*(r*r*Rho1-Rho2*Rho3*Rho3)+Rho1*Rho3*Rho3*(r*r-2.*rsinphi*rsinphi\
            +Rho1*Rho2))*exponential+(rsinphi*rsinphi +\
            Rho1*Rho3)*(Eiz3*Rho2-Eiz2x*Rho3)\
            *(Rho1-Rho3)*(Rho2-Rho3)*exponential*exponential))

def B(Rho1,Rho2,Rho3,Rhotil1,Eiz1x,Eiz2x,Eiz3,exponential): 
    return ((Eiz3*Rho1+Eiz1x*Rho3)*(Rho1+Rho3)*(Rho2+Rho3)+2.*(Eiz1x-Eiz3)*(r*r*Rho2\
            -Rho1*Rho3*Rho3-2.*Rho2*Rho3*Rho3)*exponential+(Eiz3*Rho1-Eiz1x*Rho3)*(Rho1-Rho3)\
            *(Rho2-Rho3)*exponential*exponential)+((Eiz1x*(Rhotil1-Rho1))/(Rho1*Rho1-rsinphi*rsinphi))*(-(rsinphi*rsinphi-Rho1*Rho3)\
            *(Rho1+Rho3)*(Rho2+Rho3)+2.*(rsinphi*rsinphi*(Rho1*Rho3+Rho3*Rho3)-Rho1*Rho1*Rho3*Rho3\
            +Rho1*Rho2*Rho3*Rho3)*exponential-(rsinphi*rsinphi\
            +Rho1*Rho3)*(Rho1-Rho3)*(Rho2\
            -Rho3)*exponential*exponential)


def C(Rho1,Rho2,Rho3,Rhotil1,Eiz1x,Eiz2x,Eiz3,exponential):
    return Rho1*Rho3*(-(Eiz3*Rho1+ Eiz1x*Rho3)*Rho13*Rho23+2.*Rho3*(Eiz1x-Eiz3)\
            (r*r+Rho1*Rho2)*exponential+(Eiz3*Rho1-Eiz1x*Rho3)*(Rho1-Rho3)*(Rho2\
            -Rho3)*exponential*eponential)\
            +Eiz1x*(Rhotil1-Rho1)/(Rho1*Rho1-rsinphi*rsinphi)*Rho2*Rho3*((rsinphi*rsinphi-Rho1*Rho3)\
            (Rho1+Rho3)*(Rho2+Rho3)+2.*Rho3*(Rho1*Rho1*Rho2+Rho1*Rho3*Rho3+rsinphi*rsinphi*(Rho1\
            -Rho2))*exponential-(rsinphi*rsinphi+Rho1*Rho3)*(Rho1-Rho3)*(Rho2\
            -Rho3)*exponential*exponential)

def E(Rho1,Rho2,Rho3,Rhotil1,Eiz1x,exponential):
    return 4.*Rho1*Rho2*Rho3*Rho3*Eiz1x*(Rhotil1-Rho1)/(Rho1*Rho1-rsinphi*rsinphi)*exponential

#rho13 = rho1 + rho3
#rho23 = rho2 + rho3
#rho31 = rho1 - rho3
#rho21 = rho1 - rho2
#rho32 = rho2 - rho3

p = np.zeros(shape = (len(Ls),len(ns)))
A0 = np.zeros(shape = (len(Ls),len(ns)))
A2 = np.zeros(shape = (len(Ls),len(ns)))
A2_org = np.zeros(shape = (len(Ls),len(ns)))
A2_new = np.zeros(shape = (len(Ls),len(ns)))
G = np.zeros(len(Ls))
G_new = np.zeros(len(Ls))

a_1 =   Aiz(eiz_x_051,eiz_z_051,eiz_w)
a_2 =   Aiz(eiz_x_290,eiz_z_290,eiz_w)
delta_1 = Delta(eiz_z_051,eiz_w)
delta_2 = Delta(eiz_z_290,eiz_w)
delta_prp1 = Delta_perp(eiz_x_051,eiz_w)
delta_prp2 = Delta_perp(eiz_x_290,eiz_w)
D   = delta_prp2*delta_1*(delta_1 - 2.*delta_prp1)*(delta_2-2.*delta_prp2)
DD  = DDDg2(eiz_x_051,eiz_x_290,eiz_z_051,eiz_z_290,eiz_w)
DDD = DDDg2(eiz_x_051,eiz_x_290,eiz_z_051,eiz_z_290,eiz_w)

## Integrand, A0(n=0) and A2(n=0) terms 
#f0_term0 = U*U*U * np.exp(-2.* U)\
#        *2.*(1.+3.*a_1[0])*(1.+3.*a_2[0])
#        #*((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))\
#f2_term0 =  U*U*U * np.exp(-2.* U)\
#        *(1.-a_1[0])*(1.-a_2[0])
#        #*((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))\
#Ft0_term0 = romb(f0_term0)
##Ft0_term0 = np.sum(f0_term0)
#Ft2_term0 = romb(f2_term0)
##Ft2_term0 = np.sum(f2_term0)

# Calculate 1 \geq n \leq 500 terms
for i,L in enumerate(Ls):
    print 'Computing A for separation number %3.0f of %3.0f'%(i, len(Ls))
    for j,n in enumerate(ns):
        for k,theta in enumerate(thetas):

       D = 
       (1./gamma(rho1[zs[j],eiz1x[j]],rho2[zs[j],eiz1x[j]],rho3[zs[j],eiz1x[j]],rhotilda1[zs[j],eiz1x[j],eiz1z[j]],rhotilda2[zs[j],eiz2x,eiz2z,thetas],eiz1x[j],eiz2x[j],eiz3[j],thetas[k]))
       *(
            A(rho1[zs[j],eiz1x[j]],rho2[zs[j],eiz1x[j]],rho3[zs[j],eiz1x[j]],rhotilda2[zs[j],eiz1x[j],eiz1z[j]],rhotilda2[zs[j],eiz2x,eiz2z,thetas],eiz1x[j],eiz2x[j],eiz3[j],exponential(rho3[zs[j],eiz1x[j]],Ls[i]))
            
            - eiz2x*(rhotilda2-rho2)/(rho2*rho2-r*np.sin(phi+theta)*r*np.sin(phi+theta))
            
            *(
            B(rho1[zs[j],eiz1x[j]],rho2[zs[j],eiz1x[j]],rho3[zs[j],eiz1x[j]],rhotilda1[zs[j],eiz1x[j],eiz1z[j]],rhotilda2[zs[j],eiz2x,eiz2z,thetas],eiz1x[j],eiz2x[j],eiz3[j],exponential(rho3[zs[j],eiz1x[j]],Ls[i]))
                
            *r*np.sin(phi+theta)*r*np.sin(phi+theta)
            
                -E(rho1[zs[j],eiz1x[j]],rho2[zs[j],eiz1x[j]],rho3[zs[j],eiz1x[j]],rhotilda1[zs[j],eiz1x[j],eiz1z[j]],eiz1x[j],exponential(rho3[zs[j],eiz1x[j]],Ls[i]))
                
                *(
                    2.*r*r*np.sin(phi)*np.cos(theta)*np.sin(phi+theta)+rho3*np.sin(theta)*rho3*np.sin(theta)   
                )
                +C(rho1[zs[j],eiz1x[j]],rho2[zs[j],eiz1x[j]],rho3[zs[j],eiz1x[j]],rhotilda1[zs[j],eiz1x[j],eiz1z[j]],rhotilda2[zs[j],eiz2x,eiz2z,thetas],eiz1x[j],eiz2x[j],eiz3[j],exponential(rho3[zs[j],eiz1x[j]],Ls[i]))
            )
        ) 
        
        
        
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand A0
        f0 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)\
                *(2.*(1.+3.*a_1[j])*(1.+3.*a_2[j])*T*T*T*T\
                + 4.*(1.+2.*a_1[j]+2.*a_2[j]+3.*a_1[j]*a_2[j])*T*T \
                + 4.*(1.+a_1[j])*(1.+a_2[j]))
                #*((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))\
        # Integrand A2                
        f2 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
                /(T*T+1.)\
                *((T*T*T*T +4.*(T*T)+4.)*(1.-a_1[j])*(1.-a_2[j]))
                #*((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))\
        f2_org = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)*\
                delta_1[j]*delta_2[j]\
                *((T*T*T*T +4.*(T*T)+4.)*(1.-a_1[j])*(1.-a_2[j]))
                #*((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))\
        f2_new =T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)*\
                 DDDg2(eiz_x_051[j],eiz_x_290[j],eiz_z_051[j],eiz_z_290[j],eiz_w[j])*(T*T*T*T +4.*(T*T)+4.)
        #Ft0 = np.sum(f0)
        #Ft2 = np.sum(f2)
        Ft0 = romb(f0)
        Ft2 = romb(f2)
        Ft2_org = romb(f2_org)
        Ft2_new = romb(f2_new)
        #Ft = romb(f , axis = 1)
        #Fty =romb(Ft, axis = 0)

        A0[i,j] = delta_1[j]*delta_2[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft0
        A0[i,0] = (1./2) * delta_1[0]*delta_2[0]*Ft0_term0
        #A0[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft0_term0
        
        A2[i,j] = delta_1[j]*delta_2[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
        A2_org[i,j] = p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2_org
        A2_new[i,j] = p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2_new
        A2_org[i,0] = (1./2) * delta_1[0]*delta_2[0]*Ft2_term0
        A2_new[i,0] = (1./2) * delta_1[0]*delta_2[0]*Ft2_term0
        A2[i,0] = (1./2) * delta_1[0]*delta_2[0]*Ft2_term0
        #A2_new[i,0] = (1./2) * delta_1[0]*delta_2[0]*Ft2_term0
        #A2[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft2_term0

        #A0[A0>1e6]= np.nan #nOTE: remove me later
        #A2[A2>1e6]= np.nan #nOTE: remove me later
    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
    sum_A2_org = (kbT/(32.)) * np.sum(A2_org, axis = 1)
    sum_A2_new = (kbT/(32.)) * np.sum(A2_new, axis = 1)
    #G = (((5e-10)**4)/(L**4))*(sum_A0 + sum_A2)
    G = (-np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*np.sin(theta)*L*L*L*L)*(sum_A0+sum_A2*np.cos(2.*theta))
    G_org=(-np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*np.sin(theta)*L*L*L*L)*(sum_A0+sum_A2_org*np.cos(2.*theta))
    G_new=(-np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*np.sin(theta)*L*L*L*L)*(sum_A0+sum_A2_new*np.cos(2.*theta))
#    pl.plot(ns,A2_org[i,:],'-')
#    pl.plot(ns,-A2_new[i,:],'--')
#pl.show()
np.savetxt('data/A0_51w290_perpendicular_ret.txt',sum_A0)
np.savetxt('data/A2_51w290_perpendicular_ret.txt',sum_A2)
np.savetxt('data/A0_51w290_n_l.txt',A0)
np.savetxt('data/A2_51w290_n_l.txt', A2)

pl.figure()
pl.plot(1e9*Ls,1e21*G,'b--')
#pl.plot(1e9*Ls,1e21*G_org,'g:|')
#pl.plot(1e9*Ls,1e21*G_new,'r:')
#pl.plot(1e9*Ls,1e21*G_new)
pl.show()

pl.figure()
pl.loglog(1e9*Ls,G/G[0],'b--')
#pl.plot(1e9*Ls,G_org/G_org[0],'g:|')
#pl.plot(1e9*Ls,G_new/G_new[0],'r:')
#pl.plot(1e9*Ls,1e21*G_new)
pl.show()

pl.figure()
pl.loglog(ns,(kbT/(32.)) * A2[0,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,(kbT/(32.)) * A2[1,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.loglog(ns,(kbT/(32.)) * A2[2,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.loglog(ns,(kbT/(32.)) * A2[3,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))

#pl.loglog(ns,(kbT/(32.)) * A2_org[0,:], 'g--')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.loglog(ns,(kbT/(32.)) * A2_org[1,:], 'g--')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
#pl.loglog(ns,(kbT/(32.)) * A2_org[2,:], 'g--')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
#pl.loglog(ns,(kbT/(32.)) * A2_org[3,:], 'g--')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#
#pl.loglog(ns,(kbT/(32.)) * A2_new[0,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.loglog(ns,(kbT/(32.)) * A2_new[1,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
#pl.loglog(ns,(kbT/(32.)) * A2_new[2,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
#pl.loglog(ns,(kbT/(32.)) * A2_new[3,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.legend(loc = 'best')
#pl.title(r'65w65 Matsubara terms')
#pl.title(r'90w90 Matsubara terms')
#pl.title(r'91w91 Matsubara terms')
#pl.title(r'93w93 Matsubara terms')
pl.title(r'290w290 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
#pl.savefig('plots/65_A_vs_n.pdf')
#pl.savefig('plots/90_A_vs_n.pdf')
#pl.savefig('plots/91_A_vs_n.pdf')
#pl.savefig('plots/93_A_vs_n.pdf')
pl.savefig('plots/51w290_diff_delta2g2__A_vs_n.pdf')
pl.show()
#
pl.figure()
pl.plot(ns,(kbT/(32.)) * A2[0,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.plot(ns,(kbT/(32.)) * A2[1,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.plot(ns,(kbT/(32.)) * A2[2,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.plot(ns,(kbT/(32.)) * A2[3,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.plot(ns,(kbT/(32.)) * A2_org[0,:], 'g--')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.plot(ns,(kbT/(32.)) * A2_org[1,:], 'g--')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
#pl.plot(ns,(kbT/(32.)) * A2_org[2,:], 'g--')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
#pl.plot(ns,(kbT/(32.)) * A2_org[3,:], 'g--')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.plot(ns,(kbT/(32.)) * A2_new[0,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.plot(ns,(kbT/(32.)) * A2_new[1,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
#pl.plot(ns,(kbT/(32.)) * A2_new[2,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
#pl.plot(ns,(kbT/(32.)) * A2_new[3,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.legend(loc = 'best')
#pl.title(r'65w65 Matsubara terms')
#pl.title(r'90w90 Matsubara terms')
#pl.title(r'91w91 Matsubara terms')
#pl.title(r'93w93 Matsubara terms')
pl.title(r'290w290 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
#pl.savefig('plots/65_A_vs_n.pdf')
#pl.savefig('plots/90_A_vs_n.pdf')
#pl.savefig('plots/91_A_vs_n.pdf')
#pl.savefig('plots/93_A_vs_n.pdf')
pl.savefig('plots/51w290_nonlog_compare_delta2g2__A_terms_vs_n.pdf')
pl.show()
#
pl.figure()
pl.plot(ns,(kbT/(32.)) * A0[0,:], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.plot(ns,(kbT/(32.)) * A0[1,:], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.plot(ns,(kbT/(32.)) * A0[2,:], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.plot(ns,(kbT/(32.)) * A0[3,:], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.plot(ns,(kbT/(32.)) * A2_org[0,:], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.plot(ns,(kbT/(32.)) * A2_org[1,:], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.plot(ns,(kbT/(32.)) * A2_org[2,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.plot(ns,(kbT/(32.)) * A2_org[3,:], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.legend(loc = 'best')
#plplot.title(r'65w65 Matsubara terms')
#plplot.title(r'90w90 Matsubara terms')
#plplot.title(r'91w91 Matsubara terms')
#pl.title(r'93w93 Matsubara terms')
pl.title(r'290w290 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
#pl.savefig('plots/65_A_vs_n.pdf')
#pl.savefig('plots/90_A_vs_n.pdf')
#pl.savefig('plots/91_A_vs_n.pdf')
#pl.savefig('plots/93_A_vs_n.pdf')
pl.savefig('plots/51w290_nonlog_A_terms_vs_n.pdf')
pl.show()

pl.figure()
pl.loglog(ns,(kbT/(32.)) * A0[0,:], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,(kbT/(32.)) * A0[1,:], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.loglog(ns,(kbT/(32.)) * A0[2,:], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.loglog(ns,(kbT/(32.)) * A0[3,:], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.loglog(ns,(kbT/(32.)) * A2[0,:], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,(kbT/(32.)) * A2[1,:], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.loglog(ns,(kbT/(32.)) * A2[2,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.loglog(ns,(kbT/(32.)) * A2[3,:], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.legend(loc = 'best')
#pl.title(r'65w65 Matsubara terms')
#pl.title(r'90w90 Matsubara terms')
#pl.title(r'91w91 Matsubara terms')
#pl.title(r'93w93 Matsubara terms')
pl.title(r'290w290 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
#pl.savefig('plots/65_A_vs_n.pdf')
#pl.savefig('plots/90_A_vs_n.pdf')
#pl.savefig('plots/91_A_vs_n.pdf')
#pl.savefig('plots/93_A_vs_n.pdf')
pl.savefig('plots/51w290_A_vs_n.pdf')
pl.show()
#
pl.figure()
pl.loglog(ns,A0[0,:]/sum_A0[0], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,A0[1,:]/sum_A0[1], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.loglog(ns,A0[2,:]/sum_A0[2], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.loglog(ns,A0[3,:]/sum_A0[3], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.loglog(ns,A2[0,:]/sum_A2[0], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,A2[1,:]/sum_A2[1], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.loglog(ns,A2[2,:]/sum_A2[2], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.loglog(ns,A2[3,:]/sum_A2[3], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.legend(loc = 'best')
#pl.title(r'65w65 Matsubara terms')
#pl.title(r'90w90 Matsubara terms')
#pl.title(r'91w91 Matsubara terms')
#pl.title(r'93w93 Matsubara terms')
pl.title(r'290w290 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
#pl.savefig('plots/65_A_vs_n.pdf')
#pl.savefig('plots/90_A_vs_n.pdf')
#pl.savefig('plots/91_A_vs_n.pdf')
#pl.savefig('plots/93_A_vs_n.pdf')
pl.savefig('plots/51w290_relative_A_vs_n.pdf')
pl.show()

pl.figure()
pl.plot(ns,(A0[0,:]+A2[0,:])/(sum_A0[0]+sum_A2[0]), 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.plot(ns,(A0[3,:]+A2[3,:])/(sum_A0[3]+sum_A2[3]), 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.plot(ns,(A0[6,:]+A2[6,:])/(sum_A0[6]+sum_A2[6]), 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[6]))
pl.plot(ns,(A0[9,:]+A2[9,:])/(sum_A0[9]+sum_A2[9]), 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[9]))
pl.legend(loc = 'best')
pl.title(r'290w290 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
pl.show()

DA2_l = np.diff(A2, axis = 0)

DA2_n = np.diff(A2, axis = 1)

pl.figure()
pl.semilogx(ns,DA2_l[0,:],  'b-' ,linewidth = 1.0,label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[0]))
pl.semilogx(ns,DA2_l[9,:],  'b--',linewidth = 1.0,label =r'$\ell=%2.1f\,nm$' %(1e9*Ls[9]))
pl.semilogx(ns,DA2_l[24,:], 'b-.',linewidth = 1.0,label =r'$\ell=%2.1f\,nm$'%(1e9*Ls[24]))
pl.semilogx(ns,DA2_l[49,:], 'b:' ,linewidth = 1.0,label =r'$\ell=%2.1f\,nm$'%(1e9*Ls[49]))
pl.semilogx(ns,DA2_l[3,:]-DA2_l[3,:], 'k:')# , linewidth = 1.0, label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[3]))
pl.legend(loc = 'best')
pl.title(r'dA2/dl contribution to Matsubara sum for [5,1] and [26,0]')
pl.ylabel(r'$d\mathcal{A}^{(2)}_{n}/dl$')
pl.xlabel(r'$n$')
pl.show()

pl.figure()
pl.plot(ns,(A0[0,:]+A2[0,:]), 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.plot(ns,(A0[3,:]+A2[3,:]), 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.plot(ns,(A0[6,:]+A2[6,:]), 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[6]))
pl.plot(ns,(A0[9,:]+A2[9,:]), 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[9]))
pl.legend(loc = 'best')
pl.title(r'290w290 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
pl.show()

pl.figure()
pl.loglog(ns,(A0[0,:]+A2[0,:])/(sum_A0[0]+sum_A2[0]), 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,(A0[3,:]+A2[3,:])/(sum_A0[3]+sum_A2[3]), 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.loglog(ns,(A0[6,:]+A2[6,:])/(sum_A0[6]+sum_A2[6]), 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[6]))
pl.loglog(ns,(A0[9,:]+A2[9,:])/(sum_A0[9]+sum_A2[9]), 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[9]))
pl.legend(loc = 'best')
pl.title(r'290w290 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
pl.show()
#pl.figure()
#pl.loglog(ns,(kbT/(32.)) * A0[ 0,:], 'b-',  label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.loglog(ns,(kbT/(32.)) * A0[ 4,:], 'g-',  label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[4]))
#pl.loglog(ns,(kbT/(32.)) * A0[14,:], 'r-',  label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[14]))
#pl.loglog(ns,(kbT/(32.)) * A0[24,:], 'c-',  label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[24]))
#pl.loglog(ns,(kbT/(32.)) * A0[34,:], 'k-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[34]))
#pl.loglog(ns,(kbT/(32.)) * A0[44,:], 'm-',  label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[44]))
#pl.loglog(ns,(kbT/(32.)) * A0[54,:], 'y-',  label =r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[54]))
#pl.loglog(ns,(kbT/(32.)) * A2[0,:],  'b:')
#pl.loglog(ns,(kbT/(32.)) * A2[4,:],  'g:')
#pl.loglog(ns,(kbT/(32.)) * A2[14,:], 'r:')
#pl.loglog(ns,(kbT/(32.)) * A2[24,:], 'c:')
##pl.loglog(ns,(kbT/(32.)) * A2[32,:], 'm:')
##pl.loglog(ns,(kbT/(32.)) * A2[33,:], 'm:' )
#pl.loglog(ns,(kbT/(32.)) * A2[34,:], 'k:')
##pl.loglog(ns,(kbT/(32.)) * A2[35,:], 'm:' )
##pl.loglog(ns,(kbT/(32.)) * A2[36,:], 'm:' )
##pl.loglog(ns,(kbT/(32.)) * A2[39,:], 'y:' )
#pl.loglog(ns,(kbT/(32.)) * A2[44,:], 'm:' )
#pl.loglog(ns,(kbT/(32.)) * A2[54,:], 'y:' )
#pl.legend(loc = 'best')
##pl.title(r'65w65 Matsubara terms')
##pl.title(r'90w90 Matsubara terms')
##pl.title(r'91w91 Matsubara terms')
##pl.title(r'93w93 Matsubara terms')
#pl.title(r'Matsubara terms for [5,1] and [29,0] in water')
#pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
#pl.xlabel(r'$n$')
##pl.savefig('plots/65_A_vs_n.pdf')
##pl.savefig('plots/90_A_vs_n.pdf')
##pl.savefig('plots/91_A_vs_n.pdf')
##pl.savefig('plots/93_A_vs_n.pdf')
#pl.savefig('plots/140522_Hopkins_51w290_A0A2_vs_n.pdf')
#pl.show()
#
#pl.figure()
#pl.loglog(ns,(kbT/(32.)) * A0[ 0,:]+(kbT/(32.)) * A2[ 0,:],'b-', label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[ 0]))
#pl.loglog(ns,(kbT/(32.)) * A0[ 4,:]+(kbT/(32.)) * A2[ 4,:],'g-', label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[ 4]))
#pl.loglog(ns,(kbT/(32.)) * A0[14,:]+(kbT/(32.)) * A2[14,:],'r-', label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[14]))
#pl.loglog(ns,(kbT/(32.)) * A0[24,:]+(kbT/(32.)) * A2[24,:],'c-', label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[24]))
#pl.loglog(ns,(kbT/(32.)) * A0[34,:]+(kbT/(32.)) * A2[34,:],'k-', label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[34]))
#pl.loglog(ns,(kbT/(32.)) * A0[44,:]+(kbT/(32.)) * A2[44,:],'m-', label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[44]))
#pl.loglog(ns,(kbT/(32.)) * A0[54,:]+(kbT/(32.)) * A2[54,:],'y-', label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[54]))
##pl.loglog(ns,(kbT/(32.)) * A2[0,:], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
##pl.loglog(ns,(kbT/(32.)) * A2[1,:], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
##pl.loglog(ns,(kbT/(32.)) * A2[2,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
##pl.loglog(ns,(kbT/(32.)) * A2[3,:], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.legend(loc = 'best')
##pl.title(r'65w65 Matsubara terms')
##pl.title(r'90w90 Matsubara terms')
##pl.title(r'91w91 Matsubara terms')
##pl.title(r'93w93 Matsubara terms')
#pl.title(r'$\mathcal{A}^{(0)}_{n}+\mathcal{A}^{(2)}_{n}$ Matsubara terms for [5,1] and [29,0] in water')
#pl.ylabel(r'$\mathcal{A}^{(0)}_{n}+\mathcal{A}^{(2)}_{n}$')
#pl.xlabel(r'$n$')
##pl.savefig('plots/65_A_vs_n.pdf')
##pl.savefig('plots/90_A_vs_n.pdf')
##pl.savefig('plots/91_A_vs_n.pdf')
##pl.savefig('plots/93_A_vs_n.pdf')
#pl.savefig('plots/140522_Hopkins_51w290_A0plusA2_vs_n.pdf')
#pl.show()
#
#pl.figure()
#pl.plot(ns[:499],1e21*np.diff((kbT/(32.)) * A0[ 0,:]+(kbT/(32.)) * A2[ 0,:]),'b-', label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[ 0]))
#pl.plot(ns[:499],1e21*np.diff((kbT/(32.)) * A0[ 4,:]+(kbT/(32.)) * A2[ 4,:]),'g-', label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[ 4]))
#pl.plot(ns[:499],1e21*np.diff((kbT/(32.)) * A0[14,:]+(kbT/(32.)) * A2[14,:]),'r-', label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[14]))
#pl.plot(ns[:499],1e21*np.diff((kbT/(32.)) * A0[24,:]+(kbT/(32.)) * A2[24,:]),'c-', label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[24]))
#pl.plot(ns[:499],1e21*np.diff((kbT/(32.)) * A0[34,:]+(kbT/(32.)) * A2[34,:]),'k-', label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[34]))
#pl.plot(ns[:499],1e21*np.diff((kbT/(32.)) * A0[44,:]+(kbT/(32.)) * A2[44,:]),'m-', label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[44]))
#pl.plot(ns[:499],1e21*np.diff((kbT/(32.)) * A0[54,:]+(kbT/(32.)) * A2[54,:]),'y-', label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[54]))
#pl.plot(ns[:499],ns[:499]-ns[:499], 'k:')
#pl.legend(loc = 'best')
#pl.title(r'$d(\mathcal{A}^{(0)}_{n}+\mathcal{A}^{(2)}_{n})/dn$ Matsubara terms for [5,1] and [29,0] in water')
#pl.ylabel(r'$d(\mathcal{A}^{(0)}_{n}+\mathcal{A}^{(2)}_{n})/dn\,\,\,\,[zJ]$')
#pl.xlabel(r'$n$')
#pl.savefig('plots/140522_Hopkins_51w290_derivative_A0plusA2_vs_n.pdf')
#pl.show()
#
#pl.figure()
#pl.plot(1e9*Ls[:98],1e21*np.diff((kbT/(32.))*A0[:, 0]+(kbT/(32.))*A2[:, 0],axis = 0),'b-', label = r'$n=%2.1f$'%(ns[ 0]))
#pl.plot(1e9*Ls[:98],1e21*np.diff((kbT/(32.))*A0[:, 4]+(kbT/(32.))*A2[:, 4],axis = 0),'g-', label = r'$n=%2.1f$'%(ns[ 4]))
#pl.plot(1e9*Ls[:98],1e21*np.diff((kbT/(32.))*A0[:,14]+(kbT/(32.))*A2[:,14],axis = 0),'r-', label = r'$n=%2.1f$'%(ns[14]))
#pl.plot(1e9*Ls[:98],1e21*np.diff((kbT/(32.))*A0[:,24]+(kbT/(32.))*A2[:,24],axis = 0),'c-', label = r'$n=%2.1f$'%(ns[24]))
#pl.plot(1e9*Ls[:98],1e21*np.diff((kbT/(32.))*A0[:,34]+(kbT/(32.))*A2[:,34],axis = 0),'k-', label = r'$n=%2.1f$'%(ns[34]))
#pl.plot(1e9*Ls[:98],1e21*np.diff((kbT/(32.))*A0[:,44]+(kbT/(32.))*A2[:,44],axis = 0),'m-', label = r'$n=%2.1f$'%(ns[44]))
#pl.plot(1e9*Ls[:98],1e21*np.diff((kbT/(32.))*A0[:,54]+(kbT/(32.))*A2[:,54],axis = 0),'y-', label = r'$n=%2.1f$'%(ns[54]))
#pl.plot(Ls[:98],Ls[:98]-Ls[:98], 'k:')
#pl.legend(loc = 'best')
#pl.title(r'$d(\mathcal{A}^{(0)}_{n}+\mathcal{A}^{(2)}_{n})/d\ell$ Matsubara terms for [5,1] and [29,0] in water')
#pl.ylabel(r'$d(\mathcal{A}^{(0)}_{n}+\mathcal{A}^{(2)}_{n})/d\ell\,\,\,\,[zJ]$')
#pl.xlabel(r'$\ell$  [nm]')
#pl.savefig('plots/140522_Hopkins_51w290_derivative_A0plusA2_vs_l.pdf')
#pl.show()
#
#print 'A0(separation) = ',sum_A0
#print 'A2(separation) = ',sum_A2
#print 'Contribution to A0 from n=0 term = ', (kbT/(12.*np.pi))*A0[:,0]
#print 'Contribution to A2 from n=0 term = ', (kbT/(12.*np.pi))*A2[:,0]
#
##np.savetxt('data/A0_65_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_65_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_65_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_90_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_90_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_90_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_91_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_91_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_91_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_93_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_93_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_93_perpendicular_ret.txt',Ls)
##
#np.savetxt('data/A0_51w290_perpendicular_ret.txt',sum_A0)
#np.savetxt('data/A2_51w290_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_290_perpendicular_ret.txt',Ls)
#
#A_py_par  = r'$\mathcal{A}_{\parallel}\sf{[python]}$'
#A0_py_per = r'$\mathcal{A}_{perpend}\sf{[python]}$'
#A2_py_per = r'$\mathcal{A}_{perpend}\sf{[python]}$'
#A_GH_par  = r'$\mathcal{A}_{\parallel}\sf{[ G.H. ]}$'
#
#x_ax = r'$\,\ell\,\,\,\rm{[nm]}$'
#y_ax_par = r'$\mathrm{\mathcal{A}_\parallel(\ell)}\,\,\,\rm{[zJ]}$'
#y_ax_per = r'$\mathrm{\mathcal{A}_\perp (\ell)}\,\,\,\rm{[zJ]}$'
#
#def title(cnt1,cnt2,orientation):
#	return r'$\mathrm{[%s,%s]\,\,Hamaker\,coeff:\, %s \,in\,water,\,retarded}$'%(cnt1,cnt2,orientation)
#
#def svfig(cnt1,cnt2,orientation):
#	return 'plots/140322_%sw%s_HCs_%s_ret.pdf'%(cnt1,cnt2,orientation)
#
#pl.figure()
#pl.loglog(1e9*Ls,1e21*sum_A0,      label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[2],1e21*sum_A0[2]))
#pl.loglog(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[2],1e21*sum_A2[2]))
#pl.xlabel(x_ax)
#pl.ylabel(y_ax_per)
#pl.legend(loc = 'best')
#pl.title(r'$\mathrm{[5,1]\,and\,[29,0]\,\,Hamaker\,coeff:\, %s \,in\,water,\,retarded}$')
#pl.savefig(svfig('51','290','perpendicular'))
#pl.show()
#
#pl.figure()
#pl.plot(1e9*Ls[:98],np.diff(1e21*sum_A0),      label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[34],1e21*sum_A0[34]))
#pl.plot(1e9*Ls[:98],np.diff(1e21*sum_A2),label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[34],1e21*sum_A2[34]))
#pl.plot(1e9*Ls[:98],Ls[:98]-Ls[:98], 'k:')
#pl.legend(loc = 'best')
#pl.title(r'$d\mathcal{A}^{(0)}/d\ell,\,\,d\mathcal{A}^{(2)}/d\ell$ for [5,1] and [29,0] in water')
#pl.ylabel(r'$d\mathcal{A}^{(0)}/d\ell,\,\mathcal{A}^{(2)}/d\ell\,\,\,\,[zJ]$')
#pl.xlabel(r'$\ell\,\,\,[nm]$')
#pl.savefig('plots/140522_Hopkins_51w290_dAdl_vs_l.pdf')
#pl.show()
#
#pl.figure()
#pl.plot(1e9*Ls[:98],np.diff(1e21*(sum_A0+sum_A2)),
#        label=r'$\mathcal{A}^{(0)}+\mathcal{A}^{(2)}(\ell=%1.1f nm)=%3.2f, \,\,\,\mathcal{A}^{(0)}+\mathcal{A}^{(2)}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*(sum_A0[0]+sum_A2[0]),1e9*Ls[34],1e21*(sum_A0[34]+sum_A2[34])))
#pl.plot(1e9*Ls[:98],Ls[:98]-Ls[:98], 'k:')
#pl.legend(loc = 'best')
#pl.title(r'$d(\mathcal{A}^{(0)}+\mathcal{A}^{(2)})/d\ell$ for [5,1] and [29,0] in water')
#pl.ylabel(r'$d(\mathcal{A}^{(0)}+\mathcal{A}^{(2)})/d\ell\,\,\,\,[zJ]$')
#pl.xlabel(r'$\ell\,\,\,[nm]$')
#pl.savefig('plots/140522_Hopkins_51w290_dA0plusA2dl_vs_l.pdf')
#pl.show()
#
#
#
#
#
#pl.figure()
#pl.semilogy(1e9*Ls,1e21*sum_A0,label=      r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f,\,\,\, \mathcal{A^{(0)}}(\ell=%1.1fnm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[2],1e21*sum_A0[2]))
#pl.semilogy(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1fnm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1fnm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[2],1e21*sum_A2[2]))
#pl.xlabel(x_ax)
#pl.ylabel(y_ax_per)
##pl.title(title('6','5','perpendicular'))
#pl.legend(loc = 'best')
#pl.minorticks_on()
#pl.ticklabel_format(axis = 'both')
#pl.grid(which = 'both')
#pl.tick_params(which = 'both',labelright = True)
#pl.legend(loc = 'best')
##pl.title(title('6','5','Semi-log, perpendicular'))
##pl.title(title('9','0','perpendicular'))
##pl.title(title('9','1','perpendicular'))
##pl.title(title('9','3','perpendicular'))
##pl.title(title('29','0','perpendicular'))
##
#pl.savefig(svfig('51','290','semilog_perpendicular'))
##pl.savefig(svfig('90','90','perpendicular'))
##pl.savefig(svfig('91','91','perpendicular'))
##pl.savefig(svfig('93','93','perpendicular'))
##pl.savefig(svfig('290','290','perpendicular'))
#pl.show()
#
#


## Define functions
#def Aiz(perp, par,med):
#	return (2.0*(perp-med)*med)/((perp+med)*(par-med))
#def Delta(par,med):
#	return (par - med)/med
#def Delta_perp(perp,med):
#	return (perp - med)/(perp+med)
#def DDDg2(perp1, perp2, par1,par2, med):
#    return (2./med)*((perp1-med)*(perp2-med)/(2.*med) -((perp1+med)*(par1-med)*(perp2-med)-(perp2+med)*(par2-med)*(perp1-med)+2.*med*(perp1-med)*(perp2-med))/((perp2+med)*(perp1+med)))
#def DDg2(perp1, perp2, par1,par2, med):
#    return (1./(med*med))*(par2-med)*(par1-med)*\
#            ((par1-med)/med-2.*(perp1-med)/(perp1+med))*\
#            ((par2-med)/med -2.*(perp2-med)/(perp2+med))
#def Pn(e,zn,l):
#	return np.sqrt(e)*zn*l*(1./c)
#
##Integration vars
#T  = np.linspace(0.,2.**17, 1.+2.**17)
#U  = np.linspace(0.,2.**17, 1.+2.**17)
#
## Define functions
#def Aiz(perp, par,med):
#	return (2.0*(perp-med)*med)/((perp+med)*(par-med))
#def Delta(par,med):
#	return (par - med)/med
#def Delta_perp(perp,med):
#	return (perp - med)/(perp+med)
#def DDDg2(perp1, perp2, par1,par2, med):
#    return (2./med)*((perp1-med)*(perp2-med)/(2.*med) -((perp1+med)*(par1-med)*(perp2-med)-(perp2+med)*(par2-med)*(perp1-med)+2.*med*(perp1-med)*(perp2-med))/((perp2+med)*(perp1+med)))
#def DDg2(perp1, perp2, par1,par2, med):
#    return (1./(med*med))*(par2-med)*(par1-med)*\
#            ((par1-med)/med-2.*(perp1-med)/(perp1+med))*\
#            ((par2-med)/med -2.*(perp2-med)/(perp2+med))
#def Pn(e,zn,l):
#	return np.sqrt(e)*zn*l*(1./c)
