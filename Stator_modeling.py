# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 17:46:41 2018

Build Mass matrix
"""
import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from sympy import *
import sys
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
init_printing(use_unicode=False, wrap_line=False, no_global=True)

pi = np.pi
eps0 = 8.85e-12
itom = 2.54e-2

# Material Properties
"rhos: Circular plate density [kg/m3]"
rhos = 7860
"Es: Circular plate Stiffness [N/m2]"
Es = 210e9
"nus: Circular plate poisson's ratio []"
nus = 0.35
"Qs: Compliance matrix for circular plate"
c11s = Es/(1-nus**2)
c12s = nus*Es/(1-nus**2)
c66s = Es/(1+nus)/2
Qs = np.zeros((3,3))
Qs[0,0]=c11s;   Qs[1,1]=c11s;   Qs[0,1]=c12s;   Qs[1,0]=c12s;   Qs[2,2]=c66s


"rhop: Piezo plate density [kg/m3]"
rhop = 7550
"Ep: Piezo plate Stiffness [N/m2]"
Ep = 80e9
"Piezo plate poisson's ratio []"
nup = 0.35
"Qs: Compliance matrix for circular plate"
c11p = Ep/(1-nup**2)
c12p = nup*Ep/(1-nup**2)
c66p = Ep/(1+nup)/2
Qp = np.zeros((3,3))
Qp[0,0]=c11p;   Qp[1,1]=c11p;   Qp[0,1]=c12p;   Qp[1,0]=c12p;   Qp[2,2]=c66p
"d31: Electromechanical coupling factor [m/v]"
d31 = -127e-12
"d: Electromechanical Matrix"
d = np.zeros((1,3))
d[0,0] = d31; d[0,1] = d31
"epst: Dielctric contant under constant stress"
epst = 1550*eps0
"e: Dielectric contant matrix"
e = np.matmul(d,Qp)
"epss: Dielectric contant under constant strain"
epss = epst - np.linalg.norm(np.matmul(e,np.transpose(d)))

# Geometrym
"a: radius of hole in circular plate [m]"
a = 0.5*15e-3
"b: radius of circular plate [m]"
b = 0.5*60e-3
"ap: radius of hole in piezo plate [m]"
ap = 0.5*45e-3
"bp: radius of piezo plate [m]"
bp = 0.5*60e-3
"ts: thickness of circular plate [m]"
ts = 1.5e-3
"tp: thickness of piezo plate [m]"
tp = 0.5e-3
"h: height of composite plate [m]"
h = tp + ts
"ht: tooth height [m]"
ht = 0.5e-3
"dt: distance between teeth [m]"
dt = 0.5e-3
"nt: number of teeth"
nt = 72
"r1: inner edge of teeth [m]"
r1 = 0.5*50e-3
"r2: outer edge of teeth [m]"
r2 = b
"rt: middle of teeth [m]"
rt = (r1+r2)/2
"zt: tooth center of mass height [m]                        Check"
zt = ts/2 + ht*3/4
"mt: tooth mass [kg]                                        Check"
mt = (pi*(r2**2 - r1**2)-nt*(r2-r1)*dt)*ht*rhos/nt
"Vibration wave number:"
N = 9
'''
# Driving Singal
"Amp: Voltage amplitude [V]"
Amp = 50
"drivefreq: Driving Signal Frequency [Hz]"
drivefreq = 1000
"tau: period of driving signal [s]"
tau = 1/drivefreq
"omega: Driving Signal Frequency [rad/s]"
omega = 2*pi*drivefreq
"Omega: Driving Frequency of Traveling Wave [rad/s]"
Omega = omega/N
"Tau: period of traveling wave [s]"
Tau = N*tau
"Va: Signal to phase A [V]"
Va = Amp*np.cos(omega*t)
"Vb: Signal to phase B [V]"
Vb = Amp*np.sin(omega*t)
"Ea: Electric Field to phase A [V/m]"
Ea = Va/tp
"Eb: Electric Field to phase B [V/m]"
Eb = Vb/tp
'''
"Capacitance Matrix"
con = pi*(bp**2-ap**2)/tp*epss;
Cap = np.array([[con,0],[0,con]])

"Vibration Mode shapes: "

#Centerline Displacement
r = Symbol('r')
theta=Symbol('theta')
z = Symbol('z')
U0 = Matrix(1,2,[0,0])
V0 = Matrix(1,2,[0,0])
W0 =  Matrix(1,2,[((r-a)/(b-a))**2*cos(N*theta), ((r-a)/(b-a))**2*sin(N*theta)])

"Stiffness Matrix: "
S0_rad = diff(U0, r)
S0_tan = U0/r + diff(V0,theta)/r
S0_shear = diff(U0,theta)/r + diff(V0,r) - V0/r

K_rad = -diff(W0,r,r)
K_tan = -diff(W0,r)/r - diff(W0,theta,theta)/r**2
K_shear = -2*diff(W0,r,theta)/r + 2*diff(W0,theta)/r**2

S0 = BlockMatrix(3,1,[S0_rad,S0_tan,S0_shear])
k = BlockMatrix(3,1,[K_rad,K_tan,K_shear])

Epsilon = BlockMatrix(2,1,[S0,k])

#Stress Strain Relations: 
#Qp = Plane stress piezo stiffness Matrix
Qp = Matrix(3,3,[c11p,c12p,0,c12p,c11p,0,0,0,c66p])
#Qs = Plane stress substrate stifffness matrix
Qs = Matrix(3,3,[c11s,c12s,0,c12s,c11s,0,0,0,c66s])

#A = Extensional Stiffness Matrix
Ap = integrate(Qp, (z, -h/2, -h/4))     # limita are based on the plate thickness assuming neutral axis at h/2
As = integrate(Qs, (z, -h/4, h/2))
A = MatAdd(Ap,As)

#B = Coupling Stiffness Matrix
Bp = integrate(z*Qp, (z, -h/2, -h/4))
Bs = integrate(z*Qs, (z, -h/4, h/2))
B = MatAdd(Bp,Bs)

#D = Bending Stiffness Matrix
Dp = integrate(z**2*Qp, (z, -h/2, -h/4))
Ds = integrate(z**2*Qs, (z, -h/4, h/2))
DD = MatAdd(Dp,Ds)

Kappa = BlockMatrix(2,2,[A,B,B,DD])
Kappam = BlockMatrix(2,2,[Ap,Bp,Bp,Dp])

Kkk = MatMul(Epsilon.transpose(),Kappa,Epsilon)
Kkkm = MatMul(Epsilon.transpose(),Kappam,Epsilon)

Kk00 = integrate(nsimplify(r*Kkk[0,0]), (r, ap, bp))
Kkm00 = integrate(nsimplify(r*Kkkm[0,0]), (r, a, ap)) + integrate(r*Kkkm[0,0], (r, bp, b))
Kk01 = integrate(nsimplify(r*Kkk[0,1]), (r, ap, bp))
Kkm01 = integrate(nsimplify(r*Kkkm[0,1]), (r, a, ap)) + integrate(r*Kkkm[0,1], (r, bp, b))
Kk11 = integrate(nsimplify(r*Kkk[1,1]), (r, ap, bp))
Kkm11 = integrate(nsimplify(r*Kkkm[1,1]), (r, a, ap)) + integrate(r*Kkkm[1,1], (r, bp, b))
K = integrate(Matrix(2,2,[Kk00+Kkm00,Kk01+Kkm01,Kk01+Kkm01,Kk11+Kkm11]), (theta, 0, 2*pi)).evalf()
K = np.array(K).astype(np.float64)




"Mass Matrix: "
W = W0.T*W0

#Thickness
Ms = ts*rhos*W
Mp = tp*rhop*W
#Radial
Mmm = integrate(r*Ms, (r, a, b)) + integrate(r*Mp, (r, ap, bp))
#tangential
Mm = integrate(Mmm, (theta, 0, 2*pi)).evalf()

#Add Teeth inertia and mass
thetat=pi/2/N
MW1 = nt*mt*W.subs([(r,rt),(theta,0)])
MW2 = nt*mt*W.subs([(r,rt),(theta,thetat)])
Mtw = MW1 + MW2

wor = diff(W0, r)
WOR = wor.T*wor
wot = diff(W0/r, theta)
WOT = wot.T*wot

Mu1 = nt*mt*zt**2*WOR.subs([(r,rt),(theta,0)])
Mu2 = nt*mt*zt**2*WOR.subs([(r,rt),(theta,thetat)])
Mtu = Mu1 +Mu2

Mv1 = nt*mt*zt**2*WOT.subs([(r,rt),(theta,0)])
Mv2 = nt*mt*zt**2*WOT.subs([(r,rt),(theta,thetat)])
Mtv = Mv1 + Mv2

Mt = Mtw + Mtu + Mtv

M = (Mm + Mt).evalf()
M = np.array(M).astype(np.float64)
Minv = np.linalg.inv(M)

evals, evecs = eigh(K,M)
Nat_freqs = np.sqrt(evals)/2/pi
Mode_shapes = evecs

"Damping Matrix"
zeta = 0.0064
C = 4*pi*Nat_freqs[0]*M
"Electromechanical Coupling Factor Matrix: "
coeff = [0.5,0.5]
c = coeff[0]*((r-ap)/(bp-ap))**2
d = coeff[1]*((r-ap)/(bp-ap))**3

R = c + d
phir = Matrix(1,2,[R*cos(N*theta),R*sin(N*theta)])
s_rad = -z*diff(phir,r,r)
s_tan = -z/r*(diff(phir,r)+1/r*diff(phir,theta,theta))
s_sh = -2*z/r*(diff(phir,r,theta)-1/r*diff(phir,theta))

Nr = BlockMatrix(3,1,[s_rad,s_tan,s_sh])
d = Matrix(3,3,[0,0,0,0,0,0,d31,d31,0])
et = (d*Qp).T

Nv = Matrix(3,2,[0,0,0,0,N/tp,N/tp])

arg0 = Nr.T*et*Nv
arg1a = integrate(nsimplify(arg0[0,0]),(z, -h/2, -h/4))
arg2a = integrate(nsimplify(arg1a),(r,ap,bp))
arga = integrate(arg2a,(theta,-pi/2/N,pi/2/N))

#arg1b = integrate(nsimplify(arg0[1,1]),(z, -h/2, -h/4))
#arg2b = integrate(nsimplify(arg1b),(r,ap,bp))
#argb = integrate(arg2b,(theta,0,pi/N))

Theta = Matrix(2,2,[arga,0,0,arga])
Theta = np.array(Theta).astype(np.float64)
print(Theta)
print(Nat_freqs)

'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Create the mesh in polar coordinates and compute corresponding vibration amplitude.
r = np.linspace(a, b, 100)
p = np.linspace(0, 2*np.pi, 10000)
R, P = np.meshgrid(r, p)
Z = ((R-a)/(b-a))**3*np.cos(4*P)

# Express the mesh in the cartesian system.
X, Y = R*np.cos(P), R*np.sin(P)

# Plot the surface.
ax.plot_surface(X, Y, Z, cmap=plt.cm.YlGnBu_r)

# Tweak the limits and add latex math labels.
#ax.set_zlim(0, 1)
ax.set_xlabel(r'$\phi_\mathrm{real}$')
ax.set_ylabel(r'$\phi_\mathrm{im}$')
ax.set_zlabel(r'$V(\phi)$')

plt.show()
'''
height = h/2 + ht 

coef = [0.5,0.5]
defl1 = coef[0]*((rt-a)/(b-a))**3
defl2 = coef[1]*((rt-a)/(b-a))**3
defl = defl1 + defl2

"Interface:"
#Friction Coefficienct
mu = 0.2

"Rotor: "
#Inertia
Ir = 6.4e-8
#Inertial Damping
Idamp = 2e-5
#Mass
Mr=30e-3
#Mass Damping
Mdamp = 500