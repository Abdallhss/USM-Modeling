# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 17:46:41 2018

Build Mass matrix
"""
import numpy as np
from scipy.integrate import solve_ivp, odeint, RK23
import matplotlib.pyplot as plt

def eom(t,x):
    "Displacement and Velocity Profiles"
    w = defl*np.matmul(np.array([np.cos(N*theta),np.sin(N*theta)]).T,np.array([[x[0]],[x[1]]]).reshape((2,)))
    # wave number is 4
    vel = -N/rt*height*defl*np.matmul(np.array([np.sin(N*theta),-np.cos(N*theta)]).T,np.array([[x[2]],[x[3]]]).reshape((2,1)))
    "Rotor Stiffness (Can replace with Hertzian): "
    #kr = -2.59e21*wf^2 + 5.38e16*wf + 1.23e11;
    kr = 1e8
    "Rotor speed and Flexure height"
    wr = x[5];
    vr = rt*wr;
    wf = x[6];
    
    # Contact Points
    j = [i for i in range(len(w)) if w[i] >= (wf-height)]
    "Normal Forces"
    force = np.zeros(len(theta))
    if len(j)>0:
        one1=np.ones(len(j))
        force[j] = kr*(w[j] + one1*(height - wf))
        normwork1 = -N*dt*np.sum(defl*np.multiply(np.cos(N*theta),force))
        normwork2 = -N*dt*np.sum(defl*np.multiply(np.sin(N*theta),force))
    else:
        normwork1 = 0
        normwork2 = 0
    
    "Speed Differential: "
    #Friction Region (Dynamic Friction)
    k = [i for i in range(len(vel)) if vel[i] > vr]
    l = [i for i in range(len(vel)) if vel[i] < vr]

    
    "Tangential Forces"
    tanforce = np.zeros(len(theta))
    if (len(k) > 0) or (len(l) > 0):
        tanforce[k] = mu*force[k]
        tanforce[l] = -mu*force[l]
        tanworkl = +N*dt*np.sum(N*height/rt*defl*np.multiply(np.sin(N*theta),tanforce))
        tanwork2 = -N*dt*np.sum(N*height/rt*defl*np.multiply(np.cos(N*theta),tanforce))
    else:
        tanworkl = 0;
        tanwork2 = 0;
    
    "Interface Force and torque"
    normforce = N*dt*np.sum(force)
    torque = N*rt*dt*np.sum(tanforce);
    "Electric Drive: "
    #Driving Voltage
    v1 = Amp*np.cos(omega*t);
    v2 = Amp*np.sin(omega*t);
    
    "Interface Forces: "

    fv1 = Theta[0][0]*v1
    fv2 = Theta[1][1]*v2
    
    fnl = normwork1
    ft1 = tanworkl
    
    fn2 = normwork2
    ft2 = tanwork2
    
    "Total Forces"
    f1 = fv1 + fnl + ft1
    
    f2 = fv2 + fn2 + ft2
    
    f3 = torque - load_torque
    
    f4 = normforce - load_force
    
    f = np.array([f1,f2,f3,f4]).T
    f = f.reshape((4,))
    x = x.reshape((8,))
    xdot = np.add(np.matmul(A,x),np.matmul(B,f))
    return xdot


pi=np.pi

"Create State Matrices: "
I = np.eye(2)
O = 0*I

A1 = np.concatenate((O,I,O,O),axis=1)                                    # P's
A2 = np.concatenate((-np.matmul(Minv,K),-np.matmul(Minv,C),O,O),axis=1)  # pdot's            
A3 = np.array([[0,0,0,0,0,1,0,0],                                        # alpha
              [0,0,0,0,0,-Idamp/Ir,0,0],                                 # alphadot
              [0,0,0,0,0,0,0,1],                                         # wf
              [0,0,0,0,0,0,0,-Mdamp/Mr]])                                # wfdot              
A = np.concatenate((A1,A2,A3))

B = np.concatenate((np.concatenate((O,O),axis=1),
                    np.concatenate((Minv,O),axis=1),                                 # Piezos and interface forcing enter here
                    np.concatenate((O,np.array([[0,0],[1/Ir,0]])),axis=1),           # Torque here
                    np.concatenate((O,np.array([[0,0],[0,1/Mr]])),axis=1)),axis=0)   # Force here

numpoints = 250;
t1 = 0
t2= 2*pi/N
#width of element
dt = (t2-t1)/numpoints
theta = np.linspace(t1,t2-dt,numpoints)


"Driving Signal"
Amp = 200
# Driving Voltage Frequency
freq = 42000
tau = 1/freq
omega = 2*pi*freq

"Loads"
load_force = 200
load_torque = 0.1

# Stator Intial Conitions
'''
amp = 2.5e-6

ampdot = amp*omega;
phaseloss = 90*pi / 180    # % At resonance, get 90 deg phase loss

ampl = amp*np.cos(phaseloss)
amp2 = amp*np.sin(phaseloss)

ampdotl = -ampdot*np.sin(phaseloss)
ampdot2 = ampdot*np.cos(phaseloss)
'''
ampl = 0
amp2 = 0
ampdotl = 0
ampdot2 = 0

# Rotor angle I.C.'s
alpha = 0
alphadot = 0

# Flexure height I.C.'s
wf = height
wfdot = 0

"Simulation"
#Time
t0 = 0
tf = 10*tau        # Length of time simulation

#Initial Conditions
x0 = np.array([ampl,amp2,ampdotl,ampdot2,alpha,alphadot,wf,wfdot]).T
#Start ODE
#sol = odeint(eom, x0, [t0,tf])
sol = solve_ivp(eom, [t0,tf], x0,method='RK23',vectorized=True)
#sol = RK23(eom,t0, x0,tf)
t=sol.t
speed=sol.y[5]*60/2/pi

plt.plot(t,speed)


