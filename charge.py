# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 15:03:52 2019

@author: MoritaLab_Student
"""
import numpy as np
import matplotlib.pyplot as plt
import operator

wr = sol.y[5]

# Input Voltage and its derivatives
Va = Amp*np.cos(omega*t)
Vadot = -omega*Amp*np.sin(omega*t)

Vb = Amp*np.sin(omega*t)
Vbdot = omega*Amp*np.cos(omega*t)

#Charge on electrodes
qa = Theta[0,0]*sol.y[0] + Cap[0,0]*Va;
qb = Theta[1,1]*sol.y[1] + Cap[1,1]*Vb;

#Current across electrodes
Ia = Theta[0,0]*sol.y[2] + Cap[0,0]*Vadot
Ib = Theta[1,1]*sol.y[3] + Cap[1,1]*Vbdot

#Instantaneous input power of both channels
Pa = np.multiply(Va,Ia)
Pb = np.multiply(Vb,Ib)

#Total instantaneous input power
pin = np.add(Pa,Pb)

#Total instantaneous output power
pout = load_torque*sol.y[5]

#Efficiency (assuming steady state at end of simulation)
eff = 100*pout[-1]/pin[-1]
print("Efficiency: ",eff)

I_a = Ia[int(len(Ia)/2):]
V_a = Vb[int(len(Vb)/2):]

I_A = np.fft.fft(I_a)
V_A = np.fft.fft(V_a)

idx_I, mag_I = max(enumerate(np.abs(I_A)), key=operator.itemgetter(1))
idx_V, mag_V = max(enumerate(np.abs(V_A)), key=operator.itemgetter(1))

pI = np.angle(I_A[idx_I])
pV = np.angle(V_A[idx_V])
phase_diff = (pI - pV)*180/pi

current = max(I_a)




