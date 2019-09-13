# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 15:03:52 2019

@author: MoritaLab_Student
"""
import numpy as np
import matplotlib.pyplot as plt
wf = sol.y[6][-1]
w = defl*np.matmul(np.array([np.cos(N*theta),np.sin(N*theta)]).T,np.array([[sol.y[0][-1]],[sol.y[1][-1]]]).reshape((2,)))
vel = -N/rt*height*defl*np.matmul(np.array([np.sin(N*theta),-np.cos(N*theta)]).T,np.array([[sol.y[2][-1]],[sol.y[3][-1]]]).reshape((2,1)))
j = [i for i in range(len(w)) if w[i] >= (wf-height)]


force = np.zeros(len(theta))
kr = 1.5e7
mu = 0.15
if len(j)>0:
    one1=np.ones(len(j))
    force[j] = kr*(w[j] + one1*(height - wf))
        
# equivalent but more general
plt.figure(1)
ax1=plt.subplot(2, 1, 1)
ax1.plot(theta,w)
ax1.grid(True)
ax1.hold(True)
rotor_height = (wf-height)*np.ones((len(theta)))
ax1.plot(theta,rotor_height,'y--')
ax1.plot(theta[j],w[j],'ro')
ax1.set_title('Interference Region')
ax1.set_xlabel('theta [rad]')
ax1.set_ylabel('Deflection [um]')
ax1.hold(False)

ax2=plt.subplot(2, 1, 2)
ax2.plot(theta,force)
ax2.grid(True)
ax2.hold(True)
ax2.plot(theta[j],force[j],'ro')
#ax2.set_title('Force Distribution')
#ax2.set_xlabel('theta [rad]')
ax2.set_ylabel('Force [N]')
ax2.hold(False)

#########################################
vr = rt*sol.y[5][-1]
k = [i for i in range(len(vel)) if vel[i] > vr]
l = [i for i in range(len(vel)) if vel[i] < vr]
speed = np.zeros((len(theta),1))
tanforce = np.zeros(len(theta))
if (len(k) > 0) or (len(l) > 0):
    tanforce[k] = mu*force[k]
    tanforce[l] = -mu*force[l]
        
pos = speed.copy()
pos[k]=vel[k]
pos_bool = (pos!=0).reshape(len(pos))  
force_bool = (force!=0).reshape(len(force))  
pos_force_bool = np.logical_and(pos_bool, force_bool)
pos = [i for i in range(len(theta)) if pos_force_bool[i] == True]

neg = speed.copy()
neg[l]=vel[l]
neg_bool = (neg!=0).reshape(len(neg))  
neg_force_bool = np.logical_and(neg_bool, force_bool)
neg = [i for i in range(len(theta)) if neg_force_bool[i] == True]

plt.figure(2)
ax3=plt.subplot(2, 1, 1)
ax3.plot(theta[j],vel[j],'wo')
ax3.grid(True)
ax3.hold(True)

Vr = sol.y[5][-1]*np.ones(len(theta))
ax3.plot(theta,Vr,'--')
ax3.plot(theta[pos],vel[pos]/rt,'g+')
ax3.plot(theta[neg],vel[neg]/rt,'rx')
ax3.set_title('Speed Profile')
ax3.set_xlabel('theta [rad]')
ax3.set_ylabel('Speed [rad/s]')
ax3.hold(False)

ax4=plt.subplot(2, 1, 2)
ax4.plot(theta,tanforce,'y--')
ax4.grid(True)
ax4.hold(True)
ax4.plot(theta[j],tanforce[j],'wo')
ax4.plot(theta[pos],tanforce[pos],'g+')
ax4.plot(theta[neg],tanforce[neg],'rx')
#ax4.set_title('Speed Profile')
ax4.set_xlabel('theta [rad]')
ax4.set_ylabel('Driving Force [N]')
ax4.hold(False)




