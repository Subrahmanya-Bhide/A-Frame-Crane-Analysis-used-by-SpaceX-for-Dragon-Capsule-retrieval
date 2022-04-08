# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 20:24:49 2020

@author: Lenovo
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

d_rope = 44e-3 #diameter of the rope
A = 0.25*(np.pi)*d_rope**2 #cross-section are of the rope
E = 5e+9 #Modulus of Elasticity of Nylon
tstop = 100 #stop_time
tinc = 0.1 #time_increment
Widthby2 = 3.81 #Half the width of the A-Frame
Length = 7.84 # Length of the pendulum string
Mass= 2394.12 # Mass of Dragon Capsule
g= 9.81 #Acceleration due to gravity
B=g/Length 
theta0 = np.radians(5.0) #initial angular displacement 
omega0 = 1.4 #initial angular displacement
# function for solving pair of 1st order differential equations
def solve(tstop,tinc,initialconditions,params,func,E):
    A,B    = params
    ATOL,RTOL = E
    t = np.arange(0., tstop, tinc)
    
    solution =odeint(func, initialconditions, t, args=(params,),rtol=RTOL,atol=ATOL)
    return solution

#Function for getting the position vector from angular displacement
def vectorize(Angles,Length,unit_vector):
    pos = np.zeros((len(Angles),3)) 
    if unit_vector == 0 :
        for n in range(0,len(Angles)):
            pos[n,0] = Length*np.sin(Angles[n])
            pos[n,1] = 0.0
            pos[n,2] = np.negative(Length*(np.cos(Angles[n])))
    else :
        for m in range(0,len(Angles)):
            pos[m,0] = 0.0
            pos[m,1] = Length*np.sin(Angles[m])
            pos[m,2] = np.negative(Length*(np.cos(Angles[m])))
    return pos

t = np.arange(0., tstop, tinc)

#Code for finding variation of maximum angular displacement for different initial angular velocities 
'''Omegas = np.linspace(0,3,800)
maxthetas = np.zeros(len(Omegas))

for i in range(0,len(Omegas)):
    def f(y, t, params):
        theta, omega = y      # unpack current values of y
        A , B = params        # unpack parameters
        derivs = [omega,      # list of dy/dt=f functions
             -omega*A - B*(np.sin(theta))]
        return derivs
    
    Controlparams = [1.49012e-8,1.49012e-8]
    parameters = [0.0,B] #[A,B]
    y0 = [theta0,Omegas[i]] #initial conditions
    psoln = solve(tstop,tinc,y0,parameters,f,Controlparams)
    maxthetas[i] = max(psolnx[:,0])
'''    
# Code for solving for 3 different initial angular velocities 
Omegas = np.array([1.4,1.9,2.4])
thetas = np.zeros([len(t),len(Omegas)])

for i in range(0,len(Omegas)):
    def f(y, t, params):
        theta, omega = y      # unpack current values of y
        A , B = params        # unpack parameters
        derivs = [omega,      # list of dy/dt=f functions
             -omega*A - B*(np.sin(theta))]
        return derivs
    
    Controlparams = [1.49012e-8,1.49012e-8]
    parameters = [0.0,B] #[A,B]
    y0 = [theta0,Omegas[i]] #[theta0x,omega0x]
    psoln = solve(tstop,tinc,y0,parameters,f,Controlparams)
    for j in range(len(t)):
        thetas[j,i] = psoln[j,0]


#Plotting relevant graphs
#fig = plt.figure(1, figsize=(8,8))

'''ax1 = fig.add_subplot(311)
ax1.plot(Omegas[0:595], maxthetas[0:595])
ax1.set_xlabel('Initial Angular Velocity',size='16')
ax1.set_ylabel('Max Angular Displacement',size='16')

ax2 = fig.add_subplot(312)
ax2.plot(Omegas[596:800], maxthetas[596:800])
ax2.set_xlabel('Initial Angular Velocity',size='16')
ax2.set_ylabel('Max Angular Displacement',size='16')

ax3 = fig.add_subplot(313)
ax3.plot(Omegas, maxthetas)
ax3.set_xlabel('Initial Angular Velocity',size='16')
ax3.set_ylabel('Max Angular Displacement',size='16')'''

'''ax4 = fig.add_subplot(311)
ax4.plot(t, thetas[:,0],'r-',label='$\omega0 = 1.4 rad/s$ ')
ax4.plot(t, thetas[:,1],'b-',label='$\omega0 = 1.9 rad/s$ ')
ax4.plot(t, thetas[:,2],'g-',label='$\omega0 = 2.4 rad/s$ ')
ax4.set_xlabel('time (s)',size='16')
ax4.set_ylabel('Angular Displacement',size='16')
ax4.legend(loc="best")
plt.ylim(-2.2,4.0)'''


'''plt.tight_layout()
plt.show()'''