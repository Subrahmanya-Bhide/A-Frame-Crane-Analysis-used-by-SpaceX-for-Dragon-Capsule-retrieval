# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 16:35:18 2020

@author: Lenovo
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

d_rope = 44e-3 #diameter of the rope
Area = 0.25*(np.pi)*d_rope**2 #cross-section are of the rope
E = 5e+9 #Modulus of Elasticity of Nylon
tstop = 1000 #stop_time
tinc = 0.1 #time_increment
Widthby2 = 3.81 #Half the width of the A-Frame
Length = 7.84 # Length of the pendulum string
Mass= 2394.12 # Mass of Dragon Capsule
g= 9.81 #Acceleration due to gravity
B=g/Length 
theta0x = np.radians(5.0) #initial angular displacement in x plane
theta0y = np.radians(5.0) #initial angular displacement in y plane

# function for solving set of 1st order differential equations
def solve(tstop,tinc,initialconditions,params,func,E):
    A,B    = params
    ATOL,RTOL = E
    t = np.arange(0., tstop, tinc)
    solution =odeint(func, initialconditions, t, args=(params,),rtol=RTOL,atol=ATOL)
    return solution

def f(yx, t, params):
    thetax, omegax,u, L= yx      # unpack current values of yx
    A , B = params          # unpack parameters
    derivsx = [omegax,       # list of dy/dt=f functions
             (-omegax*A - B*(np.sin(thetax))),
              (Mass*g*np.cos(thetax)  +  Mass*(omegax**2)*L)*(1/Area*E),
              (Mass*g*np.cos(thetax)  +  Mass*(omegax**2)*L)*(1/Area*E)
              ]
    return derivsx

Controlparams = [1.49012e-8,1.49012e-8]
parameters = [0.0,B] #[A,B]
y0x = [theta0x,0.49,0,Length] #[theta0x,omega0x]
psolnx = solve(tstop,tinc,y0x,parameters,f,Controlparams)

t = np.arange(0., tstop, tinc)

fig = plt.figure(1, figsize=(8,8))
ax3 = fig.add_subplot(311)
ax3.plot(t, psolnx[:,2],'r-')
ax3.set_xlabel('time (s)',size='16')
ax3.set_ylabel('Energy Difference',size='16')