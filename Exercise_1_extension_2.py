# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 22:10:15 2020

@author: Lenovo
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import Exercise_1_Basic_Problem as p 

d_rope = 44e-3 #diameter of the rope
A = 0.25*(np.pi)*d_rope**2 #cross-section are of the rope
E = 5e+9 #Modulus of Elasticity of Nylon
tstop = 100 #stop_time
tinc = 0.1 #time_increment
Widthby2 = 3.81 #Half the width of the A-Frame
Length = 8.624 # Initial Length of the pendulum string
Mass= 2394.12 # Mass of Dragon Capsule
g= 9.81 #Acceleration due to gravity
Alpha = 100 #Initial angular acceleration
tlift = 40 #Time taken to lift the crane
B= 0.784/tlift # rate of change of length with time
theta0x = np.radians(5.0) #initial angular displacement in x plane
theta0y = np.radians(5.0) #initial angular displacement in y plane
phi0 = np.radians(10.0)  #initial Angle of crane
phidot0 = np.radians(80/tlift) #initial angular velocity of the crane in case of constant angular velocity
phidot0 = Alpha*tlift #initial angular velocity of the crane in case of constant angular acceleration

# function for solving pair of 1st order differential equations
def solve(tstop,tinc,initialconditions,params,func,E):
    A,B,Alpha    = params
    ATOL,RTOL = E
    t = np.arange(0., tstop, tinc)
    
    solution =odeint(func, initialconditions, t, args=(params,),rtol=RTOL,atol=ATOL)
    Alist=[]
    for i in range(0,int(tstop/tinc)):
        if (1.57 - solution[i,3]) > 0 :
            Alist.append(i)
        else:
            break
    n = len(Alist)
    n0 = Alist[n-1]
    return [solution,n0]  #returns solution array and the index for which the crane angle phi is just less than 90

#defining array of functions for solving motion in x direction
def f(yx, t, params):
    thetax, omegax, phidot, phi, Length = yx     # unpack current values of yx
    A, B, Alpha= params                          # unpack parameters
    derivsx = [omegax,                           # list of dy/dt=f functions
             -omegax*A - (9.81/Length)*(np.sin(thetax)),
               Alpha,
               phidot,
               -B]
    return derivsx
Controlparams = [1.49012e-8,1.49012e-8]
parameters = [0.0,B,Alpha] #[A,B,Alpha]
y0x = [theta0x,0.49,phidot0,phi0,Length] #initial conditions
psolnx = solve(tstop,tinc,y0x,parameters,f,Controlparams)[0]
n0 = solve(tstop,tinc,y0x,parameters,f,Controlparams)[1]


def energy(W,R,Mass,Length):
    z=R[:,2]
    pe = np.zeros(len(z))
    for i in range(0,len(z)):
        pe[i] = (Mass*g*(Length + z[i]))      # Potential Energy
    
    ke = np.zeros(len(z))
    for i in range(0,len(z)):
        ke[i] = 0.5*Mass*(Length**2)*(W[i]**2)  #Kinetic Energy
        
    energy = np.array([ke,pe])
    return energy

def vectorize(Angles,Length,unit_vector):
    pos = np.zeros((len(Angles),3)) 
    if unit_vector == 0 :                    # for x plane
        for n in range(0,len(Angles)):
            pos[n,0] = Length*np.sin(Angles[n])
            pos[n,1] = 0.0
            pos[n,2] = np.negative(Length*(np.cos(Angles[n])))
    else :                                   #for y plane
        for m in range(0,len(Angles)):
            pos[m,0] = 0.0
            pos[m,1] = Length*np.sin(Angles[m])
            pos[m,2] = np.negative(Length*(np.cos(Angles[m])))
    return pos



#defining array of functions for solving motion in y direction
def f(yy, t, params):
    thetay, omegay, phidot, phi, Length = yy      # unpack current values of yy
    A, B, Alpha= params                           # unpack parameters
    derivsy = [omegay,                            # list of dy/dt=f functions
             -omegay*A - (9.81/Length)*(np.sin(thetay)),
               Alpha,
               phidot,
               -B]
    return derivsy
Controlparams = [1.49012e-8,1.49012e-8]
parameters = [0.0,B,Alpha] #[A,B,Alpha]
y0y = [theta0y,0.49-0.035,phidot0,phi0,Length] #initial conditions
psolny = solve(tstop,tinc,y0y,parameters,f,Controlparams)[0]

t = np.arange(0., tstop, tinc)

#A = p.psolnx[:,0] #angular displacements in x plane without change in length of string 
#E = p.psolny[:,0] #angular displacements in y plane without change in length of string 
C = psolnx[:,0] #angular displacements in x plane for changing length of string
D = psolny[:,0] #angular displacements in y plane without changing length of string

posvectX = vectorize(psolnx[:,0],Length,0) #Position Vectors at different times for X-Plane oscillations
posvectY = vectorize(psolny[:,0],Length,1)

EnergyX = energy(psolnx[:,1],posvectX,Mass,Length)  #Energy array for X-plane oscillations
EnergyY = energy(psolny[:,1],posvectY,Mass,Length)  #Energy array for Y-plane oscillations

TotalEX=np.array(EnergyX[0,:] + EnergyX[1,:] ) #Total Energy for X Plane oscillations at different times
TotalEY=np.array(EnergyY[0,:] + EnergyY[1,:] )


#Plotting relevant graphs

fig = plt.figure(1, figsize=(8,8))
ax1 = fig.add_subplot(311)
ax1.plot(t,TotalEX,'r',label = 'Constant Length X')
ax1.set_xlabel('time (s)',size='16')
ax1.set_ylabel('Theta X',size='16')
ax1.legend(loc='best')

'''ax5 = fig.add_subplot(311)
ax5.plot(t,E,'b',label = 'Constant Length Y')
ax5.plot(t,D,'m',label = 'Changing Length Y')
ax5.set_xlabel('time (s)',size='16')
ax5.set_ylabel('Theta Y ',size='16')
ax5.legend(loc='best')
plt.ylim(-0.6,1.0)'''





