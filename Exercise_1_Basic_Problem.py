# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 12:17:51 2020

@author: Subbu
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

d_rope = 44e-3 #diameter of the rope
A = 0.25*(np.pi)*d_rope**2 #cross-section are of the rope
E = 5e+9 #Modulus of Elasticity of Nylon
tstop = 20 #stop_time
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
    solution = odeint(func, initialconditions, t, args=(params,),rtol=RTOL,atol=ATOL)
    return solution

#Function for getting the position vector from angular displacement.
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

#Function for Finding out Energy of the Pendulum
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
  
#Function for finding out Tension in the rope    
def tension(pos,Mass,Length):
    theta = pos[:,0]
    w     = pos[:,1]
    tension = np.zeros(len(w))
    for i in range(0,len(w)):
        tension[i] = Mass*((9.81*np.cos(theta[i])) - Length*(w[i])**2 )
    return tension

#Function to find extension in the rope due to tension
def extension(T,A,Y):
    extend = np.zeros(len(T))
    for i in range(0,len(T)):
        extend[i] = T[i] * (1/(A*Y))*100
    return extend

#defining array of functions for solving motion in x direction and find critical angular velocity
def f_initial(yx, t, params):
    thetax, omegax = yx         # unpack current values of yx
    A , B = params              # unpack parameters
    derivsx = [omegax,          # list of dy/dt=f functions
             -omegax*A - B*(np.sin(thetax))]
    return derivsx

Omegas=np.arange(0.0,12.0,0.01)    #list of initial velocities to iterate through
Lengths=np.zeros(len(Omegas))
for i in range(0,len(Omegas)):
    Controlparams = [1.49012e-8,1.49012e-8]
    parameters = [0.0,B] #[A,B]
    y0x = [theta0x,Omegas[i]] #[theta0x,omega0x]
    psolnx_initial = solve(tstop,tinc,y0x,parameters,f_initial,Controlparams)
    a=vectorize(psolnx_initial[:,0],Length,0)
    a1=a[:,0]
    amax = np.max(a1)
    Lengths[i] = amax
Alist=[]
for i in range(0,len(Lengths)):        #condition for displacement being 
    if (Widthby2 - Lengths[i]) > 0 :   #lesser than half the width of the A-Frame
        Alist.append(i)
    else:
        break
n = len(Alist)
n0 = Alist[n-1]  #index for the critical angular velocity
 

#defining array of functions for solving motion in x direction
def f(yx, t, params):
    thetax, omegax = yx      # unpack current values of yx
    A , B = params           # unpack parameters
    derivsx = [omegax,       # list of dy/dt=f functions
             -omegax*A - B*(np.sin(thetax))]
    return derivsx

Controlparams = [1.49012e-8,1.49012e-8]
parameters = [0.0,B] #[A,B]
y0x = [theta0x,Omegas[n0]] #[theta0x,omega0x]
psolnx = solve(tstop,tinc,y0x,parameters,f,Controlparams)  #solution array in x-plane

#defining array of functions for solving motion in y direction
def f(yy, t, params):
    thetay, omegay = yy     # unpack current values of yy
    A , B = params          # unpack parameters
    derivsy = [omegay,      # list of dy/dt=f functions
             -omegay*A - B*(np.sin(thetay))]
    return derivsy

Controlparams = [1.49012e-8,1.49012e-8]
parameters = [0.0,B] #[A,B]
y0y = [theta0y,Omegas[n0] - 0.035] #[theta0y,omega0y]
psolny = solve(tstop,tinc,y0y,parameters,f,Controlparams) #solution array in y-plane

t = np.arange(0., tstop, tinc)

posvectX = vectorize(psolnx[:,0],Length,0) #Position Vectors at different times for X-Plane oscillations
posvectY = vectorize(psolny[:,0],Length,1) #Position Vectors at different times for Y-Plane oscillations

tensionX = tension(psolnx,Mass,Length)   #Tension in rope at different times for X-Plane oscillations
tensionY = tension(psolny,Mass,Length)   #Tension in rope at different times for Y-Plane oscillations
meanx =0.5*(max(tensionX) + min(tensionX)) #mean values
meany =0.5*(max(tensionY) + min(tensionY))
diffx = max(tensionX) - min(tensionX) #Ranges
diffy = max(tensionY) - min(tensionY)

extensionx=extension(tensionX,A,E)#Extension in rope at different times for X-Plane oscillations
extensiony=extension(tensionY,A,E)#Extension in rope at different times for Y-Plane oscillations

EnergyX = energy(psolnx[:,1],posvectX,Mass,Length)  #Energy array for X-plane oscillations
EnergyY = energy(psolny[:,1],posvectY,Mass,Length)  #Energy array for Y-plane oscillations

TotalEX=np.array(EnergyX[0,:] + EnergyX[1,:] ) #Total Energy for X Plane oscillations at different times
TotalEY=np.array(EnergyY[0,:] + EnergyY[1,:] ) #Total Energy for Y Plane oscillations at different times



teto = np.zeros(len(t))
for i in range(0,len(t)):
    teto[i] = TotalEX[0]

difftex = TotalEX -teto
#Plotting Relevant Graphs

#fig = plt.figure(1, figsize=(8,8))
'''ax3 = fig.add_subplot(311)
ax3.plot(t, difftex,'r-')
ax3.set_xlabel('time (s)',size='16')
ax3.set_ylabel('Energy Difference',size='16')
#plt.ylim(-0.6,1.0)'''

#ax3 = fig.add_subplot(311)
#ax3.plot(t, tensionX,'m-',label='Tension X')
#ax3.plot(t,tensionY,'g-',label='Tension Y')
#ax3.set_xlabel('time (s)',size='16')
#ax3.set_ylabel('Tension',size='16')
#ax3.legend(loc="best")
#plt.ylim(20000,33000)

'''ax4 = fig.add_subplot(311)
ax4.plot(t, EnergyY[0,:],'g-',label='Kinetic Energy')
ax4.plot(t,EnergyY[1,:],'y-',label='Potetial Energy')
ax4.set_xlabel('time (s)',size='16')
ax4.set_ylabel('Energy',size='16')
ax4.legend(loc="best")
plt.ylim(0,35000)'''

'''ax4 = fig.add_subplot(311)
ax4.plot(t, extensionx,'r-',label='Percent Extension in X-Plane')
ax4.plot(t,extensiony,'b-',label='Percent Extension in Y-Plane')
ax4.set_xlabel('time (s)',size='16')
ax4.set_ylabel('Percent Extension',size='16')
ax4.legend(loc="best")
plt.ylim(0.26,0.44)'''

'''plt.tight_layout()
plt.show()'''














    





    
    
    

