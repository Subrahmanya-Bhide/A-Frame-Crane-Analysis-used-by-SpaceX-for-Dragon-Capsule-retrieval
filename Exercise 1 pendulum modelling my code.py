# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 12:17:51 2020

@author: Subbu
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

d_rope = 44e-3
A = 0.25*(np.pi)*d_rope**2
E = 5e+9
tstop = 40
tinc = 0.1
Widthby2 = 3.81
Length = 7.84
Mass= 2394.12
g= 9.81
B=g/Length
theta0x = np.radians(5.0)
theta0y = np.radians(5.0)
# function for solving pair of 1st order differential equations
def solve1(tstop,tinc,initialconditions,params,func,E,g1y,g2y):
    A,B    = params
    ATOL,RTOL = E
    t = np.arange(0., tstop, tinc)
    
    solution =odeint(func, initialconditions, t, args=(params,),rtol=RTOL,atol=ATOL)
    return solution

def solve(tstop,tinc,initialconditions,params,func,E,g1y,g2y):
    A,B    = params
    ATOL,RTOL = E
    t = np.arange(0., tstop, tinc)
    solution =odeint(func, initialconditions, t, args=(params,),rtol=RTOL,atol=ATOL)
    #fig = plt.figure(1, figsize=(8,8))

    # Plot theta as a function of time
    '''ax1 = fig.add_subplot(311)
    ax1.plot(t, solution[:,0])
    ax1.set_xlabel('time')
    ax1.set_ylabel('theta')

    # Plot omega as a function of time
    ax2 = fig.add_subplot(312)
    ax2.plot(t, solution[:,1])
    ax2.set_xlabel('time')
    ax2.set_ylabel('omega')

    plt.tight_layout()
    plt.show()'''
    
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

#Function for Finding out Energy of the Pendulum
def energy(W,R,Mass,Length):
    z=R[:,2]
    pe = np.zeros(len(z))
    for i in range(0,len(z)):
        pe[i] = (Mass*g*(Length + z[i]))
    
    ke = np.zeros(len(z))
    for i in range(0,len(z)):
        ke[i] = 0.5*Mass*(Length**2)*(W[i]**2)
        
    energy = np.array([ke,pe])
    return energy
  
def tension(pos,Mass,Length):
    theta = pos[:,0]
    w     = pos[:,1]
    tension = np.zeros(len(w))
    for i in range(0,len(w)):
        tension[i] = Mass*((9.81*np.cos(theta[i])) + Length*(w[i])**2 )
    return tension

def extension(T,A,Y):
    extend = np.zeros(len(T))
    for i in range(0,len(T)):
        extend[i] = T[i] * (1/(A*Y))*100
    return extend

#define array of functions for solving motion in x direction
def f_initial(yx, t, params):
    thetax, omegax = yx      # unpack current values of y
    A , B = params        # unpack parameters
    derivsx = [omegax,      # list of dy/dt=f functions
             -omegax*A - B*(np.sin(thetax))]
    return derivsx

Omegas=np.arange(0.,10.0,0.01)
Lengths=np.zeros(len(Omegas))
for i in range(0,len(Omegas)):
    Controlparams = [1.49012e-8,1.49012e-8]
    parameters = [0.0,B] #[A,B]
    y0x = [theta0x,Omegas[i]] #[theta0x,omega0x]
    psolnx_initial = solve1(tstop,tinc,y0x,parameters,f_initial,Controlparams,'theta','omega')
    a=vectorize(psolnx_initial[:,0],Length,0)
    a1=a[:,0]
    amax = np.max(a1)
    Lengths[i] = amax

Alist=[]
for i in range(0,len(Lengths)):
    if (Widthby2 - Lengths[i]) > 0 :
        Alist.append(i)
    else:
        break
n = len(Alist)
n0 = Alist[n-1]
wo=Omegas[n0]

#define array of functions for solving motion in x direction
def f(yx, t, params):
    thetax, omegax = yx      # unpack current values of y
    A , B = params        # unpack parameters
    derivsx = [omegax,      # list of dy/dt=f functions
             -omegax*A - B*(np.sin(thetax))]
    return derivsx

Controlparams = [1.49012e-8,1.49012e-8]
parameters = [0.0,B] #[A,B]
y0x = [theta0x,Omegas[n0]] #[theta0x,omega0x]
psolnx = solve(tstop,tinc,y0x,parameters,f,Controlparams,'theta','omega')



#define array of functions for solving motion in y direction
def f(yy, t, params):
    thetay, omegay = yy     # unpack current values of y
    A , B = params        # unpack parameters
    derivsy = [omegay,      # list of dy/dt=f functions
             -omegay*A - B*(np.sin(thetay))]
    return derivsy

Controlparams = [1.49012e-8,1.49012e-8]
parameters = [0.0,B] #[A,B]
y0y = [theta0y,Omegas[n0] - 0.035] #[theta0y,omega0y]
psolny = solve(tstop,tinc,y0y,parameters,f,Controlparams,'theta','omega')


posvectX = vectorize(psolnx[:,0],Length,0)
posvectY = vectorize(psolny[:,0],Length,1)

tensionX = tension(psolnx,Mass,Length)
meanx =0.5*(max(tensionX) + min(tensionX))
diffx = max(tensionX) - min(tensionX)
tensionY = tension(psolny,Mass,Length)
meany =0.5*(max(tensionY) + min(tensionY))
diffy = max(tensionY) - min(tensionY)

extensionx=extension(tensionX,A,E)
extensiony=extension(tensionY,A,E)

EnergyX = energy(psolnx[:,1],posvectX,Mass,Length)
EnergyY = energy(psolny[:,1],posvectY,Mass,Length)
t = np.arange(0., tstop, tinc)
minkex=min(EnergyX[0,:])
minpex=min(EnergyX[1,:])
minkey=min(EnergyY[0,:])
minpey=min(EnergyY[1,:])
TotalEX=np.array(EnergyX[0,:] + EnergyX[1,:] )
TotalEY=np.array(EnergyY[0,:] + EnergyY[1,:] )

fig = plt.figure(1, figsize=(8,8))

#plt.plot(posvectX[:,0],posvectY[:,1])   
# Plot g1y as a function of time
'''ax3.plot(t, psolnx[:,0],'r-',label='theta X')
ax3.plot(t,psolny[:,0],'b-',label='theta Y')
ax3.set_xlabel('time (s)',size='16')
ax3.set_ylabel('theta (in rad)',size='16')
ax3.legend(loc="best")
plt.ylim(-0.6,1.0)'''

'''ax3 = fig.add_subplot(311)
ax3.plot(t, tensionX,'m-',label='Tension X')
ax3.plot(t,tensionY,'g-',label='Tension Y')
ax3.set_xlabel('time (s)',size='16')
ax3.set_ylabel('Tension',size='16')
ax3.legend(loc="best")
plt.ylim(20000,33000)'''


'''ax4 = fig.add_subplot(311)
ax4.plot(t, EnergyY[0,:],'g-',label='Kinetic Energy')
ax4.plot(t,EnergyY[1,:],'y-',label='Potetial Energy')
ax4.set_xlabel('time (s)',size='16')
ax4.set_ylabel('Energy',size='16')
ax4.legend(loc="best")
plt.ylim(0,35000)'''

ax4 = fig.add_subplot(311)
ax4.plot(t, extensionx,'r-',label='Percent Extension in X-Plane')
ax4.plot(t,extensiony,'b-',label='Percent Extension in Y-Plane')
ax4.set_xlabel('time (s)',size='16')
ax4.set_ylabel('Percent Extension',size='16')
ax4.legend(loc="best")
plt.ylim(0.26,0.44)


plt.tight_layout()
plt.show()














    





    
    
    

