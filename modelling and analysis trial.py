import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
#from sundials import cvode
#from scipytoolkits.odes.sundials.cvode import cvode

d_rope = 44e-3
A = 0.25*(np.pi)*d_rope**2
E = 5e+9
tstop = 100
tinc = 0.1
Widthby2 = 3.81
Length = 7.84
Mass= 2394.12
g= 9.81
B=g/Length
theta0 = np.radians(30)
phi0 = np.radians(30)
thetadot0 = -1.58
phidot0 = 0.916


t = np.arange(0., tstop, tinc)
#thetadot0 = (0.116*xdot[i]) + (0.0029)
#phidot0 =  (-1.57*xdot[i]) + (0.39)


def solve(tstop,tinc,initialconditions,params,func,E,g1y,g2y):
    A,B    = params
    ATOL,RTOL = E
    t = np.arange(0., tstop, tinc)
    
    solution =odeint(func, initialconditions, t, args=(params,),rtol=RTOL,atol=ATOL)
    return solution

def f(yy, t, params):
    theta, thetadot, phi, phidot = yy     # unpack current values of y
    A , B = params        # unpack parameters
    derivsy = [thetadot,      # list of dy/dt=f functions
             -B*np.sin(theta) + phidot**2 *(np.cos(theta))*np.sin(phi),
              phidot,
              -2*thetadot*phidot*(1/np.tan(theta))]
    return derivsy
Controlparams = [1.49012e-8,1.49012e-8]
parameters = [0.0,B] 
y0y = [theta0, thetadot0, phi0, phidot0] 
psolny = solve(tstop,tinc,y0y,parameters,f,Controlparams,'theta','omega')
X= Length*np.sin(psolny[:,0])*np.cos(psolny[:,2])
Y = -Length*np.sin(psolny[:,0])*np.sin(psolny[:,2])

fig = plt.figure(1, figsize=(8,8))
ax1 = fig.add_subplot(311)
ax1.plot(X,Y,label='theta X')
#ax1.plot(t,psolny[:,0])
#ax1.set_xlabel('time (s)',size='16')
#ax1.set_ylabel('theta (in rad)',size='16')


#Axes3D.plot_trisurf(psolny[:,0],psolny[:,2],t,cmap = cm.jet)
'''fig = plt.figure(figsize=(6,6))

ax = Axes3D(fig)
ax.scatter(psolny[:,0],psolny[:,2],0, marker='o')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')'''

#plt.show()


#ax.plot_trisurf(psolny[:,0],psolny[:,2],t,cmap = cm.jet)



