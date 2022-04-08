# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 15:00:30 2020
@author: Subbu
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#defining a function to solve the matrix equation for position velocity and acceleration.
def pos_vel_acc(r2,r3,r4,d,lf,D,th2,th2d,th2dd):
    A = d - (r2*np.cos(th2))
    B = -r2*np.sin(th2)
    C = np.arccos((-r3**2 - r4**2 + A**2 + B**2)/(2*r3*r4))
    a = r3 + r4*np.cos(C) 
    b = r4*np.sin(C)
    E = np.arctan(b/a)
    th3 = np.arccos(A/np.sqrt(a**2 + b**2)) - E
    th4 = C +th3 
    thf = np.arctan((lf*np.sin(th2))/(D - lf*np.cos(th2)))
    l5 = (D - lf*np.cos(th2))/np.cos(thf)
    velv = np.array([r2*th2d*np.sin(th2) , -r2*th2d*np.cos(th2)])
    velm = np.array([[r3*np.cos(th3) , r4*np.sin(th4)] , [r3*np.cos(th3) , r4*np.cos(th4)]])
    vel = np.linalg.solve(velm, velv)
    th3d = vel[0]
    th4d = vel[1]
    pistondm = np.array([[np.cos(thf) , -l5*np.sin(thf)] , [ np.sin(thf) , l5*np.cos(thf)]])
    pistondv = np.array([lf*th2d*np.sin(th2)  , lf*th2d*np.cos(th2)])
    pistond = np.linalg.solve(pistondm,pistondv)
    l5d = pistond[0]
    thfd = pistond[1]
    f1 = -(r4*np.cos(th4)*th4d**2) - (r3*np.cos(th3)*th3d**2)
    - (r2*np.cos(th2)*th2d**2) - (r2*np.sin(th2)*th2dd)
    f2 = (r4*np.sin(th4)*th4d**2) + (r3*np.sin(th3)*th3d**2) + (r2*np.sin(th2)*th2d**2) - (r2*np.cos(th2)*th2dd)
    accv = np.array([f1 , f2])
    accm = np.array([[r3*np.sin(th3) , +r4*np.sin(th4)] , [r3*np.cos(th3) , r4*np.cos(th4)]]) 
    acc = np.linalg.solve(accm, accv)
    
    pistonddm = np.array([[ np.cos(thf) , -l5*np.sin(thf)] , [ np.sin(thf) , l5*np.cos(thf) ]])
    p1 = lf*(th2dd*np.sin(th2) + th2d**2*np.cos(th2)) + 2*l5d*thfd*np.sin(thf)  + l5*thfd**2*np.cos(thf)
    p2 = lf*(th2dd*np.cos(th2) - th2d**2*np.sin(th2)) - 2*l5d*thfd*np.cos(thf)  + l5*thfd**2*np.sin(thf)
    pistonddv = np.array([p1,p2])
    
    pistondd = np.linalg.solve(pistonddm,pistonddv)
    l5dd = pistondd[0]
    thfdd = pistondd[1]
    
    return [th3,th4,vel,acc,l5d,thfd,l5dd,thfdd,l5,thf]
 
#defining constants 
tstop = 40
tinc = 0.1
t = np.arange(tinc,tstop,tinc)

d= 5 #AD
l2 = 5 #AB
l3 = 5 #BC
l4 = 5 #CD
lf = 3 #AF
D = 2 #AE

# setting initial conditions
th2dd = 0.0 #theta 2 double dot (Angular acceleration)
th2d0 = np.radians(80)/tstop
th2d =np.zeros(len(t))
th20 = np.radians(10)
th2 = np.zeros(len(t))

for i in range(0,len(t)):
    th2d[i] = th2d0 + th2dd*t[i]
    th2[i]  =  th2d0*t[i]  +  0.5*th2dd*t[i]**2

still = []
for i in range(0,len(t)):
    ans = final = pos_vel_acc(l2,l3,l4,d,lf,D,th2[i],th2d[i],th2dd)
    still.append(ans)
    
Th3 = []
Th4 = []
Th3d = []
Th4d = []
Th3dd = []
Th4dd = []
l5d = []
thfd = []
l5dd = []
thfdd =[]
Bx = []
By = []
Cx = []
Cy = []
l5 = []
thf = []

#creating data arrays
for i in range(0,len(t)):
    Th3.append(still[i][0])
    if 2*3.14 - still[i][1] > 1.57:
        Th4.append(abs(-3.14 + still[i][1]))
    elif 2*3.14 + still[i][1] < -1.57:
        Th4.append(abs(3.14 + still[i][1]))
    else :
        Th4.append(abs(2*3.14 - still[i][1]))
    Th3d.append(still[i][2][0])
    Th4d.append(still[i][2][1])
    Th3dd.append(still[i][3][0])
    Th4dd.append(still[i][3][1])   
    l5d.append(still[i][4])
    thfd.append(still[i][5])
    l5dd.append(still[i][6])
    thfdd.append(still[i][7])
    Bx.append(l2*np.cos(th2[i]))
    By.append(l2*np.sin(th2[i]))
    Cx.append(d + l4*np.cos(Th4[i]))
    Cy.append(l4*np.sin(Th4[i]))
    l5.append(still[i][8])
    thf.append(still[i][9])

#Creating animation
'''fig = plt.figure()
fig.canvas.set_window_title('Matplotlib Animation')
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-20,20), ylim=(-20,20))
ax.grid()
ax.set_title('Four-Bar motion')
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
line, = ax.plot([], [], 'o-', lw=2, color='#de2d26')


# initialization function
def init():
    line.set_data([], [])
    return line,

# animation function
def animate(i):
    x_points = [0, Bx[i], Cx[i],d]
    y_points = [0, By[i], Cy[i],0]

    line.set_data(x_points, y_points)
    return line,


ani = animation.FuncAnimation(fig, animate, init_func=init, frames=len(t), interval=180, blit=True, repeat=False)
ani.save('Case 10.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
##plt.show()
'''
#Plotting relevant graphs
fig = plt.figure(1, figsize=(8,8))
ax3 = fig.add_subplot(311)
plt.plot(t, th2, 'b-', label = 'Theta 2')
plt.plot(t, Th3, 'g-', label = 'Theta 3')
plt.plot(t, Th4, 'r-', label = 'Theta 4')
ax3.set_xlabel('time (s)',size='16')
ax3.set_ylabel('Theta',size='16')
ax3.legend(loc="best")







