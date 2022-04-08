# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 16:05:30 2020

@author: subbu
"""

import numpy as np
from numpy import sin, cos 
import matplotlib.pyplot as plt

def modulo(A):
    for i in range(0,len(A)):
        A[i] = abs(A[i])


#defining a function to solve the matrix equation for position velocity and acceleration.
def pos_vel_acc(r2,r3,r4,d,lf,D,th4,th4d,th4dd):
    A = d + (r4*np.cos(th4))
    B = r4*np.sin(th4)
    C = np.arccos((-r3**2 - r2**2 + A**2 + B**2)/(2*r2*r3))
    a = r3 + r2*np.cos(C) 
    b = r2*np.sin(C)
    E = np.arctan(b/a)
    th3 = np.arccos(A/np.sqrt(a**2 + b**2)) - E
    th2 = C +th3 
    thf = np.arctan((lf*np.sin(th2))/(D - lf*np.cos(th2)))
    l5 = (D - lf*np.cos(th2))/np.cos(thf)
    velv = np.array([r4*th4d*np.sin(th4) , r4*th4d*np.cos(th4)])
    velm = np.array([[r2*np.sin(th2) , r3*np.sin(th3)] , [r2*np.cos(th2) , r3*np.cos(th3)]])
    vel = np.linalg.solve(velm, velv)
    th2d = vel[0]
    th3d = vel[1]
    pistondm = np.array([[np.cos(thf) , -l5*np.sin(thf)] , [ np.sin(thf) , l5*np.cos(thf)]])
    pistondv = np.array([lf*th2d*np.sin(th2)  , lf*th2d*np.cos(th2)])
    pistond = np.linalg.solve(pistondm,pistondv)
    l5d = pistond[0]
    thfd = pistond[1]
    f1 = -(r4*np.cos(th4)*th4d**2) - (r3*np.cos(th3)*th3d**2)
    - (r2*np.cos(th2)*th2d**2) + (r4*np.sin(th4)*th4dd)
    f2 = (r4*np.sin(th4)*th4d**2) + (r3*np.sin(th3)*th3d**2) + (r2*np.sin(th2)*th2d**2) + (r4*np.cos(th4)*th4dd)
    accv = np.array([f1 , f2])
    accm = np.array([[r2*np.sin(th2) , +r3*np.sin(th3)] , [r2*np.cos(th2) , r3*np.cos(th3)]]) 
    acc = np.linalg.solve(accm, accv)
    th2dd = acc[0]
    pistonddm = np.array([[ np.cos(thf) , -l5*np.sin(thf)] , [ np.sin(thf) , l5*np.cos(thf) ]])
    p1 = lf*(th2dd*np.sin(th2) + th2d**2*np.cos(th2)) + 2*l5d*thfd*np.sin(thf)  + l5*thfd**2*np.cos(thf)
    p2 = lf*(th2dd*np.cos(th2) - th2d**2*np.sin(th2)) - 2*l5d*thfd*np.cos(thf)  + l5*thfd**2*np.sin(thf)
    pistonddv = np.array([p1,p2])
    
    pistondd = np.linalg.solve(pistonddm,pistonddv)
    l5dd = pistondd[0]
    thfdd = pistondd[1]
    
    return [th2,th3,vel,acc,l5d,thfd,l5dd,thfdd,l5,thf,thfd,thfdd,l5d,l5dd]
 
#function to solve for forces
def fun_solve(th2,th3,th4,thf,th2dd,th3dd,th4dd,thfdd,l5,x2dd,y2dd,x3dd,y3dd,x4dd,y4dd,x5dd,y5dd,x6dd,y6dd):
    x2=r2*cos(th2+ph2)
    x3=l2*cos(th2) + r3*cos(th3+ph3)
    x4=d - r4*cos(th4)
    x5=lf*cos(th2) + r5*cos(thf)
    x6=d - r6*cos(thf)
    y2=r2*sin(th2+ph2)
    y3=l2*sin(th2) + r3*sin(th3+ph3)
    y4=-r4*sin(th4)
    y5=lf*sin(th2) - r5*sin(thf)
    y6=r6*sin(thf)
    xb=l2*cos(th2)
    yb=l2*sin(th2)
    xc=d - l4*cos(th4)
    yc=-l4*sin(th4)
    xa=0
    ya=0
    xd=d
    yd=0
    xe=D
    ye=0
    xf=lf*cos(th2)
    yf=lf*sin(th2)
    l=l5-l51-l52
    Matrix = np.array([[1,0,1,0,0,0,0,0,0,0,0,0,1,0,0],
                       [0,1,0,1,0,0,0,0,0,0,0,0,0,1,0],
                       [y2-ya,xa-x2,y2-yb,xb-x2,0,0,0,0,0,0,0,0,y2-yf,xf-x2,0],
                       [0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0],
                       [0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0],
                       [0,0,yb-y3,x3-xb,y3-yc,xc-x3,0,0,0,0,0,0,0,0,0],
                       [0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0],
                       [0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0],
                       [0,0,0,0,yc-y4,x4-xc,y4-yd,xd-x4,0,0,0,0,0,0,0],
                       [0,0,0,0,0,0,0,0,-cos(thf),0,0,sin(thf),-1,0,0],
                       [0,0,0,0,0,0,0,0,sin(thf),0,0,cos(thf),0,-1,0],
                       [0,0,0,0,0,0,0,0,0,0,0,l51-delta,yf-y5,x5-xf,-1],
                       [0,0,0,0,0,0,0,0,cos(thf),1,0,-sin(thf),0,0,0],
                       [0,0,0,0,0,0,0,0,-sin(thf),0,1,-cos(thf),0,0,0],
                       [0,0,0,0,0,0,0,0,0,y6-ye,xe-x6,l+l52,0,0,1]
                       ])
    Vector = np.array([[m2*x2dd],
                       [m2*y2dd + m2*g],
                       [I2*th2dd],
                       [m3*x3dd],
                       [m3*(y3dd + g)],
                       [I3*th3dd],
                       [m4*x4dd],
                       [m4*(y4dd+g) + N],
                       [I4*th4dd + N*(xn-x4)],
                       [m5*x5dd],
                       [m5*(y5dd+g)],
                       [I5*(thfdd)],
                       [m6*x6dd],
                       [m6*(y6dd+g)],
                       [I6*(thfdd)]])
    solution = np.linalg.solve(Matrix,Vector)
    return solution
    
    

#defining constants 
tstop = 40
tinc = 0.1
t = np.arange(tinc,tstop,tinc)
w = 0.5
density = 2000
d= 3 #AD
l2 = 6 #AB
l3 = 3 #BC
l4 = 6 #CD
lf = 4 #AF
r2 = l2/2
r3 = l3/2
r4 = l4/2
D = 4 #AE
ph2 = ph3 = ph4 = 0
l51 = 0.5
l52 = 0.25
r5 = l51/2
r6 = l52/2
delta = 0.1
m2= density * l2
m3=density * l3
m4=density * l4
m5=density * l51
m6=density * l52
g=9.81
I2=(1/12)*m2*(w**2 + l2**2) + (1/4)*m2*l2**2
I3=(1/12)*m3*(w**2 + l3**2) + (1/4)*m3*l3**2
I4=(1/12)*m4*(w**2 + l4**2) + (1/4)*m4*l4**2
I5=(1/12)*m5*(w**2 + l51**2) + (1/4)*m5*l51**2
I6=(1/12)*m6*(w**2 + l52**2) + (1/4)*m6*l52**2
xn = 8
N = 25000


# setting initial conditions
th4dd = 0.0 #theta 2 double dot (Angular acceleration)
th4d0 = np.radians(80)/tstop
th4d =np.zeros(len(t))
th40 = np.radians(10)
th4 = np.zeros(len(t))

for i in range(0,len(t)):
    th4d[i] = th4d0 + th4dd*t[i]
    th4[i]  =  th40 + th4d0*t[i]  +  0.5*th4dd*t[i]**2

still = []
for i in range(0,len(t)):
    ans = final = pos_vel_acc(l2,l3,l4,d,lf,D,th4[i],th4d[i],th4dd)
    still.append(ans)
    
Th2 = []
Th3 = []
Th2d = []
Th3d = []
Th2dd = []
Th3dd = []
l5d = []
thfd = []
l5dd = []
thfdd =[]
l5 = []
thf = []
thfd= []
thfdd = []
x2dd = []
y2dd = []
x3dd = []
y3dd = []
x4dd = []
y4dd = []
x5dd = []
y5dd = []
x6dd = []
y6dd = []

#creating data arrays
for i in range(0,len(t)):
    Th2.append(still[i][0])
    if 2*3.14 - still[i][1] > 1.57:
        Th3.append(abs(-3.14 + still[i][1]))
    elif 2*3.14 + still[i][1] < -1.57:
        Th3.append(abs(3.14 + still[i][1]))
    else :
        Th3.append(abs(2*3.14 - still[i][1]))
    Th2d.append(still[i][2][0])
    Th3d.append(still[i][2][1])
    Th2dd.append(still[i][3][0])
    Th3dd.append(still[i][3][1])   
    l5d.append(still[i][4])
    thfd.append(still[i][5])
    l5dd.append(still[i][6])
    thfdd.append(still[i][7])
    l5.append(still[i][8])
    thf.append(still[i][9])
    thfd.append(still[i][10])
    thfdd.append(still[i][11])
    
 
    
for i in range(0,len(t)):
    x2dd.append(-r2*(Th2d[i]**2*cos(Th2[i]+ph2) + Th2dd[i]*sin(Th2[i]+ph2)))
    y2dd.append(r2*(-Th2d[i]**2*sin(Th2[i]+ph2) + Th2dd[i]*cos(Th2[i]+ph2)))
    x3dd.append(-l2*((Th2d[i]**2)*cos(Th2[i]) + Th2dd[i]*sin(Th2[i])) - (r3*((Th3d[i]**2)*cos(Th3[i]+ph3) + Th3dd[i]*sin(Th3[i]+ph3))))
    y3dd.append(l2*(-Th2d[i]**2*sin(Th2[i]) + Th2dd[i]*cos(Th2[i]))  + r3*(-Th3d[i]**2*sin(Th3[i]+ph3) + Th3dd[i]*cos(Th3[i]+ph3)))
    x4dd.append(r4*(th4d[i]**2*cos(th4[i]+ph4) + th4dd*sin(th4[i]+ph4)))
    y4dd.append(-r4*(-th4d[i]**2*sin(th4[i]+ph4) + th4dd*cos(th4[i]+ph4)))
    x5dd.append(-lf*(Th2d[i]**2*cos(Th2[i]) + Th2dd[i]*sin(Th2[i])) - r5*(thfd[i]**2*cos(thf[i]) + thfdd[i]*sin(thf[i])))
    y5dd.append(lf*(-Th2d[i]**2*sin(Th2[i]) + Th2dd[i]*cos(Th2[i])) - r5*(-thfd[i]**2*sin(thf[i]) + thfdd[i]*cos(thf[i])))
    x6dd.append(r6*(thfd[i]**2*cos(thf[i]) + thfdd[i]*sin(thf[i])))
    y6dd.append(r6*(-thfd[i]**2*sin(thf[i]) + thfdd[i]*cos(thf[i])))
   

    
still1 = []
for n in range(0,len(t)):
    ans1 = fun_solve(Th2[n],Th3[n],th4[n],thf[n],Th2dd[n],Th3dd[n],th4dd,thfdd[n],l5[n],x2dd[n],y2dd[n],x3dd[n],y3dd[n],x4dd[n],y4dd[n],x5dd[n],y5dd[n],x6dd[n],y6dd[n])
    still1.append(ans1)
 
R1=[]
R2=[]
R3=[]
R4=[]
R5=[]
R6=[]
R7=[]
R8=[]
Fp=[]
R10=[]
R11=[]
R12=[]
Fx=[]
Fy=[]
M1=[]    
    
for i in range(0,len(t)):
    R1.append(abs(still1[i][0]))
    R2.append(abs(still1[i][1]))
    R3.append(abs(still1[i][2]))
    R4.append(abs(still1[i][3]))
    R5.append(abs(still1[i][4]))
    R6.append(abs(still1[i][5]))
    R7.append(abs(still1[i][6]))
    R8.append(abs(still1[i][7]))
    Fp.append(abs(still1[i][8]))
    R10.append(abs(still1[i][9]))
    R11.append(abs(still1[i][10]))
    R12.append(abs(still1[i][11]))
    Fx.append(abs(still1[i][12]))
    Fy.append(abs(still1[i][13]))
    M1.append(abs(still1[i][14]))
 
#plotting necessary graphs
plt.plot(th4*180/3.14,Fp)
plt.ylabel('Fp', size ='14')
plt.xlabel('theta', size ='14')