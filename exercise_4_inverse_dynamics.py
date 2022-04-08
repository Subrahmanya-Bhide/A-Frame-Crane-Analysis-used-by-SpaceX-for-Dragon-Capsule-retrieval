# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 09:53:53 2020

@author: Subbu
"""
import numpy as np
np.seterr(divide='raise')

def pos_vel_acc(r2,r3,r4,d,th2,th2d,th2dd):
    A = d - (r2*np.cos(th2))
    B = -r2*np.sin(th2)
    C = np.arccos((r3**2 + r4**2 - A**2 -B**2)/(2*r3*r4))
    a = r3*np.cos(C) - 1
    b = r3*np.sin(C)
    E = np.arctan(b/a)
    th4 = np.arccos(A/np.sqrt(a**2 + b**2)) - E
    th3 = th4 + C
    velv = np.array([r2*th2d*np.sin(th2) , -r2*th2d*np.cos(th2)])
    velm = np.array([[-r3*np.cos(th3) , r4*np.sin(th4)] , [r3*np.cos(th3) , -r4*np.cos(th4)]])
    vel = np.linalg.solve(velm, velv)
    th3d = vel[0]
    th4d = vel[1]

    f1 = (r4*np.cos(th4)*th4d**2) - (r3*np.cos(th3)*th3d**2)
    - (r2*np.cos(th2)*th2d**2) - (r2*np.sin(th2)*th2dd)
    
    f2 = -(r4*np.sin(th4)*th4d**2) + (r3*np.sin(th3)*th3d**2) + (r2*np.sin(th2)*th2d**2) - (r2*np.cos(th2)*th2dd)

    accv = np.array([f1 , f2])
    accm = np.array([[r3*np.sin(th3) , -r4*np.sin(th4)] , [r3*np.cos(th3) , -r4*np.cos(th4)]]) 

    acc = np.linalg.solve(accm, accv)
    return [th3,th4,vel,acc]
    

d= 5
l2 = 2
l3 = 8
l4 = 6
th2 = np.radians(30)
th2d = 1.0
th2dd = 0.0
r2 = l2/2
r3 = l3/2
r4 = l4/2
#th2arey = np.arange(10,89,0.5)
#still = []
#for i in range(0,len(th2arey)):
#    ans = final = pos_vel_acc(r2,r3,r4,d,th2arey[i],th2d,th2dd)
#    still.append(ans)
ans = pos_vel_acc(l2,l3,l4,d,th2,th2d,th2dd)
  
th3 = ans[0]
th4 = ans[1]
th3d = ans[2][0]
th4d = ans[2][1]
th3dd =ans[3][0]
th4dd =ans[3][1]

m2=10
m3=20
m4=30
ph2=0
ph3=0
ph4=0
x2=(r2)*np.cos(th2+ph2)
x3=r2*np.cos(th2)  + r3*np.cos(th3+ph3)
x4=d + r4*np.cos(th4+ph4)
y2=(r2)*np.sin(th2+ph2)
y3=r2*np.sin(th2)  + r3*np.sin(th3+ph3)
y4= r4*np.sin(th4+ph4)
xA=0
yA=0
xB=l2*np.cos(th2)
yB=l2*np.sin(th2)
xC=l2*np.cos(th2) + l3*np.cos(th3)
yC=l2*np.sin(th2) + l3*np.sin(th3)
xD=d
yD=0
I2 = 1
I3 = 1
I4 = 1
Fx = 5
Fy = 5
xF = 1
yF = 1
g = 9.81

M = np.array([[1,0,1,0,0,0,0,0,0,0,0,-m2,0,0,0,0,0],
              [0,1,0,1,0,0,0,0,0,0,0,0,-m2,0,0,0,0],
              [y2 - yA,xA-x2,y2-yB,xB-x2,0,0,0,0,-I2,0,0,0,0,0,0,0,0],
              [0,0,-1,0,1,0,0,0,0,0,0,0,0,-m3,0,0,0],
              [0,0,0,-1,0,1,0,0,0,0,0,0,0,0,-m3,0,0],
              [0,0,yB-y3,x3-xB,y3-yC,xC-x3,0,0,0,-I3,0,0,0,0,0,0,0],
              [0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,-m4,0],
              [0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,-m4],
              [0,0,0,0,yC-y4,x4-xC,y4-yD,xD-x4,0,0,-I4,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,r2*np.sin(th2+ph2),0,0,1,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,-r2*np.cos(th2+ph2),0,0,0,1,0,0,0,0],
              [0,0,0,0,0,0,0,0,r2*np.sin(th2),r3*np.sin(th3+ph3),0,0,0,1,0,0,0],
              [0,0,0,0,0,0,0,0,-r2*np.cos(th2),-r2*np.cos(th3+ph3),0,0,0,0,1,0,0],
              [0,0,0,0,0,0,0,0,0,0,r4*np.sin(th4+ph4),0,0,0,0,1,0],
              [0,0,0,0,0,0,0,0,0,0,-r4*np.cos(th4+ph4),0,0,0,0,0,1],
              [0,0,0,0,0,0,0,0,-r2*np.sin(th2),-r3*np.sin(th3),r4*np.sin(th4),0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,r2*np.cos(th2),r3*np.cos(th3),-r4*np.cos(th4),0,0,0,0,0,0],
              ])


K =np.array([-Fx,

 -Fy  + m2*g,

 Fx*(yF - y2)  - Fy*(xF - x2),

 0,
 
 0,

 0,

 0,

 0,

 0,

 -r2*th2d**2*np.cos(th2+ph2),

 -r2*th2d**2*np.sin(th2+ph2),

 -r2*th2d**2*np.cos(th2)  -  r3*th3d**2*np.cos(th3+ph3),

 -r2*th2d**2*np.sin(th2) - r3*th3d**2*np.sin(th3+ph3),

 -r4*th4d**2*np.cos(th4+ph4),

 -r4*th4d**2*np.sin(th4+ph4),
 
 r2*th2d**2*np.cos(th2) + r3*th3d**2*np.cos(th3)  - r4*th4d**2*np.cos(th4),

 r2*th2d**2*np.sin(th2) + r3*th3d**2*np.sin(th3)  - r4*th4d**2*np.sin(th4)])

    
    
final = np.linalg.solve(M,K)    




