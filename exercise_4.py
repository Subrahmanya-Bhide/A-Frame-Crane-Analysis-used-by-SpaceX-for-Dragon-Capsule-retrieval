# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 11:22:41 2020

@author: Lenovo
"""

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
                       [m4*(y4dd+g) +N],
                       [I4*th4dd + N*(xn-x4)],
                       [m5*x5dd],
                       [m5*(y5dd+g)],
                       [I5*(thfdd)],
                       [m6*x6dd],
                       [m6*(y6dd+g)],
                       [I6*(thfdd)]])
    solution = np.linalg.solve(Matrix,Vector)
    return solution