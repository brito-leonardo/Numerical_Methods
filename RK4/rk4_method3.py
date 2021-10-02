#code name: 'rk4_method3.py'
#
#description: This code aims to apply the 4th order Runge-Kutta (RK4)
#             method to solve a system of first order ordinary differential
#             equations (ODE).
#
#             problem: we have an second order ODE :
#                   y''-4y'+4y=0, y(0)=-2, y'(0)=1
#
#             exact solution: y(t)=(-2+5t)exp(2t)
#
#             we want to find the solutions y(0.2)=x1(0.6) and y(0.6)=x2(0.6)
#             we conveniently rewrite the ODE as asystem with variables y1=y and y'=y2
#                   y1'=y2, y1(0)=-2
#                   y2'=4y2-4y1, y2(0)=1
#
#
#reference: Differential Equation with Boundary-Value Problems (2018)
#                                 (Dennis G. Zill)
#           [problem 3, section 9.4, p.385]
#
#date: October 2nd, 2021
#author: Leonardo Brito
#Physics Phd Student
#Institute of Physics (IFUSP)
#University of São Paulo (USP)
#_______________________________________________________________________
#import the required libraries
import numpy as np

#----------------------------------------------------------------------
#function f: function f(t,x)
def f(t,y):
    func=np.zeros((2))
    func[0]=y[1]
    func[1]=4*y[1]-4*y[0]
    return func
#----------------------------------------------------------------------
#exact solution  x1_ex(t)
def y_ex(t):
    y=(-2.0+5*t)*np.exp(2*t)
    return y
#----------------------------------------------------------------------
#RK4 method x(dt,t)
def rk4(dt,t,y):
    #define the auxiliary variables
    k1=np.zeros((2))
    k2=np.zeros((2))
    k3=np.zeros((2))
    k4=np.zeros((2))
    #rk4 scheme
    k1=f(t,y)
    k2=f(t+0.5*dt,y+0.5*dt*k1)
    k3=f(t+0.5*dt,y+0.5*dt*k2)
    k4=f(t+dt,y+dt*k3)
    y=y+dt*(k1+2*k2+2*k3+k4)/6.0
    return y
#_______________________________________________________________________
#main code

#time instants
t0=0.0
tf=0.2
#initial values
y0=np.zeros((2))
y=np.zeros((2))
y0[0]=-2.0
y0[1]=1.0

print('exact solution x1_ex(tf): ', y_ex(tf))
print('')

#aproximate solution for dt=0.2
#see table 94.1 from reference
t=t0
y=y0
dt=0.1
Nt=int((tf-t0)/dt)
print('RK4 approximation when dt is ', dt)
print('Columns t, x1 and x2')
for i in range(1,Nt+2):
    print(t,' ',y[0])
    y=rk4(dt,t,y)
    t=t0+i*dt
print('')

#aproximate solution for dt=0.1
#see table 94.1 from reference
t=t0
y=y0
dt=0.05
Nt=int((tf-t0)/dt)
print('RK4 approximation when dt is ', dt)
print('Columns t, x1 and x2')
for i in range(1,Nt+2):
    print(t,' ',y[0])
    y=rk4(dt,t,y)
    t=t0+i*dt
print('')

#_______________________________________________________________________