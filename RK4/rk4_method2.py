#code name: 'rk4_method2.py'
#
#description: This code aims to apply the 4th order Runge-Kutta (RK4)
#             method to solve a system of first order ordinary differential
#             equations (ODE).
#
#             problem: we have an ODE system:
#                   x'=2x+4y, x(0)=-1
#                   y'=-x+6y, y(0)=6
#             exact solutions: x(t)=(26t-1)exp(4t)
#                              y(t)=(13t+6)exp(4t)
#
#             we want to find the solutions x(0.6)=x1(0.6) and y(0.6)=x2(0.6)
#             we conveniently rewrite the system with variables x1 and x2
#                   x1'=2x1+4x2, x1(0)=-1
#                   x2'=-x1+6x2, x2(0)=6
#
#
#reference: Differential Equation with Boundary-Value Problems (2018)
#                                 (Dennis G. Zill)
#           [Example 3, section 9.4, p.383-384]
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
def f(t,x):
    func=np.zeros((2))
    func[0]=2*x[0]+4*x[1]
    func[1]=-x[0]+6*x[1]
    return func
#----------------------------------------------------------------------
#exact solution  x1_ex(t)
def x1_ex(t):
    x1=(26*t-1.0)*np.exp(4*t)
    return x1
#----------------------------------------------------------------------
#exact solution  x2_ex(t)
def x2_ex(t):
    x2=(13*t+6.0)*np.exp(4*t)
    return x2

#----------------------------------------------------------------------
#RK4 method x(dt,t)
def rk4(dt,t,x):
    #define the auxiliary variables
    k1=np.zeros((2))
    k2=np.zeros((2))
    k3=np.zeros((2))
    k4=np.zeros((2))
    #rk4 scheme
    k1=f(t,x)
    k2=f(t+0.5*dt,x+0.5*dt*k1)
    k3=f(t+0.5*dt,x+0.5*dt*k2)
    k4=f(t+dt,x+dt*k3)
    x=x+dt*(k1+2*k2+2*k3+k4)/6.0
    return x
#_______________________________________________________________________
#main code

#time instants
t0=0.0
tf=0.6
#initial values
x0=np.zeros((2))
x=np.zeros((2))
x0[0]=-1.0
x0[1]=6.0

print('exact solution x1_ex(tf): ', x1_ex(tf))
print('exact solution x2_ex(tf): ', x2_ex(tf))
print('')

#aproximate solution for dt=0.2
#see table 94.1 from reference
t=t0
x=x0
dt=0.2
Nt=int((tf-t0)/dt)
print('RK4 approximation when dt is ', dt)
print('Columns t, x1 and x2')
for i in range(1,Nt+3):
    print(t,' ',x[0],' ',x[1])
    x=rk4(dt,t,x)
    t=t0+i*dt
print('')

#aproximate solution for dt=0.1
#see table 94.1 from reference
t=t0
x=x0
dt=0.1
Nt=int((tf-t0)/dt)
print('RK4 approximation when dt is ', dt)
print('Columns t, x1 and x2')
for i in range(1,Nt+3):
    print(t,' ',x[0],' ',x[1])
    x=rk4(dt,t,x)
    t=t0+i*dt
print('')

#_______________________________________________________________________
