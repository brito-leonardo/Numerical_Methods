#code name: 'rk4_method.py'
#
#description: This code aims to apply the 4th order Runge-Kutta (RK4)
#             method to solve a first order ordinary differential
#             equation (ODE).
#
#             problem: we have an second order ODE :
#                   y'=2xy, y(1)=1
#             exact solution: y(x)=exp(x**2)/e , where 'e' is the Euler's number, e=2,71828...
#
#
#reference: Differential Equation with Boundary-Value Problems (2018)
#                                 (Dennis G. Zill)
#           [Example 1, section 9.2, p.375-376]
#
#date: October 1st, 2021
#author: Leonardo Brito
#Physics Phd Student
#Institute of Physics (IFUSP)
#University of São Paulo (USP)
#_______________________________________________________________________

#import the required libraries
import numpy as np

#----------------------------------------------------------

#function f: function f(x,y)
def f(x,y):
    func=2*x*y
    return func
#-----------------------------------------------------------
#exact solution  x1_ex(x)
def yex(x):
    e=np.exp(1.0) 
    y=np.exp(x**2)/e
    return y
#________________________________________________________________
#main code

dx=0.1 #step size
x0=1.0 #initial value
x=x0 #initial value
xf=1.5 #final value
Nx=int((xf-x0)/dx) #number of steps
y=1.0 #initial value

#evolve the system by the RK4 method
for i in range(1,Nx+1):
    k1=f(x,y)
    k2=f(x+0.5*dx,y+0.5*dx*k1)
    k3=f(x+0.5*dx,y+0.5*dx*k2)
    k4=f(x+dx,y+dx*k3)
    y=y+dx*(k1+2*k2+2*k3+k4)/6.0
    x=x0+i*dx
#print the exact and aproximate solutions
print('exact solution y(xf): ',yex(xf))
print('aproximate RK4 solution y: ',y)
#________________________________________________________________