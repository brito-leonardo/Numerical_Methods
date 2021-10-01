#code name: 'euler_improved.f90'
#
#description: This code aims to apply the Improved Euler's
#             method, to solve a ordinary differential
#             equation (ODE).
#
#             problem: we want to find y(1.0)
#
#             we have the ODE:
#                   y'=x-2y, y(0)=1
#
#
#reference: Differential Equation with Boundary-Value Problems (2018)
#                                 (Dennis G.Zill)
#           [Exercise 15, section 9.1, p.373]
#
#date: September 23rd, 2021
#author: Leonardo Brito
#Physics Phd Student
#Institute of Physics (IFUSP)
#University of SÆo Paulo (USP)
#_______________________________________________________________________

#import the required libraries
import numpy as np

#----------------------------------------------------------

#function f: function f(x,y)
def f(x,y):
    func=x-2*y
    return func
#-----------------------------------------------------------
#exact solution  x1_ex(x)
def y_ex(x):
    y=(1.0/2)*x-(1.0/4)+(5.0/4)*np.exp(-2*x)
    return y
#________________________________________________________________

#main code
    
dx=0.01 #step-size
x0=0.0
x=x0
xf=1.0
y=1.0
y=1.0 #initial value
Nx=int(xf/dx)
#evolve the system by teh Improved Euler's method
for i in range(1,Nx+1):
    xold=x
    yaux=y+dx*f(x,y)
    x=x0+i*dx
    y=y+dx*(f(xold,y)+f(x,yaux))/2.0

print('Exact solution yex(1.0): ',y_ex(xf))
print('Improved Euler method y(1.0): ',y)
    
#_______________________________________________________________________
    