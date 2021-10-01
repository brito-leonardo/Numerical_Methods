#code name: 'euler_method4.py'
#
#description: This code aims to apply the Euler's method
#             to solve a second order ordinary differential
#             equation (ODE) as a system.
#
#             problem: we have an second order ODE :
#                   y''-4y'=12x, y(0)=4, y'(0)=1
#             exact solution: y(x)=exp(-x/2)*(-cos(x)+(3/4)*sin(2*x))
#             we handle this problem rewriting it as a ODE system
#                 y1'(x)=y2(x) , y1(0)=4
#                 y2'(x)=y''(x)=4y+12x=4y1+12x, y2(0)=1
#
#
#reference: Differential Equation with Boundary-Value Problems (2018)
#                                 (Dennis G.Zill)
#           [Example 2. section 4.1.1, p.120]
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
#function f: function f(t,x)
def f(x,y):
    func=np.zeros((2))
    func[0]=y[1]
    func[1]=4*y[0]+12*x
    return func
#-----------------------------------------------------------
#exact solution  x1_ex(t)
def y_ex(x):
    y=3*np.exp(2*x)+np.exp(-2*x)-3*x
    return y
#________________________________________________________________
#main code
dx=0.01 #time-step size
x0=0.0
xf=2.0 # final time
Nx=int((xf-x0)/dx) #number of time steps
y=np.zeros((2))
#initial conditions
y[0]=4.0 
y[1]=1.0 

#define de data files
data1=np.empty([Nx+1,2])
data2=np.empty([Nx+1,2])

#write the initial values
data1[0,0]=x0
data1[0,1]=y_ex(x0)
data2[0,0]=x0
data2[0,1]=y[0]

#Evolve the system by Euler's method
for i in range(1,Nx+1):
    x=x0+i*dx
    y=y+dx*f(x,y)
    data1[i,0]=x
    data1[i,1]=y_ex(x)
    data2[i,0]=x
    data2[i,1]=y[0]

#save the data files
np.savetxt('y-exact-py.txt',data1,fmt="%10.4f")
np.savetxt('y-euler-py.txt',data2,fmt="%10.4f")

#________________________________________________________________
