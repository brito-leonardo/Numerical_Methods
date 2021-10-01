#code name: 'euler_method3.py'
#
#description: This code aims to apply the Euler's method
#             to solve an ordinary differential
#             equation (ODE) system.
#
#             problem: we have 1st order ODE system:
#                   x1'=-(2/25)*x1+(1/50)*x2, x1(0)=25
#                   x2'=(2/25)*x1-(2/25)*x2, x2(0)=0
#             exact solutions: x1(t)=(25/2)*exp(-t/25)+(25/2)*exp(-3*t/25)
#                              x2(t)=25*exp(-t/25)-25*exp(-3*t/25
#
#
#reference: Differential Equation with Boundary-Value Problems (2018)
#                                 (Dennis G.Zill)
#           [Example 3. section 4.9, p.186-187]
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
def f(t,x):
    func=np.zeros((2))
    func[0]=-(2.0/25)*x[0]+(1.0/50)*x[1]
    func[1]=(2.0/25)*x[0]-(2.0/25)*x[1]
    return func
#-----------------------------------------------------------
#exact solution  x1_ex(t)
def x1_ex(t):
    x1=(25.0/2)*np.exp(-t/25)+(25.0/2)*np.exp(-3*t/25)
    return x1
#-----------------------------------------------------------
#exact solution  x2_ex(t)
def x2_ex(t):
    x2=25.0*np.exp(-t/25)-25.0*np.exp(-3*t/25)
    return x2

#________________________________________________________________
#main code
dt=0.01 #time-step size
t0=0.0
tf=100.0 # final time
Nt=int((tf-t0)/dt) #number of time steps
x=np.zeros((2))
#initial conditions
x[0]=25.0 
x[1]=0.0 

#define de data files
data1=np.empty([Nt+1,2])
data2=np.empty([Nt+1,2])
data3=np.empty([Nt+1,2])
data4=np.empty([Nt+1,2])

#write the initial values
data1[0,0]=t0
data1[0,1]=x1_ex(t0)
data2[0,0]=t0
data2[0,1]=x[0]
data3[0,0]=t0
data3[0,1]=x2_ex(t0)
data4[0,0]=t0
data4[0,1]=x[1]

#Evolve the system by Euler's method
for i in range(1,Nt+1):
    t=t0+i*dt
    x=x+dt*f(t,x)
    data1[i,0]=t
    data1[i,1]=x1_ex(t)
    data2[i,0]=t
    data2[i,1]=x[0]
    data3[i,0]=t
    data3[i,1]=x2_ex(t)
    data4[i,0]=t
    data4[i,1]=x[1]

#save the data files
np.savetxt('x1-exact-py.txt',data1,fmt="%10.4f")
np.savetxt('x1-euler-py.txt',data2,fmt="%10.4f")
np.savetxt('x2-exact-py.txt',data3,fmt="%10.4f")
np.savetxt('x2-euler-py.txt',data4,fmt="%10.4f")

#________________________________________________________________
