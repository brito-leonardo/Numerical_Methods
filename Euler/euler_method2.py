#code name: 'euler_method2.py'
#
#description: This code aims to apply the Euler's method,
#             to solve a ordinary differential
#             equation (ODE).
#
#             problem: we have the ODE:
#                   y'=f(t,y), f(t,y)=y-t**2+1 , y(0)=0.5
#             exact solution: y(t)=1+2*t+t**2-0.5*exp(t)
#
#
#
#reference: Materials from 'Computational Physics Course (2015)'
#           Prof.Dr. Leonardo Castelano
#           Physics Department (DF)
#           Federal university of São Carlos (UFSCAR)
#
#date: September 26th, 2021
#author: Leonardo Brito
#Physics Phd Student
#Institute of Physics (IFUSP)
#University of São Paulo (USP)
#_______________________________________________________________________

import numpy as np

#-----------------------------------------------------------
#function f0: exact solution
def f0(t):
    func=1.0+2*t+t**2-0.5*np.exp(t)
    return func
#-----------------------------------------------------------

#function f: function f(t,y)
def f(t,y):
    func=y-t**2+1.0
    return func
#-----------------------------------------------------------
#________________________________________________________________
#main code
dt=0.025 #time-step size
tf=0.5 # final time
Nt=int(tf/dt) #number of time steps
t=np.zeros((Nt+1))
y=np.zeros((Nt+1))
yex=np.zeros((Nt+1))
#initial conditions
t[0]=0.0
y[0]=0.5 
yex[0]=0.5
#Evolve the system by Euler's method
for i in range(1,Nt+1):
    t[i]=i*dt
    yex[i]=f0(t[i])
    y[i]=y[i-1]+dt*f(t[i],y[i-1])

#write the results in data files
#exact results
data1=np.empty([t.size,2])
data1[:,0]=t
data1[:,1]=yex
np.savetxt('exact2-py.txt',data1,fmt="%.14E")
#euler results
data2=np.empty([t.size,2])    
data2[:,0]=t
data2[:,1]=y
np.savetxt('euler2-py.txt',data2,fmt="%.14E")

#________________________________________________________________
