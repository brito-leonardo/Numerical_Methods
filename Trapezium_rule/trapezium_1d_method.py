#code name: 'trapezium_1d_method.f90'
#
#description: This code aims to apply the Trapezium rule to integrate
#             a one-dimensional function f(x)
#
#     scheme:
#               result=dx*[f(x1)+f(x2)]/2
#                        + dx*[f(x2)+f(x3)]/2
#                        + dx*[f(x3)+f(x4)]/2
#                        + ...
#                        + dx*[f(xN-1)+f(xN)]/2
#
#     problem: we want to find integral of the function
#     defined by
#              f(x)=sin(x)**2
#     in the interval x=[0,pi].
#     We know the exact answer: Result=pi/2=1.570796327...
#
#reference: https://en.wikipedia.org/wiki/Trapezoidal_rule
#
#date: July 26th, 2021
#author: Leonardo Brito
#Physics Phd Student
#Institute of Physics (IFUSP)
#University of São Paulo (USP)
#_______________________________________________________________________
import numpy as np

#----------------------------------------------
# function to define the framework, i.e,
# the spatial-step length 'dx', spatial grid 'x'
# and the fucntion 'f' in the grid 
#
# input: N0: size of the arrays
#        a0: initial point of the interval
#        b0: final point of the interval
#
# output: dx0: spatial-step length
#         x0: x-space
#         f0: fucntion denined int he x-space    
#
def func(N0,a0,b0):
    dx0=(b0-a0)/(N0-1) #define the spatial-step length
    x0=np.zeros(N0) # define the spatial array and its size
    f0=np.zeros(N0) # define the function array and its size
    x0[0]=a
    for i in range(1,N0):
        x0[i]=i*dx0 
    for i in range(0,N0):
        f0[i]=np.sin(x0[i])**2
    return (dx0,x0,f0)    

#----------------------------------------------
#
# function to run the trapezium rule and
# get the integration result
#
# input: N0: size of the arrays
#        dx0: spatial-step length
#        f0: fucntion defined in the x-space    
#
# output: r: result

def trapezium(N0,dx0,f0):
    fx0=np.zeros(N0) #define an auxiliary array of trapezium factors
    fx0[0]=1
    for i in range(1,N0-1):
        fx0[i]=2
    fx0[N0-1]=1    
    r0=0.0
    for i in range(N0):
        r0=r0+(fx0[i]*dx0*f0[i])/2.0
    return r0

#----------------------------------------------
#______________________________________________

# main code

#define the arrays size
Nx=100
#define the integration interval x=[a,b]=[0,pi]
a=0.0
b=np.pi

#return the variables
(dx,x,f)= func(Nx,a,b)
#return the integration result by trapezium rule
r= trapezium(Nx,dx,f)

print('The trapezium rule returns the integration')
print('of the fucntion f(x)=sin(x)**2 ')
print('in the interval [a,b]=[0,pi] as')
print('the result: ', r)
print('while the exact solution is: ', np.pi/2)

#________________________________________________
