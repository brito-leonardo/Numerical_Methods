#code name: 'newton_FD_method.py'
#
#description: This code aims to apply the Newton's method,
#              to find roots.
#
#             problem: we want to find x that satisfies cos(x)=x**3
#
#             we define:
#                   f(x)=cos(x)-x**3
#             suppose we don't know its first derivative,
#             we will use finite differences (FD) method
#             to calculate the first derivative
#
#reference: https://en.wikipedia.org/wiki/Newton%27s_method
#
#date: July 9th, 2021
#author: Leonardo Brito
#Physics Phd Student
#Institute of Physics (IFUSP)
#University of SÆo Paulo (USP)
#_______________________________________________________________________
import numpy as np
#----------------------------------------------------------
def f(x):
    func=np.cos(x)-x**3
    return func
#----------------------------------------------------------
# here we approximate the first derivative by finite differences method
# see: https://en.wikipedia.org/wiki/Finite_difference_method
#
def fd(x):
    #set an arbitrary step
    dx=0.0001
    x1=x
    x2=x1+dx
    func=(f(x2)-f(x1))/dx
    return func
#----------------------------------------------------------
#__________________________________________________________

# main code

#solution guess
x0=0.5
#convergence criterium
epsilon=1.e-6
#set starting variables
xold=x0
x=x0
err=1.0
#Newton's scheme
while(err>epsilon):
    x=x-f(x)/fd(x)
    err=abs(x-xold)
    xold=x

print('The solution is x: ', x)
#_____________________________________________________________