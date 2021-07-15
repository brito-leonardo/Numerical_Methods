#code name: 'cg_method.f90'
#
#description: This code aims to apply the Conjugate gradient (CG) method,
#              to solve a linear matrix problem.
#
#             problem: A*x=b
#
#     A=(4  1)  ; x=(x1) ; b=(1)
#       (1  3)      (x2)     (2)
#
#reference: https://en.wikipedia.org/wiki/Conjugate_gradient_method
#
#date: July 9th, 2021
#author: Leonardo Brito
#Physics Phd Student
#Institute of Physics (IFUSP)
#University of São Paulo (USP)
#_______________________________________________________________________

import numpy as np

#-----------------------------------------------------------
#cg method
def cg(A,x0,b,epsilon):
    #set the initial guess
    x=x0
    r0=b-np.matmul(A,x)
    p=r0
    #set an arbitrary initial error
    err=1.0
    while(err>epsilon):
        alpha=np.dot(r0,r0)/np.dot(p,np.matmul(A,p))
        x=x+alpha*p
        r=r0-alpha*np.matmul(A,p)
        #calculate the error
        err=abs(r.max())
        #update the variables
        beta=np.dot(r,r)/np.dot(r0,r0)
        p=r+beta*p
        r0=r
    return x
#-----------------------------------------------------------

print('the initial conditions are:')
print('')
#define the matrix A
A=np.array([[4.0, 1.0],
            [1.0, 3.0]])
print('The matrix A is: ')            
print(A)
print('')
#define the right side of the equation
b=np.array([1.0, 2.0])
print('the vector b is: ', b)
print('')
#enter an inital solution guess
x0=np.zeros(2)
print('the solution guess x0 is: ', x0)
print('')
#set an convergence criterion
epsilon=1.e-6
#find the solution by cg method
x=cg(A,x0,b,epsilon)

print('The solution x is: ', x)
