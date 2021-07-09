#code name: 'sor_method.py'
#
#description: This code aims to solve a linear matrix problem by
#             Sucessive-Over relaxation (SOR) method
#
#             problem:
#                      A*x=b
#
#     consider the linear system
#      4x1 - x2 - 6x3 + 0x4 = 2
#     -5x1 - 4x2 + 10x3 + 8x4 = 21
#      0x1 + 9x2 + 4x3 - 2x4 = -12
#      1x1 + 0x2 - 7x3 + 5x4 = -6
#
#     in the matrix form, it becomes
#
#       ( 4 -1  -6  0)      (x1)     (2)
#       (-5 -4  10  8)      (x2)     (21)
#    A= ( 0  9   4 -2)  ; x=(x3) ; b=(-12)
#       ( 1  0  -7  5)      (x4)     (-6)
#
#reference: https://en.wikipedia.org/wiki/Successive_over-relaxation
#
#date: July 9th, 2021
#author: Leonardo Brito
#Physics Phd Student
#Institute of Physics (IFUSP)
#University of São Paulo (USP)
#_______________________________________________________________________
import numpy as np

#-------------------------------------------------------

# Sucessive-Over relaxation (SOR) scheme based on
# the algorith provided by:
# https://en.wikipedia.org/wiki/Successive_over-relaxation
def sor_scheme(A0,x0,b0,omega0,epsilon0):
    err=1.0
    itera=0
    x=x0
    N0=len(x0)
    while(err>epsilon0):
        for i in range(N0):
            sigma=0.0
            for j in range(N0):
                if j!=i:
                    sigma=sigma+A0[i,j]*x[j]
            x[i]=(1.0-omega0)*x[i]+omega0*(b0[i]-sigma)/A0[i,i]
        err_vec=np.matmul(A0,x)-b0
        err=err_vec.max()
        itera=itera+1
    return x


#-------------------------------------------------------

#________________________________________________________

#main code


#enter the matrix 'A'
A=np.array([[4.0,-1.0, -6.0, 0.0],
            [-5.0, -4.0, 10.0, 8.0],
            [0.0, 9.0, 4.0, -2.0],
            [1.0, 0.0, -7.0, 5.0]])
#enter the vector 'b'
b=np.array([2.0, 21.0, -12.0, -6.0])
#enter a guess solution 'x'
xi=np.zeros(4)

print('The initial conditions are')
print('')
print('A: ', A)
print('')
print('b: ', A)
print('')
print('xi: ', xi)
print('')
print('')

#set SOR relaxation parameter
omega=0.5
#set the convergence criterium
epsilon=1.e-8
#find the solution
x=sor_scheme(A,xi,b,omega,epsilon)

print('The solution is:')
print('')
print('x: ', x)
#_______________________________________________________