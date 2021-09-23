#code name: 'euler_method.py'
#
#description: This code aims to apply the Euler's method,
#             to solve a ordinary differential
#             equation (ODE).
#
#             problem: we want to find y(4)
#
#             we have the ODE:
#                   y'=y, y(0)=1
#
#
#reference: Wikipedia contributors. (2021, September 19). Euler method.
#           In Wikipedia, The Free Encyclopedia. Retrieved 15:59,
#           September 22, 2021,
#           from https://en.wikipedia.org/w/index.php?title=Euler_method&oldid=1045231921
#
#date: September 23rd, 2021
#author: Leonardo Brito
#Physics Phd Student
#Institute of Physics (IFUSP)
#University of São Paulo (USP)
#_______________________________________________________________________

#------------------------------------

#function y'(t,y)=f(t,y), f(t,y)=y

def f(y):
    f=y
    return f
#------------------------------------
#_______________________________________________________________________
    
#main code
    
h=1.0 #step-size
t=0.0 #initial instant
y=1.0 #initial value
tf=4.0
Nt=int(tf/h)
#print the initial values
print(t,'',y)
#evolve the system by Euler's method
for i in range(1,Nt+1):
    t=i*h
    y=y+h*f(y)
    print(t,'',y)
    
#_______________________________________________________________________
    