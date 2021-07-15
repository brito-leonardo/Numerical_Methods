#code name: 'bisection_method.f90'
#
#description: This code aims to apply the Bisection method,
#              to find the root of a fucntion.
#
#             problem: f(x)=x**3-x-2
#
#
#reference: https://en.wikipedia.org/wiki/Bisection_method
#
#date: July 10th, 2021
#author: Leonardo Brito
#Physics Phd Student
#Institute of Physics (IFUSP)
#University of São Paulo (USP)
#_______________________________________________________________________

#------------------------------------------------------------------

#define the function

def f(x):
    func=x**3-x-2.0
    return func    

#------------------------------------------------------------------
#___________________________________________________________________

#main code

#enter an guess interval [a,b] (we should have some insight of a good interval)
a=1.0
b=2.0
#this choice can be understood if we note
#f(a=1)=-2 and f(b=2)=+4
#there is a root within the interval [1,2]
#because the function change its sign.

#set a convergence criterion
epsilon=1.e-6
#set an arbitrary error
err=1.0
#bisection algorithm (see the reference)
while(err>epsilon):
    #calculate the middle point
    c=(a+b)/2
    if(f(a)*f(c)<0.0): #if 'true' there a change of sign in the interval [a,c]
        b=c #if 'true' we make the interval closer to 'a'
    else:
        a=c #if 'false' we make the interval closer to 'b'
    err=abs(f(c))    

#the root is the last middle point 'c' found
x=c
print('The root is : ',x)

#___________________________________________________________________
