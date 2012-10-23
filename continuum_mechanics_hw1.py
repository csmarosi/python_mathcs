#!/usr/bin/python

#Ninja! comment

from sympy import *

#Define our coordinates
X1=Symbol('X1')
X2=Symbol('X2')
X3=Symbol('X3')
X=Matrix([X1,X2,X3])

x1=Symbol('x1')
x2=Symbol('x2')
x3=Symbol('x3')
x=Matrix([x1,x2,x3])

#Our motion
FI=Matrix([X1+X3**2,X2*X1,X3])
pointP={X1:4,X2:2,X3:6}

print "a) The deformation gradient:"
F=FI.jacobian(X)
F_at_P=F.evalf(subs=pointP)
print F, "\n" 

print "b) The inverse motion:"
#Note: the trivial attempt, and even this resulted an error. Fortunately, this after an error in the specification, this step was given explicitly.
#invFI=solve([ x[0]-FI[0], x[1]-FI[1], x[2]-FI[2] ], [ X[0], X[1], X[2] ])
invFI=Matrix([x1-x3**2, x2/(x1-x3**2), x3])
print invFI
print "Inverse deformation gradient:"
print invFI.jacobian(x), "\n"

print "c) right Cauchy-Green:"
C=F.transpose()*F
C_at_P=C.evalf(subs=pointP)
print C_at_P
print "Piola:"
print C.inv().evalf(subs=pointP)
print "Green-Lagrange strain:"
E=(C-eye(3))/2.0
print E.evalf(subs=pointP)
print "left Cauchy-Green deformation:"
b=F*F.transpose()
print b.evalf(subs=pointP)
print "Cauchy deformation:"
print b.inv().evalf(subs=pointP)
print "Euler-Almansi strain:"
e=(eye(3)-b.inv())/2.0
print e.evalf(subs=pointP), "\n"

print "d) material displacement field:"
print (FI-X)
print "spatial one:"
print (invFI-x)
print "displacement vector at P"
print (FI-X).evalf(subs=pointP),"\n"

print "e) volume ratio:"
J=F.det()
print J
print "Numerically:"
J_at_P=J.evalf(subs=pointP)
print J_at_P
print "... isochoric part:"
print F_at_P*(J_at_P**(1.0/3))
print "volumetric:"
print eye(3)*(J_at_P**(1.0/3))

#Here, we will work with...
from numpy import *
#note also that the subsequent questions explicitly defines that we should calculate stuff only for P

print "f) principal stretch:"
U, s, V = linalg.svd(C_at_P) 
for i in range(3):
	print "lambda_%d= " % (i+1) , sqrt(s[i])

#note: The lambda_2 value was correct!

