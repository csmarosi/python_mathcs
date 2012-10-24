#!/usr/bin/python

print 'Continuum Mechanics, solution of the 1st homework (using python)\n'

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

#Our motion and a given material point
FI=Matrix([X1+X3**2,X2*X1,X3])
pointP={X1:4,X2:2,X3:6}
#Here we calculate the deformed point
pointp={}
for i in range(3):
	pointp[x[i]]=FI[i].evalf(subs=pointP)

print "a) The deformation gradient:"
F=FI.jacobian(X)
F_at_P=F.evalf(subs=pointP)
print F, "\n" 

print "b) The inverse motion:"
#Note: the trivial attempt, and even this resulted an error. Fortunately, after an error in the specification, the result of this step was given explicitly.
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
B_at_P=C.inv().evalf(subs=pointP)
print B_at_P
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
U_disp=(FI-X)
print U_disp 
print "spatial one:"
print (x-invFI)
print "displacement vector at P:"
print (FI-X).evalf(subs=pointP),"\n"

print "e) volume ratio is:"
J=F.det()
print J
print "numerically:"
J_at_P=J.evalf(subs=pointP)
print J_at_P
print "... isochoric part:"
print F_at_P*(J_at_P**(-1.0/3))
print "volumetric part:"
print eye(3)*(J_at_P**(1.0/3)), "\n"

#Here, we will work with...
from numpy import *
#note also that the subsequent questions explicitly defines that we should calculate stuff only for P
#IMPORTANT: it seems that the numpy array multiplication is NOT matrix multiplication

print "f) principal stretch:"
U, s, V = linalg.svd(C_at_P) 
for i in range(3):
	print "\\lambda_%d= " % (i+1) , sqrt(s[i])
print ""

print "g) the stretch at the direction of M:"
M=Matrix([6,5,1])
NM=M/M.norm()
#Note: the following variable is a matrix (with singleton dimensions), and calculating non-integer power of a 
#   matrix is tricky; therefore we must explicitly define an element to use in computation 
lambdaMsqr=(NM.transpose()*(F_at_P.transpose())*F_at_P*NM)
print "\\lambda_M= ", lambdaMsqr[0]**(1.0/2), "\n"

print "h) the unit length normal of the deformed surface:"
n=F_at_P.inv().transpose()*NM/(NM.transpose()*B_at_P*NM)[0]**(1.0/2)
print n.evalf(), "\n"

print "i) the angle change:"
N1=Matrix([1,0,0])
N2=Matrix([0,1,0])
lambda1=(N1.transpose()*C_at_P*N1)[0]**(1.0/2)
lambda2=(N2.transpose()*C_at_P*N2)[0]**(1.0/2)
cosalpha=N1.transpose()*C_at_P*N2/(lambda1*lambda2)
print "cos( \\alpha )= ", cosalpha[0] 
print "that is the angle change is ", 90-(acos(cosalpha[0])/pi*180).evalf(), " [deg]\n"

print "e1) the rotational part of the deformation gradient:"
#We implement the curl here in a hackish way...
def curl_3d(A,X):
	Ax = S(A[0])
	Ay = S(A[1])
	Az = S(A[2])
	return Matrix([Az.diff(X[1])-Ay.diff(X[2]), Ax.diff(X[2])-Az.diff(X[0]), Ay.diff(X[0])-Ax.diff(X[1])])

omega=curl_3d(invFI,x)/2
print omega.evalf(subs=pointp), "\n"

print "e2) the Hencky's strain tensor:"
Htmp=Matrix(eye(3))
for i in range(3):
	Htmp[i,i]=log(s[i])/2
H=Matrix(U)*Htmp*Matrix(V)
print H

print "Finally, to check our result, we calculate the zero matrix:"
#calculate the matrix exponential with approximation
Hexp=Matrix(eye(3))
approx=35 #how many elements we calculate from the infinite sum that gives the matrix exponential
for i in range(approx):
	Hexp=Hexp*H/(approx-i)+Matrix(eye(3))
print Hexp**2-C_at_P, '\n'

print "Goodbye world!\n"

