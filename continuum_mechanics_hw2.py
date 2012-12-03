#!/usr/bin/python

print '\n\section*{Python solver output}\n'

from sympy import *

#from printer import Printer
from sympy.printing.latex import * 
settings = {'mat_str': 'bmatrix','mat_delim':''}
myprinter = LatexPrinter(settings)

#Define our coordinates
X1=Symbol('X1')
X2=Symbol('X2')
X3=Symbol('X3')
t=Symbol('t')
X=Matrix([X1,X2,X3])

x1=Symbol('x1')
x2=Symbol('x2')
x3=Symbol('x3')
x=Matrix([x1,x2,x3])

#Our specific problem (although fairly general)
print '\n \\begin{framed}'
print '\nOur velocity field is given as:'
V=Matrix([X2*t,X2,X3+1])
pointPt={X1:5,X2:7,X3:3,t:1}
vectorm={x1:2,x2:3,x3:4}
sigma=Matrix([[x2,t,0],[t,0,1],[0,1,x1-t]])
print '\[ V(X,t)= ', myprinter.doprint(V), '\]'
print 'with'
print '\[ P_0 = ', myprinter.doprint(X.evalf(subs=pointPt)), ' \]'
print '\[ m = ', myprinter.doprint(x.evalf(subs=vectorm)), ' \]'
print '\[ \sigma = ', myprinter.doprint(sigma), ' \]\n'
print '\\end{framed}'


print '\nThe motion is the integral of velocity:'
FI=integrate(V,t)+X
FIdict={x1:FI[0], x2:FI[1], x3:FI[2]} 
print '\[ \\varphi(X,t) = \int\limits_0^t V(X,\\tau) \mathrm{d} \\tau + X = ', myprinter.doprint(FI), '\]'

print '\n\\textbf{(a)} Therefore the deformation gradient is:'
F = FI.jacobian(X)
print '\[ F = ', myprinter.doprint(F), '\]\n'
print '\[ F |_{(X=P_0)} = ',myprinter.doprint(F.evalf(subs=pointPt)), '\]'

print '\nThe inverse motion is: '
invFIdict = solve([ FI[0]-x[0],FI[1]-x[1],FI[2]-x[2] ], [ X[0],X[1],X[2] ])
invFI=Matrix([ invFIdict[X1], invFIdict[X2], invFIdict[X3] ])
print '\[ \\varphi^{-1}(x,t) = ', myprinter.doprint(invFI), '\]\n'

print '\n Thus: '
pointptVec=FI.subs(pointPt)
pointpt={x1:pointptVec[0],x2:pointptVec[1],x3:pointptVec[2],t:pointPt[t]}
print '\[ p_0 = ', myprinter.doprint( x.subs(pointpt) ), '\]'

print '\n \\textbf{(b)} Which gives for the Eulerian velocity field as'
v=V.subs(invFIdict)
print '\[ v (x,t) = V(\\varphi^{-1}(x,t),t)', myprinter.doprint(v), '\]'
print '\[ v|_{(x=p_0)} = ', myprinter.doprint( v.subs(pointpt) ), '\]'

print '\n \\textbf{(c)} The acceleration fields:'
A=V.diff(t)
print '\[ A|_{P_0} = ', myprinter.doprint( A.subs(pointPt) ), '\]'
a=A.subs(invFIdict)
print '\[ a|_{p_0} = A(\\varphi^{-1}(x,t),t) = ', myprinter.doprint( a.subs(pointpt) ), '\]'


print '\n \\textbf{(d)} The Eulerian tensors are:'
l=v.jacobian(x)
print '\[ l(x,t) = \\frac{\partial v}{\partial x} = \mathrm{grad} v = ', myprinter.doprint( l ), '\]'
latp=l.subs(pointpt)
print '\[ l|_{p_0} = ', myprinter.doprint(latp), '\]'
print '\nthe symmetric and non-symmetric decomposition is:'
datp=0.5*(latp+latp.transpose())
print '\[ d = ',myprinter.doprint( datp  ),'\]'
watp=0.5*(latp-latp.transpose())
print '\[ w = ',myprinter.doprint( watp  ),'\]'


print '\n \\textbf{(e)} The stretching along (non-unit) $m$ '
streching=(x.transpose()*datp*x).subs(vectorm) * ( (x.transpose()*x).subs(vectorm)[0] )**(-1)
print '\[ \\dot{\\lambda_m}/\\lambda_m = \\frac{m^T d m}{(\sqrt{m^T m})^2} =', myprinter.doprint( streching[0].evalf() ), '\]'


print '\n \\textbf{(f)} The material time-derivative of the Euler-Almansi tensor is:'
c=(F*F.transpose()).inv()
e=( Matrix(eye(3)) - c )/2
dote=(l.transpose()*c+c*l.transpose())
print '\[ \dot{ e }= \\frac{1}{2}\left(l^T c + cl^T\\right)=', myprinter.doprint( dote.subs(pointpt).evalf() ), '\]'


print '\n \\textbf{(g)} First and second Piola-Kirchoff:'
J=F.det().subs(pointPt)
print '\[ J = \mathrm{det}F  \]'
sigmaatp=sigma.subs(pointpt)
FatP=F.subs(pointPt)

P=J*sigmaatp*(FatP.transpose()).inv()
S=J*FatP.inv()*sigmaatp*(FatP.transpose()).inv()

print '\[ P = J \sigma {F^{T}}^{-1} = ', myprinter.doprint(P), '\]'
print '\[ S = J F^{-1} \sigma \left(F^{T}\\right)^{-1} = ', myprinter.doprint(S), '\]'


print '\n \\textbf{(h)} The invariants:'
print '\[ I_S = \mathrm{tr} S = ', myprinter.doprint(S.trace().evalf()), '\]'
iis=( (S.trace())**2 - (S**2).trace() )/2
print '\[ II_S = \\frac{ I_S^2 -\mathrm{tr} S^2 }{2}= ', myprinter.doprint( iis.evalf()  ), '\]'
print '\[ III_S = \mathrm{det} S = ', myprinter.doprint(S.det().evalf()), '\]'



print '\n \\textbf{(i)} The Jaumann rate:'
Jsigma = ( sigma.subs(FIdict) ).diff(t) + sigmaatp * watp - watp * sigmaatp
print '\[ \overset{\Delta}{\sigma} = \dot{\sigma}(X,t)+\sigma w - w \sigma = '
print myprinter.doprint( Jsigma.subs(pointPt) ), '\]'







#Goodbye
print '\\end{document}'

