
import sympy as sym
from sympy import pprint as pp
from sympy import simplify as simply
from sympy import Symbol as Sy

import numpy as np
import sys

sym.init_printing(use_latex=True,wrap_line=False)


# dNami convention :
# (x,y,z)       --> Computational Domain Coordinate System (arbitrary curvilinear space)
# (xi,eta,zeta) --> Original Cartesian Coordinate System
'''
x,y,z  = sym.symbols(['x','y','z'],Real=True)
xi     = sym.Function('xi'  )(x,y,z)
eta    = sym.Function('eta' )(x,y,z)
zeta   = sym.Function('zeta')(x,y,z)

A = sym.Matrix([xi,eta,zeta])

Jm1  = sym.gcd(tuple(A.jacobian([x,y,z]).inv()))

J      =  A.jacobian([x,y,z]).inv()


x_xi   =  J[0,0]/Jm1

pp(x_xi)

x_eta  =  J[0,1]/Jm1
x_zeta =  J[0,2]/Jm1

y_xi   =  J[1,0]/Jm1
y_eta  =  J[1,1]/Jm1
y_zeta =  J[1,2]/Jm1

z_xi   =  J[2,0]/Jm1
z_eta  =  J[2,1]/Jm1
z_zeta =  J[2,2]/Jm1

'''
# 2D 

xi,eta,t  = sym.symbols(['xi','eta','t'],Real=True)
x = sym.Function('x',Real=True)(xi,eta,t)
y = sym.Function('y',Real=True)(xi,eta,t)


u = sym.Function('u')(xi,eta,t)
v = sym.Function('v')(xi,eta,t)
T = sym.Function('T')(xi,eta,t)

ucont = sym.Function('u')(x,y,t)
vcont = sym.Function('v')(x,y,t)
Tcont = sym.Function('T')(x,y,t)

mu     = sym.Symbol('mu')
mub   = sym.Symbol('mu_b')
kappa = sym.Symbol('kappa')

SymGradU = []

Vel = [u,v]
Dir = [xi,eta]

for i in range(len(Vel)):
	Line = []
	for j in range(len(Dir)):
		Line.append(0.5*(sym.Derivative(Vel[i],Dir[j])+sym.Derivative(Vel[j],Dir[i])))
	SymGradU.append(Line)	

S  = sym.Matrix(SymGradU)
I2 = sym.eye(2)
div = 0.
for  i in range(2):
	div = div + sym.Derivative(Vel[i],Dir[i])

Tau = 2.0*mu*(S-1./3.*(div)*I2) #+ mub*div*I3

Tau = Tau.subs(u,ucont)
Tau = Tau.subs(v,vcont)

Beta = sym.Matrix([kappa*T.diff(xi)   + u*Tau[0,0]+v*Tau[0,1],
		           kappa*T.diff(eta)  + u*Tau[1,0]+v*Tau[1,1]])

Beta = Beta.subs(T,Tcont)

Fldif_x = [0,
		   x.diff(xi)*Tau[0,0]+x.diff(eta)*Tau[0,1],
		   x.diff(xi)*Tau[1,0]+x.diff(eta)*Tau[1,1],
		   x.diff(xi)*Beta[0] +x.diff(eta)*Beta[1] ]

Fldif_y = [0,
		   y.diff(xi)*Tau[0,0]+y.diff(eta)*Tau[0,1],
		   y.diff(xi)*Tau[1,0]+y.diff(eta)*Tau[1,1],
		   y.diff(xi)*Beta[0] +y.diff(eta)*Beta[1] ]


x_xi,x_eta = sym.symbols(['x_xi','x_eta'])
y_xi,y_eta = sym.symbols(['y_xi','y_eta'])

dNamiFlx_x = []
dNamiFlx_y = []

from genNSBC import sympy2dNami

for eq in Fldif_x:

	tmpeq = simply(eq)
	tmpeq = tmpeq.subs({x.diff(xi):x_xi, x.diff(eta):x_eta,
			            y.diff(xi):y_xi, y.diff(eta):y_eta,}, simultaneous=True,evaluate=False)
	tmpeqSTR = str(tmpeq).replace('(xi, eta, t)','')

	tmpeqSTRSym = sym.sympify(tmpeqSTR,evaluate=False)
	dNamiFlx_x.append(sympy2dNami(tmpeqSTRSym).replace('_wp_wp','_wp').replace('[','{').replace(']','}'))

for eq in Fldif_y:

	tmpeq = simply(eq)
	tmpeq = tmpeq.subs({x.diff(xi):x_xi, x.diff(eta):x_eta,
			            y.diff(xi):y_xi, y.diff(eta):y_eta,}, simultaneous=True,evaluate=False)
	tmpeqSTR = str(tmpeq).replace('(xi, eta, t)','')

	tmpeqSTRSym = sym.sympify(tmpeqSTR,evaluate=False)
	dNamiFlx_y.append(sympy2dNami(tmpeqSTRSym).replace('_wp_wp','_wp').replace('[','{').replace(']','}'))



# init = 2.0*mu*(0.2557968721238*x_eta*vcont.diff(x)+ucont)#-0.666666666666667*x_xi*ucont.diff(x))
# init = str(init).replace('(xi, eta, zeta, t)','')
# init = sym.sympify(init,evaluate=False)
# print(sympy2dNami(init))

