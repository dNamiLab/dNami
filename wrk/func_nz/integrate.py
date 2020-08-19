import numpy as np
from scipy.interpolate 	import interp1d

def wrapper(f,xin,y0,**args):
	if 'method' in args:
		method = args['method']
	else:
		method = RK4
		
	sol = method(f,xin,y0)
	return sol.copy()


def RK4(f,xin,y0,dtype=np.float64):
	nRK = 4
	aij = np.float64([[0.,0.,0.,0.],[.5,0.,0.,0.],[0.,.5,0.,0.],[0.,0.,1.,0.]])
	bj 	= np.float64([1./6.	, 1./3. , 1./3. , 1./6.		])
	cj	= np.float64([0.	, 1./2. , 1./2. , 1.		])
	nt 	= len(xin)
	ny 	= len(y0)
	y = np.zeros((nt,ny),dtype=dtype)
	y[0,:] = y0
	for i in range(1,nt):
		kj 	= np.zeros((nRK,ny),dtype=dtype)
		h 	= xin[i]-xin[i-1]
		y[i,:] = y[i-1,:]
		for j in range(nRK):
			Ka = np.zeros(ny,dtype=dtype)
			for k in range(nRK):
				Ka += aij[j,k]*kj[k,:]
			kj[j,:] = f(xin[i-1]+cj[j]*h,y[i,:]+h*Ka)
		Kb = np.zeros(ny,dtype=dtype)
		for k in range(nRK):
				Kb[:] += bj[k]*kj[k,:]
		y[i,:] += h*Kb
	if len(y0) == 1:return y[:,0].copy()
	return y.copy()



