# import test_fftnd

import mkl_fft as mklfft
#import numpy.fft as mklfft
import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys


def rfftn(f):
	fft = mklfft.rfftn_numpy(f,f.shape)
	#fft = mklfft.rfftn(f,f.shape)
	return fft

def compute_grad3d(f,fft_f=None):


	n   = f.shape[0]
	n_half = n/2
	n_2 = int((f.shape[0])/2+1)

	k   = np.arange(n_2  ,dtype='float64')	
	kz  = np.arange(n_2  ,dtype='float64')

	kx  = np.concatenate((k,-np.flip(k)[1:-1])) 
	ky  = kx

	kx_xy = kx.reshape(n,1 ,1  ) * np.ones((1  ,n ,n_2))*2.0*m.pi
	ky_xy = ky.reshape(1,n ,1  ) * np.ones((n  ,1 ,n_2))*2.0*m.pi
	kz_xy = kz.reshape(1,1 ,n_2) * np.ones((n  ,n ,1  ))*2.0*m.pi
	
	
	
	if fft_f == None: fft_f = mklfft.rfftn_numpy(f,f.shape)

	out      =  fft_f.copy()
	out.real = -kx_xy*fft_f.imag 
	out.imag =  kx_xy*fft_f.real		
	gradf_x  =  mklfft.irfftn_numpy(out,s=f.shape)

	
	out      =  fft_f.copy()
	out.real = -ky_xy*fft_f.imag 
	out.imag =  ky_xy*fft_f.real		
	gradf_y  =  mklfft.irfftn_numpy(out,s=f.shape)

	
	out      =  fft_f.copy()
	out.real = -kz_xy*fft_f.imag 
	out.imag =  kz_xy*fft_f.real		
	gradf_z  =  mklfft.irfftn_numpy(out,s=f.shape)
	

	return [gradf_x,gradf_y,gradf_z]

def power_spectra(u,v,w,fft_u=None,fft_v=None,fft_w=None):

	n   = u.shape[0]
	n_2 = int(n/2+1)

	scale = n**2*n_2

	if fft_u == None:
		fft_u = mklfft.rfftn_numpy(u,u.shape)/scale
		#fft_u = mklfft.rfftn(u,u.shape)/scale
	if fft_v == None:
		fft_v = mklfft.rfftn_numpy(v,v.shape)/scale
		#fft_v = mklfft.rfftn(v,v.shape)/scale
	if fft_w == None:
		fft_w = mklfft.rfftn_numpy(w,w.shape)/scale
		#fft_w = mklfft.rfftn(w,w.shape)/scale

	# del u,v,w	
	
	k   = np.arange(n_2  ,dtype='float64')	
	kz  = np.arange(n_2  ,dtype='float64')

	kx  = np.concatenate((k,-np.flip(k)[1:-1])) 
	ky  = kx

	# kx_xy = kx.reshape(n,1 ,1  ) * np.ones((1  ,n ,n_2))#*2.0*m.pi
	# ky_xy = ky.reshape(1,n ,1  ) * np.ones((n  ,1 ,n_2))#*2.0*m.pi
	# kz_xy = kz.reshape(1,1 ,n_2) * np.ones((n  ,n ,1  ))#*2.0*m.pi   

	kx_xy = k.reshape(n_2,1 ,1  ) * np.ones((1  ,n_2 ,n_2))#*2.0*m.pi
	ky_xy = k.reshape(1,n_2 ,1  ) * np.ones((n_2  ,1 ,n_2))#*2.0*m.pi
	kz_xy = k.reshape(1,1 ,n_2)   * np.ones((n_2  ,n_2 ,1  ))#*2.0*m.pi   

	fft_u = fft_u[0:n_2,0:n_2,0:n_2]
	fft_v = fft_v[0:n_2,0:n_2,0:n_2]
	fft_w = fft_w[0:n_2,0:n_2,0:n_2]

	power = (fft_u*np.conj(fft_u) + fft_v*np.conj(fft_v) + fft_w*np.conj(fft_w)).real	

	# print("sum spectra",np.sum(power)) #already half of the sum over -kinf,+kinf (see scale... )

	knorm = np.sqrt(kx_xy**2 + ky_xy**2  + kz_xy**2)
	
	knorm = knorm.reshape(np.size(knorm))
	power = power.reshape(np.size(power))

	perm        = np.argsort(knorm)
	knorm       = knorm[perm]
	pow_ordered = power[perm]


	pow_spct = []
	# Ntot = n_2**3
	# knorm = knorm[0:Ntot]
	# pow_ordered = pow_ordered[0:Ntot]
	# print(knorm[0:25])
	# print(pow_ordered[0:25])

	for kind in k:
		ind = np.nonzero(knorm < kind + 0.5)
		# print(kind,np.size(pow_ordered),np.size(pow_ordered[ind]),ind[0][-1])
		# print(np.size(pow_ordered[ind]),np.sum(pow_ordered[ind]),kind)
		# pow_spct.append( np.sum(pow_ordered[ind])/np.size(pow_ordered[ind])*(m.pi*(2.0*m.pi*kind)**2)) # kind**2 for volume of the sphere in 3D...
		pow_spct.append( np.sum(pow_ordered[ind])/np.size(pow_ordered[ind])*(m.pi*(kind)**2)) # kind**2 for volume of the sphere in 3D...
	
		# clean
		knorm       = knorm[ind[0][-1]+1:]
		pow_ordered = pow_ordered[ind[0][-1]+1:]

	return [k,np.asarray(pow_spct)	]
	

# def compute_hlmtz3d(u,v,w):

# 	n   = f.shape[0]
# 	n_half = n/2
# 	n_2 = int((f.shape[0])/2+1)

# 	k   = np.arange(n_2  ,dtype='float64')	
# 	kz  = np.arange(n_2  ,dtype='float64')

# 	kx  = np.concatenate((k,-np.flip(k)[1:-1])) 
# 	ky  = kx

# 	kx_xy = kx.reshape(n,1 ,1  ) * np.ones((1  ,n ,n_2))*2.0*m.pi
# 	ky_xy = ky.reshape(1,n ,1  ) * np.ones((n  ,1 ,n_2))*2.0*m.pi
# 	kz_xy = kz.reshape(1,1 ,n_2) * np.ones((n  ,n ,1  ))*2.0*m.pi
	

# 	fft_vel = [np.ones( (n,n,n_2) )*complex(0.0,0.0),
# 			   np.ones( (n,n,n_2) )*complex(0.0,0.0),
# 			   np.ones( (n,n,n_2) )*complex(0.0,0.0)]	

# 	fft_vel[0] = mklfft.rfftn_numpy(u,u.shape)
# 	fft_vel[1] = mklfft.rfftn_numpy(v,v.shape)
# 	fft_vel[2] = mklfft.rfftn_numpy(w,v.shape)

# # u velocity 

# 	out      =  fft_vel[0].copy()

# 	out.real =  -ky_xy*fft_vel[2].imag + kz_xy*fft_vel[1].imag
# 	out.imag =   ky_xy*fft_vel[2].real - kz_xy*fft_vel[1].real
	
# 	gradf_x  =  mklfft.irfftn_numpy(out,s=f.shape)

	
# 	out      =  fft_f.copy()
# 	out.real = -ky_xy*fft_f.imag 
# 	out.imag =  ky_xy*fft_f.real		
# 	gradf_y  =  mklfft.irfftn_numpy(out,s=f.shape)

	
# 	out      =  fft_f.copy()
# 	out.real = -kz_xy*fft_f.imag 
# 	out.imag =  kz_xy*fft_f.real		
# 	gradf_z  =  mklfft.irfftn_numpy(out,s=f.shape)
	

# 	return [gradf_x,gradf_y,gradf_z]	


# if __name__ == '__main__':
#     main()
