import matplotlib.pyplot as plt
import glob
import numpy as np
from post_io import read_restart, load_ax, read_restart_wshell

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rc
fnt    =  10      
rc('text', usetex=True)
rc('text.latex',preamble=r'\usepackage{bm}')
font = {'family': 'serif',
    'weight': 'normal',
    'size': fnt+2,
    }
fontt = {'family': 'serif',
	'color':  'white',
    'weight': 'normal',
    'size': fnt+2,
    }
matplotlib.rc('font', **font)

# ==============================================================

# -- Plot initial restart field
x,nx = load_ax('./out/axes.bin')

rho0 = np.float64(1.)
eps  = np.float64(1e-2)

plate = True 
ntfinal = 10
iPurp   = True 
iRdBu   = True 

fld = './restarts/restart_*'
print('Looking for ...', fld)
from natsort import natsorted
folders = natsorted(glob.glob(fld))
print('Found nfiles:', len(folders))

# --
x0 = 20.
def sol(t):
    xll    = x -x0 
    k1     = 1.
    k2     = np.sqrt(2.) 
    eta1_0 = 0.
    eta2_0 =  2.*np.sqrt(2)
    eta1  = k1*xll[:] - k1**3*t+ eta1_0 
    eta2  = k2*xll[:] - k2**3*t+ eta2_0 
    f     = 1. + np.exp(eta1) + np.exp(eta2) + (k1-k2)**2/(k1+k2)**2 * np.exp(eta1+eta2)
    logf  = np.log(f)
    gradlogf = np.gradient(logf, xll[:])
    ua = 2.  * np.gradient(gradlogf, xll[:]) 
    return ua

#
#plt.plot(xloc, rho, ls='-', color='k')
#plt.show()
#exit()
print('Loading ...', folders[0])

# -- INITIAL SOLUTION

t,f = read_restart_wshell(folders[0])

fig = plt.figure(figsize=(4.1,4.1))
ax0 = fig.add_subplot(111) 


ax0.plot(x - x0,f, ls=' ', color='k', marker='o', ms=2., label=r'$\texttt{dNami}$')
ax0.plot(x - x0,sol(0.), ls='-', color='k', label=r'Analytical solution')

ax0.set_xlabel(r'$x/L$')
ax0.set_ylabel(r'$u$')
ax0.set_xlim([-x0,x0])
ax0.set_ylim([-0.1,1.])

plt.savefig('kdv_initial.png',dpi=250)

plt.close()

# -- FINAL SOLUTION
t,f = read_restart_wshell(folders[-1])

fig = plt.figure(figsize=(4.1,4.1))
ax0 = fig.add_subplot(111) #xt diags

x0 = 20.

ax0.plot(x - x0,f, ls=' ', color='k', marker='o', ms=2., label=r'$\texttt{dNami}$')
ax0.plot(x - x0,sol(t), ls='-', color='k', label=r'Analytical solution')

ax0.set_xlabel(r'$x/L$')
ax0.set_ylabel(r'$u$')
ax0.set_xlim([-x0,x0])
ax0.set_ylim([-0.1,1.])
ax0.legend(loc='upper left')

plt.savefig('kdv_final.png',dpi=250)

plt.close()
