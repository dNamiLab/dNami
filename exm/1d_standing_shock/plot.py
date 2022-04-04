import matplotlib.pyplot as plt
import glob
import numpy as np
from post_io import read_restart, load_ax, read_restart_wshell

import matplotlib
matplotlib.use('Agg')
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

# -- restarts
try:
    import os
    os.mkdir('./pics/')
except FileExistsError:
    pass

from natsort import natsorted

fld = './restarts/restart_*'
print('Looking for ...', fld)
folders = natsorted(glob.glob(fld))
print('Found nfiles:', len(folders))

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)

for k, ff in enumerate(folders[:50]): 
    
    print('Doing file ...', ff)

    t,f = read_restart_wshell(ff)
    rho = f[:,0]

    if k == 0:
        x = np.linspace(0.,1.,rho.shape[0])
        line, = ax.plot(x,rho, ls='-', color='black')
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$\rho/\rho_0$')
        ax.set_xlim([0.,1.])
        ax.set_ylim([0.85,2.20])
    else:
        line.set_ydata(rho)

    fname = './pics/'+ str(k).zfill(4) + '.png'

    plt.savefig(fname,dpi=200)

plt.close()
