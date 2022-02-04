import matplotlib.pyplot as plt
import glob
import numpy as np
from post_io import read_restart, load_ax, read_restart_wshell
import matplotlib
import os 

# -- Set font 
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

# ==============================================================

if not os.path.isdir('pics/'): os.mkdir('pics/')

# -- Get list of restarts
fld = './restarts/restart_*'
print('Looking for ...', fld)
folders = sorted(glob.glob(fld))
print('Found nfiles:', len(folders))

for k, ff in enumerate(folders): 

    print('Reading file ...', ff)

    t,f = read_restart_wshell(ff)
    uxmid = f[:,:,32,1] # grid quarter height velocity profile

    if k == 0:
        # -- Create figure
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_axes([0.10,0.1,0.85,0.85])

        ar = 1. # figure aspect ratio

        # -- Create figure
        im = ax.imshow(uxmid.T, cmap='Greys', origin='lower', extent=(0, 2*np.pi, 0, 2*np.pi), vmin = -1.001, vmax= 1.001)
        ax.set_aspect(ar)

        ax.set_xlabel(r'x')
        ax.set_ylabel(r'y')

        # -- Limits
        ax.set_xlim([0.,2*np.pi])
        ax.set_ylim([0.,2*np.pi])

    else:
        im.set_array(uxmid.T)

    # -- Save figure
    fname = 'pics/' + str(k).zfill(4) + '.png' 
    plt.savefig(fname,dpi=200)

plt.close()
