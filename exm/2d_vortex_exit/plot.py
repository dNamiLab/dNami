import matplotlib.pyplot as plt
import glob
import numpy as np
from post_io import read_restart, load_ax, read_restart_wshell
import matplotlib
import os 

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
    rho = f[:,:,0]  - 1.  # extract density field fluctuations
    v   = f[:,:,2]   # extract vertical velocity field 

    if k == 0:
	
        # -- Set initial contours levels
        lvls = np.linspace(v.min(), v.max(), 12)

        # -- Create figure
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_axes([0.,0.,1.,1.])

        ar = 1. # figure aspect ratio

        # -- Create figure
        im = ax.imshow(rho.T, cmap='Greys', origin='lower', vmin = -0.025, vmax= 0.025)
        cntr = ax.contour(v.T, levels=lvls)
        ax.set_aspect(ar)

    else:
        for coll in cntr.collections: 
            ax.collections.remove(coll) 
        cntr = ax.contour(v.T, levels=lvls) 
        im.set_array(rho.T)

    # -- Save figure
    fname = 'pics/' + str(k).zfill(4) + '.png' 
    plt.savefig(fname,dpi=200)

plt.close()
