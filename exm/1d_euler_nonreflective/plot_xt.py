import matplotlib.pyplot as plt
import glob
import numpy as np
from post_io import read_restart, load_ax, read_restart_wshell
import matplotlib

# -- Set font 
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

# ==============================================================

# -- Get axis information
x,nx = load_ax('./out/axes.bin')

# -- Get list of restarts
fld = './restarts/restart_*'
print('Looking for ...', fld)
folders = sorted(glob.glob(fld))
print('Found nfiles:', len(folders))

# -- Mean density value:
rho0 = np.float64(1.)
# -- Flow speed
u0 = np.float64(0.295803989135199)



for k, ff in enumerate(folders): 	


    print('Reading file ...', ff)

    t,f = read_restart_wshell(ff)

    if k ==0:
        # -- Create an array to store the density fluctuations
        rhop = np.zeros( (f[:,0].size, len(folders))  )

    rhop[:,k] = (f[:,0]-rho0)*1e3


# -- Create figure
fig = plt.figure(figsize=(4,4))
ax = fig.add_axes([0.17,0.1,0.78,0.85])


ar = (x[-1]-x[0])/t # figure aspect ratio

# -- Create figure
ax.imshow(rhop.T, cmap='Greys', origin='lower', extent=(0, 1, 0, t))
ax.set_aspect(ar)

ax.set_xlabel(r'x')
ax.set_ylabel(r't')

# -- Add speed lines
for i in range(10):
    ax.axline(xy1=(0.,-5. +i), slope=1./u0, color='tab:blue', ls='--', lw=1. )

# -- Limits
ax.set_xlim([0.,1.])
ax.set_ylim([0.,2.])

# -- Save figure
fname = 'xt.png'
plt.savefig(fname,dpi=200)

plt.close()
