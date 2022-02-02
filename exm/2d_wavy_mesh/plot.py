import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib 
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
from post_io import load_ax, read_restart

# ---------- Load and plot the initial condition

# Load the physical grid 
ksi = read_restart('./grid/ksi_00000000')[-1][:,:,0]
eta = read_restart('./grid/eta_00000000')[-1][:,:,0]

print('ksi min/max, shape', ksi.min(), ksi.max(), ksi.shape)
print('eta min/max, shape', eta.min(), eta.max(), eta.shape)

lvls = np.linspace(-0.5,2.0,8)

# Get the list of restarts 
datlist = sorted(glob.glob('./restarts/restart_*'))


# Prepare the figure
fig = plt.figure(figsize=(7,4))
ax  = fig.add_subplot(111)

# Grad the first restart and the vorticity 
datname = datlist[0] 
nt,t,q = read_restart(datname)
datname = './out/omg/omg_' + str(0).zfill(8) 
omg = read_restart(datname)[-1][:,:,0]

# -- Add the grid
nx = np.shape(ksi)[0]
ny = np.shape(ksi)[1]
for i in range(nx):
    ax.plot( ksi[i,:], eta[i,:], ls='-',color='grey', lw=0.5  )
for j in range(ny):
    ax.plot( ksi[:,j], eta[:,j], ls='-',color='grey', lw=0.5  )

# -- Grab rho
qp = q[:,:,0]
im  = ax.pcolormesh(ksi,eta,qp, cmap='Blues_r')
cnt = ax.contour(ksi,eta,omg,levels=lvls,colors='k',zorder=99)

# -- Colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
fig.colorbar(im,cax=cax,label=r'$\rho/\rho_0$')

ax.set_aspect('equal')
fig.subplots_adjust(left=0.15, right=0.90, top=0.97)

#  ----- Label
ax.set_xlabel(r'$\xi/r_v$')
ax.set_ylabel(r'$\eta/r_v$')

#  ----- Plot colormesh
figname = '2d_wavymesh_mesh.png' 
plt.savefig(figname,format='png',dpi=200)

