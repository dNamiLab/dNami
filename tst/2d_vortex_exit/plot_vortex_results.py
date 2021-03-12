# Plot the 2D results and compare to analytical solution
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
import matplotlib as mpl
from matplotlib import rc
import matplotlib.colors as colors
import glob

fnt    =  14      
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

wp = 'float64' # working precision
nvar = 4
iVDW = True 
ndim  = 2
ivar = 0

# -- EoS

delta = 0.4 

def eos_p(rho,e):
    return delta*rho*e


# load grid
with open('out/axes.bin',"rb") as fh:
        head = np.fromfile(fh,dtype=wp,count=6)
        nxgb, nygb, nzgb = int(head[0]), int(head[1]), int(head[2])
        Lx, Ly, Lz = head[3], head[4], head[5]; del head
        ndim = 2
        if nygb == 1: ndim = 1
        if nzgb >  1: ndim = 3
        x = np.fromfile(fh,dtype=wp,count=nxgb); dx = x[1]-x[0]
        if ndim>1: y = np.fromfile(fh,dtype=wp,count=nygb); dy = y[1]-y[0]
        if ndim>2: z = np.fromfile(fh,dtype=wp,count=nzgb); dz = z[1]-z[0]
        fh.closed
nx = nxgb; ny = nygb; nz = nzgb

# restart loader
def read_restart(fname = 'restart.bin'):
        with open(fname,"rb") as fh:
            headsize = int(np.fromfile(fh,dtype=wp,count=1))
            head = np.fromfile(fh,dtype=wp,count=headsize-1)
            n = int(head[0])
            t = head[1]
            nx,ny,nz,nv = int(head[2]),int(head[3]),int(head[4]),int(head[5])
            f = np.fromfile(fh,dtype=wp,count=nx*ny*nz*nv)
            if ndim == 3:
                    f = np.reshape(f,(nx,ny,nz,nvar))
            elif ndim == 2:
                    f = np.reshape(f,(nx,ny   ,nvar))
            else:
                    f = np.reshape(f,(nx      ,nvar))
        fh.closed
        return n,t,f

# Load all the stuff/
path = './'

#nn  = np.arange(0, 1000000, 100)
nn  = [0] 

# -- Debug and show plot
iShow = False 
iHide = False 
iCompStart = True

# --- PERPARE THE ACTUAL FIGURE

# -- Grab the first field
if iCompStart:
    datname = path + 'restarts/restart_' + str(0).zfill(8) 
    print('Reading the starting flow from ....', datname)
    nt,t,q = read_restart(datname)
    rs = q[:,:,0]
    urs = q[:,:,1]
    uts = q[:,:,2]
    ets = q[:,:,3]

    # -- Create base pressure
    es = ets - 0.5*(urs**2 + uts**2)
    ps = eos_p(rs,es)
    fieldsS = [rs, urs, ps]

# -- Grab n zero field
datname = path + 'restarts/restart_' + str(nn[0]).zfill(8) 
datlist = sorted(glob.glob(path+'restarts/restart_*'))
ist = -1 
horiz = False

# --------------------------------------------------
print('Creating the plot ...')
ll = 0.1 
lr = 0.1
offs = 0.95 
lw = 2
fig = plt.figure(figsize = (6,6) )
ax  = fig.add_axes([ll,lr,offs-ll,offs-lr])

for i,datname in enumerate(datlist) :
    if np.mod(i,1) == 0 and i > ist:
        print('Reading the baseflow from ....', datname)
        nt,t,q = read_restart(datname)
        r0 = q[:,:,0]
        u0 = q[:,:,1]
        v0 = q[:,:,2]
        et0 = q[:,:,3]
        
        rref = 1.
        delr = (r0 - rref)/rref * 1e2

        # -- Min/max
        vminv = [-2,-1]
        vmaxv = [2,1]
        fdname = ['rho','omg']
        cbname = [r'$\Delta \rho / \rho_0 \times 10^{2}$',r'$\Omega$']

        # -- 
        nxx = np.size(u0[:,0])
        nyy = np.size(u0[0,:])
        xx = np.asarray(range(nxx))
        yy = np.asarray(range(nyy))

        if i ==0 and horiz:
            u0min = np.amin(u0)
            u0max = np.amax(u0)
            lvls  = np.linspace(u0min,u0max,8)

        if i ==0 and not horiz:
            v0min = np.amin(v0)
            v0max = np.amax(v0)
            lvls  = np.linspace(v0min,v0max,8)

        k = 0
        field = delr

        if i == 0:
            # -- Grab min/max
            vmin= vminv[k]
            vmax= vmaxv[k]

            print('min/max del rho', delr.min(), delr.max())
            # --- CREATE THE 2D PLOT
            cmap = 'Greys'

            # --- ACTUALLY PUT TOGETHER THE PLOT
            im  = ax.pcolormesh(field.T,cmap=cmap,vmin=vmin,vmax=vmax)

            # -- Add quiver
            if horiz:
                C = ax.contour(u0.T,levels=lvls)
            else:
                C = ax.contour(v0.T,levels=lvls)

            # Limits
            ax.set(xlim=[xx[0], xx[-1]])
            ax.set(ylim=[yy[0], yy[-1]])

            # -- Colorbar plus axes
            ax2  = fig.add_axes([0.60,0.20,0.3,0.03])
            fig.colorbar(im,cax=ax2,cmap=cmap,orientation='horizontal',extend='both',ticks=[vmin,0.5*(vmin+vmax),vmax])
            ax2.text(0.30,3.0, cbname[k],horizontalalignment='center')

            # --  Title
            str_ttl = r'Min/Max $\times 10 ^2$ : {:4f} {:4f}'.format(delr.min(), delr.max())
            title_text = ax.set_title(str_ttl)

            figname = 'pics/' + fdname[k] + str(i).zfill(4) + '.png'.format(i)
            print('SAVING FIG ...', figname)

            plt.savefig(figname,format='png',dpi=200)
            print('DONE')
        else:
            # -- Update values
            im.set_array( delr[:,:].T.ravel() )
            for coll in C.collections: 
                    ax.collections.remove(coll) 
            if horiz:
                C = ax.contour(u0.T,levels=lvls)
            else:
                C = ax.contour(v0.T,levels=lvls)

            str_ttl = r'Min/Max $\times 10 ^2$ : {:4f} {:4f}'.format(delr.min(), delr.max())
            title_text.set_text(str_ttl)

            figname = 'pics/' + fdname[k] + str(i).zfill(4) + '.png'.format(i)
            print('SAVING FIG ...', figname)

            plt.savefig(figname,format='png',dpi=200)
            print('DONE')

plt.close()
