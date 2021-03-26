# Plot the 2D results and compare to analytical solution
import numpy as np
import sys
import matplotlib.pyplot as plt

wp = 'float64' # working precision
nvar = 4
iVDW = False
ndim  = 2
ivar = 0

# ksi/eta loader
def read_coords(fname):
        with open(fname,"rb") as fh:
            headsize = int(np.fromfile(fh,dtype=wp,count=1))
            head = np.fromfile(fh,dtype=wp,count=headsize-1)
            nx,ny,nz = int(head[2]),int(head[3]),int(head[4])
            print('sizes:', nx,ny,nz )
            f = np.fromfile(fh,dtype=wp,count=nx*ny*nz)
            if ndim == 3:
                    f = np.reshape(f,(nx,ny,nz))
            elif ndim == 2:
                    f = np.reshape(f,(nx,ny   ))
            else:
                    f = np.reshape(f,(nx      ))
        fh.closed
        print('loaded ', fname, 'shape:', np.shape(f) )
        return f

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

# Load all the stuff
ksi = read_coords('./grid/ksi_00000000')
eta = read_coords('./grid/eta_00000000')

vardic = ['rho','u','v','et']
lvls = np.linspace(-0.5,2.0,8)
nn  = np.arange(0, 600000, 1000)

iAnim = True 
iCut  = False

if iAnim :
    fig = plt.figure(figsize=(10,5))
    ax  = fig.add_subplot(111)
    datname = './restarts/restart_' + str(0).zfill(8) 
    nt,t,q = read_restart(datname)
    datname = './out/omg/omg_' + str(0).zfill(8) 
    omg = read_coords(datname)
    print(np.amin(omg), np.amax(omg))
    qp = q[:,:,ivar]
    im  = ax.pcolormesh(ksi,eta,qp, cmap='Oranges')
    cnt = ax.contour(ksi,eta,omg,levels=lvls,colors='k')
    #ax.axhline(y=0)
    #ax.axvline(x=0)
    # -- Add the grid
    nx = np.shape(ksi)[0]
    ny = np.shape(ksi)[1]
    for i in range(nx):
        ax.plot( ksi[i,:], eta[i,:], 'k--', lw=0.5  )
    for j in range(ny):
        ax.plot( ksi[:,j], eta[:,j], 'k--', lw=0.5  )
    fig.colorbar(im,ax=ax)

    #  ----- Plot colormesh
    for k,n in enumerate(nn):
        print('Making fig at n:', n)
        datname = './restarts/restart_' + str(n).zfill(8) 
        nt,t,q = read_restart(datname)
        qp = q[:,:,ivar]
        datname = './out/omg/omg_' + str(n).zfill(8) 
        omg = read_coords(datname)
        print(np.amin(omg), np.amax(omg), np.shape(omg))
        im  = ax.pcolormesh(ksi,eta,qp, cmap='Oranges')
        for coll in cnt.collections: 
                plt.gca().collections.remove(coll)
        cnt = ax.contour(ksi,eta,omg,levels=lvls,colors='k')
        ax.set(title='Var: ' + vardic[ivar] + ' \n min/max: {:.4} {:.4}'.format(np.amin(qp), np.amax(qp)))
        figname = 'pics/omg_' + str(k).zfill(4) + '.png'
        plt.savefig(figname,format='png',dpi=200)

# -- Comparison
if iCut:
    ivar = 2
    fig = plt.figure(figsize=(8,8))
    ax  = fig.add_subplot(111)
    datname = './restarts/restart_' + str(0).zfill(8) 
    nt,t,q = read_restart(datname)
    datname = './out/omg/omg_' + str(0).zfill(8) 
    omg0 = read_coords(datname)
    qp0 = q[:,:,ivar]

    datname = './restarts/restart_' + str(595000).zfill(8) 
    nt,t,q = read_restart(datname)
    qp1 = q[:,:,ivar]
    datname = './out/omg/omg_' + str(595000).zfill(8) 
    omg1 = read_coords(datname)

    nymid = int(np.size(qp0[0,:])/2.) 
    x = ksi[:,nymid] 
    uinf =  0.5916079783099616
    ax.plot(x, qp0[:,nymid]/uinf, ls='', marker='o', color='k', label=r'Initial condition')
    ax.plot(x, qp1[:,nymid]/uinf, ls='-', color='k', label=r'(9pt, 8th ord) - 11 travel times ')
    ax.set_xlabel(r'x/r_v')
    ax.set_ylabel(r'v/u_inf')
    ax.set_xlim([-5,5])
    ax.legend(loc='best')

    figname = 'pics/cut.png'
    plt.savefig(figname,format='png',dpi=200)

