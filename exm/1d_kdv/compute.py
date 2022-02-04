# -----------------------------------------------------------------------------
#
# CASE: 1D kdv equation 
#       -- with kernels: rhs.py genRhs.py
#
# -----------------------------------------------------------------------------

# ======================================================================== LOAD
import dnami as dn         # dNami kernel

# Third-party libraries
from dnami import np, sys  # call non-dNami libs already available in dNami
import math as m           # ... and some others
import os

# dNami pre-processing tools?
# dNami runtime tools?

# =================================================================== Create WRK

# -- restarts
try:
    os.mkdir('./restarts/')
except FileExistsError:
    pass
# -- out
try :
    os.mkdir('./out/')
except FileExistsError:
    pass

# ========================================================================= ASK

# Solve the equation ...

# ... in space ...
L = dn.cst(40.0) 
with_length = [L]  # domain length in each direction
with_grid   = [500]   # number of points in each direction

# ... and time ...
with_dt   = dn.cst(1.e-4) # time step
nitmax    = 60000 
filtr_amp = dn.cst(0.1)    # filter amplitude

# ... as fast as possible!
with_proc = [2] # mpi proc. topology

# ===================================================================== PREPARE

# create & populate dnami data tree

dtree = dn.create_tree()

# .. shortcut key handles
eqns     = dtree['eqns']
geom     = dtree['grid']['geom']
grid     = dtree['grid']['size']
mpi      = dtree['mpi']['split']
numerics = dtree['num']

# .. assign user-defined values
eqns['coeff'][0][1] = dn.cst(6.) 
eqns['coeff'][1][1] = dn.cst(1.) 

geom['Lx'] = with_length[0] 
geom['Ly'] = 0. 
geom['Lz'] = 0. 

grid['nxgb'] = with_grid[0]
grid['nygb'] = 1 
grid['nzgb'] = 1

mpi['nxpr'] = with_proc[0]
mpi['nypr'] = 1 
mpi['nzpr'] = 1 

numerics['tint']['tstep'] = with_dt
numerics['filtr']['eps']  = filtr_amp

# .. start the message passing interface
dtree = dn.start_mpi(dtree)
dMpi = dtree['mpi']['dMpi']

# .. create the computational grid and write to file
dtree = dn.create_grid(dtree)
dn.dnami_io.write_grid(dtree)

# .. broadcast start-up info
dn.dnami_io.hello_world(dtree)

# .. allocate tree
large = 10000
dtree['libs']['cache blocking'] = [large,large,large]
dtree = dn.allocate(dtree)

# define useful aliases

xloc= geom['xloc']

Lx  = geom['Lx']  

dx  = geom['dx']  

nx  = dMpi.nx     

x   = geom['x']

dt  = numerics['tint']['tstep']
hlo = numerics['hlo']

rho = dtree['eqns']['qvec']['views']['rho']

q  = dtree['eqns']['qvec']['views']['q'] 

if 'qstored' in dtree['eqns']['qvec']['views'].keys():
        qstored  = dtree['eqns']['qvec']['views']['qstored'] 

# ================================================================== INITIALISE

# initial clock
ti = dn.cst(0.0)
ni = 1
trstart = dn.cst(0.)

# -- Initial thermo. conditions
rho0  = dn.cst(1.)

if dMpi.ioproc: 
    print(' --------------------------------------------------------')
    print(' ------------ 1D kdv test case            ---------------')
    print(' --------------------------------------------------------')

x0       = 20
xll      = np.linspace(-x0,x0,2*grid['nxgb']) #xloc[:] - 20.
k1       = 1.
k2       = np.sqrt(2.) 
eta1_0   = 0.
eta2_0   =  2.*np.sqrt(2)
eta1     = k1*xll[:] + eta1_0 
eta2     = k2*xll[:] + eta2_0 
f        = 1. + np.exp(eta1) + np.exp(eta2) + (k1-k2)**2/(k1+k2)**2 * np.exp(eta1+eta2)
logf     = np.log(f)
gradlogf = np.gradient(logf, xll[:])
sol      = 2.  * np.gradient(gradlogf, xll[:]) 
from scipy.interpolate import interp1d
sol = interp1d(xll, sol)

rho[hlo:nx+hlo] = sol(xloc[:]-x0) 

dMpi.swap(q,hlo,dtree)

# ========================================================================= RUN

intparam,fltparam,data = (dtree['libs']['fort']['integers'],
                                                  dtree['libs']['fort']['floats'],
                                                  dtree['libs']['fort']['data'])

# -- Write out initial restarts 
dn.dnami_io.write_restart(0,dn.cst(0.),0,dtree)

mod_filter = 1 
mod_output = 500 
mod_rstart = 1000 

for n in range(ni,nitmax+ni):
        ti = ti + dt

        for nrk in range(1,4):
            intparam[7] = nrk
            dMpi.swap(q,hlo,dtree) 
            if 'qstored' in dtree['eqns']['qvec']['views'].keys():
                dn.dnamiF.stored(intparam,fltparam,data,0)
                dMpi.swap(qstored,hlo,dtree)
            dn.dnamiF.time_march(intparam,fltparam,data)    

        if np.mod(n,mod_filter) == 0:
            dn.dnamiF.filter(1,intparam,fltparam,data)

        if np.mod(n,mod_rstart) == 0:
            dn.dnami_io.write_restart(n,ti,0,dtree)

        if np.mod(n,mod_output) == 0:

            if dMpi.ioproc:
                    print('____________________________________________________________')
                    print('iteration',n,' with time t =',ti)
            dn.dnami_io.globalMinMax(dtree,rho[hlo:nx+hlo],'r')
# ----------------------------------------------------------------------------


# -- Grab the max value of rho-rho0 at end of run

if dMpi.iMpi:        
    maxval = np.amax(rho[:])
    MPI    = dMpi.MPIlib
    erra   = dMpi.comm_torus.reduce(maxval,op=MPI.MAX,root=0)
    if dMpi.ioproc:
        np.savetxt('out.dat',np.asarray([erra]))
else:
    erra = np.amax(rho[:])
    np.savetxt('out.dat',np.asarray([erra]))

