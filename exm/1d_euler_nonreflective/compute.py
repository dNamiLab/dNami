# -----------------------------------------------------------------------------
#
# CASE: 1D euler equations - entropy wave propagation + non-reflective BCs 
#       -- with kernels: rhs.py genRhs.py
#
# -----------------------------------------------------------------------------

# ======================================================================== LOAD
import dnami as dn         # dNami kernel

# Third-party libraries
from dnami import np, sys  # call non-dNami libs already available in dNami
import math as m           # ... and some others
import os

# ===================================================================  FOLDERS  

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
# ... for fluid ...
delta    = dn.cst(0.4) # R/Cv

# ... in space ...
L = dn.cst(1.) 
with_length = [L]    # domain length in each direction
with_grid   = [480]  # number of points in each direction

# ... and time ...
with_dt   = dn.cst(5.00e-4) # time step
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
eqns['coeff'][0][1] = delta 

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
large = 10000 #no cache blocking in this example
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
u   = dtree['eqns']['qvec']['views']['u']
et  = dtree['eqns']['qvec']['views']['et']

q  = dtree['eqns']['qvec']['views']['q'] 

# =================================================================== FUNCTIONS

gamma = dn.cst(1. + delta)

# - Internal energy EoS
def eos_e(rho,p):
        e = p/(rho*delta)
        return e

# - Speed of sound 
def eos_sos(rho,p):
        c = np.sqrt(gamma*p/rho)
        return c

# - Pressure EoS 
def eos_p(rho,e):
        p = delta*rho*e
        return p 

# ================================================================== INITIALISE

# initial clock
ti = dn.cst(0.0)
ni = 1
trstart = dn.cst(0.)
nitmax = 4000

# -- Initial thermo. conditions
rho0  = dn.cst(1.)
p0    = dn.cst(1.)
M0    = dn.cst(0.25)
c0    = eos_sos(rho0, p0)
u0    = M0*c0

# -- Entropy wave parameters
Lp    = dn.cst(0.1)  # wavelength
amp   = dn.cst(1e-3) # amplitude 
dom   = np.s_[hlo:nx+hlo] # interior domain

if dMpi.ioproc: 
    print(' --------------------------------------------------------')
    print(' ------------         1D wave propagation ---------------')
    print(' --------------------------------------------------------')

# -- Fill density and velocity fields 
rho[:] = rho0
u  [:] = u0 

# -- Add half sin-wave perturbation to the density field 
rho[dom] += amp * ( np.cos( np.pi*(xloc[:]-dn.cst(0.5)*Lx)/Lp ) ) * ( np.abs(xloc[:] - dn.cst(0.5)*Lx) <= dn.cst(0.5)*Lp   ) 

# -- Update total energy
et [:] = eos_e(rho[:],p0) + dn.cst(0.5)*u0*u0 

# ========================================================================= RUN

intparam,fltparam,data = (dtree['libs']['fort']['integers'],
                                                  dtree['libs']['fort']['floats'],
                                                  dtree['libs']['fort']['data'])
# -- Set restart output path 
rpath = './restarts/'

# -- Write out initial restart 
dn.dnami_io.write_restart(0,dn.cst(0.),0,dtree,fpath=rpath)

mod_filter = 1   # Filter frequency
mod_output = 50  # Output information frequency
mod_rstart = 50  # Restart write frequency

from timeit import default_timer as timer
t0 = timer()

for n in range(ni,nitmax+ni):
    ti = ti + dt

    for nrk in range(1,4):
        intparam[7] = nrk
        dMpi.swap(q,hlo,dtree) 
        dn.dnamiF.time_march(intparam,fltparam,data)    

    if np.mod(n,mod_filter) == 0:
        dMpi.swapX(q,hlo,dtree) 
        dn.dnamiF.filter(1,intparam,fltparam,data)

    if np.mod(n,mod_rstart) == 0:
        dn.dnami_io.write_restart(n,ti,0,dtree,fpath=rpath)

    if np.mod(n,mod_output) == 0:

        if dMpi.ioproc:
                print('____________________________________________________________')
                print('iteration',n,' with time t =',ti)
        e = et - .5*(u*u)
        p = eos_p(rho,e)
        c = eos_sos(rho[hlo:nx+hlo],p[hlo:nx+hlo])
        dn.dnami_io.globalMinMax(dtree,rho[hlo:nx+hlo],'r')
        dn.dnami_io.globalMinMax(dtree,u[hlo:nx+hlo],'u')
        dn.dnami_io.globalMinMax(dtree,et[hlo:nx+hlo],'et')
        dn.dnami_io.globalMinMax(dtree,np.abs( u[hlo:nx+hlo])/c,'M')
        if dMpi.ioproc:
                print('convective CFL numbers')
                sys.stdout.flush()
        cfl = dt*np.abs(u[hlo:nx+hlo])/dx
        dn.dnami_io.globalMax(dtree,cfl,'cfl-x')
        if dMpi.ioproc:
                print('acoustic CFL numbers')
                sys.stdout.flush()
        cfl = dt*(np.abs(u[hlo:nx+hlo])+c)/dx
        dn.dnami_io.globalMax(dtree,cfl,'cfl-x')

# ----------------------------------------------------------------------------

# -- Grab the max value of rho-rho0 at end of run

if dMpi.iMpi:        
    maxval = np.amax(rho[:]-rho0)
    MPI    = dMpi.MPIlib
    erra   = dMpi.comm_torus.reduce(maxval,op=MPI.MAX,root=0)
    if dMpi.ioproc:
        np.savetxt('out.dat',np.asarray([erra]))
else:
    erra = np.amax(rho[:]-rho0)
    np.savetxt('out.dat',np.asarray([erra]))

