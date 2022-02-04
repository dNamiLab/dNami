# -----------------------------------------------------------------------------
#
# CASE: 3D Navier-Stokes equations - Taylor Green Vortex case 
#       -- with kernels: rhs.py genRhs.py
#
# -----------------------------------------------------------------------------

# ======================================================================== LOAD
import dnami as dn # dNami kernel

from dnami import np, sys # re-import non dnami modules
import os

# ===================================================================  FOLDERS

# -- restarts 
try :
    os.mkdir('./restarts/')
except FileExistsError:
    pass

# ========================================================================= ASK

# Parameters for the case ...
alpha = dn.cst(2.5)
Ma    = dn.cst(0.45)
Re    = dn.cst(1600.0)
Pr    = dn.cst(0.71)
gamma = dn.cst(1.4)
dimPcr = False

if dimPcr:
	Cv = dn.cst(1.0/(gamma-1.0))
else:
	Cv = dn.cst(1.0/(Ma**2*gamma*(gamma-1.0)))

filtr_amp = dn.cst(0.1)    # filter amplitude

# ... in time ...
with_dt = dn.cst(1e-3)
#nitmax  = 50000 # for actual run
nitmax  = 1000  # for test case 

# ... in space ...
L = dn.cst(2.*np.pi) 
with_length = [L,L,L]      # domain length in each direction
with_grid   = [64,64,64]   # number of points in each direction

# ... as fast as possible!
with_proc     = [2,2,1] # mpi proc. topology

# ===================================================================== PREPARE

dtree = dn.create_tree()

# .. assign user-defined values
dtree['eqns']['coeff'][0][1] = dn.cst(1.0/Re)
dtree['eqns']['coeff'][1][1] = dn.cst(1.0/( (gamma-1.0)*Ma**2*Re*Pr ))
dtree['eqns']['coeff'][2][1] = dn.cst(gamma-1.)
dtree['eqns']['coeff'][3][1] = dn.cst(1.0/Cv)

# .. shortcut key handles
numerics = dtree['num']
grid     = dtree['grid']['size']
geom     = dtree['grid']['geom']
mpi      = dtree['mpi']['split']

grid['nxgb'] = with_grid[0] 
grid['nygb'] = with_grid[1]
grid['nzgb'] = with_grid[2]

geom['Lx'] = with_length[0] 
geom['Ly'] = with_length[1] 
geom['Lz'] = with_length[2] 

mpi['nxpr'] = with_proc[0] 
mpi['nypr'] = with_proc[1] 
mpi['nzpr'] = with_proc[2] 

# .. start the message passing interface
dtree = dn.start_mpi(dtree)
dMpi = dtree['mpi']['dMpi']

nx  = dMpi.nx
ny  = dMpi.ny
nz  = dMpi.nz

hlo  = numerics['hlo']

# .. create the computational grid and write to file
dtree = dn.create_grid(dtree)
dn.dnami_io.hello_world(dtree)

# define useful aliases
xloc, yloc, zloc  = geom['xloc'], geom['yloc'], geom['zloc']
Lx  , Ly  , Lz    = geom['Lx']  , geom['Ly']  , geom['Lz']
dx  , dy  , dz    = geom['dx']  , geom['dy']  , geom['dz']


numerics['tint']['tstep'] = with_dt 
dt = numerics['tint']['tstep']
numerics['filtr']['eps'] = filtr_amp 

# .. allocate tree
large = 10000 #no cache blocking in this example
dtree['libs']['cache blocking'] = [large,large,large]
dtree = dn.allocate(dtree)

# - Primitive variables
rh = dtree['eqns']['qvec']['views']['rho']
ux = dtree['eqns']['qvec']['views']['u']
uy = dtree['eqns']['qvec']['views']['v']
uz = dtree['eqns']['qvec']['views']['w']
et = dtree['eqns']['qvec']['views']['et']

q  = dtree['eqns']['qvec']['views']['q'] 

# - Store variables aliases if any
if 'qstored' in dtree['eqns']['qvec']['views'].keys():
	qstored  = dtree['eqns']['qvec']['views']['qstored'] 

# ================================================================== FUNCTIONS 

def sound_speed():
    e = et - .5*(ux*ux+uy*uy)
    T = (1./alpha)*( e*Ma*Ma )
    c = np.sqrt( T*(1.+1./alpha) )/Ma
    return c	

# ================================================================== INITIALISE

# initial clock
ti = dn.cst(0.0)
ni = 1
trstart = dn.cst(0.)

#init thermo
T0   = dn.cst(1.0)
P0   = dn.cst(1.0)/(Ma**2*gamma)
Rho0 = dn.cst(1.0)#P0/T0*Ma**2*gamma

#numpy slice refering to the core of the domain
dom = np.s_[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo]

rh[dom] = Rho0
ux[dom] = (np.sin(xloc[:, np.newaxis, np.newaxis])
			*np.cos(yloc[np.newaxis, :, np.newaxis])
			*np.cos(zloc[np.newaxis, np.newaxis, :]))

uy[dom] = (-np.cos(xloc[:, np.newaxis, np.newaxis])
			*np.sin(yloc[np.newaxis, :, np.newaxis])
			*np.cos(zloc[np.newaxis, np.newaxis, :]))

uz[dom] = dn.cst(0.0)

p =(P0 + dn.cst(Rho0/16.0)*(np.cos( dn.cst(2.0) * (xloc[:, np.newaxis, np.newaxis]) )
	 +np.cos( dn.cst(2.0) * (yloc[np.newaxis, :, np.newaxis]) ))
	*(np.cos( dn.cst(2.0) * (zloc[np.newaxis, np.newaxis, :]) ) + dn.cst(2.0)))

et[dom] = (  p/rh[dom]*dn.cst(1./(gamma-1.)) 
						 + dn.cst(0.5)*(ux[dom]*ux[dom]
						 			   +uy[dom]*uy[dom]
						 			   +uz[dom]*uz[dom]))

# -- Swap 
dMpi.swap(q,hlo,dtree)

if 'qstored' in dtree['eqns']['qvec']['views'].keys():
	dn.dnamiF.stored(intparam,fltparam,data)	
	dMpi.swap(qstored,hlo,dtree)

# -- Write the first restart
dn.dnami_io.write_restart(0,ti,0,dtree)

# ========================================================================= RUN

intparam,fltparam,data = (dtree['libs']['fort']['integers'],
						  dtree['libs']['fort']['floats'],
						  dtree['libs']['fort']['data'])

mod_filter = 1
mod_output = 1000
mod_info   = 100


for n in range(1,nitmax+1):
    ti = ti + dt

    # - RK loop
    for nrk in range(1,4):
        intparam[7] = nrk
        dMpi.swap(	q,hlo,dtree)

        if 'qstored' in dtree['eqns']['qvec']['views'].keys():
        	dn.dnamiF.stored(intparam,fltparam,data)
        	dMpi.swap(	qstored,hlo,dtree)

        dn.dnamiF.time_march(intparam,fltparam,data)	

    # - Filter
    if np.mod(n,mod_filter) == 0:

        dMpi.swapXc(q,hlo,dtree)
        dn.dnamiF.filter(1,intparam,fltparam,data)
        dMpi.swapYc(q,hlo,dtree)
        dn.dnamiF.filter(2,intparam,fltparam,data)
        dMpi.swapZc(q,hlo,dtree)
        dn.dnamiF.filter(3,intparam,fltparam,data)

    # - Output restarts 
    if np.mod(n,mod_output) == 0:
            dn.dnami_io.write_restart(n,ti,0,dtree)

    # - Output information
    if np.mod(n,mod_info) == 0:

        if dMpi.ioproc:
            print('____________________________________________________________')
            print('iteration',n,' with time t =',ti)
            sys.stdout.flush()
        dn.dnami_io.globalMinMax(dtree,rh,'r')
        dn.dnami_io.globalMinMax(dtree,ux,'u')
        dn.dnami_io.globalMinMax(dtree,uy,'v')
        dn.dnami_io.globalMinMax(dtree,uz,'w')
        dn.dnami_io.globalMinMax(dtree,et,'et')

# ----------------------------------------------------------------------------

# -- Grab the max value of rho-rho0 at end of run

if dMpi.iMpi:        
    maxval = np.amax(rh[:]-Rho0)
    MPI    = dMpi.MPIlib
    erra   = dMpi.comm_torus.reduce(maxval,op=MPI.MAX,root=0)
    if dMpi.ioproc:
        np.savetxt('out.dat',np.asarray([erra]))
else:
    erra = np.amax(rh[:]-Rho0)
    np.savetxt('out.dat',np.asarray([erra]))

