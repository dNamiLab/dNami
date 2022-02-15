# -----------------------------------------------------------------------------
#
# CASE: 1D I/O test 
#       -- with kernels: rhs.py genRhs.py
#
# -----------------------------------------------------------------------------

# ======================================================================== LOAD
import dnami as dn # dNami kernel

from dnami import np, sys # re-import non dnami modules
import os

# ========================================================================= ASK

# ... in space ...
with_length = [1.]      # domain length in each direction
with_grid   = [64]   # number of points in each direction

# ... proc topology 
with_proc     = [2] # mpi proc. topology

# ===================================================================== PREPARE

dtree = dn.create_tree()

# .. shortcut key handles
numerics = dtree['num']
grid     = dtree['grid']['size']
geom     = dtree['grid']['geom']
mpi      = dtree['mpi']['split']

grid['nxgb'] = with_grid[0] 
grid['nygb'] = 1 
grid['nzgb'] = 1 

geom['Lx'] = with_length[0] 
geom['Ly'] = 0. 
geom['Lz'] = 0. 

mpi['nxpr'] = with_proc[0] 
mpi['nypr'] = 1 
mpi['nzpr'] = 1 

# .. start the message passing interface
dtree = dn.start_mpi(dtree)
dMpi = dtree['mpi']['dMpi']

nx  = dMpi.nx

hlo  = numerics['hlo']

# .. create the computational grid and write to file
dtree = dn.create_grid(dtree)
dn.dnami_io.hello_world(dtree)

# .. allocate tree
large = 10000 #no cache blocking in this example
dtree['libs']['cache blocking'] = [large,large,large]
dtree = dn.allocate(dtree)

# - Primitive variables
a = dtree['eqns']['qvec']['views']['a']
b = dtree['eqns']['qvec']['views']['b']

q  = dtree['eqns']['qvec']['views']['q'] 

# ================================================================== INITIALISE

dom = np.s_[hlo:nx+hlo]

a[dom] = np.random.rand(nx) 
b[dom] = np.random.rand(nx)

# -- Swap 
dMpi.swap(q,hlo,dtree)

a0 = a.copy()
b0 = b.copy()

# ================================================================== I/O TEST 

# -- Write a restart 

if dMpi.ioproc:
    print('Writing restart ...')

dn.dnami_io.write_restart(0,0.,0,dtree,fpath='./')

if dMpi.ioproc:
    print('Done.')

# -- Zero the arrays
a[:] = dn.cst(0.)
b[:] = dn.cst(0.)

# -- Read it back in  
dn.dnami_io.read_restart(dtree,fname='restart_00000000')

# -- Swap 
dMpi.swap(q,hlo,dtree)

# -- Compute difference and print min/max
adiff = a[:] - a0[:]
bdiff = b[:] - b0[:]

dn.dnami_io.globalMinMax(dtree,adiff,'adiff')
dn.dnami_io.globalMinMax(dtree,bdiff,'bdiff')

# -- Grab the max value of the difference 

if dMpi.iMpi:        
    maxval = np.amax(np.abs(adiff))
    MPI    = dMpi.MPIlib
    erra   = dMpi.comm_torus.reduce(maxval,op=MPI.MAX,root=0)
    if dMpi.ioproc:
        np.savetxt('out.dat',np.asarray([erra]))
else:
    erra = np.amax(np.abs(adiff))
    np.savetxt('out.dat',np.asarray([erra]))
