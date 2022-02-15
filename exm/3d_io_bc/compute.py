# -----------------------------------------------------------------------------
#
# CASE: 3D I/O test 
#       -- with kernels: rhs.py genRhs.py
#
# -----------------------------------------------------------------------------

# ======================================================================== LOAD
import dnami as dn # dNami kernel

from dnami import np, sys # re-import non dnami modules
import os
import glob
import shutil

# ========================================================================= ASK

# ... in space ...
with_length = [1.,1.,1.]      # domain length in each direction
with_grid   = [64,64,64]   # number of points in each direction

# ... proc topology 
with_proc     = [2,1,2] # mpi proc. topology

# ===================================================================== PREPARE

dtree = dn.create_tree()

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

# .. allocate tree
large = 10000 #no cache blocking in this example
dtree['libs']['cache blocking'] = [large,large,large]
dtree = dn.allocate(dtree)

# - Primitive variables
a = dtree['eqns']['qvec']['views']['a']
b = dtree['eqns']['qvec']['views']['b']

q  = dtree['eqns']['qvec']['views']['q'] 

# ================================================================== INITIALISE

dom = np.s_[:,:,:]

a[dom] = np.random.rand(nx+2*hlo,ny+2*hlo,nz+2*hlo) 
b[dom] = np.random.rand(nx+2*hlo,ny+2*hlo,nz+2*hlo)

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
a[:,:,:] = dn.cst(0.)
b[:,:,:] = dn.cst(0.)

# -- Swap 
dMpi.swap(q,hlo,dtree)

# -- Copy files for reading in purposes
if dMpi.ioproc:
    for shell in glob.glob('restartshell*'):
        rname = shell.split('_')
        rname = rname[0] + '_' + rname[-1]
        shutil.copyfile(shell, rname)

# -- Read it back in  
dn.dnami_io.read_restart(dtree,fname='restart_00000000')

# -- Swap 
dMpi.swap(q,hlo,dtree)

# -- Compute difference and print min/max
adiff = a[:] - a0[:]
bdiff = b[:] - b0[:]

if dMpi.ioproc:
    print('Diff in core:')
dn.dnami_io.globalMinMax(dtree,adiff[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo],'adiff')
dn.dnami_io.globalMinMax(dtree,bdiff[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo],'bdiff')
if dMpi.ioproc:
    print('Diff over core + shells:')
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
