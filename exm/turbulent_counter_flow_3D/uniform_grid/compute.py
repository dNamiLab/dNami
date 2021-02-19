# -----------------------------------------------------------------------------
#
# CASE: 3D Counter-Flow case with a uniform grid 
#       
#       -- with kernels: rhs_???.py genKer_???.py
#                                                          -djl 15/02/2021
# -----------------------------------------------------------------------------

# ======================================================================== LOAD
import dnami as dn         # dNami kernel
# Third-party libraries
from dnami import np, sys  # call non-dNami libs already available in dNami
import math as m           # ... and some others
import os
import time
# =================================================================== Create WRK (if not already done)

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
try:
    os.mkdir('./out/derived/')
except FileExistsError:
    pass
# -- out/liv
try :
    os.mkdir('./out/liv')
except FileExistsError:
    pass
# -- out
try :
    os.mkdir('./pics/')
except FileExistsError:
    pass

# ========================================================================= ASK

# Solve the equation ...
iVDW = False
# ... in space ...
Nx, Ny, Nz = 128, 231, 64
with_length = [12.0, 2.0, 6.0]            #Domain length in each direction
with_grid   = [Nx,Ny,Nz]        #Number of points in each direction
# ... and time ...i
with_cfl  = dn.cst(0.1)       #Fixed CFL number (-> this will dynamically determine dt)
nitmax    = 1000000       #Number of time steps 
filtr_amp = dn.cst(0.1)        #Filter amplitude

# ... as fast as possible!
MPI_x, MPI_y, MPI_z = 4,7,4
#MPI_x, MPI_y, MPI_z = 16,31,8
with_proc = [MPI_x, MPI_y, MPI_z]              #Mpi proc. topology

# ... with parameters :
Re        = dn.cst(190.7)
M0        = dn.cst(0.5)       #Reference Mach number
gamma     = dn.cst(1.4) # Gamma Cp/Cv
delta     = gamma - 1.0
Pr        = dn.cst(0.70) # Prandtl number
Twall     = dn.cst(1.0)
chi       = dn.cst(10.0) # Bulk contribution
mu        = dn.cst(1.0) # Shear contribution
# ===================================================================== PREPARE

if chi == 0.0:
    chi_cfl = 1.0
else:
    chi_cfl = chi

restart = False
#restart = True
# create & populate dnami data tree
dtree = dn.create_tree()

# .. shortcut key handles
eqns, geom, grid, mpi, numerics = dtree['eqns'], dtree['grid']['geom'], dtree['grid']['size'], dtree['mpi']['split'], dtree['num']

# .. assign user-defined values
eqns['coeff'][0][1] = 1./Re # ReInv
eqns['coeff'][1][1] = M0**2 # M0_sq
eqns['coeff'][2][1] = delta # delta
eqns['coeff'][3][1] = dn.cst(1.0/( (gamma-1.0)*M0**2*Re*Pr )) # kappa
eqns['coeff'][4][1] = mu # mu
eqns['coeff'][5][1] = chi # mub
eqns['coeff'][6][1] = gamma # gamma
eqns['coeff'][7][1] = Twall # Tw

# Domain lengths
geom['Lx'], geom['Ly'], geom['Lz'] = with_length[0], with_length[1], with_length[2]
# Number of grid points per direction
grid['nxgb'], grid['nygb'], grid['nzgb'] = with_grid[0], with_grid[1], with_grid[2]
# Number of MPI processes per direction
mpi['nxpr'], mpi['nypr'], mpi['nzpr'] = with_proc[0], with_proc[1], with_proc[2]
# Filter strength
numerics['filtr']['eps']  = filtr_amp

# .. start the message passing interface
dtree = dn.start_mpi(dtree)
dMpi = dtree['mpi']['dMpi']

# .. create the computational grid and write to file
y1, y2 = -1.0, 1.0
dtree = dn.create_grid(dtree, y1=y1, y2=y2)
dn.dnami_io.write_grid(dtree)

# .. broadcast start-up info
dn.dnami_io.hello_world(dtree)

# .. allocate tree
large = 100000
dtree['libs']['cache blocking'] = [large,large,large]
dtree = dn.allocate(dtree)

# define useful aliases
xloc, yloc, zloc = geom['xloc'], geom['yloc'], geom['zloc']
Lx, Ly, Lz = geom['Lx'], geom['Ly'], geom['Lz']
dx, dy, dz = geom['dx'], geom['dy'], geom['dz']
x, y, z = geom['x'], geom['y'], geom['z']
nx, ny, nz = dMpi.nx, dMpi.ny, dMpi.nz 
hlo = numerics['hlo']

# -- Aliases for the fields (density, velocity, etc)
eqns = dtree['eqns']['qvec']['views']
q  = eqns['q']
rho, u, v, w, et = eqns['rho'], eqns['u'], eqns['v'], eqns['w'], eqns['et']

# -- Stored values (here the divergence field)
if 'qstored' in dtree['eqns']['qvec']['views'].keys():
        qstored  = dtree['eqns']['qvec']['views']['qstored'] 

divV = qstored[:,:,:,0]
c = qstored[:,:,:,1]
M = qstored[:,:,:,2]
phi = qstored[:,:,:,3]

# Statistics
variables = ['rho', 'u', 'v', 'w', 'M', 'T', 'et']
rho_mean = qstored[:,:,:,4]
u_mean = qstored[:,:,:,5]
v_mean = qstored[:,:,:,6]
w_mean = qstored[:,:,:,7]
M_mean = qstored[:,:,:,8]
T_mean = qstored[:,:,:,9]
et_mean = qstored[:,:,:,10]

# Set the body forcing term
for i in range(nx+2*hlo):
    for k in range(nz+2*hlo):
        phi[i,:,k] = np.tanh(100.0*yloc)
# =================================================================== FUNCTIONS

# -- Equations of State
if iVDW:
    gas_type = 1
    Zc = np.float64(3./8.)
    gas_type_str = 'vdw'
    # delta = 1./alpha

    def eos_e(rho,p):
        n_o_e = np.float64(9.)/np.float64(8.)
        thr = np.float64(3.)
        egh = np.float64(8.)
        T   = (p+thr*rho**2)/egh/rho*(thr-rho)
        e         = T / delta - n_o_e*rho
        return e

    def eos_sos(rho,p):
        one       = np.float64(1.)
        two       = np.float64(2.)
        thr       = np.float64(3.)
        n_o_e     = np.float64(9.)/np.float64(8.)
        thr       = np.float64(3.)
        egh       = np.float64(8.)
        T         = (p+thr*rho**2)/egh/rho*(thr-rho)
        c_squared = T*(one + delta)*( one / ( one - rho / thr) )**2 - two*rho*n_o_e
        c         = np.sqrt( c_squared )
        return c

    def eos_p(rho,e):
        thr = np.float64(3.)
        egh = np.float64(8.)
        thr    = np.float64(3.)
        egh    = np.float64(8.)
        n_o_e  = np.float64(9.)/np.float64(8.)
        T = delta*(e+n_o_e*rho)
        p   = egh*rho*T/(thr-rho) - thr*rho**2
        return p 

else:
    gas_type = 0
    Zc = np.float64(1.)
    gas_type_str = 'ig'
    alpha = 1.0/delta
    gamma = 1. + 1./alpha

    def eos_e(rho,p): # Internal energy
        e = p/(rho*(gamma-1.))
        return e

    def eos_sos(rho,p): # Speed of sound
        c =  np.sqrt(gamma*p/rho)
        return c

    def eos_p(rho,e): # Pressure
            p = (gamma-1.)*rho*e
            return p

def set_dt(target_cfl):
    dtMax = None

    dts = target_cfl*dx/np.abs(u)
    dtMax1 = dn.dnami_io.dtMax(dtree,dts)
    dts = target_cfl*dy/np.abs(v)
    dtMax2 = dn.dnami_io.dtMax(dtree,dts)
    dts = target_cfl*dz/np.abs(w)
    dtMax3 = dn.dnami_io.dtMax(dtree,dts)
    if dMpi.ioproc:
        dtMaxConv = np.amin([dtMax1,dtMax2,dtMax3])
        sys.stdout.flush()

    dts = target_cfl*dx/(np.abs(u)+c)
    dtMax1 = dn.dnami_io.dtMax(dtree,dts)
    dts = target_cfl*dy/(np.abs(v)+c)
    dtMax2 = dn.dnami_io.dtMax(dtree,dts)
    dts = target_cfl*dz/(np.abs(w)+c)
    dtMax3 = dn.dnami_io.dtMax(dtree,dts)
    if dMpi.ioproc:
        dtMaxAcou = np.amin([dtMax1,dtMax2,dtMax3])
        sys.stdout.flush()

        dtMaxVisc = np.amin([target_cfl*Re*dx*dx/chi, target_cfl*Re*dy*dy/chi, target_cfl*Re*dz*dz/chi])
        sys.stdout.flush()
        dtMax = np.amin([dtMaxConv,dtMaxAcou,dtMaxVisc])
    
    dn.dnami_io.set_dt(dtree,dtMax) 
    dt = numerics['tint']['tstep']
    return dt

# Initial info
if dMpi.ioproc: 
    print(' ------------ 3D Counter-Flow simulation ---------------')
    print(' Restart: {:}'.format(restart))
    print(' Grid dimensions (Lx, Ly, Lz) =  ({:.4}, {:.4}, {:.4})'.format(Lx, Ly, Lz))
    print(' Grid distribution (Nx, Ny, Nz): ({:}, {:})'.format(Nx, Ny, Nz))
    print(' MPI distribution (MPI_x, MPI_y, MPI_z): ({:}, {:}, {:})'.format(MPI_x, MPI_y, MPI_z))
    print(' VdW Gas: {:}'.format(iVDW))
    print(' Twall: {:.3f}, Mach: {:.3f}, Pr: {:.3f}'.format(Twall, M0, Pr))
    print(' ------------------------------------------------')

def initial_channel(x, y, z, rho, u, v, w, et):
    vonkar, b = 2.5, 5.0
    mask = Re*(1-np.abs(y)) < 10.0
    ubar = np.zeros_like(rho)
    ubar[mask] = Re*(1 - np.abs(y[mask]))
    ubar[~mask] = b + vonkar*np.log10(Re*(1.0 - np.abs(y[~mask])))
    A = 0.05*(vonkar*np.log10(Re) + b)

    # Trigonometric functions
    Lx, Ly, Lz = with_length[0], with_length[1], with_length[2]
    sx, sy, sz = np.sin(4.*np.pi*x / Lx), np.sin(np.pi*y), np.sin(2.*np.pi*z / Lz)
    cx, cy, cz = np.cos(4.*np.pi*x / Lx), 1.0 + np.cos(np.pi*y), np.cos(2.*np.pi*z / Lz)

    u[:] = 0.05*ubar + 0.5*A*Lx*cx*sy*cz
    v[:] = -A*sx*cy*sz
    w[:] = -A*(Lz/2.)*sx*sy*cz
    # Assume constant temperature and density
    p = 1.0 / (gamma*M0**2)
    rho[:] = 1.0
    et[:] = p/(rho*(gamma-1.0)) + 0.5*(u**2 + v**2 + w**2)
    return rho, u, v, w, et
# ================================================================== INITIALISE

# Read in the initial condition
if restart:
    dn.dnami_io.read_restart(dtree)
    ni = dtree['num']['tint']['itn']
    ti = dtree['eqns']['time']
    trstart = dtree['eqns']['time']
else:
    y_mesh = np.repeat(yloc[np.newaxis,:], nx+2*hlo, axis=0)
    y_mesh = np.repeat(y_mesh[:,:,np.newaxis], nz+2*hlo, axis=2)
    x_mesh = np.zeros((nx+2*hlo, ny+2*hlo, nz+2*hlo))
    z_mesh = np.zeros((nx+2*hlo, ny+2*hlo, nz+2*hlo))
    for j in range(ny+2*hlo):
        for k in range(nz+2*hlo):
            x_mesh[hlo:nx+hlo,j,k] = xloc[:]

    for i in range(nx+2*hlo):
        for j in range(ny+2*hlo):
            z_mesh[i,j,hlo:nz+hlo] = zloc[:]
    rho, u, v, w, et = initial_channel(x_mesh, y_mesh, z_mesh, rho, u, v, w, et)
    # Set the initial time
    ti = dn.cst(0.0)
    ni = 1
    trstart = dn.cst(0.)

# -- Swap fields in mpi
dMpi.swap(q,hlo,dtree)

if dMpi.ioproc:
    print(' -------------------------- ' )
    print('Initial values:')
dn.dnami_io.globalMinMax(dtree,rho,'rho')
dn.dnami_io.globalMinMax(dtree,u,'u')
dn.dnami_io.globalMinMax(dtree,v,'v')
dn.dnami_io.globalMinMax(dtree,w,'w')
dn.dnami_io.globalMinMax(dtree,et,'et')
dn.dnami_io.globalMinMax(dtree,phi,'phi')
dn.dnami_io.globalMinMax(dtree,yloc,'yloc')

# ========================================================================= RUN

intparam,fltparam,data = (dtree['libs']['fort']['integers'],
                                                  dtree['libs']['fort']['floats'],
                                                  dtree['libs']['fort']['data'])

# -- Compute the stored values 
if 'qstored' in dtree['eqns']['qvec']['views'].keys():
        dMpi.swap(q,hlo,dtree)
        dMpi.swap(qstored,hlo,dtree)
        dn.dnamiF.stored(intparam,fltparam,data,1)        

if not restart:
    # -- Write first restart (and first initial state in live_view)
    dn.dnami_io.write_restart(0,ti,0,dtree)
    dn.dnami_io.write_restart(0,ti,1,dtree)

# -- Pause to see initial field
#if dMpi.ioproc:
#   input("Press Enter to continue...")

# -- Set how often things should happen (filtering, printing, etc)

# - Set in multiples of time step
mod_filter = 10
mod_info   = 100
mod_cfl    = 1
mod_NaN    = 1000

# - Set in physical times 
tw_out    = 0.1
tw_rst    = 0.1
to         = ti + tw_out
tr         = ti + tw_rst

dt = 5.0e-5
# -- Initialise dt
#fixed_timestep = False
fixed_timestep = True
if fixed_timestep:
    fltparam[3] = dt
else:
    dt = set_dt(with_cfl); fltparam[3] = dt

count_stats = 0
start_time = time.time()
for n in range(ni,nitmax+ni):
        ti = ti + dt

        #RK loop
        for nrk in range(1,4):
            intparam[7] = nrk
            dMpi.swap(q,hlo,dtree) 
            if 'qstored' in dtree['eqns']['qvec']['views'].keys():
                    dMpi.swap(qstored,hlo,dtree)
                    dn.dnamiF.stored(intparam,fltparam,data)
            dn.dnamiF.time_march(intparam,fltparam,data)

        #Filtering
        if np.mod(n,mod_filter) == 0:
            dMpi.swapX(q,hlo,dtree)
            dn.dnamiF.filter(1,intparam,fltparam,data)
            dMpi.swapY(q,hlo,dtree)
            dn.dnamiF.filter(2,intparam,fltparam,data)
            dMpi.swapZ(q,hlo,dtree)
            dn.dnamiF.filter(3,intparam,fltparam,data)

        #Output to restart
        if np.abs(ti-tr) < 0.5*dt:  
            tr = tr + tw_out  
            dn.dnami_io.write_restart(n,ti,0,dtree)
            
        #Output to liveview
        if np.abs(ti-to) < 0.5*dt:
            to = to + tw_out  
            dn.dnami_io.write_restart(n,ti,1,dtree)
            dn.dnami_io.write_data(['divV','M'],n,ti,dtree,'./out/derived/','output')
        #Adjusting CFL
        if fixed_timestep:
            pass
        else:
            if np.mod(n,mod_cfl) == 0:
                dt = set_dt(with_cfl); fltparam[3] = dt

        #Print info CFL, time, etc
        if np.mod(n,mod_info) == 0:
            if dMpi.ioproc:
                    print('____________________________________________________________')
                    print('iteration',n,' with time t =',ti ,' with time-step dt =',dt)
            dn.dnami_io.globalMinMax(dtree,u,'u')
            dn.dnami_io.globalMax(dtree,(et - 0.5*(u**2 + v**2 + w**2)) * M0**2 * gamma*(gamma - 1.0),'T')
            dn.dnami_io.globalMax(dtree,M,'M')
            sys.stdout.flush()
        count_stats += 1

end_time = time.time()

if dMpi.ioproc:
    print("Simulation time was: {:.3f}".format(end_time-start_time))
