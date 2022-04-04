# -----------------------------------------------------------------------------
#
# CASE: 1D standing shock 
#       -- with kernels: rhs_???.py genKer_???.py
#                                                          -     sw 01-SEP-21
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
# -- grid  
try :
    os.mkdir('./grid/')
except FileExistsError:
    pass

# ========================================================================= ASK

# Solve the equation ...
# ... for fluid ...
alpha    = dn.cst(1./0.4) # Cv/R

# ... in space ...
L = dn.cst(1.) 
with_length = [L]    # domain length in each direction
with_grid   = [480]  # number of points in each direction

# ... and time ...
with_dt   = dn.cst(9.00e-4) # time step
filtr_amp = dn.cst(0.1)    # filter amplitude

# ... as fast as possible!
with_proc = [2] # mpi proc. topology

# ... forcing values for the acoustic wave
amp   = 10e-2
cin   = np.sqrt(1.4)
omega = 80*2.*np.pi*cin/(4.)

# -- First few cycles
T = 4./cin
nitmax = int(4*T/with_dt) 

# -- Shock capturing parameters
mub = dn.cst(2e-3)
Pr  = dn.cst(0.01) 

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
eqns['coeff'][0][1] = dn.cst(1.0/alpha) #set delta
eqns['coeff'][1][1] = amp  
eqns['coeff'][2][1] = omega 
eqns['coeff'][4][1] = mub 
eqns['coeff'][5][1] = Pr*mub 

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
u   = dtree['eqns']['qvec']['views']['u']
et  = dtree['eqns']['qvec']['views']['et']

q  = dtree['eqns']['qvec']['views']['q'] 

# =================================================================== FUNCTIONS

# -- Equations of State needed 
one   = dn.cst(1.)
gamma = one + one/alpha

def eos_e(rho,p):
        e = p/(rho*(gamma-one))
        return e

def eos_p(rho,e):
        p = (gamma-one)*rho*e
        return p 

def eos_sos(rho,p):
        c =  np.sqrt(gamma*p/rho)
        return c

# ================================================================== INITIALISE

# initial clock
ti = dn.cst(0.0)
ni = 1
trstart = dn.cst(0.)

# -- Initial thermo. conditions
gamma = dn.cst(1. + 1./alpha)

if dMpi.ioproc: 
    print(' --------------------------------------------------------------------------')
    print(' ------------      1D standing shock initial conditions ---------------')
    print(' --------------------------------------------------------------------------')

xl     = x[dMpi.ibeg-1:dMpi.iend+2*hlo] 

one    = dn.cst(1.)
two    = dn.cst(2.)
half   = dn.cst(.5)

M0   = dn.cst(1.5)
rho0 = dn.cst(1.) 
p0   = dn.cst(1.) 
u0   = M0*np.sqrt(gamma*p0/rho0) 
c0   = u0/M0
p    = np.zeros_like(rho)

rho1 = (gamma + one)*M0**2/(two + (gamma-one)*M0**2)* rho0
p1   = (two*gamma*M0**2 - (gamma-one))/(gamma+one)  * p0 
M1   = np.sqrt ( (two + (gamma-one) * M0**2) / (two*gamma*M0**2 - (gamma-one)) ) 
u1   = M1*np.sqrt(gamma*p1/rho1)
c1   = u1/M1

if dMpi.ioproc:
    print(' --------------------------')
    print('Upstream values')
    print('rho0, p0, M0', rho0, p0, M0)
    print('Downstream values')
    print('rho1, p1, M1', rho1, p1, M1)
    print(' --------------------------')

# -- Fill initial field
rho[:] = rho0*( xl < half )  + rho1*(xl >= half) 
u[:]   = u0*( xl < half )    + u1*(xl >= half) 
p[:]   = p0*( xl < half )    + p1*(xl >= half) 
et [:] = eos_e(rho[:],p[:]) + dn.cst(0.5)*u[:]*u[:]

dMpi.swap(q,hlo,dtree)

# ========================================================================= RUN

intparam,fltparam,data = (dtree['libs']['fort']['integers'],
                                                  dtree['libs']['fort']['floats'],
                                                  dtree['libs']['fort']['data'])

rpath = './restarts/'

# -- grid  
try :
    os.mkdir(rpath)
except FileExistsError:
    pass

# -- Write out initial restarts and live_view files
dn.dnami_io.write_restart(0,dn.cst(0.),0,dtree,fpath=rpath)

mod_filter = 1 #1 # 1 #10 #25 #50 #1000 #25
mod_output = 50 #200 # 500 #500

from timeit import default_timer as timer
t0 = timer()

for n in range(ni,nitmax+ni):
    ti = ti + dt

    fltparam[8] = 1*ti

    for nrk in range(1,4):
        intparam[7] = nrk
        dMpi.swap(q,hlo,dtree) 
        dn.dnamiF.time_march(intparam,fltparam,data)    

    if np.mod(n,mod_filter) == 0:
        dMpi.swapX(q,hlo,dtree) 
        dn.dnamiF.filter(1,intparam,fltparam,data)

    if np.mod(n,mod_output) == 0:
        dn.dnami_io.write_restart(n,ti,0,dtree,fpath=rpath)

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

