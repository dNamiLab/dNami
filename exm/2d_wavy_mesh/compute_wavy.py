# -----------------------------------------------------------------------------
#
# CASE: 2D wavy grid validation case 
#       -- with kernels: rhs_???.py genKer_???.py
#                                                                 -et 23-JAN-20
#                                                          -(modif)sw 21-FEB-20
# -----------------------------------------------------------------------------

# ======================================================================== LOAD
import dnami as dn         # dNami kernel

# Third-party libraries
from dnami import np, sys  # call non-dNami libs already available in dNami
import math as m           # ... and some others
import os

# dNami pre-processing tools?
# dNami runtime tools?

# ========================================================================= ASK

# -- restarts 
try :
    os.mkdir('./restarts/')
except FileExistsError:
    pass

# -- grid  
try :
    os.mkdir('./grid/')
except FileExistsError:
    pass

# --out  
try :
    os.mkdir('./out/')
except FileExistsError:
    pass

# -- out/omg
try :
    os.mkdir('./out/omg')
except FileExistsError:
    pass
# ========================================================================= ASK

iRestart = False 
#iRestart = True 

# Solve the equation ...
# ... with EoS model ...
iVDW = False

# ... for fluid ...
alpha    = dn.cst(1./0.4) # Cv/R

# ... in space ...
L = dn.cst(1.0) 
with_length = [L,L]         # domain length in each direction
#with_grid   = [96,48]   # number of points in each direction
#with_grid   = [81,46]   # number of points in each direction
with_grid   = [100,50]   # number of points in each direction

# ... and time ...
#with_dt   = dn.cst(7.5e-4) # time step
with_dt   = dn.cst(1.0e-3) # time step
#nitmax    = 10000  #10000#400000 # total number of time steps
nitmax = int(720*0.5/0.5916/with_dt)
filtr_amp = dn.cst(0.1)    # filter amplitude

# ... as fast as possible!
with_proc     = [4,5] # mpi proc. topology

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

geom['Lx'] = with_length[0] 
geom['Ly'] = with_length[1] 
geom['Lz'] = 0. 

grid['nxgb'] = with_grid[0]
grid['nygb'] = with_grid[1]
grid['nzgb'] = 1

mpi['nxpr'] = with_proc[0]
mpi['nypr'] = with_proc[1]
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
yloc= geom['yloc']

Lx  = geom['Lx']  
Ly  = geom['Ly']  

dx  = geom['dx']  
dy  = geom['dy']  

nx  = dMpi.nx     
ny  = dMpi.ny     

x   = geom['x']
y   = geom['y']

dt  = numerics['tint']['tstep']
hlo = numerics['hlo']

rho = dtree['eqns']['qvec']['views']['rho']
u   = dtree['eqns']['qvec']['views']['u']
v   = dtree['eqns']['qvec']['views']['v']
et  = dtree['eqns']['qvec']['views']['et']

ksi = dtree['eqns']['qvec']['views']['ksi']
eta = dtree['eqns']['qvec']['views']['eta']

dksidx = dtree['eqns']['qvec']['views']['dksidx']
detady = dtree['eqns']['qvec']['views']['detady']
dksidy = dtree['eqns']['qvec']['views']['dksidy']
detadx = dtree['eqns']['qvec']['views']['detadx']
Jm1    = dtree['eqns']['qvec']['views']['Jm1']

U      = dtree['eqns']['qvec']['views']['U']
V      = dtree['eqns']['qvec']['views']['V']

omg = dtree['eqns']['qvec']['views']['omg']

q  = dtree['eqns']['qvec']['views']['q'] 

if 'qstored' in dtree['eqns']['qvec']['views'].keys():
        qstored  = dtree['eqns']['qvec']['views']['qstored'] 

# =================================================================== FUNCTIONS

gamma = 1. + 1./alpha

def eos_e(rho,p):
    e = p/(rho*(gamma-1.))
    return e

def eos_sos(rho,p):
    c =  np.sqrt(gamma*p/rho)
    return c

def eos_p(rho,e):
        p = (gamma-1.)*rho*e
        return p 

# ================================================================== INITIALISE

# initial clock
ti = dn.cst(0.0)
ni = 1
trstart = dn.cst(0.)

# =====================================================================
# ============================== GRID =================================
# =====================================================================

# -- Grid parameters 
rv      = dn.cst(1.)
ksi_min = -dn.cst(12.)*rv 
eta_min = -dn.cst(6.)*rv
A_x = dn.cst(0.4)*rv 
A_y = dn.cst(1.6)*rv
llx = dn.cst(24.)
lly = dn.cst(12.)
Mi = dn.cst(0.5)
Mv = dn.cst(0.5)
gam = 1*gamma
x0 = dn.cst(0.0)
y0 = dn.cst(0.0)

# -- WARNING: the whole of ksi and eta need to be filled for the metrics to be correct  

xx = np.arange(-hlo*dx,1.+hlo*dy, dx)
yy = np.arange(-hlo*dy,1.+hlo*dy, dy)

ksi[:,:] = ksi_min + xx[dMpi.ibeg-1:dMpi.iend+2*hlo,np.newaxis]*llx + A_x*np.sin(dn.cst(2.)*np.pi*yy[np.newaxis,dMpi.jbeg-1:dMpi.jend+2*hlo]) 
eta[:,:] = eta_min + yy[np.newaxis,dMpi.jbeg-1:dMpi.jend+2*hlo]*lly + A_y*np.sin(dn.cst(4.)*np.pi*xx[dMpi.ibeg-1:dMpi.iend+2*hlo,np.newaxis]) 

if dMpi.ioproc:
    print(' -------------------------- ' )
    print('Done reading and filling.')

dn.dnami_io.write_data(['ksi'],0,0,dtree,'./grid/','ksi')
dn.dnami_io.write_data(['eta'],0,0,dtree,'./grid/','eta')

# =====================================================================
# ============================== INIT =================================
# =====================================================================

one = dn.cst(1.)
two = dn.cst(2.)
rhoinf = 1.*one
pinf   = 1.*one
uinf = Mi*eos_sos(rhoinf, pinf)

# -- Fill with the actual initial solution
for i in range(hlo,nx+hlo):
    for j in range(hlo,ny+hlo):
        rij = np.sqrt( (ksi[i,j] - x0)**2 + (eta[i,j]-y0)**2   )/rv
        u[i,j] = uinf * (one - Mv/Mi*(eta[i,j]-y0)/rv*np.exp( (one - rij**2)/two ) )
        v[i,j] = uinf * (Mv/Mi*(ksi[i,j]-x0)/rv*np.exp( (one-rij**2)/two  ))
        rho[i,j] = rhoinf * ( one - (gam-one)/two*Mv**2*np.exp(one-rij**2) )**( one/(gam-one)  )
        pij = pinf * ( one - (gam-one)/two*Mv**2*np.exp(one-rij**2) )**( gam/(gam-one)  )
        et[i,j] = eos_e(rho[i,j], pij) + dn.cst(0.5)* ( u[i,j]**2 + v[i,j]**2 )

dMpi.swap(q,hlo,dtree)

# ========================================================================= RUN

intparam,fltparam,data = (dtree['libs']['fort']['integers'],
                                                  dtree['libs']['fort']['floats'],
                                                  dtree['libs']['fort']['data'])

# -- Compute the stored values i.e. the metrics using whatever scheme is used here
if 'qstored' in dtree['eqns']['qvec']['views'].keys():
        dn.stored(intparam,fltparam,data,1) # compute static
        dn.stored(intparam,fltparam,data,0)        

dMpi.swap(qstored,hlo,dtree)

dn.dnami_io.write_data(['omg'],0,0,dtree,'./out/omg/','omg')
dn.dnami_io.write_restart(0,dn.cst(0.),0,dtree)

if dMpi.ioproc:
    print(' -------------------------- ' )
    print('About to start, printing some min/max:')
dn.dnami_io.globalMinMax(dtree,rho,'r')
dn.dnami_io.globalMinMax(dtree,u,'u')
dn.dnami_io.globalMinMax(dtree,v,'v')
dn.dnami_io.globalMinMax(dtree,et,'et')
if dMpi.ioproc:
    print(' -------------------------- ' )
    print('Metrics:')
dn.dnami_io.globalMinMax(dtree,dksidx,'dksidx')
dn.dnami_io.globalMinMax(dtree,dksidy,'dksidy')
dn.dnami_io.globalMinMax(dtree,detadx,'detadx')
dn.dnami_io.globalMinMax(dtree,detady,'detady')
dn.dnami_io.globalMinMax(dtree,Jm1,'Jm1')
dn.dnami_io.globalMinMax(dtree,omg,'omg')
if dMpi.ioproc:
    print(' -------------------------- ' )
    print('Needs nitmax:', nitmax)

mod_filter = 1 
mod_output = 1000 
mod_rstart = 1000 

#scheme = '11_10'
#scheme = '9_8'
#scheme = '7_6'
#scheme = '5_4'
#scheme = '3_2'

#rpath = f'./restarts_{scheme}/'
rpath = f'./restarts/'

try:
    import os
    os.mkdir(rpath)
except:
    pass

from timeit import default_timer as timer
t0 = timer()

for n in range(ni,nitmax+ni):
        ti = ti + dt

        for nrk in range(1,4):
            intparam[7] = nrk
            dMpi.swap(q,hlo,dtree) 
            if 'qstored' in dtree['eqns']['qvec']['views'].keys():
                dMpi.swap(q,hlo,dtree) 
                dn.stored(intparam,fltparam,data,0)
            dn.time_march(intparam,fltparam,data)

        if np.mod(n,mod_filter) == 0:
            dMpi.swapX(q,hlo,dtree)
            dn.filter(1,intparam,fltparam,data)
            dMpi.swapY(q,hlo,dtree)
            dn.filter(2,intparam,fltparam,data)

        if np.mod(n,mod_rstart) == 0 or n == nitmax + ni -1:
            dn.dnami_io.write_restart(n,ti,0,dtree,fpath=rpath)

        if np.mod(n,mod_output) == 0  :
            dn.dnami_io.write_data(['omg'],n,ti,dtree,'./out/omg/','omg')

            if dMpi.ioproc:
                    print('____________________________________________________________')
                    print('iteration',n,' with time t =',ti)
            #e = et - .5*(u*u+v*v)
            #p = eos_p(rho,e)
            #c = eos_sos(rho[hlo:nx+hlo,hlo:ny+hlo],p[hlo:nx+hlo,hlo:ny+hlo])
            dn.dnami_io.globalMinMax(dtree,rho[hlo:nx+hlo,hlo:ny+hlo],'r')
            dn.dnami_io.globalMinMax(dtree,u[hlo:nx+hlo,hlo:ny+hlo],'u')
            dn.dnami_io.globalMinMax(dtree,v[hlo:nx+hlo,hlo:ny+hlo],'v')
            dn.dnami_io.globalMinMax(dtree,et[hlo:nx+hlo,hlo:ny+hlo],'et')
            #dn.dnami_io.globalMinMax(dtree,np.abs(u[hlo:nx+hlo,hlo:ny+hlo]/c),'M')
            #dn.dnami_io.globalMinMax(dtree,U,'U')
            #dn.dnami_io.globalMinMax(dtree,V,'V')
            #dn.dnami_io.globalMinMax(dtree,Jm1,'Jm1')
            #dn.dnami_io.globalMinMax(dtree,omg,'omg')
            #if dMpi.ioproc:
            #        print('convective CFL numbers')
            sys.stdout.flush()
            #cfl = dt*np.abs(u[hlo:nx+hlo,hlo:ny+hlo])/dx
            #dn.dnami_io.globalMax(dtree,cfl,'cfl-x')
            #if dMpi.ioproc:
            #        print('acoustic CFL numbers')
            #        sys.stdout.flush()
            #cfl = dt*(np.abs(u[hlo:nx+hlo,hlo:ny+hlo])+c)/dx
            #dn.dnami_io.globalMax(dtree,cfl,'cfl-x')

