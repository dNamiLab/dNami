# -----------------------------------------------------------------------------
#
# CASE: 2D vortex test case 
#       -- with kernels: rhs_???.py genKer_???.py
#                                                          -sw 15-JUN-20
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

iRestart = False  

# Solve the equation ...
# ... with EoS model ...
iVDW = True

# ... for fluid ...
alpha    = dn.cst(1./0.4) # Cv/R

# ... in space ...
L = dn.cst(1.0) 
with_length = [L,L]  # domain length in each direction

with_grid   = [576,576]   # number of points in each direction

# ... and time ...
with_dt   = dn.cst(5.0e-4) # time step
nitmax    = 1000
filtr_amp = dn.cst(0.1)    # filter amplitude

# ... as fast as possible!
with_proc = [1,1] # mpi proc. topology

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
eqns['coeff'][1][1] = dn.cst(0.)      #set mub 

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
u  = dtree['eqns']['qvec']['views']['u']
v  = dtree['eqns']['qvec']['views']['v']
et  = dtree['eqns']['qvec']['views']['et']

q  = dtree['eqns']['qvec']['views']['q'] 

if 'qstored' in dtree['eqns']['qvec']['views'].keys():
        qstored  = dtree['eqns']['qvec']['views']['qstored'] 

# =================================================================== FUNCTIONS

# -- Grab functions to fill out nozzle baseflow 
sys.path.append('./func_nz/')

if iVDW:
    gas_type = 1
    Zc = np.float64(3./8.)
    gas_type_str = 'vdw'
else:
    gas_type = 0
    Zc = np.float64(1.)
    gas_type_str = 'ig'

import func_nz.values as val
val.fill_values(Zc, 1./alpha, gas_type )
import func_nz.eos as eos
from func_nz.read_geo import read_geo 

# ================================================================== INITIALISE

# initial clock
ti = dn.cst(0.0)
ni = 1
trstart = dn.cst(0.)

# BC to test  (0: right, 1: top, 2: left, 3: bot)
bc_test = 0

# -- Parameters
r0 = dn.cst(1.)
p0 = dn.cst(1.)
c0 = eos.eos_sos(r0,p0)
Rc = dn.cst(0.05)*L
Ma = dn.cst(0.5)

if bc_test ==0:
    Uinf = Ma*c0
    Vinf = dn.cst(0.)
elif bc_test == 1:
    Uinf = dn.cst(0.) 
    Vinf = Ma*c0
elif bc_test == 2:
    Uinf = -Ma*c0
    Vinf = dn.cst(0.)
elif bc_test == 3:
    Uinf = dn.cst(0.) 
    Vinf = -Ma*c0
else:
    if dMpi.ioproc: 
        print('PLEASE SPECIFIC BC')
        exit()

x0 = dn.cst(0.5)*Lx   			#Vortex center
y0 = dn.cst(0.5)*Ly   			#Vortex center
Uvort = dn.cst(0.1)   			#Maximum induced |velocity| -- specify direction of baseflow above
Gam = Uvort*Rc*np.sqrt(np.exp(1.))      #Vortex strength

#-- INIT VORTEX

for i in range(nx+2*hlo):
    for j in range(ny+2*hlo):

        xll  = xloc[0] + (i-hlo)*dx
        yll  = yloc[0] + (j-hlo)*dy

        xx = xll - x0 
        yy = yll - y0

        pf       = - r0*Gam**2/(2*Rc**2)*np.exp( -(xx**2+yy**2)/(Rc**2)  ) 
        dpsidy   = - yy/Rc**2*Gam*np.exp( -(xx**2+yy**2)/(2*Rc**2)  ) 
        dpsidx    = - xx/Rc**2*Gam*np.exp( -(xx**2+yy**2)/(2*Rc**2)  ) 

        rho[i,j] = 1*r0
        u[i,j]   = Uinf - dpsidy 
        v[i,j]   = Vinf + dpsidx
        pij      = p0 + pf 
        et[i,j]  = eos.eos_e(r0,pij) + dn.cst(0.5)*(u[i,j]**2 + v[i,j]**2)

# --------------------------------------

if dMpi.ioproc: 
    print('Done filling field ....')

dMpi.swap(q,hlo,dtree)

# ========================================================================= RUN

intparam,fltparam,data = (dtree['libs']['fort']['integers'],
                                                  dtree['libs']['fort']['floats'],
                                                  dtree['libs']['fort']['data'])

# -- Compute the stored values i.e. the metrics using whatever scheme is used here
if 'qstored' in dtree['eqns']['qvec']['views'].keys():
        dMpi.swap(q,hlo,dtree)
        dMpi.swap(qstored,hlo,dtree)
        dn.dnamiF.stored(intparam,fltparam,data,1)        

if dMpi.ioproc:
    print(' -------------------------- ' )
    print('About to start, printing some min/max:')
dn.dnami_io.globalMinMax(dtree,rho,'r')
dn.dnami_io.globalMinMax(dtree,u,'u')
dn.dnami_io.globalMinMax(dtree,v,'v')
dn.dnami_io.globalMinMax(dtree,et,'et')


# -- Write first restart
dn.dnami_io.write_restart(0,ti,0,dtree)
dn.dnami_io.write_restart(0,ti,1,dtree)

mod_filter = 1 # 1 #10 #25 #50 #1000 #25
mod_output = 50 # 500 #500
mod_rstart = 50 

from timeit import default_timer as timer
t0 = timer()
tini = ti

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

        if np.mod(n,mod_filter) == 0:
            dn.dnamiF.filter(1,intparam,fltparam,data)
            dn.dnamiF.filter(2,intparam,fltparam,data)

        if np.mod(n,mod_rstart) == 0:
            dn.dnami_io.write_restart(n,ti,0,dtree)

        if np.mod(n,mod_output) == 0:
            dn.dnami_io.write_restart(n,ti,1,dtree)

            if dMpi.ioproc:
                    print('____________________________________________________________')
                    print('iteration',n,' with time t =',ti)
            e = et - .5*(u*u+v*v)
            p = eos.eos_p(rho,e)
            c = eos.eos_sos(rho[hlo:nx+hlo,hlo:ny+hlo],p[hlo:nx+hlo,hlo:ny+hlo])
            U = np.sqrt(u*u + v*v)
            dn.dnami_io.globalMinMax(dtree,rho[hlo:nx+hlo,hlo:ny+hlo],'r')
            dn.dnami_io.globalMinMax(dtree,u[hlo:nx+hlo,hlo:ny+hlo],'u')
            dn.dnami_io.globalMinMax(dtree,v[hlo:nx+hlo,hlo:ny+hlo],'v')
            dn.dnami_io.globalMinMax(dtree,et[hlo:nx+hlo,hlo:ny+hlo],'et')

            if dMpi.ioproc:
                    print('convective CFL numbers')
                    sys.stdout.flush()
            cfl = dt*np.abs(U[hlo:nx+hlo,hlo:ny+hlo])/dx
            dn.dnami_io.globalMax(dtree,cfl,'cfl-x')
            if dMpi.ioproc:
                    print('acoustic CFL numbers')
                    sys.stdout.flush()
            cfl = dt*(np.abs(U[hlo:nx+hlo,hlo:ny+hlo])+c)/dx
            dn.dnami_io.globalMax(dtree,cfl,'cfl-x')

t1  = timer()
print('Total comp time:', t1-t0)
