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
import os

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

# ... for fluid ...
alpha       = dn.cst(1./0.4) # define ideal gas via Cv/R

# ... in space ...
L           = dn.cst(1.0) 
with_length = [L,L]       # domain length in each direction
with_grid   = [128,128]   # number of points in each direction

# ... and time ...
nitmax    = 5000          # number of time steps
with_dt   = dn.cst(1.0e-3) # time step
filtr_amp = dn.cst(0.1)    # filter amplitude

# ... as fast as possible!
with_proc = [2,2]          # mpi proc. topology

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
eqns['coeff'][0][1] = dn.cst(1./alpha) #set R/cv 

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
dMpi  = dtree['mpi']['dMpi']

# .. create the computational grid and write to file
dtree = dn.create_grid(dtree)
dn.dnami_io.write_grid(dtree)

# .. broadcast start-up info
dn.dnami_io.hello_world(dtree)

# .. allocate tree
large = 100000 #no cache blocking in this example
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

q   = dtree['eqns']['qvec']['views']['q'] 

# ================================================================== INITIALISE

# initial clock
ti    = dn.cst(0.)
ni    = 1
delta = dn.cst(1./alpha)

# -------------------------------------------------------------------
# BC to test //  specify the components of the velocity field (e.g. [1,0] for right boundary, [-1,-1] for bottom left) 
#
bc_test = np.asarray([0.,1.],dtype='float64') # TOP
#bc_test = np.asarray([1.,0.],dtype='float64') # RIGHT
#
# -------------------------------------------------------------------

# -- Parameters
r0 = dn.cst(1.)
p0 = dn.cst(1.)
c0 = np.sqrt( dn.cst(1. + delta)*p0/r0  ) 
Rc = dn.cst(0.05)*L
Ma = dn.cst(0.5)

Uinf, Vinf = Ma*c0*bc_test/np.sqrt( bc_test[0]**2 + bc_test[1]**2  )

if dMpi.ioproc:
    print('=========================================')
    print('Testing  with', bc_test , ' (x,y)'        )
    print('=========================================')
    
#-- Initialise vortex at the center of the domain
x0     = dn.cst(0.5*Lx)   		#Vortex center x
y0     = dn.cst(0.5*Ly)  		#Vortex center y 
Uvort  = dn.cst(0.25)   		#Maximum induced Delta |velocity| 
Gam    = Uvort*Rc*np.sqrt(np.exp(1.))   #Vortex strength
Rcsq   = Rc*Rc

XX, YY = xloc[:,np.newaxis], yloc[np.newaxis,:]
RR     = (XX-x0)*(XX-x0)+(YY-y0)*(YY-y0) 
two    = dn.cst(2.)

pf     = - r0*Gam*Gam/(two*Rcsq)*np.exp( -RR/(Rcsq)      ) 
dpsidy = - (YY-y0)/Rcsq*Gam     *np.exp( -RR/(two*Rcsq)  ) 
dpsidx = - (XX-x0)/Rcsq*Gam     *np.exp( -RR/(two*Rcsq)  ) 

# -- Define domain slice
domain = np.s_[0:nx+2*hlo,0:ny+2*hlo]

# -- Fill initial condition
rho[domain] = r0
u[domain]   = Uinf - dpsidy 
v[domain]   = Vinf + dpsidx
et[domain]  = (p0 + pf)/ delta/ r0 + dn.cst(0.5)*(u[domain]**2. + v[domain]**2.)

dMpi.swap(q,hlo,dtree)

# ========================================================================= RUN

intparam,fltparam,data = (dtree['libs']['fort']['integers'],
                                                  dtree['libs']['fort']['floats'],
                                                  dtree['libs']['fort']['data'])

# -- Write first restart
dn.dnami_io.write_restart(0,ti,0,dtree)

# -- Frequency of info
mod_filter = 1 
mod_output = 1000 
mod_rstart = 1000 

for n in range(ni,nitmax+ni):
        ti = ti + dt

        #RK loop
        for nrk in range(1,4):
            intparam[7] = nrk
            dMpi.swap(q,hlo,dtree) 
            dn.time_march(intparam,fltparam,data)    

        if np.mod(n,mod_filter) == 0:
            dMpi.swapX(q,hlo,dtree) 
            dn.filter(1,intparam,fltparam,data)
            dMpi.swapY(q,hlo,dtree) 
            dn.filter(2,intparam,fltparam,data)

        if np.mod(n,mod_rstart) == 0:
            dn.dnami_io.write_restart(n,ti,0,dtree)

        if np.mod(n,mod_output) == 0:

            if dMpi.ioproc:
                    print('____________________________________________________________')
                    print('iteration',n,' with time t =',ti)
            e = et - dn.cst(.5)*(u*u+v*v)
            c = np.sqrt( dn.cst(1. + delta)* (delta * rho [domain]* e[domain]) / rho[domain] )  
            dn.dnami_io.globalMinMax(dtree,rho[domain],'r')
            dn.dnami_io.globalMinMax(dtree,u[domain],'u')
            dn.dnami_io.globalMinMax(dtree,v[domain],'v')
            dn.dnami_io.globalMinMax(dtree,et[domain],'et')

            if dMpi.ioproc:
                    print('convective CFL numbers')
                    sys.stdout.flush()
            cfl = dt*np.abs(u[domain])/dx
            dn.dnami_io.globalMax(dtree,cfl,'cfl-x')
            cfl = dt*np.abs(v[domain])/dy
            dn.dnami_io.globalMax(dtree,cfl,'cfl-y')
            if dMpi.ioproc:
                    print('acoustic CFL numbers')
                    sys.stdout.flush()
            cfl = dt*(np.abs(u[domain])+c)/dx
            dn.dnami_io.globalMax(dtree,cfl,'cfl-x')
            cfl = dt*(np.abs(v[domain])+c)/dy
            dn.dnami_io.globalMax(dtree,cfl,'cfl-y')


# ========================================================================= END 

# -- Get global maximum of pressure as verification 
dMpi.swap(q,hlo,dtree) 
p    = delta * rho * (et - dn.cst(.5)*(u*u+v*v))  
if dMpi.ioproc:
    print(' ----------------------------------------------- ' )
    print('    Maximum of presure field at last iteration   ' )
pmax = dn.dnami_io.globalMax(dtree,p[domain],'')
if dMpi.ioproc:
    print(' ----------------------------------------------- ' )
# -- Write out
np.savetxt('out.dat',np.asarray([pmax]))
