# -----------------------------------------------------------------------------
#
# CASE: Two-way coupled shallow water equations for the Tonga explosion event 
#       -- with kernels: rhs.py genKer.py
#
# Required files:
#  o baseflow values 
#       - inputs/thermo.bin
#       - inputs/bathy.bin
#       - inputs/velocity.bin
#  o reference values:
#       - values.py
#  o pre-processing utility 
#       - misc.py
#  o positions of sensor to interpol fields to 
#       - inputs/dart.csv
#       - inputs/psens.pos
#
# All the inputs can be found in the ./data/dnami/inputs folder of the 
#  Zenodo submission: 10.5281/zenodo.7197431 
#  
# /!\ Prior to running the compute.py, users should run 'python3 misc.py NX NY'
# where NX and NY are the desired grid size to interpolate the baseflow fields  
# 
# Due to MPI communications being used to enforce spherical boundary conditions
# this case must be run with at least 2 processors (with an even split in y)
#
# Please keep an eye out for updates on the repository  
#                       on GitHub  : https://github.com/dNamiLab/dNami
#
#
#                                                          -sw 15-JUN-20
# -----------------------------------------------------------------------------

# ======================================================================== LOAD
import dnami as dn         # dNami kernel

# Third-party libraries
from dnami import np, sys  # call non-dNami libs already available in dNami
import math as m           # ... and some others
import os

# Custom functions 
import misc
import time
import glob 
import pandas as pd
from scipy.ndimage import map_coordinates as MC

# ========================================================================= ASK

iRestart   = False

# Grab atmospheric property class and check if 1 or 2 layer
from equations import iAtmos
from values import vals
val = vals()

# Solve the equation ...

# ... in space ...
Lt = dn.cst(np.pi)
Lp = dn.cst(2*np.pi)
with_length = [Lt,Lp]  # domain length in each direction
with_grid   = [800,1600]   # number of points in each direction

# ... and time ...
filtr_amp = dn.cst(0.1)    # filter amplitude
with_cfl  = dn.cst(0.10) # CFL constraint for time step

# ... as fast as possible!
with_proc = [1,2] # mpi proc. topology

# -- 
Re      = 6371e3 #Earth radius

# -- Ratio of earth radius to average depth 
Rad = Re/val.lref

# =================================================================== Create WRK

if iAtmos:
    resf = './restarts_2layer/'
    outf = './out_2layer/'
else:
    resf = './restarts_1layer/'
    outf = './out_1layer/'

# -- Create folders
if not iRestart:
    for fldr in [resf, outf, outf+'historical', outf+'force', outf+'stations', outf+'dart']:
        try :
            os.mkdir(fldr)
        except FileExistsError:
            pass

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
eqns['coeff'][0][1] = dn.cst(0.) #perturbation time dependence 
eqns['coeff'][1][1] = dn.cst(1./Rad) 
eqns['coeff'][2][1] = dn.cst(1.)  
eqns['coeff'][3][1] = dn.cst(0.) 

if iAtmos:
    eqns['coeff'][4][1] = val.gamma 
    eqns['coeff'][5][1] = val.hERA5/val.lref 

geom['Lx'] = with_length[0] 
geom['Ly'] = with_length[1] 
geom['Lz'] = 0. 

grid['nxgb'] = with_grid[0]
grid['nygb'] = with_grid[1]
grid['nzgb'] = 1

mpi['nxpr'] = with_proc[0]
mpi['nypr'] = with_proc[1]
mpi['nzpr'] = 1 

numerics['tint']['tstep'] = dn.cst(0.) 
numerics['filtr']['eps']  = filtr_amp

# .. start the message passing interface
dtree = dn.start_mpi(dtree)
dMpi = dtree['mpi']['dMpi']

# .. create the computational grid and write to file
dtree = dn.create_grid(dtree)
dn.dnami_io.write_grid(dtree,fpath=outf)

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

# -- Water related 
eta1 = dtree['eqns']['qvec']['views']['eta1']
vt   = dtree['eqns']['qvec']['views']['vt']
vp   = dtree['eqns']['qvec']['views']['vp']

if iAtmos:
    # -- Air related 
    rhoa = dtree['eqns']['qvec']['views']['rhoa']
    ut   = dtree['eqns']['qvec']['views']['ut']
    up   = dtree['eqns']['qvec']['views']['up']
    pa   = dtree['eqns']['qvec']['views']['pa']

# --  Compute/fill at start
sint = dtree['eqns']['qvec']['views']['sint'] 
sintinv = dtree['eqns']['qvec']['views']['sintinv'] 
cost = dtree['eqns']['qvec']['views']['cost'] 
eta0 = dtree['eqns']['qvec']['views']['eta0'] 

# -- Historical values
hmax = dtree['eqns']['qvec']['views']['hmax']
tmax = dtree['eqns']['qvec']['views']['tmax']

# -- Forcing terms
fss = {} 
for var in dtree['eqns']['qvec']['solved']:
    varstr = var[0]
    fss[varstr] = dtree['eqns']['qvec']['views']['rhs_'+varstr]
fsskeys  = fss.keys()

# -- Initial condition
Seta1 = dtree['eqns']['qvec']['views']['Seta1']
Svt   = dtree['eqns']['qvec']['views']['Svt']
Svp   = dtree['eqns']['qvec']['views']['Svp']

if iAtmos:
    Srhoa = dtree['eqns']['qvec']['views']['Srhoa']
    Sut   = dtree['eqns']['qvec']['views']['Sut']
    Sup   = dtree['eqns']['qvec']['views']['Sup']
    Spa   = dtree['eqns']['qvec']['views']['Spa']

q  = dtree['eqns']['qvec']['views']['q'] 

if 'qstored' in dtree['eqns']['qvec']['views'].keys():
    qstored  = dtree['eqns']['qvec']['views']['qstored'] 

# --- Modify MPI topology for spherical exchanges

# -- Require an even split
if with_proc[1] % 2 != 0:
    if dMpi.ioproc: print('[ERROR] y-direction split must be even ...')
    exit()

# -- Define a new position/negative direction
pid = dMpi.procid
mod = np.mod(pid, with_proc[1])

if mod <= int(with_proc[1]/2) - 1 :
    #Copy 
    nxm = 1*dMpi.neighxm  
    nxp = 1*dMpi.neighxp  
    #Swap 
    dMpi.neighxm = nxp
    dMpi.neighxp = nxm
    dMpi.flipped = True
    #Flip ud and lr 
    if dMpi.ibeg-1 == 0 or dMpi.iend == grid['nxgb']:
        dMpi.reverse = True

# -- South pole procs
if dMpi.ibeg-1 == 0:
    mid = int(with_proc[1]/2) - 1  
    if pid <= mid:
        dp  = abs(mid-pid)
        dMpi.neighxp = mid + dp + 1 
        dMpi.edge = True
    else:
        dp  = abs(mid+1-pid)
        dMpi.neighxm = mid - dp 
        #Flip ud and lr 
        dMpi.reverse = True

# -- North pole procs
if dMpi.iend == grid['nxgb']:
    mid = with_proc[0]*with_proc[1] -1 - int(with_proc[1]/2) 
    if pid <= mid:
        dp  = abs(mid-pid)
        dMpi.neighxm = mid + dp + 1 
    else:
        dp  = abs(mid+1-pid)
        dMpi.neighxp = mid - dp 
        #Flip ud and lr 
        dMpi.reverse = True
        dMpi.edge = True

## -- Velocity flipping - lat component must be opposite (in this case x component)
dMpi.sym    = dMpi.reverse and dMpi.edge 
if iAtmos:
    dMpi.velidx = [1,4] #vt, ut
else:
    dMpi.velidx = [1] #vt

# -- Display updated "torus"
dMpi.showTorus()


# ================================================================= FUNCTIONS

def gauss(x0,y0,amp,rad0):
    # - extra modules
    from haversine import haversine, haversine_vector, Unit
    # - convert back to degrees 
    rad0m  = rad0*Re
    x0d    = x0*180./np.pi - 90.
    y0d    = y0*180./np.pi - 180.
    xld    = xloc[:,np.newaxis]*180./np.pi - 90.
    yld    = yloc[np.newaxis,:]*180./np.pi - 180.
    XX,YY  = np.meshgrid(xld, yld, indexing='ij')
    darray = haversine_vector(list(zip(x0d*np.ones_like(XX).ravel(),y0d*np.ones_like(YY).ravel())),list(zip(XX.ravel(), YY.ravel())), Unit.METERS).reshape((nx,ny))
    return amp*np.exp(- (darray**2)/(2.*rad0m**2)) 

def set_dt(target_cfl):
    dtMax = None      
    epsi  = 1e-18
    # --- Grid coeff 
    coeff_c = 10.
    coeff = np.clip(sint[dom], np.amin(np.sin(x[:])) * coeff_c , a_max=None)
    # -- CFL CONV 
    dts = target_cfl*(dx*Rad)/np.abs(vt[dom] + epsi)
    dtMax1 = dn.dnami_io.dtMax(dtree,dts)             
    dts = target_cfl*(dy*Rad*coeff)/np.abs(vp[dom]+ epsi)
    dtMax2 = dn.dnami_io.dtMax(dtree,dts)              
    if dMpi.ioproc:                                   
        dtMaxConv = np.amin([dtMax1,dtMax2])     
        sys.stdout.flush()                      
    # -- CFL SWE 
    c_swe = np.sqrt(np.abs(eta0)) 
    c_swe[above_water] = dn.cst(0.)
    dts = target_cfl*(dx*Rad)/np.abs(vt[dom] + c_swe[dom] + epsi)
    dtMax1 = dn.dnami_io.dtMax(dtree,dts)             
    dts = target_cfl*(dy*Rad*coeff)/np.abs(vp[dom] + c_swe[dom] + epsi)
    dtMax2 = dn.dnami_io.dtMax(dtree,dts)              
    if dMpi.ioproc:                                   
        dtMaxSWE = np.amin([dtMax1,dtMax2])     
        sys.stdout.flush()                      
    # -- CFL SWE 
    if dMpi.ioproc:  
        dtMax = np.amin([dtMaxConv, dtMaxSWE])

    if iAtmos:
        # -- CFL CONV 
        dts = target_cfl*(dx*Rad)/np.abs(ut[dom] + epsi)
        dtMax1 = dn.dnami_io.dtMax(dtree,dts)             
        dts = target_cfl*(dy*Rad*coeff)/np.abs(up[dom]+ epsi)
        dtMax2 = dn.dnami_io.dtMax(dtree,dts)              
        if dMpi.ioproc:                                   
            dtMaxConv = np.amin([dtMax1,dtMax2])     
            sys.stdout.flush()                      
        # -- CFL SOUND 
        c_air = np.sqrt(val.gamma * pa/rhoa) 
        dts = target_cfl*(dx*Rad)/np.abs(ut[dom] + c_air[dom] + epsi)
        dtMax1 = dn.dnami_io.dtMax(dtree,dts)             
        dts = target_cfl*(dy*Rad*coeff)/np.abs(up[dom] + c_air[dom] + epsi)
        dtMax2 = dn.dnami_io.dtMax(dtree,dts)              
        if dMpi.ioproc:                                   
            dtMaxAir = np.amin([dtMax1,dtMax2])     
            sys.stdout.flush()                      
        # -- CFL SWE 
        if dMpi.ioproc:  
            dtMax = np.amin([dtMaxConv, dtMaxAir, dtMax])

    dn.dnami_io.set_dt(dtree,dtMax)                      
    dt = numerics['tint']['tstep']
    return dt          


# ================================================================== INITIALISE

# initial clock
ti = dn.cst(0.0)
ni = 1
trstart = dn.cst(0.)

# -- Parameters
r0 = dn.cst(1.)

# -- Tonga volcano coordinates - will recenter to closest grid point
lat, lon = -20.5428391431001, -175.38900969908855
lat = lat + 90 
lon = lon + 180 

# -- Fill stored:
dom = np.s_[hlo:nx+hlo,hlo:ny+hlo]

sintinv[dom] = 1./np.sin(xloc[:,np.newaxis])
sint[dom]    = np.sin(xloc[:,np.newaxis])
cost[dom]    = np.cos(xloc[:,np.newaxis])

xT = lat/180*np.pi
yT = lon/180*np.pi

# ================== WATER 
# -- Load bathymetry (Warning here values are in meters)
dn.dnami_io.read_data(dtree, ['eta0'], fname='./bathy/bathy_{}_{}.bin'.format(with_grid[1], with_grid[0]))

# -- Normalise bathymetry
eta0[dom] = eta0[dom]/val.lref
above_water = eta0 >= 0.

# -- Initial condition function
def init(info=False):

    # =======================  Water
    eta1[dom] = dn.cst(0.) 
    vt  [dom] = dn.cst(0.) 
    vp  [dom] = dn.cst(0.)

    hmax[dom] = dn.cst(0.) 
    tmax[dom] = dn.cst(0.) 

    # ======================= Air 
    if iAtmos:
        # - Load dimensional fields (created by misc.py)
        dn.dnami_io.read_data(dtree, ['pa'], fname='./thermo/press_{}_{}.bin'.format(with_grid[1], with_grid[0]))
        dn.dnami_io.read_data(dtree, ['rhoa'], fname='./thermo/dens_{}_{}.bin'.format(with_grid[1], with_grid[0]))
        dn.dnami_io.read_data(dtree, ['ut'], fname='./velocity/utheta_{}_{}.bin'.format(with_grid[1], with_grid[0]))
        dn.dnami_io.read_data(dtree, ['up'], fname='./velocity/uphi_{}_{}.bin'.format(with_grid[1], with_grid[0]))

        # -- Rescale with reference values
        pa[:,:] = pa [:,:]/val.pref
        rhoa[:,:] = rhoa [:,:]/val.rref
        ut[:,:] = ut[:,:]/val.uref
        up[:,:]= up[:,:]/val.uref

    return 

# -- Initial the field for the first time (used to generate the forcing vectors)
init(info=True)


if dMpi.ioproc: 
    print('Done filling field ....')
    ('----------------------------')

dMpi.swap(q,hlo,dtree)
dMpi.swap(qstored,hlo,dtree,spherical=False)


# ========================================================================= RUN

intparam,fltparam,data = (dtree['libs']['fort']['integers'],
                                                  dtree['libs']['fort']['floats'],
                                                  dtree['libs']['fort']['data'])


## -- Compute the stored values if any 
if 'qstored' in dtree['eqns']['qvec']['views'].keys():
    dMpi.swap(q,hlo,dtree)
    dMpi.swap(qstored,hlo,dtree,spherical=False)
    dn.dnamiF.stored(intparam, fltparam, data, 0)
    dn.dnamiF.stored(intparam, fltparam, data, 1)
    dMpi.swap(qstored,hlo,dtree,spherical=False)

## -- Print some values for verification 
if dMpi.ioproc:
    print(' -------------------------- ' )
    print('About to start, printing some min/max:')
    print(' -------------------------- ' )
    print('Water')
dn.dnami_io.globalMinMax(dtree,eta1[dom],'eta1 ')
dn.dnami_io.globalMinMax(dtree,vt  [dom],  'vt   ')
dn.dnami_io.globalMinMax(dtree,vp  [dom],  'vp   ')
if iAtmos:
    if dMpi.ioproc:
        print(' -------------------------- ' )
        print('Air')
    dn.dnami_io.globalMinMax(dtree,rhoa[dom], 'rhoa ')
    dn.dnami_io.globalMinMax(dtree,ut  [dom], 'ut   ')
    dn.dnami_io.globalMinMax(dtree,up  [dom], 'up   ')
    dn.dnami_io.globalMinMax(dtree,pa  [dom], 'pa   ')
if dMpi.ioproc:
    print(' -------------------------- ' )
    print('Misc')
dn.dnami_io.globalMinMax(dtree,sintinv[dom],'sintinv ')
dn.dnami_io.globalMinMax(dtree,sint[dom],'sint ')
dn.dnami_io.globalMinMax(dtree,cost[dom],'cost ')
dn.dnami_io.globalMinMax(dtree,eta0[dom],'eta0 ')
if dMpi.ioproc:
    print(' -------------------------- ' )


# -- Write first restart
if not iRestart:
    dn.dnami_io.write_restart(0,ti,0,dtree,fpath=resf)

# --- Add perturbation function 

def add_pert(x0,y0):

    # == Zero all the source terms 
    Seta1[:,:] = dn.cst(0.)
    Svt[:,:]   = dn.cst(0.)
    Svp[:,:]   = dn.cst(0.)
    if iAtmos:
        Srhoa[:,:] = dn.cst(0.)
        Sut[:,:]   = dn.cst(0.)
        Sup[:,:]   = dn.cst(0.)
        Spa[:,:]   = dn.cst(0.)
    
    # - Center on closest point to target
    nxp = int( x0*grid['nxgb']/Lx )    
    nyp = int( y0*grid['nygb']/Ly )    
    
    x0 = x[nxp]
    y0 = y[nyp]
    
    # ------ Compute source parameters 
    diam      = 50e3               # Volcano size estimate [m]
    max_d     = 5*max(dx,dy) 
    tonga_rad = max(max_d,diam/ Re) # m to angle           [rad] 
    res_ratio = max_d/(diam/Re) 

    # ------ Pressure perturbation
    ppertdim  = 5.2e2 # [Pa] 
    sigmax    = tonga_rad*Re/val.lref
    amp       = ppertdim/val.pref/np.sqrt(2.*sigmax**2)

    # ------ Gaussian spatial support 
    gaussf = gauss(x0,y0,1.,tonga_rad) 

    # =================================================
    H0    = np.abs(eta0[dom]) * (eta0[dom] < 0. )
    h0    = val.hERA5/val.lref - np.maximum(eta0, np.zeros_like(eta0))[dom]

    if iAtmos:
        Ca2   = h0
        C02   = val.gamma*pa[dom]/rhoa[dom] 
        beta  = rhoa[dom] 
        phi   = beta*Ca2/C02
        Cw2   = H0
        CT2   = pa[dom]/rhoa[dom] 

        lam2  = 0.5*(Cw2 + C02 + np.sqrt((Cw2-C02)**2 + 4.*beta*Cw2*(C02-CT2+Ca2)))

        etahat = phi*H0
        rhohat = rhoa[dom]    *(lam2/C02 + Cw2/C02*(beta-1.))
        phat   = rhoa[dom]*C02*(lam2/C02 + Cw2/C02*(beta-1.))

        #                         factor two from (A+) + (A-)
        Seta1[dom] += etahat/phat *2.*amp*gaussf
        Srhoa[dom] += rhohat/phat *2.*amp*gaussf
        Spa[dom]   += 1.          *2.*amp*gaussf
    else:
        Ca2   = val.hERA5/val.lref 
        C02   = val.gamma*(val.pERA5/val.pref)/(val.rERA5/val.rref)
        beta  = (val.rERA5/val.rref) 
        phi   = beta*Ca2/C02
        Cw2   = H0
        CT2   = (val.pERA5/val.pref)/(val.rERA5/val.rref)

        lam2  = 0.5*(Cw2 + C02 + np.sqrt((Cw2-C02)**2 + 4.*beta*Cw2*(C02-CT2+Ca2)))

        etahat = phi*H0
        rhohat = val.rERA5/val.rref    *(lam2/C02 + Cw2/C02*(beta-1.))
        phat   = val.rERA5/val.rref*C02*(lam2/C02 + Cw2/C02*(beta-1.))

        #                         factor two from (A+) + (A-)
        Seta1[dom] += etahat/phat *2.*amp*gaussf

    # =================================================
    dMpi.swap(qstored,hlo,dtree,spherical=False)

    dn.dnami_io.globalMinMax(dtree,eta1[dom],'eta1 ')
    dn.dnami_io.globalMinMax(dtree,vt[dom],  'vt   ')
    dn.dnami_io.globalMinMax(dtree,vp[dom],  'vp   ')
    if iAtmos:
        if dMpi.ioproc:
            print(' -------------------------- ' )
            print('Air')
        dn.dnami_io.globalMinMax(dtree,rhoa[dom], 'rhoa ')
        dn.dnami_io.globalMinMax(dtree,ut[dom]  , 'ut   ')
        dn.dnami_io.globalMinMax(dtree,up[dom]  , 'up   ')
        dn.dnami_io.globalMinMax(dtree,pa[dom]  , 'pa   ')
    if dMpi.ioproc:
        print(' -------------------------- ' )
        print('Perturbation')
    dn.dnami_io.globalMinMax(dtree,Seta1[dom],'Seta1 ')
    dn.dnami_io.globalMinMax(dtree,Svt[dom],  'Svt   ')
    dn.dnami_io.globalMinMax(dtree,Svp[dom],  'Svp   ')
    if iAtmos:
        dn.dnami_io.globalMinMax(dtree,Srhoa[dom],'Srhoa ')
        dn.dnami_io.globalMinMax(dtree,Sut[dom]  ,'Sut   ')
        dn.dnami_io.globalMinMax(dtree,Sup[dom]  ,'Sup   ')
        dn.dnami_io.globalMinMax(dtree,Spa[dom]  ,'Spa   ')
    sys.stdout.flush()

    return

# -- Temporal forcing signal shape
def update_ft(ta):

    taustar = 0.75*45*60    #seconds
    tau     = taustar/val.tref  
    ft      = 0.
    pp      = 1. 
    pm      = -0.10

    a = (512*pm - 512*pp)/(9*tau**5)
    b = (-1152*pm + 1408*pp)/(9*tau**4) 
    c =  (768*pm - 1280*pp)/(9*tau**3) 
    d = (-128*pm + 384*pp)/(9*tau**2)

    # - Force if time less than tau 
    if ta <= tau:
        ft= a*ta**5 + b*ta**4 + c*ta**3 + d*ta**2 

    # -- Update coefficient
    fltparam[5] = ft
    dtree['eqns']['coeff'][0][1] = ft
    return 

# ============== Load station coordinates for on-the-fly output ========================

# Create dictionnary per proc
coords_atm = {}
coords_dart = {}

if iAtmos:
    # -- ATMOSPHERIC PRESSURE STATIONS
    
    atms = pd.read_csv('inputs/psens.pos')
    lats = atms['lat']
    lons = atms['lon']
    codes = atms['code']

    for (key, lat, lon) in zip(codes,lats,lons):

        # -- Add coord to dic if point in proc 
        xp = (lat+90)/180*np.pi
        yp = (lon+180)/180*np.pi

        if (xp > xloc[0]-0.5*dx and xp<xloc[-1]+0.5*dx ) and (yp>yloc[0]-0.5*dy and yp<yloc[-1]+0.5*dy):
            coords_atm[key] = [xp, yp]

# -- DART BUOYS 
darts = pd.read_csv('inputs/dart.csv')
codes = darts['Number']
lats  = darts['Lat ']
lons  = darts['Lon']

for (key, lat, lon) in zip(codes,lats,lons):

    # -- Add coord to dic if point in proc 
    xp = (lat+90)/180*np.pi
    yp = (lon+180)/180*np.pi

    if (xp > xloc[0]-0.5*dx and xp<xloc[-1]+0.5*dx ) and (yp>yloc[0]-0.5*dy and yp<yloc[-1]+0.5*dy):
        coords_dart[key] = [xp, yp]

# =========================================================================================

mod_filter = 1 # 1 #10 #25 #50 #1000 #25
mod_output = 500 # 500 #500
mod_cfl    = 10000 
time0      = time.time() 

tini = ti
n    = 1*ni

# -- Initialise dt
dt = set_dt(with_cfl); fltparam[3] = dt

if dMpi.ioproc:
    if iRestart:
        print('Resuming ... with dt:', dt)
    else:
        print('Starting ... with dt:', dt)

tdim    = val.tref 
tdimmax = 18.5 * 3600  # X hours in seconds - total run time

tdcstep = 10.*60. # X seconds - restart output frequency
tdup    = 30.     # X seconds - hmax/tmax update frequency
tstatout= 90.     # X seconds - station interpolation frequency 

# - Zero some counters
tdc     = 0.
tdmax   = 0.
tstat   = 0.
iGetRhs = True 
 
t0 = time.time()

# -- Loop until max physical time is achieved
while ( ti * tdim ) < tdimmax:

    # - Advance time and non-dimensional time
    ti  = ti  + dt
    tdc = tdc +  dt* tdim
    tdmax = tdmax +  dt* tdim  
    tstat = tstat  + dt* tdim

    # -- Get RHS forcing terms 
    if iGetRhs and not iRestart:
        dMpi.swap(q,hlo,dtree)

        # - Compute forcing terms -- ONLY DO THIS ONCE
        if 'qstored' in dtree['eqns']['qvec']['views'].keys():
            dn.dnamiF.stored(intparam, fltparam, data)
            dMpi.swap(qstored,hlo,dtree,spherical=False)

        # - Write them out
        rhsnames = []
        for key in fsskeys:
            rhsnames.append( 'rhs_' + key  )
        if dMpi.ioproc:
            print(' > Saving ...', rhsnames)
        dn.dnami_io.write_data(rhsnames,0,0,dtree,fpath=outf+'force/',fname='rhs')

        # Fill perturbation array
        add_pert(xT, yT)

        # -- Switch off this step 
        iGetRhs = False

    # - Update temporal forcing term
    update_ft(ti)

    #RK loop
    for nrk in range(1,4):
        intparam[7] = nrk
        dMpi.swap(q,hlo,dtree) 
        dn.dnamiF.time_march(intparam,fltparam,data)

        ## -- Force zero velocity
        vt[above_water] = dn.cst(0.)
        vp[above_water] = dn.cst(0.)

    #Adjusting CFL
    if np.mod(n,mod_cfl) == 0:
        dt = set_dt(with_cfl); fltparam[3] = dt

    if np.mod(n,mod_filter) == 0:
        dMpi.swapXc(q,hlo,dtree)
        dn.dnamiF.filter(1, intparam, fltparam, data)
        dMpi.swapYc(q,hlo,dtree)
        dn.dnamiF.filter(2, intparam, fltparam, data)

    # -- Update historical hmax and tmax
    if tdmax > tdup :
        idx = np.nonzero( np.abs(eta1[dom]) > (hmax[dom]+1e-8) ) # add epsilon to remove noise (e.g. filtering) from time history
        hmax[dom][idx] = np.abs(eta1[dom][idx])
        tmax[dom][idx] = ti
        tdmax = 0.

    if tdc > tdcstep : 
        # -- Write restart
        dn.dnami_io.write_restart(n,ti,0,dtree,fpath=resf)
        # -- Write hmax/tmax 
        dn.dnami_io.write_data(['hmax', 'tmax'],n,ti,dtree,fpath=outf+'historical/',fname='hmax_tmax')
        # -- Reset time
        tdc = 0.

    if tstat > tstatout: 
        # -- If proc has a station, output value 
        #
        # -- ATM
        if iAtmos:
            if len(coords_atm)> 0 :
                for key in coords_atm.keys():
                    # 
                    xp, yp = coords_atm[key]
                    # 
                    nxp = xp*grid['nxgb']/Lx    
                    nyp = yp*grid['nygb']/Ly   
                    #
                    nxpl = nxp - (dMpi.ibeg -1)
                    nypl = nyp - (dMpi.jbeg -1)
                    # Interpolation
                    outpa = MC(pa[dom], [[nxpl], [nypl]], order=1, mode='nearest')
                    oute1 = MC(eta1[dom], [[nxpl], [nypl]], order=1, mode='nearest')
                    outra = MC(rhoa[dom], [[nxpl], [nypl]], order=1, mode='nearest')
                    arr = np.array([ti, outpa.item(), oute1.item(), outra.item()])
                    # Write out
                    with open(outf+'/stations/'+key[:-4]+ '.dat', 'a+') as fh:
                        np.savetxt(fh,arr, newline=' ', delimiter=',')
                        fh.write("\n")
                    fh.close()
        # -- DART 
        if len(coords_dart)> 0:
            for key in coords_dart.keys():
                # 
                xp, yp = coords_dart[key]
                # 
                nxp = xp*grid['nxgb']/Lx     
                nyp = yp*grid['nygb']/Ly     
                #
                nxpl = nxp - (dMpi.ibeg -1)
                nypl = nyp - (dMpi.jbeg -1)
                # Interpolation
                if iAtmos:
                    oute1 = MC(eta1[dom], [[nxpl], [nypl]], order=1, mode='nearest')
                    oute0 = MC(eta0[dom], [[nxpl], [nypl]], order=1, mode='nearest')
                    outpa = MC(pa[dom], [[nxpl], [nypl]], order=1, mode='nearest')
                    outra = MC(rhoa[dom], [[nxpl], [nypl]], order=1, mode='nearest')
                    arr = np.array([ti, oute1.item(), oute0.item(), outpa.item(), outra.item()])
                else:
                    oute1 = MC(eta1[dom], [[nxpl], [nypl]], order=1, mode='nearest')
                    arr = np.array([ti, oute1])

                # Write out
                with open(outf+'/dart/'+str(key)+ '.dat', 'a+') as fh:
                    np.savetxt(fh,arr, newline=' ', delimiter=',')
                    fh.write("\n")
                fh.close()

        # -- Reset time
        tstat = 0.

    if np.mod(n,mod_output) == 0:

        # -- Timing
        time1 = time.time()
        delnt = time1 - time0 
        time0 = 1*time1

        if dMpi.ioproc:
            print('____________________________________________________________')
            print('iteration',n,' with time t =',ti,' with dt = ', dt, ' with tdim [h] = {:.0f} [min] = {:.0f}'.format(int(ti*tdim/3600),int( ((ti*tdim/3600) - int(ti*tdim/3600))*60)))
            print('____________________________________________________________')
            print('Water')
        dn.dnami_io.globalMinMax(dtree,eta1[dom],'eta1 ')
        dn.dnami_io.globalMinMax(dtree,vt[dom],  'vt   ')
        dn.dnami_io.globalMinMax(dtree,vp[dom],  'vp   ')
        if iAtmos:
            if dMpi.ioproc:
                print(' -------------------------- ' )
                print('Air')
            dn.dnami_io.globalMinMax(dtree,rhoa[dom], 'rhoa ')
            dn.dnami_io.globalMinMax(dtree,ut[dom]  , 'ut   ')
            dn.dnami_io.globalMinMax(dtree,up[dom]  , 'up   ')
            dn.dnami_io.globalMinMax(dtree,pa[dom]  , 'pa   ')
        if dMpi.ioproc:
                print(' -------------------------- ' )
                print(' Time for ', mod_output, ' iterations' )
                print(' dphystime                      :', delnt)
                print(' Ratio (phys time to real time)):', delnt/(mod_output*dt*tdim))
                print(' ft                             :', fltparam[5])
        sys.stdout.flush()

    n += 1



