import math as m
import dnami as dn
from dnami import np, sys # re-import non dnami modules
import os

dtree = dn.create_tree()

alpha = dn.cst(2.5)
iVDW  = False

Ma    = dn.cst(0.1)
Re    = dn.cst(1600.0)
Pr    = dn.cst(0.71)
gamma = dn.cst(1.4)
dimPcr = False
if dimPcr:
	Cv = dn.cst(1.0/(gamma-1.0))
else:
	Cv = dn.cst(1.0/(Ma**2*gamma*(gamma-1.0)))

dtree['eqns']['coeff'][0][1] = dn.cst(1.0/Re)
dtree['eqns']['coeff'][1][1] = dn.cst(1.0/( (gamma-1.0)*Ma**2*Re*Pr ))
dtree['eqns']['coeff'][2][1] = dn.cst(gamma-1.)
dtree['eqns']['coeff'][3][1] = dn.cst(1.0/Cv)


numerics = dtree['num']
grid     = dtree['grid']['size']
geom     = dtree['grid']['geom']
mpi      = dtree['mpi']['split']

grid['nxgb'] = 360
grid['nygb'] = 360
grid['nzgb'] = 360


geom['Lx'] = dn.cst(2.0*m.pi) 
geom['Ly'] = dn.cst(2.0*m.pi)
geom['Lz'] = dn.cst(2.0*m.pi)


mpi['nxpr'] = 1
mpi['nypr'] = 5
mpi['nzpr'] = 8

dtree = dn.start_mpi(dtree)

dMpi = dtree['mpi']['dMpi']

nx  = dMpi.nx
ny  = dMpi.ny
nz  = dMpi.nz


dtree = dn.create_grid(dtree)
#dn.dnami_io.write_grid(dtree)
dn.dnami_io.hello_world(dtree)

xloc, yloc, zloc  = geom['xloc'], geom['yloc'], geom['zloc']
Lx  , Ly  , Lz    = geom['Lx']  , geom['Ly']  , geom['Lz']
dx  , dy  , dz    = geom['dx']  , geom['dy']  , geom['dz']

#=====================================================================================================================
#
#                                                ALLOCATE DATA
#
#=====================================================================================================================

numerics['tint']['tstep'] = dn.cst(0.0000001)
dt = numerics['tint']['tstep']

numerics['filtr']['eps'] = dn.cst(0.1)


dtree['libs']['cache blocking'] = [1024,2,6]


dtree = dn.allocate(dtree)

rh = dtree['eqns']['qvec']['views']['rho']
ux = dtree['eqns']['qvec']['views']['u']
uy = dtree['eqns']['qvec']['views']['v']
uz = dtree['eqns']['qvec']['views']['w']
et = dtree['eqns']['qvec']['views']['et']

q  = dtree['eqns']['qvec']['views']['q'] 

if 'qstored' in dtree['eqns']['qvec']['views'].keys():
	qstored  = dtree['eqns']['qvec']['views']['qstored'] 

#=====================================================================================================================
#
#                                                Initialisation 
#
#=====================================================================================================================

hlo  = numerics['hlo']

#init thermo
T0   = dn.cst(1.0)
P0   = dn.cst(1.0)/(Ma**2*gamma)
Rho0 = dn.cst(1.0)#P0/T0*Ma**2*gamma


rh[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo] = Rho0

ux[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo] = (np.sin(xloc[:         , np.newaxis, np.newaxis])
									   *np.cos(yloc[np.newaxis, :         , np.newaxis])
									   *np.cos(zloc[np.newaxis, np.newaxis, :         ]))

uy[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo] = (-np.cos(xloc[:         , np.newaxis, np.newaxis])
										*np.sin(yloc[np.newaxis, :         , np.newaxis])
										*np.cos(zloc[np.newaxis, np.newaxis, :         ]))

uz[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo] = dn.cst(0.0)

p =(P0 + dn.cst(Rho0/16.0)*(np.cos( dn.cst(2.0) * (xloc[:         , np.newaxis, np.newaxis]) )
						   +np.cos( dn.cst(2.0) * (yloc[np.newaxis, :         , np.newaxis]) ))
						  *(np.cos( dn.cst(2.0) * (zloc[np.newaxis, np.newaxis, :         ]) ) + dn.cst(2.0)))


et[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo] = (  p/rh[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo]*dn.cst(1./(gamma-1.)) 
						 + dn.cst(0.5)*(ux[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo]*ux[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo]
						 			   +uy[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo]*uy[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo]
						 			   +uz[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo]*uz[hlo:nx+hlo,hlo:ny+hlo,hlo:nz+hlo]))


def sound_speed():
	e = et - .5*(ux*ux+uy*uy)
	if iVDW:
		T = (1./alpha)*( e*Ma*Ma +1.125*rh)
		c = np.sqrt( T*(1.+1./alpha)/((1.-rh/3.)**2) - 2.*(9./8.)*rh )/Ma
	else:
		T = (1./alpha)*( e*Ma*Ma )
		c = np.sqrt( T*(1.+1./alpha) )/Ma
	return c	

#dn.dnami_io.write_restart(0,dn.cst(0.),0,dtree)

dMpi.swap(q,hlo,dtree)

if 'qstored' in dtree['eqns']['qvec']['views'].keys():
	#dn.dnamiF.stored(intparam,fltparam,data)	
	dn.stored(intparam,fltparam,data)	
	dMpi.swap(qstored,hlo,dtree)

#=====================================================================================================================
#
#                                                dNami time loop 
#
#=====================================================================================================================


ti = dn.cst(0.0)

intparam,fltparam,data = (dtree['libs']['fort']['integers'],
						  dtree['libs']['fort']['floats'],
						  dtree['libs']['fort']['data'])


# dtree = dn.dnami_io.read_restart(dtree)

nitmax = 40

mod_filter = 1
mod_output = 1000
mod_io     = 1000

from timeit import default_timer as timer

t0    = timer()

tcom = 0.0
tcmp = 0.0
tflt = 0.0

try:
	NOMP=int(os.environ['OMP_NUM_THREADS'])	
except:
	NOMP=1


for n in range(1,nitmax+1):
	ti = ti + dt
	# t00 = timer()
	for nrk in range(1,4):
		intparam[7] = nrk

		# dn.dnamiF.bc(0,intparam,fltparam,data)
		tbeg = timer()
		dMpi.swap(	q,hlo,dtree)
		tend = timer()
		tcom = tcom + tend-tbeg

		# tcmp_0 = timer()
		# dMpi.swap(	qstored,hlo,dtree)
		if 'qstored' in dtree['eqns']['qvec']['views'].keys():
			#dn.dnamiF.stored(intparam,fltparam,data)
			dn.stored(intparam,fltparam,data)
			dMpi.swap(	qstored,hlo,dtree)

		# tmarkcmp = timer()	
		tbeg = timer()
		#dn.dnamiF.time_march(intparam,fltparam,data)	
		dn.time_march(intparam,fltparam,data)	
		tend = timer()	
		tcmp = tcmp + tend-tbeg

	if np.mod(n,mod_filter) == 0:

		# tmark = timer()
		tbeg = timer()
		dMpi.swapXc(q,hlo,dtree)
		tend = timer()
		tcom = tcom + tend-tbeg

		# dn.dnamiF.bc(1,intparam,fltparam,data)
		# tmarkflt = timer()
		tbeg = timer()
		#dn.dnamiF.filter(1,intparam,fltparam,data)
		dn.filter(1,intparam,fltparam,data)
		tend = timer()
		tflt = tflt + tend-tbeg

		# tmark = timer()
		tbeg = timer()
		dMpi.swapYc(q,hlo,dtree)
		tend = timer()
		tcom = tcom + tend-tbeg

		# dn.dnamiF.bc(2,intparam,fltparam,data)
		# tmarkflt = timer()
		tbeg = timer()
		#dn.dnamiF.filter(2,intparam,fltparam,data)
		dn.filter(2,intparam,fltparam,data)
		tend = timer()
		tflt = tflt + tend-tbeg

		# tmark = timer()
		tbeg = timer()
		dMpi.swapZc(q,hlo,dtree)
		tend = timer()
		tcom = tcom + tend-tbeg

		# dn.dnamiF.bc(3,intparam,fltparam,data)
		# # tcmp_0 = timer()
		# dn.dnamiF.bc(3,intparam,fltparam,data)	
		# tmarkflt = timer()
		tbeg = timer()
		#dn.dnamiF.filter(3,intparam,fltparam,data)
		dn.filter(3,intparam,fltparam,data)
		tend = timer()
		tflt = tflt + tend-tbeg

print("Tend* :",ti*2.0*m.pi)

t1 = timer()



print('____________________________________________________________')
print('iteration',n,' with time t =',ti)
# sys.stdout.flush()
dn.dnami_io.globalMinMax(dtree,rh,'r')
dn.dnami_io.globalMinMax(dtree,ux,'u')
dn.dnami_io.globalMinMax(dtree,uy,'v')
dn.dnami_io.globalMinMax(dtree,et,'et')
if dMpi.ioproc:
	print('convective CFL numbers')
	sys.stdout.flush()
cfl = dt*np.abs(ux[hlo:nx+hlo,hlo:ny+hlo])/dx
dn.dnami_io.globalMax(dtree,cfl,'cfl-x')
cfl = dt*np.abs(uy[hlo:nx+hlo,hlo:ny+hlo])/dy
dn.dnami_io.globalMax(dtree,cfl,'cfl-y')
if dMpi.ioproc:
	print('acoustic CFL numbers')
	sys.stdout.flush()
c = sound_speed()

cfl = dt*(np.abs(ux)+c)/dx
dn.dnami_io.globalMax(dtree,cfl,'cfl-x')
cfl = dt*(np.abs(uy)+c)/dy
dn.dnami_io.globalMax(dtree,cfl,'cfl-y')
U = np.sqrt(ux*ux+uy*uy); locMach = U/c
dn.dnami_io.globalMinMax(dtree,locMach,'M')
if dMpi.ioproc:
	print('diffusive CFL numbers')
	print('dt/(Re*dx*dx)',dt/(Re*dx*dx))
	sys.stdout.flush()
if dMpi.ioproc: sys.stdout.flush()	

# NO IOS	if np.mod(n,mod_io) == 0:	
# NO IOS		dMpi.swap(	q,hlo,dtree)
# NO IOS		if 'qstored' in dtree['eqns']['qvec']['views'].keys():
# NO IOS			dn.dnamiF.stored(intparam,fltparam,data,1)
# NO IOS		field = ['u','rho','Omega']
# NO IOS		write_data(field,path,n,ti,dtree)



dn.dnami_io.write_restart(n,ti,0,dtree)



print('Nb. thread is: ',NOMP)

print("Time total  ",(t1-t0)/(grid['nxgb']*grid['nygb']*grid['nzgb']*3*nitmax)*NOMP*mpi['nxpr']*mpi['nypr']*mpi['nzpr'],t1-t0)
# print("Time compute",(tcmp)/(grid['nxgb']*grid['nygb']*grid['nzgb']*3*nitmax)*NOMP*mpi['nxpr']*mpi['nypr']*mpi['nzpr'],tcmp)

tfile = open('timeblck.dat','a')
tfile.write(str(dtree['libs']['cache blocking'][0])+'  '+
			str(dtree['libs']['cache blocking'][1])+'  '+
			str(dtree['libs']['cache blocking'][2])+'  '+
			str(t1-t0)+'\n')

if dMpi.ioproc:
	time_file = open('time_details.dat','a')
	time_file.write('{:5}   {:5}   {:5}   {:5}   {:5}   {:5}   {:5}   {:5}   {:5}   {:3.5f}   {:3.5f}   {:2.8f}   '.format(360,360,360,1,5,8,1024,2,6,tcom,tcmp,tflt)+'\n')
