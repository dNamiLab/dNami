import sys
import os
import numpy as np
import time

# =============================================================================
# input data
# =============================================================================

#nxgb, nygb, nzgb = 1024,  512,  512
#nxpr, nypr, nzpr =    4,    2,    2

#nxgb, nygb, nzgb = 128, 128, 1
#nxpr, nypr, nzpr =   4,   4, 1

#nxgb, nygb, nzgb = 4,4,1
#nxpr, nypr, nzpr = 2,1,1

#nxgb, nygb, nzgb = 2,2,4
#nxpr, nypr, nzpr = 1,1,2

#nxgb, nygb, nzgb = 64,1,1
#nxpr, nypr, nzpr = 16,1,1

nxgb, nygb, nzgb = 96,128,1 #256
nxpr, nypr, nzpr = 4,4,1

hlo  = 2
nvar = 2

iTalk = False # verbose -- print more things in console

wp = 'float64' # working precision

# =============================================================================
# some quick & dirty checks
# =============================================================================
ndim = 1
if nygb > 1:
	ndim = 2
	if nzgb >1:
		ndim = 3

# =============================================================================
# functions
# =============================================================================

def cst(x):
	if wp == 'float32': x = np.float32(x)
	if wp == 'float64': x = np.float64(x)
	return x

# =============================================================================
# START MPI (or not...)
# =============================================================================

if nxpr == 1 and nypr == 1 and nzpr == 1:
	iMpi = False # running with no mpi
else:
	iMpi = True

if iMpi:
	from mpi4py import MPI # [note] mpi_init is done in 'import mpi'
	
	# create torus
	comm_torus = MPI.COMM_WORLD.Create_cart([nxpr,nypr,nzpr],periods=True,reorder=True)
	nprocs = comm_torus.Get_size()
	procid = comm_torus.Get_rank()
	coords = comm_torus.Get_coords(procid)
	
	# assign who is in charge of console I/Os
	ioproc = False
	if procid == 0: ioproc = True # do not change ioproc rank value, 0 is assumed elsewhere
	
	# check requested number of cores makes sense
	if nxpr*nypr*nzpr != nprocs:
		if ioproc: print('[error] nxpr*nypr*nzpr does not match the number of cores requested')
		sys.exit()
	
	# check that nxgb,nygb,nzgb are multiples of nxpr,nypr,nzpr
	if (nxgb % nxpr) == 0:
		nx = int(nxgb/nxpr)
	else:
		if ioproc: print('[error] nxgb is not a multiple of nxpr')
		sys.exit()
	if (nygb % nypr) == 0:
		ny = int(nygb/nypr)
	else:
		if ioproc: print('[error] nygb is not a multiple of nypr')
		sys.exit()
	if (nzgb % nzpr) == 0:
		nz = int(nzgb/nzpr)
	else:
		if ioproc: print('[error] nzgb is not a multiple of nzpr')
		sys.exit()
	
	# create global/local index mapping
	i=1
	for itor in range(0,nxpr):
		if coords[0] == itor: ibeg = i; iend = ibeg+nx-1
		i = i+nx
	i=1
	for itor in range(0,nypr):
		if coords[1] == itor: jbeg = i; jend = jbeg+ny-1
		i = i+ny
	i=1
	for itor in range(0,nzpr):
		if coords[2] == itor: kbeg = i; kend = kbeg+nz-1
		i = i+nz
	
	# create neighbour map for swapping processes
	dummy,neighxm = comm_torus.Shift(0,-1); dummy,neighxp = comm_torus.Shift(0, 1)
	dummy,neighym = comm_torus.Shift(1,-1); dummy,neighyp = comm_torus.Shift(1, 1)
	dummy,neighzm = comm_torus.Shift(2,-1); dummy,neighzp = comm_torus.Shift(2, 1)
	
	# visualise torus?
	if iTalk:
		if ioproc:
			print('\033[1;104m'+'proc id|   with coords  | ibeg jbeg kbeg |      x-,     y-,     z-  neighbours'+'\033[0m')
			print('\033[1;104m'+'       |                | iend jend kend |      x+,     y+,     z+   (proc id)'+'\033[0m')		
			print("\033[1;30;47m{0:7d}|({1:4d},{2:4d},{3:4d})| {4:4d} {5:4d} {6:4d} | {7:7d},{8:7d},{9:7d}\033[0m".format(procid,coords[0],coords[1],coords[2],ibeg,jbeg,kbeg,neighxm,neighym,neighzm))
			print("\033[1;30;47m       |                | {0:4d} {1:4d} {2:4d} | {3:7d},{4:7d},{5:7d}\033[0m".format(iend,jend,kend,neighxp,neighyp,neighzp))
			for pro in range(1,nprocs):	
				data = np.empty(16,dtype='i')
				comm_torus.Recv([data,MPI.INT],source=pro,tag=pro)
				if (pro % 2) == 0:
					print("\033[1;30;47m{0:7d}|({1:4d},{2:4d},{3:4d})| {4:4d} {5:4d} {6:4d} | {7:7d},{8:7d},{9:7d}\033[0m".format(data[0],data[1],data[2],data[3],data[4],data[6],data[8],data[10],data[12],data[14]))
					print("\033[1;30;47m       |                | {0:4d} {1:4d} {2:4d} | {3:7d},{4:7d},{5:7d}\033[0m".format(data[5],data[7],data[9],data[11],data[13],data[15]))
				else:
					print("\033[1;37;40m{0:7d}|({1:4d},{2:4d},{3:4d})| {4:4d} {5:4d} {6:4d} | {7:7d},{8:7d},{9:7d}\033[0m".format(data[0],data[1],data[2],data[3],data[4],data[6],data[8],data[10],data[12],data[14]))
					print("\033[1;37;40m       |                | {0:4d} {1:4d} {2:4d} | {3:7d},{4:7d},{5:7d}\033[0m".format(data[5],data[7],data[9],data[11],data[13],data[15]))
		else:
			data = np.empty(16,dtype='i')
			data[:] = [procid,coords[0],coords[1],coords[2],ibeg,iend,jbeg,jend,kbeg,kend,neighxm,neighxp,neighym,neighyp,neighzm,neighzp]
			comm_torus.Send([data,MPI.INT],dest=0,tag=procid)
else:
	nx = nxgb; ny = nygb; nz = nzgb; 
	ibeg = 1; iend = nx; jbeg = 1; jend = ny; kbeg = 1; kend = nz
	ioproc = True

# =============================================================================
# CREATE DOMAIN
# =============================================================================
#                                    pt 1     2          n-1    n     period
#                                    0                             L  |
# full domain is [0,L] but carefull: |--o--|--o--| ... |--o--|--o--|--V
#                                        \___________________________/ 
Lx  = cst(1.)
dx  = Lx/cst(nxgb)
x   = np.arange(dx/cst(2.),Lx,dx,dtype=wp)
x = x**2
if nygb > 1:
	Ly  = cst(3.)
	dy  = Ly/cst(nygb)
	y   = np.arange(dy/cst(2.),Ly,dy,dtype=wp)
else:
	Ly = cst(0.); y = []
if nzgb > 1:
	Lz  = cst(5.)
	dz  = Lz/cst(nzgb)
	z   = np.arange(dz/cst(2.),Lz,dz,dtype=wp)
else:
	Lz = cst(0.); z = []
# write grid to file
if ioproc:
	head = np.empty(6,dtype=wp)
	head[:] = [nxgb,nygb,nzgb,Lx,Ly,Lz]
	with open('./out/axes.bin',"wb") as fh:
		np.concatenate((head,x,y,z)).tofile(fh)
		fh.closed

# global             iend ! ibeg                    iend  !  ibeg
# pt#          n-1    n   !   1     2          n-1    n   !   1     2
#        ... |--o--|--o--|!|--o--|--o--| ... |--o--|--o--|!|--o--|--o--| ...
# w/ hlo     [<------------------------------------------------------->]
# loc py ind    0     1      hlo                   hlo+n-1      n+2*hlo-1
# glo py ind             ibeg+hlo-1               iend+hlo-1

# dummy field
if ndim == 3:
	xloc = x[ibeg-1:iend] # without halos
	yloc = y[jbeg-1:jend]
	zloc = z[kbeg-1:kend]
	xx,yy,zz = np.meshgrid(xloc,yloc,zloc,sparse=False,indexing='ij')
	f = np.zeros((nx+2*hlo,ny+2*hlo,nz+2*hlo,nvar),dtype=wp)
	f[hlo:hlo+nx,hlo:hlo+ny,hlo:hlo+nz,0] = np.sin(cst(2.)*np.pi*xx)+np.sin(cst(2.)*np.pi*yy)+np.sin(cst(2.)*np.pi*zz)
	f[hlo:hlo+nx,hlo:hlo+ny,hlo:hlo+nz,1] = np.cos(cst(2.)*np.pi*xx)+np.cos(cst(2.)*np.pi*yy)+np.cos(cst(2.)*np.pi*zz)
elif ndim == 2:
	xloc = x[ibeg-1:iend] # without halos
	yloc = y[jbeg-1:jend]
	xx,yy = np.meshgrid(xloc,yloc,sparse=False,indexing='ij')
	f = np.zeros((nx+2*hlo,ny+2*hlo,nvar),dtype=wp)
	f[hlo:hlo+nx,hlo:hlo+ny,0] = np.sin(cst(2.)*np.pi*xx)+np.sin(cst(2.)*np.pi*yy)
	f[hlo:hlo+nx,hlo:hlo+ny,1] = np.cos(cst(2.)*np.pi*xx)+np.cos(cst(2.)*np.pi*yy)
else:
	xloc = x[ibeg-1:iend] # without halos
	f = np.zeros((nx+2*hlo,nvar),dtype=wp)
	f[hlo:hlo+nx,0] = np.sin(cst(2.)*np.pi*xloc)

# =============================================================================
# WELCOME MESSAGES
# =============================================================================
if ioproc:
	print('\033[1;40;94m'+'  /\     /\_   '+'\033[1;40;37m'+'  ||\  |  _       o'+'\033[1;40;31m'+'   /\      __         '+'\033[0m')	
	print('\033[1;40;94m'+'\/  \  _/   \  '+'\033[1;40;37m'+' _|| \ | /_\ |\/| |'+'\033[1;40;31m'+'  /  \    /  \        '+'\033[0m')
	print('\033[1;40;94m'+'     \/      \ '+'\033[1;40;37m'+'(_||  \|/   \|  | |'+'\033[1;40;31m'+' /    \__/    \__ v0.0'+'\033[0m')
	print('\033[1m'+'Geometry:'+'\033[0m')	
	print('-- grid size (nxgb,nygb,nzgb):',nxgb,nygb,nzgb)
	if iMpi:
		print('-- running with mpi')
		print('-- number of requested processes:',nxpr,'x',nypr,'x',nzpr,'=',nprocs)
	else:
		print('-- running with NO mpi')
	print('-- local grid size (nx,ny,nz):',nx,ny,nz)
	print('-- xmin/max,dx:',x[0],x[-1],dx)
	if nygb > 1: print('-- ymin/max,dy:',y[0],y[-1],dy)
	if nzgb > 1: print('-- zmin/max,dz:',z[0],z[-1],dz)
	sys.stdout.flush()

def swap1d_x(a,nv,nh):
	# 'a' is an array in x with 'nv' fields and an halo size 'nh'	
	# positive x-dir
	buf = np.empty(nh,dtype=wp)
	buf = a[nx:nh+nx,0:nv].copy().reshape((nh*nv,1))
	if iMpi: comm_torus.Sendrecv_replace(buf,neighxp,sendtag=0,source=neighxm,recvtag=0,status=None)
	a[0:nh,0:nv] = buf.reshape((nh,nv)).copy()
	# negative x-dir
	buf = np.empty(nh,dtype=wp)
	buf = a[nh:2*nh,0:nv].copy().reshape((nh*nv,1))
	if iMpi: comm_torus.Sendrecv_replace(buf,neighxm,sendtag=0,source=neighxp,recvtag=0,status=None)
	a[nx+nh:nx+2*nh,0:nv] = buf.reshape((nh,nv)).copy()
	return

def swap2d_x(a,nj,nv,nh):
	# 'a' is an array with 'nj' values in j and 'nv' fields with halo size 'nh'	
	# positive x-dir	
	buf = np.empty(nh*nj*nv,dtype=wp)
	buf = a[nx:nh+nx,nh:nh+nj,0:nv].copy().reshape((nh*nj*nv,1))
	if iMpi: comm_torus.Sendrecv_replace(buf,neighxp,sendtag=0,source=neighxm,recvtag=0,status=None)
	a[0:nh,nh:nh+nj,0:nv] = buf.reshape((nh,nj,nv)).copy()
	# negative x-dir
	buf = np.empty(nh*nj*nv,dtype=wp)
	buf = a[nh:2*nh,nh:nh+nj,0:nv].copy().reshape((nh*nj*nv,1))
	if iMpi: comm_torus.Sendrecv_replace(buf,neighxm,sendtag=0,source=neighxp,recvtag=0,status=None)
	a[nx+nh:nx+2*nh,nh:nh+nj,0:nv] = buf.reshape((nh,nj,nv)).copy()
	return

def swap3d_x(a,nj,nk,nv,nh):
	# 'a' is an array with 'nj,nk' values in j,k and 'nv' fields with halo size 'nh'	
	# positive x-dir	
	buf = np.empty(nh*nj*nk*nv,dtype=wp)
	buf = a[nx:nh+nx,nh:nh+nj,nh:nh+nk,0:nv].copy().reshape((nh*nj*nk*nv,1))
	if iMpi: comm_torus.Sendrecv_replace(buf,neighxp,sendtag=0,source=neighxm,recvtag=0,status=None)
	a[0:nh,nh:nh+nj,nh:nh+nk,0:nv] = buf.reshape((nh,nj,nk,nv)).copy()
	# negative x-dir
	buf = np.empty(nh*nj*nk*nv,dtype=wp)
	buf = a[nh:2*nh,nh:nh+nj,nh:nh+nk,0:nv].copy().reshape((nh*nj*nk*nv,1))
	if iMpi: comm_torus.Sendrecv_replace(buf,neighxm,sendtag=0,source=neighxp,recvtag=0,status=None)
	a[nx+nh:nx+2*nh,nh:nh+nj,nh:nh+nk,0:nv] = buf.reshape((nh,nj,nk,nv)).copy()
	return

def swap2d_y(a,ni,nv,nh):
	# 'a' is an array with 'ni' values in i and 'nv' fields with halo size 'nh'	
	# positive y-dir	
	buf = np.empty(nh*ni*nv,dtype=wp)
	buf = a[nh:nh+ni,ny:nh+ny,0:nv].copy().reshape((nh*ni*nv,1))
	if iMpi: comm_torus.Sendrecv_replace(buf,neighyp,sendtag=0,source=neighym,recvtag=0,status=None)
	a[nh:nh+ni,0:nh,0:nv] = buf.reshape((ni,nh,nv)).copy()
	# negative y-dir
	buf = np.empty(nh*ni*nv,dtype=wp)
	buf = a[nh:nh+ni,nh:2*nh,0:nv].copy().reshape((nh*ni*nv,1))
	if iMpi: comm_torus.Sendrecv_replace(buf,neighym,sendtag=0,source=neighyp,recvtag=0,status=None)
	a[nh:nh+ni,ny+nh:ny+2*nh,0:nv] = buf.reshape((ni,nh,nv)).copy()
	return

def swap3d_y(a,ni,nk,nv,nh):
	# 'a' is an array with 'ni,nk' values in i,k and 'nv' fields with halo size 'nh'	
	# positive y-dir
	buf = np.empty(nh*ni*nk*nv,dtype=wp)
	buf = a[nh:nh+ni,ny:nh+ny,nh:nh+nk,0:nv].copy().reshape((nh*ni*nk*nv,1))
	if iMpi: comm_torus.Sendrecv_replace(buf,neighyp,sendtag=0,source=neighym,recvtag=0,status=None)
	a[nh:nh+ni,0:nh,nh:nh+nk,0:nv] = buf.reshape((ni,nh,nk,nv)).copy()
	# negative y-dir
	buf = np.empty(nh*ni*nk*nv,dtype=wp)
	buf = a[nh:nh+ni,nh:2*nh,nh:nh+nk,0:nv].copy().reshape((nh*ni*nk*nv,1))
	if iMpi: comm_torus.Sendrecv_replace(buf,neighym,sendtag=0,source=neighyp,recvtag=0,status=None)
	a[nh:nh+ni,ny+nh:ny+2*nh,nh:nh+nk,0:nv] = buf.reshape((ni,nh,nk,nv)).copy()
	return

def swap3d_z(a,ni,nj,nv,nh):
	# 'a' is an array with 'ni,nj' values in i,j and 'nv' fields with halo size 'nh'	
	# positive z-dir	
	buf = np.empty(nh*ni*nj*nv,dtype=wp)
	buf = a[nh:nh+ni,nh:nh+nj,nz:nh+nz,0:nv].copy().reshape((nh*ni*nj*nv,1))
	if iMpi: comm_torus.Sendrecv_replace(buf,neighzp,sendtag=0,source=neighzm,recvtag=0,status=None)
	a[nh:nh+ni,nh:nh+nj,0:nh,0:nv] = buf.reshape((ni,nj,nh,nv)).copy()
	# negative z-dir
	buf = np.empty(nh*ni*nj*nv,dtype=wp)
	buf = a[nh:nh+ni,nh:nh+nj,nh:2*nh,0:nv].copy().reshape((nh*ni*nj*nv,1))
	if iMpi: comm_torus.Sendrecv_replace(buf,neighzm,sendtag=0,source=neighzp,recvtag=0,status=None)
	a[nh:nh+ni,nh:nh+nj,nz+nh:nz+2*nh,0:nv] = buf.reshape((ni,nj,nh,nv)).copy()
	return

if iMpi:
	if wp == 'float64': MPIWP = MPI.DOUBLE
	if wp == 'float32': MPIWP = MPI.FLOAT

def write_restart(n,t,flag):
	# write bare minimum to restart a run from the instant of the function call
	if flag == 0:	
		fname = './restarts/restart_' + str(n).zfill(8)
	else:
		fname = './out/liv/restart.bin'
	# header
	headsize = 7
	head    = np.empty(headsize,dtype=wp)
	head[:] = [headsize,n,t,nxgb,nygb,nzgb,nvar] # implicit type conversion
	# body
	dat = np.empty((nx,ny,nz,nvar),dtype=wp)
	if ndim == 3:
		dat = f[hlo:hlo+nx,hlo:hlo+ny,hlo:hlo+nz,np.newaxis].copy()
	elif ndim == 2:
		dat = f[hlo:hlo+nx,hlo:hlo+ny,np.newaxis,np.newaxis].copy()
	else:
		dat = f[hlo:hlo+nx,np.newaxis,np.newaxis,np.newaxis].copy()	
	if iMpi:
		header = MPIWP.Create_contiguous(headsize)
		header.Commit()	
		fh = MPI.File.Open(comm_torus,fname,MPI.MODE_WRONLY|MPI.MODE_CREATE)
		fh.Set_view(0,MPIWP,header)
		if ioproc: fh.Write_at(0,head)
		header.Free()	
		subarray = MPIWP.Create_subarray((nxgb,nygb,nzgb,nvar),(nx,ny,nz,nvar),(ibeg-1,jbeg-1,kbeg-1,0))
		subarray.Commit()
		disp = MPIWP.Get_size()*headsize
		fh.Set_view(disp,MPIWP,subarray)
		fh.Write_all(dat)
		subarray.Free()
		fh.Close()
	else:
		with open(fname,"wb") as fh:
			np.concatenate((head,np.reshape(dat,(nx*ny*nz*nvar)))).tofile(fh)
		fh.closed
	return

def read_restart():
	# read in restart
	fname = 'restart.bin'
	headsize = 7
	if iMpi:
		head = np.empty(headsize,dtype=wp)
		dat  = np.empty((nx,ny,nz,nvar),dtype=wp)
		header = MPIWP.Create_contiguous(headsize)
		header.Commit()	
		fh = MPI.File.Open(comm_torus,fname,MPI.MODE_RDONLY)
		fh.Set_view(0,MPIWP,header)
		fh.Read_at(0,head)
		header.Free()	
		subarray = MPIWP.Create_subarray((nxgb,nygb,nzgb,nvar),(nx,ny,nz,nvar),(ibeg-1,jbeg-1,kbeg-1,0))
		subarray.Commit()
		disp = MPIWP.Get_size()*headsize
		fh.Set_view(disp,MPIWP,subarray)
		fh.Read_all(dat)
		subarray.Free()
		fh.Close()
	else:
		with open(fname,"rb") as fh:
			head = np.fromfile(fh,dtype=wp,count=headsize)
			dat  = np.fromfile(fh,dtype=wp,count=nx*ny*nz*nvar)
			dat  = np.reshape(dat,(nx,ny,nz,nvar))
		fh.closed
	# run some checks & print some info to terminal:
	if ioproc:
		print('\033[1;32m'+'================================'+'\033[0m')
		print('\033[1;32m'+'RESTARTING FROM RESTART.BIN FILE'+'\033[0m')
		print('\033[1;32m'+'================================'+'\033[0m')
		print('\033[1m'+'From header:'+'\033[0m')
		print('n,time:',int(head[1]),',',head[2])
		print('nxgb,nygb,nzgb,nvar:',int(head[3]),',',int(head[4]),',',int(head[5]),',',int(head[6]))
		if nxgb != head[3]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nxgb does not match that of restart.bin'); sys.exit()
		if nygb != head[4]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nygb does not match that of restart.bin'); sys.exit()
		if nzgb != head[5]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nzgb does not match that of restart.bin'); sys.exit()
		if nvar != head[6]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nvar does not match that of restart.bin'); sys.exit()
		print('We are good to go!')
	# set iteration number and time:
	n  = int(head[1])
	ti = head[2]
	# set fields:
	if ndim == 3:
		f[hlo:hlo+nx,hlo:hlo+ny,hlo:hlo+nz,np.newaxis] = dat.copy()
	elif ndim == 2:
		f[hlo:hlo+nx,hlo:hlo+ny,np.newaxis,np.newaxis] = dat.copy()
	else:
		f[hlo:hlo+nx,np.newaxis,np.newaxis,np.newaxis] = dat.copy()
	return

if ndim == 3:
	if nvar > 1:
		swap3d_x(f[:,:,:,:],ny,nz,nvar,hlo)
		swap3d_y(f[:,:,:,:],nx,nz,nvar,hlo)
		swap3d_z(f[:,:,:,:],nx,ny,nvar,hlo)
	else:
		swap3d_x(f[:,:,:,np.newaxis],ny,nz,1,hlo)
		swap3d_y(f[:,:,:,np.newaxis],nx,nz,1,hlo)
		swap3d_z(f[:,:,:,np.newaxis],nx,ny,1,hlo)
elif ndim == 2:
	swap2d_x(f[:,:,:],ny,1,hlo)
	swap2d_y(f[:,:,:],nx,1,hlo)
else:
	swap1d_x(f[:,:],1,hlo)
#print(procid,f)
#print(f)
#write_restart()
#read_restart()

dt = cst(0.1)
nn = np.arange(1,1001,1)
ti = cst(0.)
n = 0

#liveView = True
#exec("liveViewNd = liveView" + str(ndim) + 'd') # automatically calls the 1d, 2d or 3d viewer

# time loop
def tloop():
	global n, f, ti
	for n in nn:
		ti = ti + dt
		if ndim == 1: f[hlo:hlo+nx,                      0] = np.sin(cst(2.)*np.pi*(xloc-ti))
		if ndim == 2: f[hlo:hlo+nx,hlo:hlo+ny,           0] = np.sin(cst(2.)*np.pi*(xx-ti))+np.sin(cst(2.)*np.pi*(yy-ti))
		if ndim == 3: f[hlo:hlo+nx,hlo:hlo+ny,hlo:hlo+nz,0] = np.sin(cst(2.)*np.pi*(xx-ti))+np.sin(cst(2.)*np.pi*(yy-ti))+np.sin(cst(2.)*np.pi*(zz-ti))
		if ioproc and (n%25 == 0): print('n',n,'time',ti); sys.stdout.flush()
		# --- for live view ---
		liveViewNow = False
		if ioproc:
			if os.path.exists('./out/liv/snap'):
				snap_time = os.path.getmtime('./out/liv/snap')
				if os.path.exists('./out/liv/restart.bin'):
					data_time = os.path.getmtime('./out/liv/restart.bin')
					if snap_time > data_time: liveViewNow = True
				else:
					liveViewNow = True
		if iMpi: liveViewNow = comm_torus.bcast(liveViewNow,root=0)
		if liveViewNow: write_restart(n,ti,1)
		# --- end of for live view ---
	return

tloop()
# tell the live viewer that the time loop has ended, and save last flow state
os.system('touch ./out/liv/stop')
write_restart(n,ti,1)

'''
if liveView:
	import matplotlib.pyplot as plt
	from matplotlib.widgets import TextBox
	import threading
	sys.setrecursionlimit(100000000) #;print(sys.getrecursionlimit())
	tloop_thread = threading.Thread(target=tloop, args=())
	vu = liveViewNd()
else:
	tloop()
'''
