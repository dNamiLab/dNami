import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
import psutil
sys.setrecursionlimit(10000000)

# hard-coded user inputs:
wp = 'float64' # working precision
nvar = 4
unif = True  # uniform mesh?

# do some cleaning if needed
if os.path.exists('snap'): os.system('rm snap')
if os.path.exists('stop'): os.system('rm stop')
#if os.path.exists('restart.bin'): os.system('rm restart.bin')

# load grid
with open('../axes.bin',"rb") as fh:
	head = np.fromfile(fh,dtype=wp,count=6)
	nxgb, nygb, nzgb = int(head[0]), int(head[1]), int(head[2])
	Lx, Ly, Lz = head[3], head[4], head[5]; del head
	ndim = 2
	if nygb == 1: ndim = 1
	if nzgb >  1: ndim = 3
	x = np.fromfile(fh,dtype=wp,count=nxgb); dx = x[1]-x[0]
	if ndim>1: y = np.fromfile(fh,dtype=wp,count=nygb); dy = y[1]-y[0]
	if ndim>2: z = np.fromfile(fh,dtype=wp,count=nzgb); dz = z[1]-z[0]
	fh.closed
nx = nxgb; ny = nygb; nz = nzgb

if not unif and ndim>1: xmesh,ymesh = np.meshgrid(np.append(x-.5*dx,Lx)/Lx,np.append(y-.5*dy,Ly)/Ly,sparse=False,indexing='ij')

# ///////////////////////////////////////////////////////////////////////////// functions
def safe_to_read(fname):
	# check if a file is safe to read (to avoid race condition)
	for pids in psutil.process_iter():
		try:
			for fhandle in pids.open_files():
				if fname == fhandle.path:
					return True
		except Exception:
			pass
	return False

def read_restart():
	notsafe = True
	while notsafe:
		notsafe = safe_to_read('restart.bin')
	with open('restart.bin',"rb") as fh:
		headsize = int(np.fromfile(fh,dtype=wp,count=1))
		head = np.fromfile(fh,dtype=wp,count=headsize-1)
		n = int(head[0])
		t = head[1]
		nx,ny,nz,nv = int(head[2]),int(head[3]),int(head[4]),int(head[5])
		if nx != nxgb or ny != nygb or nz != nzgb:
			print('[error] grid file needs to be re-loaded (not compatible with restart.bin)')
			sys.exit()
		if nv != nvar:
			print('[error] nvar given in header (',nvar,'), does not match that of restart.bin (',nv,')')
			sys.exit()
		f = np.fromfile(fh,dtype=wp,count=nx*ny*nz*nv)
		if ndim == 3:
			f = np.reshape(f,(nx,ny,nz,nvar))
		elif ndim == 2:
			f = np.reshape(f,(nx,ny   ,nvar))
		else:
			f = np.reshape(f,(nx      ,nvar))
	fh.closed
	return n,t,f

def liveViewCollect(pcut,ncut,ivar):
	# gather field to plot in live view
	if ndim == 3:
		if pcut == 'z':
			dat = f[0:nx,0:ny,ncut-1,ivar-1]
		elif pcut == 'y':
			dat = f[0:nx,ncut-1,0:nz,ivar-1]
		else:
			dat = f[ncut-1,0:ny,0:nz,ivar-1].T
	elif ndim == 2:
		dat = f[0:nx,0:ny,ivar-1]
	else:
		dat = f[0:nx,ivar-1]
	if ndim>1 and unif: dat = dat.T
	return dat

# ///////////////////////////////////////////////////////////////////////////// viewer class
class liveView:
	# main viewer class
	def __init__(self):
		print('[info] starting live viewer')
		self.pl = None; self.ix = None; self.qt = None   # plane normal & index, variable to plot
		self.fi = None; self.ax = None; self.im = None   # figure, axes and data handles
		self.ti = None; self.xl = None; self.yl = None   # title, x/y-label handles
		self.tb = None; self.tq = None                   # plane and variable text box handles
		self.lo = None; self.cutloc = None               # slice locator handle and actual location
		self.gb = None; self.imax   = None               # index range indicator handle and max val		
		self.co = None; self.tc = None; self.isoc = None # iso-contour obj, text box handle & contour value
		self.au = None; self.auto = False                # text handle and switch for auto refresh
		
		self.rcur = 0       # limiter on recursive calls (to deal with it before sys does)
		self.wait = None    # text handle for "wait" message
		self.newCut = False # flag to handle cut changes on non-uniform meshes
		
		self.setDefault() # set view to default values
		
		dat = liveViewCollect(self.pl,self.ix,self.qt)
		
		self.fi = plt.figure(figsize=(9,10))
		self.fi.show()
		self.ax = self.fi.add_axes([.06,.03,.88,0.88*0.9])
		if ndim > 1:
			self.ax.set_xlim(0.,1.); self.ax.set_ylim(0.,1.)
			if unif:
				self.im = self.ax.imshow(dat,origin='lower',cmap="gray",zorder=1,vmin=np.amin(dat),vmax=np.amax(dat))
			else:
				self.setCmesh(dat)
		else:
			self.ax.set_xlim(0.,1.); self.ax.set_ylim(np.amin(dat),np.amax(dat))
			self.im, = self.ax.plot(x/Lx,dat,'bo-')
		self.ti = self.fi.text(0.01,0.98,'iteration 0  | time = '+str(time)+'  | min/max: '+str(np.amin(dat))+' / '+str(np.amax(dat)))
		self.xl = self.fi.text(0.9,0.02,'',ha='right',va='top')
		self.yl = self.fi.text(0.05,0.8,'',ha='right',va='top',rotation=90)
		self.cutUpdate()
		self.fi.canvas.mpl_connect('key_press_event', self.on_press)
		self.fi.text(0.055,0.95,'Active',color='k',ha='right')
		self.fi.text(0.055,0.95,' key strikes:',color='0.4',ha='left')
		if ndim == 3: 
			col = '0.0'
		else:
			col = '0.8'
		self.fi.text(0.015,0.93,'x',color=col,ha='right');self.fi.text(0.02,0.93,'| YZ plane (x slice)',color=col)
		self.fi.text(0.015,0.91,'y',color=col,ha='right');self.fi.text(0.02,0.91,'| XZ plane (y slice)',color=col)
		self.fi.text(0.015,0.89,'z',color=col,ha='right');self.fi.text(0.02,0.89,'| XY plane (z slice)',color=col)
		self.fi.text(0.215,0.93,'u',ha='right');self.fi.text(0.22,0.93,'| REFRESH view'); self.wait = self.fi.text(0.5,0.93,'',color='r')
		self.fi.text(0.215,0.91,'r',ha='right');self.fi.text(0.22,0.91,'| RESET to default')
		self.fi.text(0.215,0.89,'q',ha='right');self.fi.text(0.22,0.89,'| EXIT')
		self.fi.text(0.415,0.89,'c',ha='right');self.fi.text(0.42,0.89,'| clear (for new run)')
		self.tb = self.fi.add_axes([0.83,0.92,0.15,0.022])
		self.fi.tbox_ind = TextBox(self.tb, 'Plane index? ', initial=str(self.ix))
		self.fi.tbox_ind.on_submit(self.get_ind)
		self.lo = self.fi.text(0.83,0.905,' loc. '+str(self.cutloc),color='0.7')
		self.gb = self.fi.text(0.82,0.910,'in {1,...,'+str(self.imax)+'}',color='0.7',ha='right')
		self.tq = self.fi.add_axes([0.83,0.87,0.15,0.022])
		self.fi.tbox_qty = TextBox(self.tq, 'Variable index? ', initial=str(self.qt))
		self.fi.tbox_qty.on_submit(self.get_var)
		self.tc = self.fi.add_axes([0.83,0.96,0.15,0.022])
		self.fi.tbox_iso = TextBox(self.tc, 'Iso-contour to plot? ', initial=str(self.isoc))
		self.fi.tbox_iso.on_submit(self.get_iso)
		self.fi.text(0.015,0.85,'a',ha='right');self.fi.text(0.02,0.85,'| auto refresh -- current status: ')
		self.au = self.fi.text(0.265,0.85,'OFF',color='r')
	
	def setCmesh(self,C):
		global xmesh, ymesh
		if self.im == None:
			self.im = self.ax.pcolormesh(xmesh,ymesh,C,cmap="gray",zorder=1,vmin=np.amin(C),vmax=np.amax(C))
		else:
			if self.newCut == False:
				self.im.set_array(C.ravel())
			else:
				self.im.remove(); self.im = None; xmesh = []; ymesh = []; self.newCut = False
				if self.pl == 'x':
					xmesh,ymesh = np.meshgrid(np.append(z-.5*dz,Lz)/Lz,np.append(y-.5*dy,Ly)/Ly,sparse=False,indexing='ij')
				elif self.pl == 'y':
					xmesh,ymesh = np.meshgrid(np.append(x-.5*dx,Lx)/Lx,np.append(z-.5*dz,Lz)/Lz,sparse=False,indexing='ij')
				else:
					xmesh,ymesh = np.meshgrid(np.append(x-.5*dx,Lx)/Lx,np.append(y-.5*dy,Ly)/Ly,sparse=False,indexing='ij')
				self.setCmesh(C)
	
	def setDefault(self):
		self.pl = 'z'
		self.ix = 1
		self.qt = 1
		if ndim == 3: 
			self.cutloc = z[self.ix-1]
		else:
			self.cutloc = 0.
		self.imax = nzgb
		self.isoc = None
		
	def timeUpdate(self):
		global n,time,f
		n,time,f = read_restart()
		dat = liveViewCollect(self.pl,self.ix,self.qt)
		if ndim>1:
			if unif:
				self.im.set_data(dat)
			else:
				self.setCmesh(dat)
			self.im.set_clim(np.amin(dat),np.amax(dat))
		else:
			self.im.set_data(x/Lx,dat)
		self.ti.set_text('iteration '+str(n)+'  | time = '+str(time)+'  | min/max: '+str(np.amin(dat))+' / '+str(np.amax(dat)))
		self.getIsoLine(dat)
		self.fi.canvas.draw()
	
	def cutUpdate(self):
		dat = liveViewCollect(self.pl,self.ix,self.qt)
		if ndim>1:
			if unif:
				self.im.set_data(dat)
			else:
				self.setCmesh(dat)
			self.im.set_clim(np.amin(dat),np.amax(dat))
		else:
			self.im.set_data(x/Lx,dat)
		self.ti.set_text('iteration '+str(n)+'  | time = '+str(time)+'  | min/max: '+str(np.amin(dat))+' / '+str(np.amax(dat)))
		if self.pl == 'z':
			self.xl.set_text('$x/\ell_x$')
			if ndim>1:
				self.yl.set_text('$y/\ell_y$')
				if unif: self.im.set_extent([(x[0]-.5*dx)/Lx,(x[-1]+.5*dx)/Lx,(y[0]-.5*dy)/Ly,(y[-1]+.5*dy)/Ly])
				if ndim == 3: self.cutloc = z[self.ix-1]
			else:
				self.yl.set_text('$f$')
			self.imax = nzgb
		if self.pl == 'y':
			self.xl.set_text('$x/\ell_x$'); self.yl.set_text('$z/\ell_z$')
			if unif: self.im.set_extent([(x[0]-.5*dx)/Lx,(x[-1]+.5*dx)/Lx,(z[0]-.5*dz)/Lz,(z[-1]+.5*dz)/Lz])
			self.cutloc = y[self.ix-1]
			self.imax = nygb
		if self.pl == 'x':
			self.xl.set_text('$z/\ell_z$'); self.yl.set_text('$y/\ell_y$')
			if unif: self.im.set_extent([(z[0]-.5*dz)/Lz,(z[-1]+.5*dz)/Lz,(y[0]-.5*dy)/Ly,(y[-1]+.5*dy)/Ly])
			self.cutloc = x[self.ix-1]
			self.imax = nxgb
		if self.lo != None: self.lo.set_text(' loc. '+str(self.cutloc))
		if self.gb != None: self.gb.set_text('in {1,...,'+str(self.imax)+'}')
		self.getIsoLine(dat)
		self.fi.canvas.draw()

	def getIsoLine(self,dat):
		if self.co != None:
			for coll in self.co.collections:
				coll.remove()
			self.co = None
		if self.isoc != None:
			if not unif: dat = dat.T
			if self.pl == 'z': self.co = self.ax.contour(x/Lx,y/Ly,dat,[self.isoc,],colors='r')
			if self.pl == 'y': self.co = self.ax.contour(x/Lx,z/Lz,dat,[self.isoc,],colors='r')
			if self.pl == 'x': self.co = self.ax.contour(z/Lz,y/Ly,dat,[self.isoc,],colors='r')
		
	def reqUpdate(self):
		self.wait.set_text('WAIT...'); self.fi.canvas.draw()		
		os.system('touch snap')
		snap_time = os.path.getmtime('snap')
		# wait for a data file if none present
		hello = os.path.exists('restart.bin')
		while not hello:
			hello = os.path.exists('restart.bin')
		# wait for a new restart file (if dNami is still running)
		if not os.path.exists('stop'):
			data_time = os.path.getmtime('restart.bin')
			while snap_time > data_time:
				data_time = os.path.getmtime('restart.bin')
				plt.pause(0.001)
		# load & show new data
		self.timeUpdate()
		self.wait.set_text(''); self.fi.canvas.draw()

	def autoSwitch(self):
		if self.auto:
			self.auto = False
			self.au.set_text('OFF'); self.au.set_color('r')
		else:
			self.auto = True
			self.au.set_text('ON'); self.au.set_color('g')
		self.fi.canvas.draw()
	
	def autoUpdate(self):
		self.rcur += 1
		finished = os.path.exists('stop')
		if self.auto and self.rcur < 1000000 and not finished:
			self.reqUpdate()
			plt.pause(0.001)
			self.autoUpdate()
		else:
			self.auto = False; self.au.set_text('OFF'); self.au.set_color('r'); self.fi.canvas.draw()
			self.rcur = 0
	
	def on_press(self,event):
		if event.key == 'a': self.autoSwitch(); self.autoUpdate()
		if event.key == 'c':
			if os.path.exists('snap'): os.system('rm snap')
			if os.path.exists('stop'): os.system('rm stop')
		if event.key == 'q': print('[info] closing live viewer'); sys.exit()
		if event.key == 'r': 
			self.setDefault(); self.fi.tbox_ind.set_val(str(self.ix)); self.fi.tbox_qty.set_val(str(self.qt))
			self.fi.tbox_iso.set_val(str(self.isoc))
			self.newCut = True; self.cutUpdate()
		if event.key == 'u': self.reqUpdate()
		if event.key == 'x' and ndim == 3: self.pl = 'x'; self.newCut = True; self.cutUpdate()
		if event.key == 'y' and ndim == 3: self.pl = 'y'; self.newCut = True; self.cutUpdate()
		if event.key == 'z' and ndim == 3: self.pl = 'z'; self.newCut = True; self.cutUpdate()
	
	def get_ind(self,text):
		self.ix = int(text)
		self.cutUpdate()
			
	def get_var(self,text):
		self.qt = int(text)
		self.cutUpdate()
			
	def get_iso(self,text):
		if text == 'None':
			self.isoc = None
		else:
			self.isoc = float(text)
		self.cutUpdate()

# ============================================================================= fire liveView interface up!
# initialise
if os.path.exists('restart.bin'):
	n,time,f = read_restart()
else:
	n = -1; time = -1.
	if ndim == 3: f = np.zeros((nx,ny,nz,nvar),dtype=wp)
	if ndim == 2: f = np.zeros((nx,ny   ,nvar),dtype=wp)
	if ndim == 1: f = np.zeros((nx      ,nvar),dtype=wp)
# visualise
vu = liveView()
plt.show()

