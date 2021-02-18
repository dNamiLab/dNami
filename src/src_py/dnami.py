import numpy as np
import dnamiF
import sys
from rhsinfo import wp
import dnami_io

from dnami_mpi import type_mpi

# =============================================================================
# TOOL BOX
# =============================================================================

def cst(x):
	from rhsinfo import wp
	if wp == 'float32': x = np.float32(x)
	if wp == 'float64': x = np.float64(x)
	return x

# =============================================================================
# INIT + MEM ALLOC
# =============================================================================

def create_tree():
	dtree = {}

	dtree['stats'] = None
	dtree['libs']  = {'fort':
							   {'integers': None, 'floats': None, 'data': None}
					 ,'cache blocking': None}
	dtree['grid']  = {'size' :
							   {'nxgb': None, 'nygb': None, 'nzgb': None}
					 ,'geom' :
							   {'Lx'  : None, 'Ly'  : None, 'Lz'  : None
							   ,'dx'  : None, 'dy'  : None, 'dz'  : None
							   ,'x'   : None, 'y'   : None, 'z'   : None
							   ,'xloc': None, 'yloc': None, 'zloc': None}}							   					 		   
	dtree['eqns']  = {'qvec' : 
							   {'nvars': None, 'solved': None, 'stored': None, 'views': None}
					 ,'coeff': None
					 ,'time' : None
					 ,'ndim' : None}

	dtree['misc']  = {'verbose': None, 'working precision': None}
	dtree['mpi']   = {'split': 
							   {'nxpr': None, 'nypr': None, 'nzpr': None}
					 ,'dMpi' : None}		   	
	dtree['num']   = {'hlo'  : None
					 ,'deriv': 
							   {'order': None, 'stencil': None, 'hlo': None} 
					 ,'filtr': 
				     	       {'order': None, 'stencil': None, 'hlo': None,'eps': None}
					 ,'tint' : 
					           {'tstep': None, 'cfl': None, 'itn': None}}

	dtree['bc']    = {'wall': 
							   {'isoT': None, 'zeroQ': None, 'slip': None}}	

	dtree['usr']   = None
	dtree['ios']   = None
	
	from rhsinfo import dim, stencil, order, coefficients, varname, varsolved, varstored, varbc, wp,hlo_rhs

	dtree['eqns']['qvec']['solved'] = []
	dtree['eqns']['qvec']['stored'] = []
	dtree['eqns']['qvec']['bcs']    = {'face':{'i' :[],'j' :[],'k' :[]},
									   'edge':{'ij':[],'jk':[],'ik':[]}}

	for v in varsolved:
		dtree['eqns']['qvec']['solved'].append([v,varname[v]])	

	for v in varstored:
		dtree['eqns']['qvec']['stored'].append([v,varstored[v]['ind']])

	for v in varbc:
		for bcloc in ['face','edge']:
			if bcloc in varbc[v]:
				dtree['eqns']['qvec']['bcs'][bcloc][varbc[v][bcloc]].append([v,varbc[v]['ind']])


	dtree['eqns']['coeff'] = []
	for v in coefficients:
		dtree['eqns']['coeff'].append([v,coefficients[v]])	

	dtree['eqns']['qvec']['nvars']   = len(varname)#+len(dtree['eqns']['qvec']['stored'])	
	dtree['num']['deriv']['stencil'] = stencil
	dtree['num']['deriv']['hlo']     = hlo_rhs #int((stencil-1)/2)
	dtree['num']['deriv']['order']   = order
	
	# if dtree['num']['filtr']['hlo'] != None:
	# 	dtree['num']['hlo'] = max(dtree['num']['deriv']['hlo'],dtree['num']['filtr']['hlo'])
	# else:
	# 	dtree['num']['hlo'] = dtree['num']['deriv']['hlo']	
                
	dtree['num']['hlo'] = hlo_rhs			
	
	dtree['eqns']['ndim'] = dim
	dtree['misc']['working precision'] = wp
	dtree['misc']['verbose'] = True

	# dtree['libs']['cache blocking'] = [256,2,6] # good for 11 pts 3D, div. forme of N.S.
	
	dtree['libs']['cache blocking'] = [2560000,2,6]

	# recover BCs info:
	
	try:
	 from rhsinfo import bc_info
	except: 
		bc_info = [{},{}]

	dtree['bc']	 = {'allbc':bc_info[1],'mybc':[]} # OVERWRITE predefined 'bc' key.

	return dtree

def start_mpi(tree):
	dMpi = type_mpi(tree)
	tree['mpi']['dMpi'] = dMpi
	return tree

def unpack_bcs(tree):
	dMpi  = tree['mpi']['dMpi'] 
	iSwap = dMpi.iSwap
	bcs   = tree['bc']
	nxgb, nygb, nzgb = tree['grid']['size']['nxgb'], tree['grid']['size']['nygb'], tree['grid']['size']['nzgb']
	
	for dirBC in bcs['allbc']:

		if dirBC == 'i1':
			if dMpi.ibeg == 1:
				iSwap['i1'] = False
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'imax':
			if dMpi.iend == nxgb:
				iSwap['imax'] = False
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'j1':
			if dMpi.jbeg == 1:
				iSwap['j1'] = False
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'jmax':
			if dMpi.jend == nygb:
				iSwap['jmax'] = False
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'k1':
			if dMpi.kbeg == 1:
				iSwap['k1'] = False
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'kmax':
			if dMpi.kend == nzgb:
				iSwap['kmax'] = False
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]

# edges	

		if dirBC == 'i1j1':
			if dMpi.ibeg == 1 and dMpi.jbeg == 1:
				iSwap['i1'] = False
				iSwap['j1'] = False				
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'i1jmax':
			if dMpi.ibeg == 1 and dMpi.jend == nygb:
				iSwap['i1'] = False
				iSwap['jmax'] = False				
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'imaxj1':
			if dMpi.iend == nxgb and dMpi.jbeg == 1:
				iSwap['imax'] = False
				iSwap['j1'] = False				
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'imaxjmax':
			if dMpi.iend == nxgb and dMpi.jend == nygb:
				iSwap['imax'] = False
				iSwap['jmax'] = False				
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]			

		if dirBC == 'i1k1':
			if dMpi.ibeg == 1 and dMpi.kbeg == 1:
				iSwap['i1'] = False
				iSwap['k1'] = False				
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'i1kmax':
			if dMpi.ibeg == 1 and dMpi.kend == nzgb:
				iSwap['i1']   = False
				iSwap['kmax'] = False				
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'imaxk1':
			if dMpi.iend == nxgb and dMpi.kbeg == 1:
				iSwap['imax'] = False
				iSwap['k1']   = False				
				bcs['mybc']   = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'imaxkmax':
			if dMpi.iend == nxgb and dMpi.kend == nzgb:
				iSwap['imax'] = False
				iSwap['kmax'] = False				
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'j1k1':
			if dMpi.jbeg == 1 and dMpi.kbeg == 1:
				iSwap['j1'] = False
				iSwap['k1'] = False				
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'j1kmax':
			if dMpi.jbeg == 1 and dMpi.kend == nzgb:
				iSwap['j1']   = False
				iSwap['kmax'] = False				
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'jmaxk1':
			if dMpi.jend == nygb and dMpi.kbeg == 1:
				iSwap['imax'] = False
				iSwap['k1']   = False				
				bcs['mybc']   = bcs['mybc'] + bcs['allbc'][dirBC]

		if dirBC == 'jmaxkmax':
			if dMpi.iend == nygb and dMpi.kend == nzgb:
				iSwap['jmax'] = False
				iSwap['kmax'] = False				
				bcs['mybc'] = bcs['mybc'] + bcs['allbc'][dirBC]			


	dMpi.iSwap = iSwap
	tree['bc'] = bcs

	return tree



def allocate(tree):	

	wp     = tree['misc']['working precision']	
	dim    = tree['eqns']['ndim']
	nvar   = tree['eqns']['qvec']['nvars']
	nvarst = len(tree['eqns']['qvec']['stored'])	
	hlod   = tree['num']['deriv']['hlo']
	hlof   = tree['num']['filtr']['hlo']
  
	hlo    = tree['num']['hlo']  
	nx     = tree['mpi']['dMpi'].nx 
	ny     = tree['mpi']['dMpi'].ny 
	nz     = tree['mpi']['dMpi'].nz 

	variables = []
	if tree['eqns']['qvec']['solved']:
		for v in tree['eqns']['qvec']['solved']:
			variables.append(v[1])

	if tree['eqns']['qvec']['stored']:
		for v in tree['eqns']['qvec']['stored']:
			variables.append(v[1])	

	variables_face = {}
	for dir in ['i','j','k']:
		variables_face[dir] = []
		for v in tree['eqns']['qvec']['bcs']['face'][dir]:
				variables_face[dir].append(v[1])			

	variables_edge = {}
	for dir in ['ij','jk','ik']:
		variables_edge[dir] = []
		for v in tree['eqns']['qvec']['bcs']['edge'][dir]:
				variables_edge[dir].append(v[1])		

	nvar_face = {'i': len(variables_face['i']),
				 'j': len(variables_face['j']),
				 'k': len(variables_face['k'])}

	nvar_edge = {'ij': len(variables_edge['ij']),
				 'jk': len(variables_edge['jk']),
				 'ik': len(variables_edge['ik'])}			 
				
	coeff = []
	if tree['eqns']['coeff']:
		for v in tree['eqns']['coeff']:
			coeff.append(v[1])

	sizex = nx + 2*hlo
	
	if(ny == 1):
		sizey = 1
	else:		
		sizey = ny + 2*hlo
	
	if(nz == 1) : 
		sizez = 1
	else:	
		sizez = nz + 2*hlo
	

	ndimpt  = sizex*sizey*sizez
	ndimtot = ndimpt*nvar


	# alloc bcs fields:
	ndimptbcs = {}
	ndimptbcs['i']  =sizey*sizez
	ndimptbcs['j']  =sizex*sizez
	ndimptbcs['k']  =sizex*sizey
	ndimptbcs['ij'] =sizez 
	ndimptbcs['jk'] =sizex 
	ndimptbcs['ik'] =sizey 

	sizebcs = {}
	sizebcs['i']  =(sizey,sizez)
	sizebcs['j']  =(sizex,sizez)
	sizebcs['k']  =(sizex,sizey)
	sizebcs['ij'] =(sizez) 
	sizebcs['jk'] =(sizex) 
	sizebcs['ik'] =(sizey) 

	# faces:
	nfacei = max(1,nvar_face['i']*ndimptbcs['i'])
	nfacej = max(1,nvar_face['j']*ndimptbcs['j'])
	nfacek = max(1,nvar_face['k']*ndimptbcs['k'])
	# edges:
	nedgeij = max(1,nvar_edge['ij']*ndimptbcs['ij'])
	nedgejk = max(1,nvar_edge['jk']*ndimptbcs['jk'])
	nedgeik = max(1,nvar_edge['ik']*ndimptbcs['ik'])


	# unpack bc info:
	tree = unpack_bcs(tree)
	
	mybc = tree['bc']['mybc']
	nbc = len(mybc)

	# Integers parameters to be passed to the fortran layer
	param_intDim = 12 + nvar + nvarst + 1 + nbc + 6 + 6
	param_int = np.empty(shape=param_intDim, dtype=np.int32,   order='F')     
	param_int[0]              = hlo
	param_int[1]              = nx
	param_int[2]              = ny
	param_int[3]              = nz
	param_int[4]              = nvar
	param_int[5]              = nvarst
	param_int[6]              = ndimtot
	param_int[7]              = 3                        # for RK sub steps in Python
	param_int[8]              = ndimpt
	param_int[9:9+3]          = tree['libs']['cache blocking']
	adr = 9+3
	param_int[adr:adr+nvar]   = variables[0:0+nvar]  # variables location (in q).
	adr = adr+nvar
	param_int[adr:adr+nvarst] = variables[nvar:nvar+nvarst]  # variables location (in qst).
	adr = adr+nvarst
	param_int[adr]            = nbc
	adr = adr + 1
	param_int[adr:adr+nbc]    = mybc
	adr = adr + nbc
	param_int[adr:adr+6]      = list(nvar_face.values()) + list(nvar_edge.values())
	adr = adr + 6
	param_int[adr:adr+6]      = [nfacei,nfacej,nfacek,nedgeij,nedgejk,nedgeik]

	if tree['eqns']['coeff']: 
		ncoef    = len(tree['eqns']['coeff'])
	else:
		ncoef = 0	

	# Floating point parameters to be passed to the Fortran layer (3 additional floats for the metrics, uniforme grid + 1 for dt  +1 for eps filter)
	param_float = np.zeros(shape=ncoef+5, dtype=wp, order='F')   

	param_float[0] = cst(1.0)/tree['grid']['geom']['dx']
	if(dim>=2)   : param_float[1] = cst(1.0)/tree['grid']['geom']['dy']
	if(dim == 3) : param_float[2] = cst(1.0)/tree['grid']['geom']['dz']

	param_float[3] = tree['num']['tint']['tstep']
	param_float[4] = tree['num']['filtr']['eps']

	for i,v in enumerate(tree['eqns']['coeff']):
		param_float[i+5] = v[1]

	# Floating point array (contiguous, aligned, actually NOT aligned yet...)
	
	nfieldbcs = sum([nfacei,nfacej,nfacek,nedgeij,nedgejk,nedgeik])

	if nvarst != 0 :
		data     = np.empty(shape=ndimtot*4 + ndimpt*nvarst + nfieldbcs, dtype=wp,order='F') # 4 -->  q,q1,q2 + rhs + nvarstored
	else:	
		data     = np.empty(shape=ndimtot*4 + 1             + nfieldbcs, dtype=wp,order='F') # 4 -->  q,q1,q2, rhs, + 1 (address for qst in fortran layer)

	# Explicit view of data (only references here, no copy) 
	views = {}
	nvsolved = len(tree['eqns']['qvec']['solved'])
	
	# WARNING assume contiguous addresses of stored variables in data_float...
	if nvarst != 0:
		addrstored_beg = ndimtot*4 
		addrstored_end = addrstored_beg + ndimpt*nvarst
	else:
		addrstored_beg = ndimtot*4	
		addrstored_end = addrstored_beg + 1

	addrbcfields_beg      = addrstored_end	
	addrbcfields_edge_beg = addrbcfields_beg + sum([nfacei,nfacej,nfacek])

	if dim == 3:
		views['q']       = data[0:ndimpt*nvsolved].view().reshape(sizex,sizey,sizez,nvsolved, order='F')
		if nvarst !=0:
			views['qstored'] = data[addrstored_beg:addrstored_end].view().reshape(sizex,sizey,sizez,nvarst, order='F')
		
		for v in tree['eqns']['qvec']['solved']:
			views[v[0]] = data[(v[1]-1)*(ndimpt):(v[1])*ndimpt].view().reshape(sizex,sizey,sizez, order='F')
		for v in tree['eqns']['qvec']['stored']:
			views[v[0]] = data[(v[1]-1)*(ndimpt)+addrstored_beg:(v[1])*ndimpt+addrstored_beg].view().reshape(sizex,sizey,sizez, order='F')			

		# bc faces:	
		shift = addrbcfields_beg
		for dir in ['i','j','k']:
			for v in tree['eqns']['qvec']['bcs']['face'][dir]:				
				views[v[0]] = data[shift:ndimptbcs[dir]+shift].view().reshape(sizebcs[dir], order='F')			
				shift       = shift + ndimptbcs[dir]	

		# bc edges:	
		shift = addrbcfields_edge_beg
		for dir in ['ij','jk','ik']:
			for v in tree['eqns']['qvec']['bcs']['edge'][dir]:	
				views[v[0]] = data[shift:ndimptbcs[dir]+shift].view().reshape(sizebcs[dir], order='F')			
				shift       = shift + ndimptbcs[dir]
	elif dim == 2:
		views['q']       = data[0:ndimpt*nvsolved].view().reshape(sizex,sizey,nvsolved, order='F')
		if nvarst !=0:
			views['qstored'] = data[addrstored_beg:addrstored_end].view().reshape(sizex,sizey,nvarst, order='F')
		
		for v in tree['eqns']['qvec']['solved']:
			views[v[0]] = data[(v[1]-1)*(ndimpt):(v[1])*ndimpt].view().reshape(sizex,sizey, order='F')	
		for v in tree['eqns']['qvec']['stored']:
			views[v[0]] = data[(v[1]-1)*(ndimpt)+addrstored_beg:(v[1])*ndimpt+addrstored_beg].view().reshape(sizex,sizey, order='F')	

		# bc faces:	
		shift = addrbcfields_beg
		for dir in ['i','j']:
			for v in tree['eqns']['qvec']['bcs']['face'][dir]:	
				views[v[0]] = data[shift:ndimptbcs[dir]+shift].view().reshape(sizebcs[dir], order='F')			
				shift       = shift + ndimptbcs[dir]

						
	else:
		views['q'] =  data[0:ndimpt*nvsolved].view().reshape(sizex,nvsolved, order='F')
		if nvarst !=0:
			views['qstored'] = data[addrstored_beg:addrstored_end].view().reshape(sizex,nvarst, order='F')
		
		for v in tree['eqns']['qvec']['solved']:
			views[v[0]] = data[(v[1]-1)*(ndimpt):(v[1])*ndimpt].view().reshape(sizex, order='F')		
		for v in tree['eqns']['qvec']['stored']:
			views[v[0]] = data[(v[1]-1)*(ndimpt)+addrstored_beg:(v[1])*ndimpt+addrstored_beg].view().reshape(sizex, order='F')									
	
	tree['libs']['fort']['integers'] = param_int
	tree['libs']['fort']['floats']   = param_float
	tree['libs']['fort']['data']     = data
	tree['eqns']['qvec']['views']    = views

	dnamiF.init(param_int,param_float,data)	
	
	return tree

# =============================================================================
# GEOMETRY
# =============================================================================

def local_coordinates(tree, var, beg, end, idx, ngb, nloc):
	""" Returns the coordinate values local to an MPI process."""
	hlo, bc = tree['num']['hlo'], tree['bc']['allbc']
	# Each local array is padded by halos on either side, to be consistent with the shape of the q vector arrays
	out = np.zeros(nloc + 2*hlo, dtype=wp) # Explicitly define the working precision
	if ('%s1' % idx in bc) or ('%smax' % idx in bc):
		if beg == 1 and end == ngb: # No MPI => local coordinate has size n+2*hlo
			out = var[:]
		elif beg == 1: # Bottom of the domain, fill from the first halo to local n+hlo
			out[0:nloc+hlo] = var[beg-1:end+hlo]
			# Overlap with the layer above
			out[nloc+hlo:] = var[end+hlo:end+2*hlo]
		elif end == ngb: # Top of the domain, fill from hlo to n+2*hlo
			out[hlo:nloc+2*hlo] = var[beg-1+hlo:end+2*hlo]
			# Overlap with the layer below
			out[0:hlo] = var[beg-1:beg-1+hlo]
		else: # Away from boundaries, place the coordinates padded by hlo on either side
			out[hlo:nloc+hlo] = var[beg-1+hlo:end+hlo]
			# Overlap with the layers above/below
			out[0:hlo] = var[beg-1:beg-1+hlo]
			out[nloc+hlo:] = var[end+hlo:end+2*hlo]
	else: # Periodic boundary condition
		out = np.zeros(nloc)
		out[:] = var[beg-1:end]
	return out

def create_grid(tree, x1=None, x2=None, y1=None, y2=None, z1=None, z2=None, x_coordinates=None, y_coordinates=None, z_coordinates=None):
	
	wp   = tree['misc']['working precision'] 
	ndim = tree['eqns']['ndim']
	grid = tree['grid']['size']
	geom = tree['grid']['geom']

	nxgb, nygb, nzgb  = grid['nxgb'], grid['nygb'], grid['nzgb'] 
	Lx  , Ly  , Lz    = geom['Lx']  , geom['Ly']  , geom['Lz'] 

	dmpi = tree['mpi']['dMpi']
	ibeg,jbeg,kbeg = dmpi.ibeg,dmpi.jbeg,dmpi.kbeg
	iend,jend,kend = dmpi.iend,dmpi.jend,dmpi.kend
	
	hlo = tree['num']['hlo']

	if not x1:
		x1 = cst(0.0)
	if not x2:
		x2 = Lx
	if not y1:
		y1 = cst(0.0)
	if not y2:
		y2 = Ly
	if not z1:
		z1 = cst(0.0)
	if not z2:
		z2 = Lz

	# create domain
	#                                   pt 1     2          n-1    n     period
	#                                   0                             L  |
	# full domain is [0,L] but careful: |--o--|--o--| ... |--o--|--o--|--V
	#                                       \___________________________/ 
	
	bc = tree['bc']['allbc']
	if x_coordinates is None:
		if ('i1' in bc) or ('imax' in bc):
			dx = Lx/cst(nxgb+2*hlo-1)
			x  = np.linspace(x1,x2,nxgb+2*hlo,dtype=wp)
		else:	
			dx = Lx/cst(nxgb)
			x  = np.arange(dx/cst(2.),Lx,dx,dtype=wp)
	else:
		x = x_coordinates[:]
		dx = Lx/cst(nxgb+2*hlo-1) # Only added dx for non-periodic

	if nygb > 1:
		if y_coordinates is None:
			if ('j1' in bc) or ('jmax' in bc):
				dy = Ly/cst(nygb+2*hlo-1)
				# y  = np.linspace(cst(0.0),Ly,nygb+2*hlo,dtype=wp)
				y = np.linspace(y1, y2, nygb+2*hlo,dtype=wp)
			else:		
				dy = Ly/cst(nygb)
				y  = np.arange(dy/cst(2.),Ly,dy,dtype=wp)
		else:
			y = y_coordinates[:]
			dy = Ly/cst(nygb+2*hlo-1) # Only added dy for non-periodic
	else:
		Ly = cst(0.); y = []; dy = cst(0.)
		
	if nzgb > 1:
		if z_coordinates is None:
			if ('k1' in bc) or ('kmax' in bc):
				dz = Lz/cst(nzgb+2*hlo-1)
				z  = np.linspace(x1,Lz,nzgb+2*hlo,dtype=wp)
			else:			
				dz = Lz/cst(nzgb)
				z  = np.arange(dz/cst(2.),Lz,dz,dtype=wp)
		else:
			z = z_coordinates[:]
			dz = Lz/cst(nzgb+2*hlo-1) # Only added dz for non-periodic
	else:
		Lz = cst(0.); z = []; dz = cst(0.)

	geom['dx'], geom['dy'], geom['dz'] = dx, dy, dz
	geom['x'] , geom['y'] , geom['z']  = x , y , z

	# global             iend ! ibeg                    iend  !  ibeg
	# pt#          n-1    n   !   1     2          n-1    n   !   1     2
	#        ... |--o--|--o--|!|--o--|--o--| ... |--o--|--o--|!|--o--|--o--| ...
	# w/ hlo     [<------------------------------------------------------->]
	# loc py ind    0     1      hlo                   hlo+n-1      n+2*hlo-1
	# glo py ind             ibeg+hlo-1               iend+hlo-1

	# Create local coordinate arrays on each MPI process depending on whether the boundary is periodic or not
	# Local MPI sizes
	nx = tree['mpi']['dMpi'].nx 
	ny = tree['mpi']['dMpi'].ny
	nz = tree['mpi']['dMpi'].nz
	geom['xloc'] = local_coordinates(tree, x, ibeg, iend, 'i', nxgb, nx)
	if ndim == 2:
		geom['yloc'] = local_coordinates(tree, y, jbeg, jend, 'j', nygb, ny)
		# xx,yy = np.meshgrid(xloc,yloc,sparse=False,indexing='ij')
	if ndim == 3:
		geom['yloc'] = local_coordinates(tree, y, jbeg, jend, 'j', nygb, ny)
		geom['zloc'] = local_coordinates(tree, z, kbeg, kend, 'k', nzgb, nz)
		# xx,yy,zz = np.meshgrid(xloc,yloc,zloc,sparse=False,indexing='ij')
	if ndim > 3:
		raise NotImplementedError("Local coordinates are only implemented up to ndim=3.")
	return tree
