"""
This module contains functions for generating 
the FORTRAN code
"""

import re
import os
import numpy as np
import sys

dir_path = os.getcwd()

try:  
	os.mkdir(dir_path+'/pymod')
except:
  print('\033[1;40;94m'+'[build info]'+'\033[0m'+' \'/pymod\'  directory already exists') 
  
try:  
	os.mkdir(dir_path+'/src_for/includes/gen')
except:
  print('\033[1;40;94m'+'[build info]'+'\033[0m'+' \'/src_for/includes/gen\' directory already exists') 
 
class rhs_info:
	def __init__(self):
		from genRhs	import dim, wp

		try:
			from genRhs import consvar
		except:
			consvar = []
		try:
			from genRhs import varbc
		except:
			varbc = {}
		try:
			from genRhs import varstored
		except:
			varstored = {}	
		try:
			from genRhs import varloc
		except:
			varloc = {}		
		try:
			from genRhs import coefficients
		except:
			coefficients = {}		
		try:
			from genRhs import varsolved
		except:
			exception('varsolved not defined',message='error')

		try:
			from genRhs import varname
		except:
			exception('varname not defined',message='error')

		self.wp     	  = wp     	 
		self.dim     	  = dim     	 
		self.stencil 	  = 1 	 
		self.order        = 1      
		self.consvar      = consvar 
		self.coefficients = coefficients
		self.varname      = varname  
		self.varstored    = varstored  
		self.varsolved    = varsolved
		self.varbc        = varbc
		self.hlo_rhs      = 1
		self.bc_info      = [{},{}]



#################################################################################
#
# Finite differences coefficients database (could be automatically computed)
#
#################################################################################

# NOTE: CONSIDER USING SYMPY
# sympy.finite_diff_weights(2,[-2,-1,0,1,2],0)[-1][-1]
# where sympyp.finite_diff_weights(order of derivative, stencil,0)[-1][-1]


fdd5o4 = [-0.08333333333333333,
		   1.3333333333333333,
		  -2.5,
		   1.3333333333333333,
		  -0.08333333333333333]

fdd3o2 = [1.0,
		 -2.0,
		  1.0]

fddrs4o2 = [ 2., # right-sided 2nd order 2nd derivative
			-5.,
			+4.,
			-1.]		  

fd3o2 = [-0.5,0.5]

fdrs3o2 = [-1.5, # right-sided 2nd order first derivative
		  2.0,
		   -0.5]

fd5o4 = [0.08333333333333333,
        -0.666666666666667,     
         0.666666666666667,
        -0.08333333333333333]

fd11o10 = [-0.0007936507936507937,
		    0.009920634920634920,
		   -0.05952380952380952,
		    0.238095238095238,
		   -0.833333333333333,
  			0.833333333333333,
  		   -0.238095238095238,
  		    0.05952380952380952,
  		   -0.009920634920634920,
  		    0.0007936507936507937]

fd9o8 = [0.003571428571428571,
 		-0.03809523809523810 ,  
 		 0.200000000000000   ,
 		-0.800000000000000   ,    
 		 0.800000000000000   ,
 		-0.200000000000000   ,  
 		 0.03809523809523810 , 
 		-0.003571428571428571]


fd7o6 = [-0.01666666666666667,
 		  0.150000000000000  , 
 		 -0.750000000000000  ,
 		  0.750000000000000  , 
 		 -0.150000000000000  ,
 		  0.01666666666666667]

 

fd13o4 = [0.0007292761780923,
		 -0.0068802952724707,
		  0.0333125681563120,
		 -0.1126517314083857,
		  0.3128958166199386,
		 -0.8910608923461686,
		  0.8910608923461686,
		 -0.3128958166199386,
		  0.1126517314083857,
		 -0.0333125681563120,
		  0.0068802952724707,
		 -0.0007292761780923]


flt13o8 = [0.0006479157828976,
           -0.0057720247824029,
            0.0242842514706709,
           -0.0648197120657933,
            0.1238033248972423,
           -0.1794082631518038,
            0.2025290156983785,
           -0.1794082631518038,
            0.1238033248972423,
           -0.0648197120657933,
            0.0242842514706709,
           -0.0057720247824029,
            0.0006479157828976]       

flt3o2 = [-0.2500000000000000,
           0.5000000000000000,
          -0.2500000000000000]

flt5o4 = [6.250000000000000E-002,
		 -0.250000000000000,     
  		  0.375000000000000,
  		 -0.250000000000000,
  		  6.250000000000000E-002]          

flt11o2 = [ -0.002999540835,
		    0.018721609157,
		   -0.059227575576, 
		    0.123755948787,
		   -0.187772883589,
		    0.215044884112,
		   -0.187772883589,
		    0.123755948787, 
		   -0.059227575576,
		    0.018721609157,
		   -0.002999540835]    

flt11o10 = [-0.0009765625000000000,
				    0.009765625000000000,
				   -0.04394531250000000,
				    0.117187500000000,
				   -0.205078125000000,     
				    0.246093750000000,
				   -0.205078125000000,
				    0.117187500000000 ,
				   -0.04394531250000000,
				    0.009765625000000000,
				   -0.0009765625000000000]

flt_10_4 = [ 0.008391235145,   # 4 - 6
			-0.047402506444,
			 0.121438547725, 
			-0.200063042812,
			 0.240069047836,
			-0.207269200140,
			 0.122263107844,
			-0.047121062819, 
			 0.009014891495,
			 0.001855812216,
			-0.001176830044]

flt_10_3 = [ -0.000054596010, # 3 - 7
		      0.042124772446,
		     -0.173103107841,
		      0.299615871352,
		     -0.276543612935,
		      0.131223506571,
		     -0.023424966418,
		      0.013937561779,
		     -0.024565095706,
		      0.013098287852,
		     -0.002308621090]

flt_10_2 = [ 0.052523901012, # 2 - 8
			-0.206299133811,
			 0.353527998250,
			-0.348142394842,
			 0.181481803619,
			 0.009440804370,
			-0.077675100452,
			 0.044887364863,
			-0.009971961849,
			 0.000113359420,
			 0.000113359420]


flt_7_1 = [-0.085777408970,
		    0.277628171524,
		   -0.356848072173,
		    0.223119093072,
		   -0.057347064865,
		   -0.000747264596,
		   -0.000027453993]

flt_4_0 = [0.320882352941,
		  -0.465000000000,
		   0.179117647059,
		  -0.035000000000]             


derDic = {};der2Dic = {}; fltDic = {}

# ['stencil']['order']
derDic[13] = {4 :fd13o4 }
derDic[11] = {10:fd11o10}
derDic[9]  = {8 :fd9o8  }
derDic[7]  = {6 :fd7o6  }
derDic[5]  = {4 :fd5o4  }
derDic[3]  = {2 :fd3o2  }

der2Dic[5] = {4:fdd5o4}
der2Dic[3] = {2:fdd3o2}

fltDic[13] = {8 :flt13o8}
fltDic[11] = {10:flt11o10,2 :flt11o2}
fltDic[5]  = {4 :flt5o4}
fltDic[3]  = {2 :flt3o2}

derDicBc = {};der2DicBc = {}; fltDicBc = {}

der2DicBc[4] = {2: fddrs4o2}
derDicBc[3]  = {2: fdrs3o2 }
# ['layer'] 
#fltDicBc[4] = flt_10_4
fltDicBc[4] = flt_10_4 
fltDicBc[3] = flt_10_3
fltDicBc[2] = flt_10_2
fltDicBc[1] = [i*0.5 for i in flt_7_1]
fltDicBc[0] = [i*0.5 for i in flt_4_0]

operators = ('\+', '\*','\-','\/','\^','\(','\)','\,') 

opsplit = re.compile("|".join(operators)).split
oplist  = re.compile('[\^+\-*^/()\,]')

Variables  = lambda input: [s.strip() for s in opsplit(input)]
Operations = lambda input: re.findall(oplist,input)


def dNamiVar(var,rangei,rangej,rangek):
	
	from genRhs import dim

	try:
		from genRhs import varbc
	except:
		varbc = {}

	bcbydir   = {'face':{'i' :[],'j' :[],'k' :[]},
	    	     'edge':{'ij':[],'jk':[],'ik':[]}}

	for v in varbc:
		if 'face' in varbc[v]:	
			loctype = ''.join(sorted(varbc[v]['face'].replace('1','').replace('max','')))
			# print(v,loctype,varbc[v])
			bcbydir['face'][loctype].append(varbc[v]['ind'])
		elif 'edge' in varbc[v]:
			loctype = ''.join(sorted(varbc[v]['edge'].replace('1','').replace('max','')))
			bcbydir['edge'][loctype].append(varbc[v]['ind'])
		else:
			exception('Missing bc location for variable '+v+' in varbc (must contain a \'face\' or \'edge\' key)',message='error')
					
		
	
	for dirloc in bcbydir:
		for dir in bcbydir[dirloc]:
			knownbcs = bcbydir[dirloc][dir]
			
			from collections import Counter	
			if 	not (Counter(Counter(knownbcs).values())[1] == len(knownbcs)):
				exception('Several occurrences of identical indices for a given bc location (i.e. face or edge) in "varbc" ('+dirloc+', '+dir+','+str(knownbcs)+'). This ambiguity is not managed yet — generation aborted.',message='error')

	try:
		from genRhs import varstored
	except:
		varstored = {}	
	try:
		from genRhs import varloc
	except:
		varloc = {}		
	try:
		from genRhs import coefficients
	except:
		coefficients = {}		
	try:
		from genRhs import varsolved
	except:
		exception('varsolved not defined',message='error')

	try:
		from genRhs import varname
	except:
		exception('varname not defined',message='error')	


	sizebc   = {'face':{'i'  :{3:'('+rangej+','+rangek,
							   2:'('+rangej,
							   1:'(1'},
						'j'  :{3:'('+rangei+','+rangek,
							   2:'('+rangei,
							   1:'(1'},
						'k'  :{3:'('+rangei+','+rangek,
							   2:'('+rangei,
							   1:'(1'}},
		        'edge':{'ij' :{3:'('+rangek,
							   2:'(1',
							   1:'(1'},
						'jk' :{3:'('+rangei,
							   2:'(1',
							   1:'(1'},
						'ik' :{3:'('+rangej,
							   2:'(1',
							   1:'(1'}}} 

	try:
		from genRhs import hlo_glob
	except:
		exception('The number of hlo cells, hlo_glob, needs to be defined in genRhs. It is not automatically computed yet.',message='error')								   

	domainBorder = {'face':{'i1'  : 1-hlo_glob,
							  'imax': 'nx +'+str(hlo_glob),
						   	'j1'  : 1-hlo_glob,
						   	'jmax': 'ny +'+str(hlo_glob),
						   	'k1'  : 1-hlo_glob,
						   	'kmax': 'nz +'+str(hlo_glob)}}
	

	dvar  = var

	if var[0:3] == 'd1_' or var[0:3] == 'd2_':
		dvar = genVname(var,rangei,rangej,rangek)
	elif var in varname:
		if dim == 3 :
			dvar = 'q('+rangei+','+rangej+','+rangek+',indvars('+str(varname[var])+'))'
		elif dim == 2 :
			dvar = 'q('+rangei+','+rangej+',indvars('+str(varname[var])+'))'	
		elif dim == 1 :	
			dvar = 'q('+rangei+',indvars('+str(varname[var])+'))'
	elif var in varstored:
		if dim == 3 :
			dvar = 'qst('+rangei+','+rangej+','+rangek+',indvarsst('+str(varstored[var]['ind'])+'))'
		elif dim == 2 :
			dvar = 'qst('+rangei+','+rangej+',indvarsst('+str(varstored[var]['ind'])+'))'	
		elif dim == 1 :	
			dvar = 'qst('+rangei+',indvarsst('+str(varstored[var]['ind'])+'))'	
	elif var in varbc:
		for loc in ['face','edge']:	
			if loc in varbc[var]:	
				
				import sympy as sym

				rangeisym = sym.sympify(rangei)
				rangejsym = sym.sympify(rangej)
				rangeksym = sym.sympify(rangek)
				
				rangeloc  = {'i':rangeisym,
							 'j':rangejsym,
							 'k':rangeksym}

				loctype = ''.join(sorted(varbc[var][loc].replace('1','').replace('max','')))
				if rangeloc[loctype] == sym.sympify(domainBorder[loc][varbc[var][loc]]):
					dvar = 'q'+loc+'_'+loctype+sizebc[loc][loctype][dim]+','+str(varbc[var]['ind'])+')'
				elif var in coefficients:
					dvar = 'param_float('+str(coefficients[var])+' + 5)' 				
				elif var in varloc:
					dvar = op_to_dNami(varloc[var],rangei,rangej,rangek)		

	elif var in coefficients:
		dvar = 'param_float('+str(coefficients[var])+' + 5)' 				
	elif var in varloc:
		dvar = op_to_dNami(varloc[var],rangei,rangej,rangek)
		
	return dvar

def op_to_dNami(source,i='i',j='j',k='k'):
	rangei = i
	rangej = j
	rangek = k
	rhs_line = ''
	for v, op in zip(Variables(source)[0:-1],Operations(source)):
		if op=='^' : op = '**'
		newterm = (dNamiVar(v,rangei,rangej,rangek) + op)
		
		if (len(rhs_line)>1) and ( (rhs_line[-1] == '+' ) or (rhs_line[-1] == '-' ) ):

			rhs_line = rhs_line + '&'+'\n'
			rhsadd = '{:>20}{:}'.format('',newterm)
			rhs_line = rhs_line + rhsadd
		else:	
			rhs_line = rhs_line + dNamiVar(v,rangei,rangej,rangek) + op 
	rhs_line = rhs_line +  dNamiVar(Variables(source)[-1],rangei,rangej,rangek)
	return rhs_line

def comment(message):
	msg = '\n\n'
	msg = msg +'!'.ljust(60,'*')+'\n'
	msg = msg +'! '.ljust(60,' ')+'\n'
	msg = msg +('! '+message+' ').ljust(60,'*')+'\n'
	msg = msg +'! '.ljust(60,' ')+'\n'
	msg = msg +'!'.ljust(60,'*')+'\n'
	msg = msg +'\n\n'
	return msg

def genNbg(expr, dir , stencil, i='i',j='j',k='k',vname='v',dirBc=None,indbc='',der='first'):

	#***********************************************************
	#                                                           
	# prepare neighborhood for derivatives 
	#
	#     inputs :
	#             expr    = arithmetical expression to be derived (user frindly synthax)
	#             dir     = 'x','y','z' (direction used for the derivation)
	#             stencil = stencil size (center included)
	#             i,j,k   = location of the desired derivative 
	#             vname   = input name for the result (local variable)
	#
	#     outputs :
  #             exprnbg : list of 'expr' expressed for each points in the neighborhood
  #             vnamebg : corresponding temprary variable name
  # 
  #  Exemple :  
  #              genNbg('rho*u','x',3,vname='derhou_') 
  #
  #             --> exprnbg = ['q(i-1,j,k,nvrh)*q(i-1,j,k,nvux)', 
  #                            'q(i+0,j,k,nvrh)*q(i+0,j,k,nvux)', 
  #                            'q(i+1,j,k,nvrh)*q(i+1,j,k,nvux)']
  #
  #             --> vnamebg = ['derhou_im1jk', 'derhou_ip0jk', 'derhou_ip1jk']
	#           
	#                                                
	#***********************************************************

	exprnbg = []
	vnamebg = []

	hlo = int((stencil-1)/2)

	if (np.mod((stencil-1),2)==1) :
		import sys 
		exception('[error](genNbg) Wrong stencil length for centrered finite differences'+'\n                given value = '+str(stencil)+' (set in genRhs.py)',message='error')


	if(dirBc != None):
		if ( dir == 'x') and (dirBc[0] == 'i'):
			hlo = min(indbc,hlo)
		if ( dir == 'y') and (dirBc[0] == 'j'):
			hlo = min(indbc,hlo)
		if ( dir == 'z') and (dirBc[0] == 'k'):
			hlo = min(indbc,hlo)		
		
	if (hlo != 0):	
		for shift in range(-hlo,hlo+1,1):
			if   dir == 'x' :
				vnamebg.append(genVname(vname,  i=i+'{:+d}'.format(shift),j=j,k=k))
				exprnbg.append(op_to_dNami(expr,i=i+'{:+d}'.format(shift),j=j,k=k))
			elif dir == 'y' :
				vnamebg.append(genVname(vname,  i=i,j=j+'{:+d}'.format(shift),k=k))
				exprnbg.append(op_to_dNami(expr,i=i,j=j+'{:+d}'.format(shift),k=k))
			elif dir == 'z':
				vnamebg.append(genVname(vname,  i=i,j=j,k=k+'{:+d}'.format(shift)))
				exprnbg.append(op_to_dNami(expr,i=i,j=j,k=k+'{:+d}'.format(shift)))
	else:
		lenOneSidedBC = 3
		if (der == 'second'): lenOneSidedBC = 4

		offsetdir = 1
		if (dirBc[1:] == 'max'):
			offsetdir = -1
		for shift in range(0,lenOneSidedBC,1):
			if   dir == 'x' :
				vnamebg.append(genVname(vname,  i=i+'{:+d}'.format(offsetdir*shift),j=j,k=k))
				exprnbg.append(op_to_dNami(expr,i=i+'{:+d}'.format(offsetdir*shift),j=j,k=k))
			elif dir == 'y' :
				vnamebg.append(genVname(vname,  i=i,j=j+'{:+d}'.format(offsetdir*shift),k=k))
				exprnbg.append(op_to_dNami(expr,i=i,j=j+'{:+d}'.format(offsetdir*shift),k=k))
			elif dir == 'z':
				vnamebg.append(genVname(vname,  i=i,j=j,k=k+'{:+d}'.format(offsetdir*shift)))
				exprnbg.append(op_to_dNami(expr,i=i,j=j,k=k+'{:+d}'.format(offsetdir*shift)))			

	return [exprnbg, vnamebg]

def genVname(vname,i='i',j='j',k='k'):

	indi = re.sub(r'\-','m',i)
	indi = re.sub(r'\+','p',indi)

	indj = re.sub(r'\-','m',j)
	indj = re.sub(r'\+','p',indj)

	indk = re.sub(r'\-','m',k)
	indk = re.sub(r'\+','p',indk)

	vname = vname.strip() +indi+indj+indk
	return vname.strip()

def der(vnamebg,order,stencil,varname='derv',dirBC=None,indbc=None,bc=None):
	if stencil not in derDic:
		exception('Wrong stencil length for centrered finite differences'+'\n        given value = '+str(stencil)+' (set in genRhs.py)',message='error')

	if not(bc):
		if(len(vnamebg)!= stencil):
			raise ValueError('Wrong neighborhood ! ',len(vnamebg),stencil)
	
		fd = derDic[stencil][order]		
		der = varname +' ='
		hlo = int((stencil-1)/2)
	else:
		if len(vnamebg) == stencil:
			fd = derDic[len(vnamebg)][order]
		else:	
			fd = derDic[len(vnamebg)][len(vnamebg)-1]		
		der = varname +' ='
		hlo = int((len(vnamebg)-1)/2)	

	for shift in range(0,hlo):
		
		if der[-1]== '=':
			der = der + ' '+str(fd[shift])+'_wp*'+vnamebg[shift]	
		else:
			der = der + '{:+}'.format(fd[shift])+'_wp*'+vnamebg[shift]

	for shift in range(hlo,2*hlo):

		if der[-1]== '=':
			der = der + str(fd[shift])+'_wp*'+vnamebg[shift]	
		else:
			der = der + '{:+}'.format(fd[shift])+'_wp*'+vnamebg[shift+1]

	# Special case. First layer needs off-sided derivatives (WARNING enforce a 2nd order off-sided derivative !!!):

	if (indbc == 0 and bc):
		fd = derDicBc[len(vnamebg)][2]		
		der = varname +' ='
		hlo = int((len(vnamebg)-1)/2)

		offsetdir = 1
		if(dirBC[1:] == 'max'):
			offsetdir = -1
		for shift in range(0,len(vnamebg)):
			
			if der[-1]== '=':
				der = der + ' '+str(offsetdir*fd[shift])+'_wp*'+vnamebg[shift]	
			else:
				der = der + '{:+}'.format(offsetdir*fd[shift])+'_wp*'+vnamebg[shift]		


	der1 = re.sub('\-','-&\n          ' ,der )
	der  = re.sub('\+','+&\n          ',der1)		
	return der+'\n\n'		

def dder(vnamebg,order,stencil,varname='derv2',dirBC=None,indbc=None,bc=None):

	if not(bc):
		if(len(vnamebg)!= stencil):
			raise ValueError('Wrong neighborhood ! ',len(vnamebg),stencil)
	
		fd = der2Dic[stencil][order]		
		der = varname +' ='
		hlo = int((stencil-1)/2)
	elif(indbc != 0):
		if len(vnamebg) == stencil:
			fd = der2Dic[len(vnamebg)][order]
		else:	
			fd = der2Dic[len(vnamebg)][len(vnamebg)-1]		
		der = varname +' ='
		hlo = int((len(vnamebg)-1)/2)	

	# Special case. First layer needs off-sided derivatives (WARNING enforce a 2nd order off-sided derivative !!!):	
	if (indbc == 0 and bc):
		fd = der2DicBc[len(vnamebg)][2]		
		der = varname +' ='
		hlo = int((len(vnamebg)-1)/2)

		for shift in range(0,len(vnamebg)):
			
			if der[-1]== '=':
				der = der + ' '+str(fd[shift])+'_wp*'+vnamebg[shift]	
			else:
				der = der + '{:+}'.format(fd[shift])+'_wp*'+vnamebg[shift]	
	else:			
		for shift in range(0,hlo+1):
			
			if der[-1]== '=':
				der = der + ' '+str(fd[shift])+'_wp*'+vnamebg[shift]	
			else:
				der = der + '{:+}'.format(fd[shift])+'_wp*'+vnamebg[shift]
	
		for shift in range(hlo,2*hlo):
	
			if der[-1]== '=':
				der = der + str(fd[shift])+'_wp*'+vnamebg[shift]	
			else:
				der = der + '{:+}'.format(fd[shift+1])+'_wp*'+vnamebg[shift+1]
		
	der1 = re.sub('\-','-&\n          ' ,der )
	der  = re.sub('\+','+&\n          ',der1)		
	return der+'\n\n'		

def createNbg(vnamebg,exprnbg,file,der='first',indbc=None,bc=False):

	stencil = len(vnamebg)
	
	hlo = int((stencil-1)/2)
	lbound = hlo

	if (der == 'second' or (indbc==0 and bc)):
		lbound = lbound + 1

	if( np.mod(stencil,2) == 1 ):

		for shift in range(0,lbound):	
			file.write(vnamebg[shift  ]+' = '+exprnbg[shift  ])
			file.write('\n\n')
		for shift in range(hlo,2*hlo):
			file.write(vnamebg[shift+1]+' = '+exprnbg[shift+1])
			file.write('\n\n')		
	else:

		for shift in range(0,stencil):	
			file.write(vnamebg[shift  ]+' = '+exprnbg[shift  ])
			file.write('\n\n')

def updateRHS(vname,expr,i='i',j='j',k='k',update=False,updateq=False):
	
	from genRhs import varname, dim

	if dim == 3:
		rhs = 'rhs('+i+','+j+','+k+',indvars('+str(varname[vname])+'))'
	elif dim == 2:	
		rhs = 'rhs('+i+','+j+',indvars('+str(varname[vname])+'))'
	elif dim == 1:	
		rhs = 'rhs('+i+',indvars('+str(varname[vname])+'))'

	if not update:	
		out = rhs + ' = ' + ' ' + expr
	else:
		out = rhs + ' = ' + rhs + ' ' + expr	

	if updateq:	
		if dim == 3:
			q = 'q('+i+','+j+','+k+',indvars('+str(varname[vname])+'))'
		elif dim == 2:	
			q = 'q('+i+','+j+',indvars('+str(varname[vname])+'))'
		elif dim == 1:	
			q = 'q('+i+',indvars('+str(varname[vname])+'))'
	
		if not update:	
			out = q + ' = ' + ' ' + expr
		else:
			out = q + ' = ' + q + ' ' + expr		

	return out

def updateStored(vname,expr,i='i',j='j',k='k',update=False):
	
	from genRhs import varstored, dim

	if dim == 3:
		stored = 'qst('+i+','+j+','+k+',indvarsst('+str(varstored[vname]['ind'])+'))'
	elif dim == 2:	
		stored = 'qst('+i+','+j+',indvarsst('+str(varstored[vname]['ind'])+'))'
	elif dim == 1:	
		stored = 'qst('+i+',indvarsst('+str(varstored[vname]['ind'])+'))'

	if not update:	
		storedout = stored + ' = ' + ' ' + expr
	else:
		storedout = stored + ' = ' + stored + ' ' + expr	
	
	return storedout
	

def updateVarbc(vname,expr,i='i',j='j',k='k',update=False):
	
	from genRhs import varbc, dim

	sizebc   = {'face':{'i'  :{3:'('+j+','+k,
							   2:'('+j,
							   1:'('},
						'j'  :{3:'('+i+','+k,
							   2:'('+i,
							   1:'('},
						'k'  :{3:'('+i+','+k,
							   2:'('+i,
							   1:'('}},
				'edge':{'ij' :{3:'('+k,
							   2:'(',
							   1:'('},
						'jk' :{3:'('+i,
							   2:'(',
							   1:'('},
						'ik' :{3:'('+j,
							   2:'(',
							   1:'('}}} 

	for loc in ['face','edge']:	
		if loc in varbc[vname]:
			loctype = ''.join(sorted(varbc[vname][loc].replace('1','').replace('max','')))
			qbc = 'q'+loc+'_'+loctype+sizebc[loc][loctype][dim]+','+str(varbc[vname]['ind'])+')'	

	if not update:	
		qbcout = qbc + ' = ' + ' ' + expr
	else:
		qbcout = qbc + ' = ' + qbc + ' ' + expr	
	
	return qbcout
	

def compute_stored(StoredVar,Stencil,Order,output,localvar,update=False,rhs=None):		

		if update and rhs != None:
			hlo_rhs = rhs.hlo_rhs
		else:	
			hlo_rhs =int((Stencil-1)/2)
	
		coefficients = rhs.coefficients
		varname      = rhs.varname
		wp           = rhs.wp
		dim          = rhs.dim

		locvar = []
	
		output.write(comment('Start computing stored variables'))
		
		history={'x':{'var':[],'symb':[]},'y':{'var':[],'symb':[]},'z':{'var':[],'symb':[]}}
		history2={'x':[[],[]],'y':[[],[]],'z':[[],[]]}
		dhistory={'x':[[],[],[],[]],'y':[[],[],[],[]],'z':[[],[],[],[]]}

		loop_create('begin',output)

		for vstored in  StoredVar.keys():

			output.write(comment('building stored variable '+ vstored))
		
			op  = StoredVar[vstored]['symb'].replace(" ", "")+'\n'
			
			output.write('!'.ljust(60,'~')+'\n')
			output.write('!'+'\n')
			output.write('! '+op)
			output.write('!'+'\n')
			output.write('!'.ljust(60,'~')+'\n\n')	
			
			# for vn in dirrhs:				

			[Out,locvar,history]  = genSymbDer1(StoredVar[vstored]['symb'],output,locvar,order=Order,stencil=Stencil,indi='i',indj='j',indk='k',vname=vstored,history=history,dhistory=dhistory)			
			[Out,locvar,history2] = genSymbDer2(Out,output,locvar,order=Order,stencil=Stencil,indi='i',indj='j',indk='k',vname=vstored,history=history2)

			for v in ['x','y','z']:
				if len(history[v]['symb']) != 0:
					for s in history[v]['symb']:	
						if(re.findall('_d'+v+'d'+v+'_',s) != []):
							hlo_rhs = max(hlo_rhs, 2*int((Stencil-1)/2))

			if(history2 != {'x':[[],[]],'y':[[],[]],'z':[[],[]]} ): hlo_rhs = max(hlo_rhs, int((Stencil-1)/2))

			output.write(comment('Update stored variables '+StoredVar[vstored]['symb']))
			
			exp_stored = op_to_dNami(Out) 
		
			output.write(updateStored(vstored,exp_stored,update=update)+'\n\n')
		
		tmpvar = ''
		for var1 in locvar: 
			for var2 in var1:
				if not tmpvar: 
					tmpvar = tmpvar + ' ' +var2
				else:
					tmpvar = tmpvar + ',' + var2
					
			tmpvar = tmpvar + ' &\n            '
		
		tmpvar = tmpvar.rstrip()
			
		if(tmpvar !=''): 
			localvar.write('\n\n real(wp) :: ')
			localvar.write(tmpvar[:-1])

		rhs.hlo_rhs = hlo_rhs	
		rhs.stencil = max(rhs.stencil,Stencil)
		rhs.order   = max(rhs.order,Order)

		loop_create('end',output)

def compute_storedbc(StoredVar,Stencil,Order,output,localvar,dirBC,update=False,updateqbc=False,rhs=None):		

		# if update and rhs != None:
		hlo_rhs      = rhs.hlo_rhs
		dim          = rhs.dim
		coefficients = rhs.coefficients
		varname      = rhs.varname
		wp           = rhs.wp

		# else:	
		# 	hlo_rhs =int((Stencil-1)/2)

	
		indi = 'i'
		indj = 'j'
		indk = 'k'

		indiri = indi
		indirj = indj
		indirk = indk

		if dirBC == 'i1': 
			indi = '1-'+str(hlo_rhs)
		if dirBC == 'imax': 
			indi = 'nx+'+str(hlo_rhs)

		if dirBC == 'j1': 
			indj = '1-'+str(hlo_rhs)
		if dirBC == 'jmax': 
			indj = 'ny+'+str(hlo_rhs)


		if dirBC == 'k1': 
			indk = '1-'+str(hlo_rhs)
		if dirBC == 'kmax': 
			indk = 'nz+'+str(hlo_rhs)


		DirDic = {'i':{'dir':None,'indbc':None},
			      'j':{'dir':None,'indbc':None},
			      'k':{'dir':None,'indbc':None}}

		DirDic[dirBC[0]] = {'dir': dirBC , 'indbc' : None}
				

		rhsname = {}

		for v in StoredVar:
			rhsname[v] = v

		locvar = []

		# output.write(comment('Start building RHS layers for BC (3D): '+dirBC))
		
		updatest = True
		layerend = hlo_rhs
		if updateqbc: 
		   updatest = False
		   layerend = 1

		for layer in range(0,layerend): #BC layers 

			DirDic[dirBC[0]]['indbc'] = layer

			gen_eqns_bc(StoredVar,output,localvar,
			            rhsname,Order=Order,Stencil=Stencil,
			            indi     =indi,indj =indj,indk =indk,
			            DirDic   = DirDic,
			            vname    = rhsname,
			            update   = update,
			            updatest = updatest,
			            updateqbc= updateqbc)

def append_Rhs(Flx,Stencil,Order,rhsname,vname,update=False,rhs=None,stored=False):		

		from genRhs import incPATH

		dim = rhs.dim
		wp  = rhs.wp
		coefficients = rhs.coefficients
		varname      = rhs.varname
		varstored    = rhs.varstored
		varbc        = rhs.varbc
		varsolved    = rhs.varsolved
		consvar      = rhs.consvar



		if update and rhs != None:
			hlo_rhs = rhs.hlo_rhs
		else:	
			hlo_rhs =int((Stencil-1)/2)

		localvar       = open(incPATH+'includeRHS_locVar.f90','a+')
		output         = open(incPATH+'include_RHSLoops.f90','a+')

		from pathlib import Path

		Path(incPATH+'include_StoredVar.f90').touch()
		Path(incPATH+'includeRHS_locVarStored.f90').touch()
		
		Path(incPATH+'include_StoredVarStatic.f90').touch()
		Path(incPATH+'includeRHS_locVarStoredStatic.f90').touch()

		Path(incPATH+'include_qbcVar.f90').touch()
		Path(incPATH+'includeRHS_locVarqbc.f90').touch()

		Path(incPATH+'include_qbcVarStatic.f90').touch()
		Path(incPATH+'includeRHS_locVarqbcStatic.f90').touch()

		Path(incPATH+'include_bccmpstored.f90').touch()
		Path(incPATH+'include_bccmpstoredstatic.f90').touch()

		if stored:
			stincludenpath = incPATH+'include_StoredVar.f90'
			stlocalvarpath = incPATH+'includeRHS_locVarStored.f90'

			stincludenpathstatic = incPATH+'include_StoredVarStatic.f90'
			stlocalvarpathstatic = incPATH+'includeRHS_locVarStoredStatic.f90'			

			if os.stat(stincludenpath).st_size == 0:
				storedinclude  = open(stincludenpath,'w')
				storedlocalvar = open(stlocalvarpath,'w')
				storedincludeStatic  = open(stincludenpathstatic,'w')
				storedlocalvarStatic = open(stlocalvarpathstatic,'w')				
			else:
				exception('Stored variables already generated. \n        '
					      'dNami can\'t generate stored variables with multiple numerical parameters (yet!)...',message='error')
			
			static  = {}
			dynamic = {}

			for var in varstored:
				if varstored[var]['static']:
					static[var] = varstored[var]
				else:
					dynamic[var] = varstored[var]

			compute_stored(dynamic,Stencil,Order,storedinclude,storedlocalvar, update=False,rhs=rhs)
			if static:
				compute_stored(static,Stencil,Order,storedincludeStatic,storedlocalvarStatic, update=False,rhs=rhs)
			exception('varstored generated with parameters (stencil, order): ('+str(Stencil)+','+str(Order)+')',message='com')	
			# print(color('varstored generated with parameters (stencil, order): ('+str(Stencil)+','+str(Order)+')'))

		elif varstored != {}:
			exception('varstored not generated',message='com')	


		if update == False:
			rhs_for        = open(incPATH[:-13]+'rhs.for','w')
			rhs_template   = open(incPATH[:-13]+'rhs_template.for','r')
			lines = rhs_template.readlines()
			for l in lines:
				rhs_for.write(l)
			
			lbeg    	   = open(incPATH+'LOOP_BEGIN','w')
			lend           = open(incPATH+'LOOP_END','w')
			init           = open(incPATH+'init_numa.f90','w')
			cns            = open(incPATH+'primitive_to_conservative.f90','w')

			#init includes BCs:
			dummy  = open(incPATH+'select_phybc_rhs.f90' ,'a+')
			dummy  = open(incPATH+'select_phybc_q.f90'   ,'a+')			
			dummy  = open(incPATH+'index_filterbc.f90'   ,'a+')

			dummy  = open(incPATH+'selectstoredbc.f90'   ,'a+')
			dummy  = open(incPATH+'selectstoredstaticbc.f90'   ,'a+')

			dummy  = open(incPATH+'selectvarbcbc.f90'   ,'a+')
			dummy  = open(incPATH+'selectvarbcstaticbc.f90'   ,'a+')

			for d in ['x','y','z']:
				dummy  = open(incPATH+'selectfilterbc_'+d+'.f90' ,'a+')
				dummy  = open(incPATH+'selectupdate_filterbc_'+d+'.f90' ,'a+')			
				dummy  = open(incPATH+'update_filterbc_'+d+'.f90' ,'a+')
				dummy  = open(incPATH+'filterbc_'+d+'.f90' ,'a+')	

				dummy  = open(incPATH+'updateFilter_'+d,'w') 
				dummy  = open(incPATH+'Filter_'+d,'w')       								

			gendtype()
			globvar()
			loop(lbeg,lend)
			geninit(init,len(varsolved),rhs=rhs)
			if consvar != []:
				genP2C(cns,'conservative',rhs=rhs)
			else:
				genP2C(cns,'standard',rhs=rhs)	
			genbcsrc(len(varsolved),rhs=rhs)

		locvar = []
	
		output.write(comment('Start building RHS with source terms ('+str(dim)+'D)'))
		
		history={'x':{'var':[],'symb':[]},'y':{'var':[],'symb':[]},'z':{'var':[],'symb':[]}}
		history2={'x':[[],[]],'y':[[],[]],'z':[[],[]]}
		dhistory={'x':[[],[],[],[]],'y':[[],[],[],[]],'z':[[],[],[],[]]}

		loop_create('begin',output)

		for rhsvar in  Flx.keys():

			output.write(comment('building source terms in RHS for '+rhsname[rhsvar]))
		
			op  = Flx[rhsvar].replace(" ", "")+'\n'
			
			output.write('!'.ljust(60,'~')+'\n')
			output.write('!'+'\n')
			output.write('! '+op)
			output.write('!'+'\n')
			output.write('!'.ljust(60,'~')+'\n\n')	
			
			# for vn in dirrhs:				

			[Flx[rhsvar],locvar,history]  = genSymbDer1(Flx[rhsvar],output,locvar,order=Order,stencil=Stencil,indi='i',indj='j',indk='k',vname=vname[rhsvar],history=history,dhistory=dhistory)			
			[Flx[rhsvar],locvar,history2] = genSymbDer2(Flx[rhsvar],output,locvar,order=Order,stencil=Stencil,indi='i',indj='j',indk='k',vname=vname[rhsvar],history=history2)

			for v in ['x','y','z']:
				if len(history[v]['symb']) != 0:
					for s in history[v]['symb']:	
						if(re.findall('_d'+v+'d'+v+'_',s) != []):
							hlo_rhs = max(hlo_rhs, 2*int((Stencil-1)/2))

			if(history2 != {'x':[[],[]],'y':[[],[]],'z':[[],[]]} ): hlo_rhs = max(hlo_rhs, int((Stencil-1)/2))

			output.write(comment('Update RHS terms for '+rhsname[rhsvar]))
			
			exp_rhs = ' - '
			exp_rhs = exp_rhs + ' ( ' + op_to_dNami(Flx[rhsvar]) + ' ) '
		
			output.write(updateRHS(rhsvar,exp_rhs,update=update)+'\n\n')
		
		tmpvar = ''
		for var1 in locvar: 
			for var2 in var1:
				if not tmpvar: 
					tmpvar = tmpvar + ' ' +var2
				else:
					tmpvar = tmpvar + ',' + var2
					
			tmpvar = tmpvar + ' &\n            '
		
		tmpvar = tmpvar.rstrip()

		if(tmpvar !=''): 
			localvar.write('\n\n real(wp) :: ')
			localvar.write(tmpvar[:-1])
	
		rhs.hlo_rhs = hlo_rhs	
		rhs.stencil = max(rhs.stencil,Stencil)
		rhs.order   = max(rhs.order,Order)

		loop_create('end',output)

def genBC(Eqns,Stencil,Order,rhsname,vname,setbc=[False,{'bcname':{'i1':['rhs']}}],update=False,rhs=None,stored=False):

		from genRhs import incPATH

		dim          = rhs.dim
		wp           = rhs.wp
		coefficients = rhs.coefficients
		varname      = rhs.varname
		varstored    = rhs.varstored
		varsolved    = rhs.varsolved
		varbc        = rhs.varbc

		if rhs.stencil == 1:
			exception('RHS must be defined before BCs',message='error')

		else:	
			hlo_rhs     = rhs.hlo_rhs
			stencil_rhs = rhs.stencil
			order_rhs   = rhs.order 


		bctype = 'rhs'

		bc_info = rhs.bc_info[0]
			
		dir2gen = ['i1','imax','j1','jmax','k1','kmax']	

		cornerperm = {}
		for d1 in dir2gen:
			cornerperm[d1] = d1
			for d2 in dir2gen:
				cornerperm[d2] = d2
				if d1 != d2 and d1[0] != d2[0]:
						 cornerperm[d1+d2] = d2+d1

		if len(rhs.bc_info[1].values()) != 0:
			# print(rhs.bc_info[1].values())
			bcnum = max(max(rhs.bc_info[1].values()))
		else:
			bcnum = 0	

		newtype = False

		from itertools import permutations

		if setbc[0]: 			
			for bckey in setbc[1].keys():									
				bcname = bckey
				if bcname in bc_info: # physical BC name already generated

					for bcdir in setbc[1][bcname].keys():
						
						perms = [''.join(p) for p in permutations(bcdir)]

						if bcdir in bc_info[bcname] or cornerperm[bcdir] in bc_info : # physical BC direction already generated
									
							
							bcdirinf = list(set([cornerperm[bcdir],bcdir])&set(bc_info[bcname]))[0]

							for bckey2 in setbc[1][bcname][bcdir]:
								bctype = bckey2
								if bctype in bc_info[bcname][bcdirinf]:
									exception('bc already generated for: '+bcname+' with type '+bctype+' in dir '+bcdir,message='error')

								else:                 # new type for direction bcdir of bcname 
									# print('NEW TYPE')
									newtype = True
									bc_info[bcname][bcdirinf].append(bctype)

						else:                        # new direction for bcname
							bc_info[bcname][bcdir] = []
							for bckey3 in setbc[1][bcname][bcdir]:
								bctype = bckey3
								bc_info[bcname][bcdir].append(bctype)		
							
							if bcdir == cornerperm[bcdir]:

								# generate static/dynamic stored variables:			
	
								staticstored  = {}
								dynamicstored = {}
								staticvarbc   = {}
								dynamicvarbc  = {}

								for var in varstored:
									if varstored[var]['static']:
										staticstored[var] = varstored[var]['symb']
									else:
										dynamicstored[var] = varstored[var]['symb']

								addvarbc = False

								for var in varbc:
									if 'face' in varbc[var]:
										if varbc[var]['face'] == bcdir: addvarbc = True
									elif 'edge' in varbc[var]:
										if varbc[var]['edge'] == bcdir: addvarbc = True
									else:
										addvarbc = False

									if addvarbc:	
										if varbc[var]['static']:
											staticvarbc[var] = varbc[var]['symb']
										else:
											dynamicvarbc[var] = varbc[var]['symb']

									addvarbc = False				

	
								var2process = {'storedstatic':staticstored, 
											   'stored'      :dynamicstored,
											   'varbcstatic' :staticvarbc,
											   'varbc'       :dynamicvarbc}
	
								for k in var2process:
									if var2process[k] != {}:
	
										st       = open(incPATH+'bcsrc'+k+'_'+bcdir+'.for','w') # set to empty													
										locst    = open(incPATH+'include_BClocVar_'+k+'_'+bcdir+'.f90','w') # set to empty							
										loopst   = open(incPATH+'include_BCLoops_'+k+'_'+bcdir+'.f90','w')  # set to empty
			
										indrangei = 'indbc(1)=1\nindbc(2)=1\n'
										indrangej = 'indbc(3)=1\nindbc(4)=1\n'
										indrangek = 'indbc(5)=1\nindbc(6)=1\n'
							
										if dim == 2:
											if bcdir[0] == 'i':
												indrangej = 'indbc(3)=1\nindbc(4)=ny\n'
											elif bcdir[0] == 'j':	
												indrangei = 'indbc(1)=1\nindbc(2)=nx\n'
										elif dim == 3:
											if   bcdir[0] == 'i':
												indrangej = 'indbc(3)=1\nindbc(4)=ny\n'
												indrangek = 'indbc(5)=1\nindbc(6)=nz\n'
											elif bcdir[0] == 'j':	
												indrangei = 'indbc(1)=1\nindbc(2)=nx\n'	
												indrangek = 'indbc(5)=1\nindbc(6)=nz\n'
											elif bcdir[0] == 'k':	
												indrangei = 'indbc(1)=1\nindbc(2)=nx\n'	
												indrangej = 'indbc(3)=1\nindbc(4)=ny\n'
		
		
										indrange = {'i':indrangei,
										            'j':indrangej,
										            'k':indrangek}	
										
										updateqbc = False
										
										if k[0:5] == 'varbc': updateqbc = True
										
										compute_storedbc(var2process[k],Stencil,Order,loopst,locst,bcdir,update=update,updateqbc=updateqbc, rhs=rhs)
			
										create_bcsubroutine(st,k+'_faces_'+bcdir,
														 'include_BClocVar_'+k+'_'+bcdir+'.f90',
														 'include_BCLoops_'+k+'_' +bcdir+'.f90',indrange)


								

				else: # new bc name
					bc_info[bcname] = {}
					for bcdir in setbc[1][bcname].keys():
						bc_info[bcname][bcdir] = []
						for bctype in setbc[1][bcname][bcdir]:
							bc_info[bcname][bcdir].append(bctype)	   	


						if bcdir == cornerperm[bcdir]:

							# generate static/dynamic stored variables:			
								
							# from genRhs import varstored
							
							# generate static/dynamic stored variables:			
	
							staticstored  = {}
							dynamicstored = {}
							staticvarbc   = {}
							dynamicvarbc  = {}

							for var in varstored:
								if varstored[var]['static']:
									staticstored[var] = varstored[var]['symb']
								else:
									dynamicstored[var] = varstored[var]['symb']

							addvarbc = False

							for var in varbc:
								if 'face' in varbc[var]:
									if varbc[var]['face'] == bcdir: addvarbc = True
								elif 'edge' in varbc[var]:
									if varbc[var]['edge'] == bcdir: addvarbc = True
								else:
									addvarbc = False

								if addvarbc:	
									if varbc[var]['static']:
										staticvarbc[var] = varbc[var]['symb']
									else:
										dynamicvarbc[var] = varbc[var]['symb']		

								addvarbc = False		

								

							var2process = {'storedstatic':staticstored, 
										   'stored'      :dynamicstored,
										   'varbcstatic' :staticvarbc,
										   'varbc'       :dynamicvarbc}
	
							for k in var2process:
								if var2process[k] != {}:
									st       = open(incPATH+'bcsrc'+k+'_'+bcdir+'.for','w') # set to empty													
									locst    = open(incPATH+'include_BClocVar_'+k+'_'+bcdir+'.f90','w') # set to empty							
									loopst   = open(incPATH+'include_BCLoops_'+k+'_' +bcdir+'.f90','w')  # set to empty
		
									indrangei = 'indbc(1)=1\nindbc(2)=1\n'
									indrangej = 'indbc(3)=1\nindbc(4)=1\n'
									indrangek = 'indbc(5)=1\nindbc(6)=1\n'
						
									if dim == 2:
										if bcdir[0] == 'i':
											indrangej = 'indbc(3)=1\nindbc(4)=ny\n'
										elif bcdir[0] == 'j':	
											indrangei = 'indbc(1)=1\nindbc(2)=nx\n'
									elif dim == 3:
										if   bcdir[0] == 'i':
											indrangej = 'indbc(3)=1\nindbc(4)=ny\n'
											indrangek = 'indbc(5)=1\nindbc(6)=nz\n'
										elif bcdir[0] == 'j':	
											indrangei = 'indbc(1)=1\nindbc(2)=nx\n'	
											indrangek = 'indbc(5)=1\nindbc(6)=nz\n'
										elif bcdir[0] == 'k':	
											indrangei = 'indbc(1)=1\nindbc(2)=nx\n'	
											indrangej = 'indbc(3)=1\nindbc(4)=ny\n'
	
	
									indrange = {'i':indrangei,
									            'j':indrangej,
									            'k':indrangek}	

									updateqbc = False
										
									if k[0:5] == 'varbc': updateqbc = True

									compute_storedbc(var2process[k],Stencil,Order,loopst,locst,bcdir,update=update,updateqbc=updateqbc,rhs=rhs)
		
									create_bcsubroutine(st,k+'_faces_'+bcdir,
													 'include_BClocVar_'+k+'_'+bcdir+'.f90',
													 'include_BCLoops_'+k+'_' +bcdir+'.f90',indrange)	
									
								
							
			dir2gen = []
			dir2gen.append(bcdir)					

		rhs.bc_info[0]	= bc_info 
#
#  recover bc info for generation and update bc_info with new data...
#  ... then start generating bc

		if dim == 3:	
				dir2gen2 = ['i1','imax','j1','jmax','k1','kmax']
		elif dim == 2:
				dir2gen2 = ['i1','imax','j1','jmax']
		elif dim == 1:
				dir2gen2 = ['i1','imax']
		
		edgeBCs = []
		for d1 in dir2gen2:
			for d2 in dir2gen2:
				if d1 != d2:
					edgeBCs.append(d1+d2)

		
		for dirBC in dir2gen:
			if dirBC not in edgeBCs:
				DirDic = {'i':{'dir':None,'indbc':None},
				          'j':{'dir':None,'indbc':None},
				          'k':{'dir':None,'indbc':None}}
	
	
				DirDic[dirBC[0]] = {'dir': dirBC , 'indbc' : None}
	
				bcnum = bcnum + 1
	
				indrangei = 'indbc(1)=1\nindbc(2)=1\n'
				indrangej = 'indbc(3)=1\nindbc(4)=1\n'
				indrangek = 'indbc(5)=1\nindbc(6)=1\n'
	
				if dim == 2:
					if dirBC[0] == 'i':
						indrangej = 'indbc(3)=1\nindbc(4)=ny\n'
					elif dirBC[0] == 'j':	
						indrangei = 'indbc(1)=1\nindbc(2)=nx\n'
				elif dim == 3:
					if   dirBC[0] == 'i':
						indrangej = 'indbc(3)=1\nindbc(4)=ny\n'
						indrangek = 'indbc(5)=1\nindbc(6)=nz\n'
					elif dirBC[0] == 'j':	
						indrangei = 'indbc(1)=1\nindbc(2)=nx\n'	
						indrangek = 'indbc(5)=1\nindbc(6)=nz\n'
					elif dirBC[0] == 'k':	
						indrangei = 'indbc(1)=1\nindbc(2)=nx\n'	
						indrangej = 'indbc(3)=1\nindbc(4)=ny\n'
	
				indrange = {'i':indrangei,
							'j':indrangej,
							'k':indrangek}			
	
				idrhs = {}		
				idflt = {}	
	
				idrhs['i1']   = 'idrhs(1)=idarray(1)\n'
				idrhs['imax'] = 'idrhs(2)=idarray(2)\n'
				idrhs['j1']   = 'idrhs(3)=idarray(3)\n'
				idrhs['jmax'] = 'idrhs(4)=idarray(4)\n'
				idrhs['k1']   = 'idrhs(5)=idarray(5)\n'
				idrhs['kmax'] = 'idrhs(6)=idarray(6)\n'
	
				idflt['i1']   = 'idflt(1)=idarray(1)\n'
				idflt['imax'] = 'idflt(2)=idarray(2)\n'
				idflt['j1']   = 'idflt(3)=idarray(3)\n'	
				idflt['jmax'] = 'idflt(4)=idarray(4)\n'						
				idflt['k1']   = 'idflt(5)=idarray(5)\n'	
				idflt['kmax'] = 'idflt(6)=idarray(6)\n'	
				
				if update == False:
					# bcnum = bcnum + 1
					if setbc[0]: 
	
						if dirBC in rhs.bc_info[1]:
							rhs.bc_info[1][dirBC].append(bcnum)
						else:
							rhs.bc_info[1][dirBC] = []
							rhs.bc_info[1][dirBC].append(bcnum)	
						
						localvar = open(incPATH+'include_PhyBClocVar_'+bcname+'_'+dirBC+'_'+bctype+'.f90','w')
						output   = open(incPATH+'include_PhyBCLoops_'+bcname+'_'+dirBC+'_'+bctype+'.f90' ,'w')
						rhs_for  = open(incPATH+'PhyBC'+bcname+'_'+dirBC+'_'+bctype+'.for','a+')	
	
						create_bcsubroutine(rhs_for,bcname+'bc'+'_faces_'+str(dirBC)+'_'+str(bcnum),
										 'include_PhyBClocVar_'+bcname+'_'+dirBC+'_'+bctype+'.f90',
										 'include_PhyBCLoops_'+bcname+'_'+dirBC+'_'+bctype+'.f90' ,indrange, phy=bctype)
	
					else:
						
						outlayer   = []
						locallayer = []
	
						for layer in range(0,hlo_rhs): #BC layers  
	
							localvar = open(incPATH+'include_BClocVar_'+dirBC+'_'+str(layer)+'.f90','w')
							output   = open(incPATH+'include_BCLoops_'+dirBC+'_'+str(layer)+'.f90' ,'w')	
							rhs_for  = open(incPATH+'bcsrc'+dirBC+'_'+str(layer)+'.for','a+')
	
							outlayer.append(output)
							locallayer.append(localvar)
	
							create_bcsubroutine(rhs_for,'_faces_'+str(dirBC)+'_'+str(layer),
										 'include_BClocVar_'+dirBC+'_'+str(layer)+'.f90',
										 'include_BCLoops_'+dirBC+'_'+str(layer)+'.f90',indrange)
	
	
				else:	
	
					if setbc[0]:
						exception('Update not supported for physical BCs',message='error')
		
					
					else:
	
						outlayer   = []
						locallayer = []
	
						for layer in range(0,hlo_rhs): #BC layers  
	
							localvar       = open(incPATH+'include_BClocVar_'+dirBC+'_'+str(layer)+'.f90','a+')
							output         = open(incPATH+'include_BCLoops_'+dirBC+'_'+str(layer)+'.f90' ,'a+')	
	
							outlayer.append(output)
							locallayer.append(localvar)
	
				
				history={'x':{'var':[],'symb':[]},'y':{'var':[],'symb':[]},'z':{'var':[],'symb':[]}}
				history2={'x':[[],[]],'y':[[],[]],'z':[[],[]]}
				dhistory={'x':[[],[],[],[]],'y':[[],[],[],[]],'z':[[],[],[],[]]}
		
				indi = 'i'
				indj = 'j'
				indk = 'k'
		
				indiri = indi
				indirj = indj
				indirk = indk
		
				if dirBC == 'i1': 
					indi = '1-'+str(hlo_rhs)
				if dirBC == 'imax': 
					indi = 'nx+'+str(hlo_rhs)
		
				if dirBC == 'j1': 
					indj = '1-'+str(hlo_rhs)
				if dirBC == 'jmax': 
					indj = 'ny+'+str(hlo_rhs)
		
		
				if dirBC == 'k1': 
					indk = '1-'+str(hlo_rhs)
				if dirBC == 'kmax': 
					indk = 'nz+'+str(hlo_rhs)
		
			    # generate sources from equations:
	
				if setbc[0]:
		
					locvar = []
					
					layer = 0
	
					DirDic[dirBC[0]]['indbc'] = layer
	
					if bctype == 'q':
						updateq = True
					else:
						updateq = False
	
					gen_eqns_bc(Eqns,output,localvar,
					            rhsname,Order=Order,Stencil=Stencil,
					            indi    =indi,indj =indj,indk =indk,
					            DirDic  = DirDic,
					            vname   = vname,
					            update  = update,
					            updateq = updateq)				
					
				else:	
	
					for layer in range(0,hlo_rhs): #BC layers 
	
						DirDic[dirBC[0]]['indbc'] = layer
	
						localvar = locallayer[layer]
						output   = outlayer[layer]
		
						gen_eqns_bc(Eqns,output,localvar,
						            rhsname,Order=Order,Stencil=Stencil,
						            indi    =indi,indj =indj,indk =indk,
						            DirDic  = DirDic,
						            vname   = vname,
						            update  = update)
#
#		Generate edges	
#
		if not setbc[0]: 

			if stored:



				staticstored  = {}
				dynamicstored = {}
				staticvarbc   = {}
				dynamicvarbc  = {}

				for var in varstored:
					if varstored[var]['static']:
						staticstored[var] = varstored[var]['symb']
					else:
						dynamicstored[var] = varstored[var]['symb']	

	
				var2process = {'storedstatic':staticstored, 
							   'stored'      :dynamicstored,
							   'varbcstatic' :staticvarbc,
							   'varbc'       :dynamicvarbc}

				slcbc_stored = {} 


			if len(rhs.bc_info[1].values()) != 0:
				bcnum = max(max(rhs.bc_info[1].values()))
			else:
				bcnum = 0	

			if dim == 3:	
				dir2gen = ['i1','imax','j1','jmax','k1','kmax']
			elif dim == 2:
				dir2gen = ['i1','imax','j1','jmax']
			elif dim == 1:
				dir2gen = ['i1','imax']

			DirDic = {'i':{'dir':None,'indbc':None},
			          'j':{'dir':None,'indbc':None},
			          'k':{'dir':None,'indbc':None}}

			idrhs = {}		

			idrhs['i1']   = 'idrhs(1)=idarray(1)\n'
			idrhs['imax'] = 'idrhs(2)=idarray(2)\n'
			idrhs['j1']   = 'idrhs(3)=idarray(3)\n'	
			idrhs['jmax'] = 'idrhs(4)=idarray(4)\n'						
			idrhs['k1']   = 'idrhs(5)=idarray(5)\n'	
			idrhs['kmax'] = 'idrhs(6)=idarray(6)\n'

			indi = 'i'
			indj = 'j'
			indk = 'k'

			indiri = indi
			indirj = indj
			indirk = indk         
		
			edone = []

			for dir1 in dir2gen:
				DirDic[dir1[0]]['dir'] = dir1	
				for dir2 in dir2gen:
						if dir2[0] != dir1[0]:

							if ((dir2+dir1 not in edone) and (dir1+dir2 not in edone)):
								if ((dir1+dir2 not in rhs.bc_info[1]) and (dir2+dir1 not in rhs.bc_info[1])) or update:

									if not update:

										bcnum = bcnum + 1
			 
										rhs.bc_info[1][dir1+dir2] = []
										rhs.bc_info[1][dir1+dir2].append(bcnum)
				
									edone.append(dir1+dir2)

									DirDic[dir2[0]]['dir'] = dir2									
									for layer1 in range(0,hlo_rhs): #BC layers dir1

										DirDic[dir1[0]]['indbc'] = layer1	

#										
#										Create first layer of subroutines (for layer1), that includes calls to second layer
#										
										
										if not update:
											edgecallname = '_edges_'+dir1+'_'+dir2+'_'+str(layer1)
											efname       = open(incPATH+'bcsrc_edgescall_'+dir1+'_'+dir2+'_'+str(layer1)+'.for','w')

										# stored call name + open bcsrc static/dyna
										if stored:	
											edgecallname_stored = {}
											efname_stored       = {}
											bcall_stored        = {}      

											var2process['varbcstatic'] = {}
											var2process['varbc']       = {}	
											staticvarbc = {}
											dynamicvarbc= {}

											addvarbc = False
											for var in varbc:
												if 'face' in varbc[var]:
													if varbc[var]['face'] == dir1 : addvarbc = True
												elif 'edge' in varbc[var]:
													if varbc[var]['edge'] == dir1+dir2: addvarbc = True
												else:
													addvarbc = False
								
												if addvarbc:	
													if varbc[var]['static']:
														staticvarbc[var] = varbc[var]['symb']
													else:
														dynamicvarbc[var] = varbc[var]['symb']

												addvarbc = False			
	
											var2process['varbcstatic'] = staticvarbc
											var2process['varbc']       = dynamicvarbc
												
	
											for k in var2process:
												if var2process[k] != {}:
													edgecallname_stored[k] = k+'_edges_'+dir1+'_'+dir2+'_'+str(layer1)
													efname_stored[k]       = open(incPATH+'bcsrc_'+k+'_edgescall_'+dir1+'_'+dir2+'_'+str(layer1)+'.for','w')
													bcall_stored[k]        = '\n'

										bcall = '\n'
										for layer2 in range(0,hlo_rhs): #BC layers dir2

#											
#											Start generating boundary equations for points (layer 1, layer 2)
#											
			 	
											DirDic[dir2[0]]['indbc'] = layer2
											localvar = open(incPATH+'include_edges_'+dir1+'_'+dir2+'_'+str(layer1)+'_'+str(layer2)+'_BClocVar.f90','a+')
											output   = open(incPATH+'include_edges_'+dir1+'_'+dir2+'_'+str(layer1)+'_'+str(layer2)+'_BCLoops.f90' ,'a+')

											# stored open include edges static/dynamic
											if stored:	
												localvar_stored = {}
												output_stored       = {} 

												var2process['varbcstatic'] = {}
												var2process['varbc']       = {}	
												staticvarbc = {}
												dynamicvarbc= {}
	
												if layer1*layer2 == 0:
													addvarbc = False
													for var in varbc:
														if 'face' in varbc[var]:
															if varbc[var]['face'] == dir1 and layer1==0: addvarbc = True
															if varbc[var]['face'] == dir2 and layer2==0: addvarbc = True

														elif 'edge' in varbc[var]:
															if varbc[var]['edge'] == dir1+dir2: addvarbc = True
														else:
															addvarbc = False
									
														if addvarbc:	
															if varbc[var]['static']:
																staticvarbc[var] = varbc[var]['symb']
															else:
																dynamicvarbc[var] = varbc[var]['symb']
														addvarbc = False			
		
													var2process['varbcstatic'] = staticvarbc
													var2process['varbc']       = dynamicvarbc
												else:
													var2process['varbcstatic'] = {}
													var2process['varbc']       = {}	
	
	
												for k in var2process:
													if var2process[k] != {}:
														localvar_stored[k] = open(incPATH+'include_'+k+'_edges_'+dir1+'_'+dir2+'_'+str(layer1)+'_'+str(layer2)+'_BClocVar.f90','a+')
														output_stored[k]   = open(incPATH+'include_'+k+'_edges_'+dir1+'_'+dir2+'_'+str(layer1)+'_'+str(layer2)+'_BCLoops.f90' ,'a+')
			 
											if DirDic['i']['dir'] == 'i1': 
												indi = '1-'+str(hlo_rhs)
											if DirDic['i']['dir'] == 'imax': 
												indi = 'nx+'+str(hlo_rhs)
									 
											if DirDic['j']['dir'] == 'j1': 
												indj = '1-'+str(hlo_rhs)
											if DirDic['j']['dir'] == 'jmax': 
												indj = 'ny+'+str(hlo_rhs)
									 
											if DirDic['k']['dir'] == 'k1': 
												indk = '1-'+str(hlo_rhs)
											if DirDic['k']['dir'] == 'kmax': 
												indk = 'nz+'+str(hlo_rhs) 
				 
											gen_eqns_bc(Eqns,output,localvar,
									            rhsname,Order=Order,Stencil=Stencil,
									            indi    = indi,indj = indj,indk = indk,
									            DirDic  = DirDic,
									            vname   = vname,
									            update  = update)

											# stored generate fortran src static/dyna
											if stored:
												for k in var2process:
													if var2process[k] != {}:
														vname_st = {}
														for v in var2process[k]:
															vname_st[v] = v

														updateqbc = False
														updatest  = True
												
														if k[0:5] == 'varbc': 
															updateqbc = True
															updatest  = False

														gen_eqns_bc(var2process[k],output_stored[k],localvar_stored[k],
												            vname_st,Order=Order,Stencil=Stencil,
												            indi     = indi,indj = indj,indk = indk,
												            DirDic   = DirDic,
												            vname    = vname_st,
												            update   = False,
												            updateqbc= updateqbc,
												            updatest = updatest)

#														
#											Create the second layer of subroutines with the newly generated eqns (for layer 2)	            											
#											

											if not update:
			 
												indrangei = 'indbc(1)=1\nindbc(2)=1\n'
												indrangej = 'indbc(3)=1\nindbc(4)=1\n'
												indrangek = 'indbc(5)=1\nindbc(6)=1\n'
				 
												edir = 'ijk'
												edir = edir.replace(dir1[0],'').replace(dir2[0],'')
				 
												if dim == 3:
													if   edir == 'i':
														indrangei = 'indbc(1)=1\nindbc(2)=nx\n'
													elif edir == 'j':
														indrangej = 'indbc(1)=1\nindbc(2)=ny\n'
													elif edir == 'k':
														indrangek = 'indbc(1)=1\nindbc(2)=nz\n'	
															# 
				 
												indrange = {'i':indrangei,
															'j':indrangej,
															'k':indrangek}
	 
		 										
												bcedges  = open(incPATH+'bcsrc_edges_'+dir1+'_'+dir2+'_'+str(layer1)+'_'+str(layer2)+'.for','w')
												create_bcsubroutine(bcedges,'_edges_'+dir1+'_'+dir2+'_'+str(layer1)+'_'+str(layer2),
																	'include_edges_'+dir1+'_'+dir2+'_'+str(layer1)+'_'+str(layer2)+'_BClocVar.f90',
																	'include_edges_'+dir1+'_'+dir2+'_'+str(layer1)+'_'+str(layer2)+'_BCLoops.f90',indrange)
												bcedges.close()
												bcedges  = open(incPATH+'bcsrc_edges_'  +dir1+'_'+dir2+'_'+str(layer1)+'_'+str(layer2)+'.for','r')
												cname = bcedges.readlines()[8][10:]
												bcall    = bcall + '      call '+ cname

											# stored create subroutines static/dynamic
											if stored:
												bcedges_stored = {}
												for k in var2process:
													if var2process[k] != {}:

														bcedges_stored[k]  = open(incPATH+'bcsrc_'+k+'_edges_'+dir1+'_'+dir2+'_'+str(layer1)+'_'+str(layer2)+'.for','w')														
														create_bcsubroutine(bcedges_stored[k],k+'_edges_'+dir1+'_'+dir2+'_'+str(layer1)+'_'+str(layer2),
																'include_'+k+'_edges_'+dir1+'_'+dir2+'_'+str(layer1)+'_'+str(layer2)+'_BClocVar.f90',
																'include_'+k+'_edges_'+dir1+'_'+dir2+'_'+str(layer1)+'_'+str(layer2)+'_BCLoops.f90' ,indrange)
														bcedges_stored[k].close()
														bcedges_stored[k]  = open(incPATH+'bcsrc_'+k+'_edges_'+dir1+'_'+dir2+'_'+str(layer1)+'_'+str(layer2)+'.for','r')
														cname = bcedges_stored[k].readlines()[8][10:]
														bcall_stored[k]    = bcall_stored[k] + '      call '+ cname

#										
#										Create calls of the second layer of subroutines
#														

										if not update:		
											create_bccalls(efname,edgecallname,bcall)
											efname.close()

										# stored add calls static/dynamic
										if stored:
											addvarbc = False
											for var in varbc:
												if 'face' in varbc[var]:
													if varbc[var]['face'] == dir1: addvarbc = True
												else:
													addvarbc = False
							
												if addvarbc:	
													if varbc[var]['static']:
														staticvarbc[var] = varbc[var]['symb']
													else:
														dynamicvarbc[var] = varbc[var]['symb']
												addvarbc = False			

											var2process['varbcstatic'] = staticvarbc
											var2process['varbc']       = dynamicvarbc
											for k in var2process:
												if var2process[k] != {}:
													create_bccalls(efname_stored[k],edgecallname_stored[k],bcall_stored[k])
													efname_stored[k].close()
													
#       
#       Generates physical BC edges
#
		else: # setbc[0] True (Physical BC)

			# if dim == 3:	
			# 	dir2gen2 = ['i1','imax','j1','jmax','k1','kmax']
			# elif dim == 2:
			# 	dir2gen2 = ['i1','imax','j1','jmax']
			# elif dim == 1:
			# 	dir2gen2 = ['i1','imax']

			DirDic = {'i':{'dir':None,'indbc':None},
			          'j':{'dir':None,'indbc':None},
			          'k':{'dir':None,'indbc':None}}
			indi = 'i'
			indj = 'j'
			indk = 'k'

			indiri = indi
			indirj = indj
			indirk = indk         
		
			edone = []

			# test ='i1j1'

			# edgeBCs = []
			# for d1 in dir2gen2:
			# 	for d2 in dir2gen2:
			# 		if d1 != d2:
			# 			edgeBCs.append(d1+d2)

#
#           Generate faces BCs (i1,j1,k1,imax,...) and subsequent edges combinations (i1j1,...). 
#           Physical BC at layer 0 0 for (i1j1,...) can be imposed separately from the faces eqns (see below)
#			
			for dir1 in dir2gen:
				if dir1 not in edgeBCs:
					DirDic[dir1[0]]['dir'] = dir1	
					for dir2 in dir2gen2:
							if dir2[0] != dir1[0]:
	
								if ((dir2+dir1 not in edone) and (dir1+dir2 not in edone)):
	
										edone.append(dir1+dir2)
										DirDic[dir2[0]]['dir']   = dir2
										DirDic[dir1[0]]['indbc'] = 0
#
#										Create first layer1 subroutine, that includes calls to second layer
#
										edgecallname =                       '_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_'+dir2+'_0'
										efname       = open(incPATH+'bcsrc_edgescall_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_'+dir2+'_0.for','w')
	
										bcall = '\n'
	
										for layer2 in range(0,hlo_rhs): #BC layers dir2
				 
											DirDic[dir2[0]]['indbc'] = layer2
	
#
#											Generate boundary eqns for point (layer1,layer2)	
#
											localvar = open(incPATH+'include_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_'+dir2+'_0'+'_'+str(layer2)+'_BClocVar.f90','a+')
											outputPhy   = open(incPATH+'include_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_'+dir2+'_0'+'_'+str(layer2)+'_BCLoops.f90' ,'a+')
	
											if DirDic['i']['dir'] == 'i1': 
												indi = '1-'+str(hlo_rhs)
											if DirDic['i']['dir'] == 'imax': 
												indi = 'nx+'+str(hlo_rhs)
										
											if DirDic['j']['dir'] == 'j1': 
												indj = '1-'+str(hlo_rhs)
											if DirDic['j']['dir'] == 'jmax': 
												indj = 'ny+'+str(hlo_rhs)
										
											if DirDic['k']['dir'] == 'k1': 
												indk = '1-'+str(hlo_rhs)
											if DirDic['k']['dir'] == 'kmax': 
												indk = 'nz+'+str(hlo_rhs) 
	
											if bctype == 'q':
												updateq = True
											else:
												updateq = False	
					 
											gen_eqns_bc(Eqns,outputPhy,localvar,
										        rhsname,Order=Order,Stencil=Stencil,
										        indi    = indi,indj = indj,indk = indk,
										        DirDic  = DirDic,
										        vname   = vname,
										        update  = update,
										        updateq = updateq)
	
#
#											Create second layer subroutine with the newly generates eqns (and accumulate calls)	
#			 
											indrangei = 'indbc(1)=1\nindbc(2)=1\n'
											indrangej = 'indbc(3)=1\nindbc(4)=1\n'
											indrangek = 'indbc(5)=1\nindbc(6)=1\n'
					 
											edir = 'ijk'
											edir = edir.replace(dir1[0],'').replace(dir2[0],'')
					 
											if dim == 3:
												if   edir == 'i':
													indrangei = 'indbc(1)=1\nindbc(2)=nx\n'
												elif edir == 'j':
													indrangej = 'indbc(1)=1\nindbc(2)=ny\n'
												elif edir == 'k':
													indrangek = 'indbc(1)=1\nindbc(2)=nz\n'	
													
					 						
											indrange = {'i':indrangei,
														'j':indrangej,
														'k':indrangek}
		 
			 								
											bcedges  = open(incPATH+'bcsrc_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_'+dir2+'_0'+'_'+str(layer2)+'.for','w')
											
											create_bcsubroutine(bcedges,'_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_'+dir2+'_0'+'_'+str(layer2),
																'include_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_'+dir2+'_0'+'_'+str(layer2)+'_BClocVar.f90',
																'include_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_'+dir2+'_0'+'_'+str(layer2)+'_BCLoops.f90',indrange)
											bcedges.close()
											bcedges  = open(incPATH+'bcsrc_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_'+dir2+'_0'+'_'+str(layer2)+'.for','r')
											cname = bcedges.readlines()[8][10:]
											bcall    = bcall + '      call '+ cname
	
											
										create_bccalls(efname,edgecallname,bcall)
										efname.close()

				else: # dir1 is a new physical BC related to a 0 0 layer (i1j1,...)

					# recover direction from dirname (another brute force approach)

					for a in dir2gen2:
						for b in dir2gen2:
							if (a+b == dir1) or (b+a == dir1):
								d1 = a
								d2 = b

					DirDic[d1[0]]['dir']   = d1	
					DirDic[d2[0]]['dir']   = d2
					DirDic[d1[0]]['indbc'] = 0
					DirDic[d2[0]]['indbc'] = 0

					localvar = open(incPATH+'include_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_0_0'+'_BClocVar.f90','a+')
					outputPhy   = open(incPATH+'include_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_0_0'+'_BCLoops.f90' ,'a+')
					if DirDic['i']['dir'] == 'i1': 
						indi = '1-'+str(hlo_rhs)
					if DirDic['i']['dir'] == 'imax': 
						indi = 'nx+'+str(hlo_rhs)
				
					if DirDic['j']['dir'] == 'j1': 
						indj = '1-'+str(hlo_rhs)
					if DirDic['j']['dir'] == 'jmax': 
						indj = 'ny+'+str(hlo_rhs)
				
					if DirDic['k']['dir'] == 'k1': 
						indk = '1-'+str(hlo_rhs)
					if DirDic['k']['dir'] == 'kmax': 
						indk = 'nz+'+str(hlo_rhs) 
					if bctype == 'q':
						updateq = True
					else:
						updateq = False	

					gen_eqns_bc(Eqns,outputPhy,localvar,
				        rhsname,Order=Order,Stencil=Stencil,
				        indi    = indi,indj = indj,indk = indk,
				        DirDic  = DirDic,
				        vname   = vname,
				        update  = update,
				        updateq = updateq)
#
#						Create second layer subroutine with the newly generates eqns (and accumulate calls)	
#
					indrangei = 'indbc(1)=1\nindbc(2)=1\n'
					indrangej = 'indbc(3)=1\nindbc(4)=1\n'
					indrangek = 'indbc(5)=1\nindbc(6)=1\n'

					edir = 'ijk'
					edir = edir.replace(d1[0],'').replace(d2[0],'')

					# if dim == 3:
					# 	if   edir == 'i':
					# 		indrangei = 'indbc(1)=1\nindbc(2)=nx\n'
					# 	elif edir == 'j':
					# 		indrangej = 'indbc(1)=1\nindbc(2)=ny\n'
					# 	elif edir == 'k':
					# 		indrangek = 'indbc(1)=1\nindbc(2)=nz\n'	
							

					indrange = {'i':indrangei,
								'j':indrangej,
								'k':indrangek}
						
					bcedges  = open(incPATH+'bcsrc_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_0_0'+'.for','w')
					
					create_bcsubroutine(bcedges,'_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_0_0',
										'include_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_0_0'+'_BClocVar.f90',
										'include_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_0_0'+'_BCLoops.f90',indrange)
					bcedges.close()				



def genBC_calls(rhs):

	from genRhs import incPATH
	
	dim = rhs.dim
	wp  = rhs.wp
	varname = rhs.varname
	varstored = rhs.varstored
	varsolved = rhs.varsolved
	varbc     = rhs.varbc
	
	hlo_rhs = rhs.hlo_rhs
	stencil_rhs = rhs.stencil
	order_rhs   = rhs.order

	bc_info = rhs.bc_info[0]

# Extract bcdir:

	bcdir_all = []
	for bcname in bc_info:
		bcdir_all = bcdir_all + list(bc_info[bcname].keys())

	bcdir_all = sorted(bcdir_all)

# Extract phybc details:

	bcphy_all = {}

	for bcdir in bcdir_all:
		bcphy_all[bcdir] = {}
		for bcname in bc_info:
			if bcdir in list(bc_info[bcname].keys()):
				bcphy_all[bcdir][bcname] = bc_info[bcname][bcdir] 

	staticstored  = {}
	dynamicstored = {}
	staticvarbc   = {}
	dynamicvarbc  = {}

	for var in varstored:
		if varstored[var]['static']:
			staticstored[var] = varstored[var]['symb']
		else:
			dynamicstored[var] = varstored[var]['symb']	

	var2process = {'storedstatic':staticstored, 
				   'stored'      :dynamicstored,
				   'varbcstatic' :staticvarbc,
				   'varbc'       :dynamicvarbc}

	slcbc_stored  = {} 
	efname_stored = {}

	DirDic = {'i':{'dir':None,'indbc':None},
	          'j':{'dir':None,'indbc':None},
	          'k':{'dir':None,'indbc':None}}

	idrhs = {}		

	idrhs['i1']   = 'idrhs(1)=idarray(1)\n'
	idrhs['imax'] = 'idrhs(2)=idarray(2)\n'
	idrhs['j1']   = 'idrhs(3)=idarray(3)\n'	
	idrhs['jmax'] = 'idrhs(4)=idarray(4)\n'						
	idrhs['k1']   = 'idrhs(5)=idarray(5)\n'	
	idrhs['kmax'] = 'idrhs(6)=idarray(6)\n'
	indi = 'i'
	indj = 'j'
	indk = 'k'

	indiri = indi
	indirj = indj
	indirk = indk         
	

	bcdone = {}
	bcnum  = 0

# # ADD EDGES :

	if dim == 3:	
		dirlist = ['i1','imax','j1','jmax','k1','kmax']
	elif dim == 2:
		dirlist = ['i1','imax','j1','jmax']
	elif dim == 1:
		dirlist = ['i1','imax']


	edgeBCs = []
	for d1 in dirlist:
		for d2 in dirlist:
			if d1 != d2:
				edgeBCs.append(d1+d2)

	from itertools import permutations

	edone = []
	for dir1 in bcdir_all:
		if dir1 not in edgeBCs:
			for dir2 in bcdir_all:
				if dir2 not in edgeBCs:
					if dir2[0] != dir1[0]:
	
						if ((dir2+dir1 not in edone) and (dir1+dir2 not in edone)):
							if ((dir1+dir2 not in bcdone) and (dir2+dir1 not in bcdone)):
	
								bcnum = bcnum + 1
	
								bcdone[dir1+dir2] = []
								bcdone[dir1+dir2].append(bcnum)		
	
								slcbc = {'rhs': open(incPATH+'select_phybc_rhs.f90' ,'a+')}
	
								typebc = []
								for bcname in bcphy_all[dir1]:
									for bctype in bcphy_all[dir1][bcname]:
										typebc.append(bctype)	
	
								if 'rhs' not in typebc:  # By default we extend interiror eqns for rhs at the edges
										slcbc['rhs'].write('CASE ('+str(bcnum)+')\n')	
	
								for bctype in typebc:
									slcbc[bctype] = open(incPATH+'select_phybc_'+bctype+'.f90' ,'a+')
									slcbc[bctype].write('CASE ('+str(bcnum)+')\n')
	
								# rhs bc edge:
	
								for layer1 in range(0,hlo_rhs): #BC layers dir1
									efname       = open(incPATH+'bcsrc_edgescall_'+dir1+'_'+dir2+'_'+str(layer1)+'.for','r')
									slcbc['rhs'].write('      call '+efname.readlines()[8][10:])
								slcbc['rhs'].write('      '+idrhs[dir1])
								slcbc['rhs'].write('      '+idrhs[dir2])
	
	
								fephy = {}
	
	
								if dir1 in bcphy_all:
									# slcbc[bctype] = open(incPATH+'select_phybc_'+bctype+'.f90' ,'a+')
									# slcbc[bctype].write('CASE ('+str(bcnum)+')\n')
									for bcname in bcphy_all[dir1]:
										for bctype in bcphy_all[dir1][bcname]:
	
											fephy[bctype] = open(incPATH+'bcsrc_edgescall_PhyBC_'+bcname+'_'+bctype+'_'+dir1+'_'+dir2+'_0.for','r')
											slcbc[bctype].write('      call '+fephy[bctype].readlines()[8][10:])
	
								if dir2 in bcphy_all:
									for bcname in bcphy_all[dir2]:
										for bctype in bcphy_all[dir2][bcname]:
											if bctype not in slcbc:
												slcbc[bctype] = open(incPATH+'select_phybc_'+bctype+'.f90' ,'a+')
												slcbc[bctype].write('CASE ('+str(bcnum)+')\n')
											fephy[bctype] = open(incPATH+'bcsrc_edgescall_PhyBC_'+bcname+'_'+bctype+'_'+dir2+'_'+dir1+'_0.for','r')
											slcbc[bctype].write('      call '+fephy[bctype].readlines()[8][10:])

								perms = [''.join(p) for p in permutations(dir1+dir2)]
										
								if (dir1+dir2) in bcphy_all or (list(set(perms)&set(bcphy_all)) != [] ):
									
									dir1dir2 = list(set(perms)&set(bcphy_all))[0]

									for bcname in bcphy_all[dir1dir2]:
										for bctype in bcphy_all[dir1dir2][bcname]:
											if bctype not in slcbc:
												slcbc[bctype] = open(incPATH+'select_phybc_'+bctype+'.f90' ,'a+')
												slcbc[bctype].write('CASE ('+str(bcnum)+')\n')
											fephy[bctype] = open(incPATH+'bcsrc_edges_PhyBC_'+bcname+'_'+bctype+'_'+dir1dir2+'_0_0.for','r')
											slcbc[bctype].write('      call '+fephy[bctype].readlines()[8][10:])
	
								edone.append(dir1+dir2)
								

# # Stored OPEN select and static/dyn extraction	
							var2process['varbcstatic'] = {}
							var2process['varbc']       = {}	
							staticvarbc = {}
							dynamicvarbc= {}
							addvarbc = False
							for var in varbc:
								if 'face' in varbc[var]:
									if varbc[var]['face'] == dir1: addvarbc = True
								elif 'edge' in varbc[var]:
									if varbc[var]['edge'] == dir1+dir2: addvarbc = True
								else:
									addvarbc = False
							
								if addvarbc:	
									if varbc[var]['static']:
										staticvarbc[var] = varbc[var]['symb']
									else:
										dynamicvarbc[var] = varbc[var]['symb']
								addvarbc = False			
	
							var2process['varbcstatic'] = staticvarbc
							var2process['varbc']       = dynamicvarbc

							for k in var2process:
								if var2process[k] != {}:	

									slcbc_stored[k] = open(incPATH+'select'+k+'bc.f90','a+')
									slcbc_stored[k].write('CASE ('+str(bcnum)+')\n')
									# if k[0:5] == 'varbc': 										
									#    layerend = 1
									# else:
									   # layerend = hlo_rhs   
									layerend = hlo_rhs   
									for layer1 in range(0,layerend): #BC layers dir1
										efname_stored[k] = open(incPATH+'bcsrc_'+k+'_edgescall_'+dir1+'_'+dir2+'_'+str(layer1)+'.for','r')
										slcbc_stored[k].write('      call '+efname_stored[k].readlines()[8][10:])								

# #		extends bc filtering along the edges:		

							fltbnd = {'i1'  :{'dir': 'x','bound': 'idloop(1) = idarray(1)\n'},
									  'imax':{'dir': 'x','bound': 'idloop(2) = idarray(2)\n'},
									  'j1'  :{'dir': 'y','bound': 'idloop(3) = idarray(3)\n'},
									  'jmax':{'dir': 'y','bound': 'idloop(4) = idarray(4)\n'},
									  'k1'  :{'dir': 'z','bound': 'idloop(5) = idarray(5)\n'},
									  'kmax':{'dir': 'z','bound': 'idloop(6) = idarray(6)\n'}}
							try:	
								slcbd.close()	  	  
							except:
								None

							slcbd = open(incPATH+'selectfilterbc_'+fltbnd[dir1]['dir']+'.f90','a+')
							slcbd.write('\n CASE ('+str(bcnum)+')\n\n')
							slcbd.write('   '+fltbnd[dir2]['bound'])
							slcbd.close()

							slcbd = open(incPATH+'selectfilterbc_'+fltbnd[dir2]['dir']+'.f90','a+')
							slcbd.write('\n CASE ('+str(bcnum)+')\n\n')
							slcbd.write('   '+fltbnd[dir1]['bound'])
							slcbd.close()										

							slcbd = open(incPATH+'selectupdate_filterbc_'+fltbnd[dir1]['dir']+'.f90','a+')
							slcbd.write('\n CASE ('+str(bcnum)+')\n\n')
							slcbd.write('   '+fltbnd[dir2]['bound'])
							slcbd.close()

							slcbd = open(incPATH+'selectupdate_filterbc_'+fltbnd[dir2]['dir']+'.f90','a+')
							slcbd.write('\n CASE ('+str(bcnum)+')\n\n')
							slcbd.write('   '+fltbnd[dir1]['bound'])
							slcbd.close()						


# # ADDS NORMAL-TO-BOUNDARY FILTERS FOR BC DIRECTION:
	for dir1 in bcdir_all:
		if dir1 not in edgeBCs:
		
			bcnum = bcnum + 1
	
			bcdone[dir1] = []
			bcdone[dir1].append(bcnum)
						
			axes  = {'i':'x','j':'y','k':'z'}
	
			slcbcflt   = open(incPATH+'selectfilterbc_'+axes[dir1[0]]+'.f90'        ,'a+')
			slcbcfltup = open(incPATH+'selectupdate_filterbc_'+axes[dir1[0]]+'.f90' ,'a+')							
			slcbcflt.write('\n CASE ('+str(bcnum)+')\n\n')
			slcbcfltup.write('\n CASE ('+str(bcnum)+')\n\n')
	
			for layer in range(0,hlo_rhs):
				     genFilter(stencil_rhs,order_rhs,len(varsolved),dirBC=dir1,indbc=layer,fltbeg=0,rhs=rhs)
			
			up    = open(incPATH+'update_filterbc_'+axes[dir1[0]]+'.f90','r') # set to empty
			fltbc = open(incPATH+'filterbc_'+axes[dir1[0]]+'.f90'       ,'r') # set to empty	
	
			fltbclines = fltbc.readlines()
			uplines    = up.readlines()
	
			for l in fltbclines:
				slcbcflt.write('   '+l)
			for l in uplines:
				slcbcfltup.write('   '+l)		
		
# #		extends in-plane filtering along the bc:
			
			slcbcflt.close()
			slcbcfltup.close()
	
			if dir1[1:] == 'max':
				bndx  = 'idflt(2) = idarray(2)\n'							
				bndy  = 'idflt(4) = idarray(4)\n'
				bndz  = 'idflt(6) = idarray(6)\n'
			else:	
				bndx  = 'idflt(1) = idarray(1)\n'							
				bndy  = 'idflt(3) = idarray(3)\n'
				bndz  = 'idflt(5) = idarray(5)\n'									
	
			fltbnd = {'x':[{},{'y':bndx},{'y':bndx,'z':bndx}],
					  'y':[{},{'x':bndy},{'x':bndy,'z':bndy}],
					  'z':[{},{        },{'x':bndz,'y':bndz}]}
			bounds = fltbnd[axes[dir1[0]]][dim-1]
			
			for axebd in bounds:
				slcbd = open(incPATH+'selectfilterbc_'+axebd+'.f90','a+')
				slcbd.write('\n CASE ('+str(bcnum)+')\n\n')
				slcbd.write('   '+bounds[axebd])		
							   	
# # GENERATE STATIC/DYNAMIC STORED VARIABLES:			
			var2process['varbcstatic'] = {}
			var2process['varbc']       = {}	
			staticvarbc = {}
			dynamicvarbc= {}


			addvarbc = False
			for var in varbc:
				if 'face' in varbc[var]:
					if varbc[var]['face'] == dir1: addvarbc = True
				elif 'edge' in varbc[var]:
					if varbc[var]['edge'] == dir1: addvarbc = True
				else:
					addvarbc = False
			
				if addvarbc:	
					if varbc[var]['static']:
						staticvarbc[var] = varbc[var]['symb']
					else:
						dynamicvarbc[var] = varbc[var]['symb']	
				addvarbc = False		
	
			var2process['varbcstatic'] = staticvarbc
			var2process['varbc']       = dynamicvarbc
					
			for k in var2process:
				if var2process[k] != {}:	
	
					slcbc    = open(incPATH+'select'+k+'bc.f90','a+')
					st       = open(incPATH+'bcsrc'+k+'_'+dir1+'.for','r') # set to empty	
					
					slcbc.write('CASE ('+str(bcnum)+')\n')
					callname = st.readlines()[8][10:]
					slcbc.write('      call '+ callname)
	
# # ADD CALL TO NEW PHYSICAL BC :					
			
			slcbc = {}
			# for bctype in ['rhs','q']:
			# 	slcbc[bctype] = open(incPATH+'select_phybc_'+bctype+'.f90' ,'a+')
			# 	slcbc[bctype].write('CASE ('+str(bcnum)+')\n')
	
			if dir1 in bcphy_all:
				for bcname in bcphy_all[dir1]:
					for bctype in bcphy_all[dir1][bcname]:
	
						slcbc[bctype] = open(incPATH+'select_phybc_'+bctype+'.f90' ,'a+')
						slcbc[bctype].write('CASE ('+str(bcnum)+')\n')
														
					# slcbc      = open(incPATH+'select_phybc_'+bctype+'.f90' ,'a+')
					# rhs_for.close()
					
						phyBCread  = open(incPATH+'PhyBC'+bcname+'_'+dir1+'_'+bctype+'.for','r')	
						rhshlo = []
						for layer in range(0,hlo_rhs):
							rhshlo.append(open(incPATH+'bcsrc'+dir1+'_'+str(layer)+'.for','r'))
			
						# slcbc[bctype].write('CASE ('+str(bcnum)+')\n')
						if(bctype=='rhs'):
							for layer in range(0,hlo_rhs):
									slcbc[bctype].write('      call '+rhshlo[layer].readlines()[8][10:])
							slcbc[bctype].write('      '+idrhs[dir1])	
						slcbc[bctype].write('      call '+phyBCread.readlines()[9][10:])								
				if slcbc != {}:
					if 'rhs' not in slcbc:
							slcbc['rhs'] = open(incPATH+'select_phybc_rhs.f90' ,'a+')
							slcbc['rhs'].write('CASE ('+str(bcnum)+')\n')		
							rhshlo = []
							for layer in range(0,hlo_rhs):
								rhshlo.append(open(incPATH+'bcsrc'+dir1+'_'+str(layer)+'.for','r'))	
							for layer in range(0,hlo_rhs):
									slcbc['rhs'].write('      call '+rhshlo[layer].readlines()[8][10:])
							slcbc['rhs'].write('      '+idrhs[dir1])							
	
	rhs.bc_info[1] = bcdone					

def loop(beg,end):
	loop_create("begin",beg)
	loop_create("end"  ,end)

def loop_create(type,input,bc='all',edge='all',corner='all'):

	from genRhs import dim

	if dim == 3:		
		if type == 'begin':
			if((bc != 'k') and (edge != 'k') and (corner != 'k')):input.write('   do k=idloop(5),idloop(6)'+'\n')
			if((bc != 'j') and (edge != 'j') and (corner != 'j')):input.write('      do j=idloop(3),idloop(4) '+'\n')
			# beg.write('!$omp simd '+'\n'+'      	   do i=idloop(1),idloop(2) '+'\n')
			if((bc != 'i') and (edge != 'i') and (corner != 'i')):input.write(' '+'\n'+'      	   do i=idloop(1),idloop(2) '+'\n')
		elif type=='end':	
			if((bc != 'k') and (edge != 'k') and (corner != 'k')):input.write('        enddo'+'\n')
			if((bc != 'j') and (edge != 'j') and (corner != 'j')):input.write('      enddo'+'\n')
			if((bc != 'i') and (edge != 'i') and (corner != 'i')):input.write('   enddo'+'\n')		

	elif dim == 2 :
		if type == 'begin':
			if((bc != 'j') and (edge != 'j') and (corner != 'j')):input.write('     do j=idloop(3),idloop(4) '+'\n')
			# beg.write('!$omp simd'+'\n'+'      do i=idloop(1),idloop(2) '+'\n')
			if((bc != 'i') and (edge != 'i') and (corner != 'i')):input.write(' '+'\n'+'      do i=idloop(1),idloop(2) '+'\n')
		elif type=='end':	
			if((bc != 'j') and (edge != 'j') and (corner != 'j')):input.write('     enddo'+'\n')
			if((bc != 'i') and (edge != 'i') and (corner != 'i')):input.write('   enddo'+'\n')
	else:
		if type == 'begin':	
			# beg.write('!$omp simd'+'\n'+'      do i=idloop(1),idloop(2) '+'\n')
			if((bc != 'i') and (edge != 'i') and (corner != 'i')):input.write(' '+'\n'+'      do i=idloop(1),idloop(2) '+'\n')
		elif type=='end':	
			if((bc != 'i') and (edge != 'i') and (corner != 'i')):input.write('   enddo'+'\n')		

def globvar():

	from genRhs import dim, incPATH

	try:
		from genRhs import varstored
	except:
		varstored = {}	

	try:
		from genRhs import varbc
	except:
		varbc = {}

	globalvarRhs   = open(incPATH+'includeRHS_globVar.f90','w')  
	globalvarRk3   = open(incPATH+'includeRHS_globVar_rk3.f90','w')
	globalvarFlt   = open(incPATH+'includeRHS_globVar_filter.f90','w')
	storedglobvar  = open(incPATH+'includeRHS_globVarStored.f90','w')
	storedglobvarF = open(incPATH+'includeF_globVarStored.f90','w')
	qbc            = open(incPATH+'includeqbc_var.f90','w')
	qbcrk          = open(incPATH+'includeqbc_varrk.f90','w')

	nvbc     = {'face':{'i' :0,'j' :0,'k' :0},
		        'edge':{'ij':0,'jk':0,'ik':0}}

	nvbc_out = {'face':{'i' :'nvar_f(1)','j' :'nvar_f(2)','k' :'nvar_f(3)'},
		        'edge':{'ij':'nvar_f(1)','jk':'nvar_f(2)','ik':'nvar_f(3)'}}

	nvbcrk_out = {'face':{'i' :'nfacei','j' :'nfacej','k' :'nfacek'},
		          'edge':{'ij':'nedgeij','jk':'nedgejk','ik':'nedgeik'}}


	sizebcrk   = {'face':{'i'  :{3:'(1-hlo:ny+hlo,1-hlo:nz+hlo',
							   2:'(1-hlo:ny+hlo',
							   1:'(1'},
						'j'  :{3:'(1-hlo:nx+hlo,1-hlo:nz+hlo',
							   2:'(1-hlo:nx+hlo',
							   1:'(1'},
						'k'  :{3:'(1-hlo:nx+hlo,1-hlo:ny+hlo',
							   2:'(1-hlo:nx+hlo',
							   1:'(1'}},
				'edge':{'ij' :{3:'(1-hlo:nz+hlo',
							   2:'(1',
							   1:'(1'},
						'jk' :{3:'(1-hlo:nx+hlo',
							   2:'(1',
							   1:'(1'},
						'ik' :{3:'(1-hlo:ny+hlo',
							   2:'(1',
							   1:'(1'}}}        


	sizebc   = {'face':{'i'  :{3:'(idarray(3):idarray(4),idarray(5):idarray(6)',
							   2:'(idarray(3):idarray(4)',
							   1:'(1'},
						'j'  :{3:'(idarray(1):idarray(2),idarray(5):idarray(6)',
							   2:'(idarray(1):idarray(2)',
							   1:'(1'},
						'k'  :{3:'(idarray(1):idarray(2),idarray(3):idarray(4)',
							   2:'(idarray(1):idarray(2)',
							   1:'(1'}},
		        'edge':{'ij' :{3:'(idarray(5):idarray(6)',
							   2:'(1',
							   1:'(1'},
						'jk' :{3:'(idarray(1):idarray(2)',
							   2:'(1',
							   1:'(1'},
						'ik' :{3:'(idarray(3):idarray(4)',
							   2:'(1',
							   1:'(1'}}}        


	for v in varbc:
		for bcloc in ['face','edge']:
			if bcloc in varbc[v]:
				loctype = ''.join(sorted(varbc[v][bcloc].replace('1','').replace('max','')))
				nvbc[bcloc][loctype] = nvbc[bcloc][loctype] + 1

	qbc_out = 'real(wp),intent(inout) ::'
	qbcrk_out = 'real(wp),intent(inout) ::'

	for bcloc in ['face','edge']:
			for dir in nvbc[bcloc]:
					if nvbc[bcloc][dir] > 0:
							if qbc_out[-2:] == '::':
								qbc_out = qbc_out + ' q'+bcloc+'_'+dir+sizebc[bcloc][dir][dim]+','+nvbc_out[bcloc][dir]+'),&\n'
							else:
								qbc_out = qbc_out + '                       q'+bcloc+'_'+dir+sizebc[bcloc][dir][dim]+','+nvbc_out[bcloc][dir]+'),&\n'	

							if qbcrk_out[-2:] == '::':
								qbcrk_out = qbcrk_out + ' q'+bcloc+'_'+dir+sizebcrk[bcloc][dir][dim]+','+nvbcrk_out[bcloc][dir]+'),&\n'
							else:
								qbcrk_out = qbcrk_out + '                       q'+bcloc+'_'+dir+sizebcrk[bcloc][dir][dim]+','+nvbcrk_out[bcloc][dir]+'),&\n'		
					else:
							if qbc_out[-2:] == '::':
								qbc_out = qbc_out + ' q'+bcloc+'_'+dir+'(1),&\n'
							else:
								qbc_out = qbc_out + '                       q'+bcloc+'_'+dir+'(1),&\n'	

							if qbcrk_out[-2:] == '::':
								qbcrk_out = qbcrk_out + ' q'+bcloc+'_'+dir+'(1),&\n'
							else:
								qbcrk_out = qbcrk_out + '                       q'+bcloc+'_'+dir+'(1),&\n'	


	qbc.write(qbc_out[:-3])
	qbcrk.write(qbcrk_out[:-3])
	

	if dim == 3:
		storedglobvar.write('real(wp),intent(inout) :: q(idarray(1):idarray(2),idarray(3):idarray(4),idarray(5):idarray(6),neq),&'+'\n')
		if varstored!= {}:
			storedglobvar.write('                        qst(idarray(1):idarray(2),idarray(3):idarray(4),idarray(5):idarray(6),neqst)'+'\n')	
		else:	
			storedglobvar.write('                        qst(1)'+'\n')	
	elif dim == 2 :
		storedglobvar.write('real(wp),intent(inout) :: q(idarray(1):idarray(2),idarray(3):idarray(4),neq),&'+'\n')
		if varstored!= {}:
			storedglobvar.write('                        qst(idarray(1):idarray(2),idarray(3):idarray(4),neqst)'+'\n')	
		else:	
			storedglobvar.write('                        qst(1)'+'\n')	
	else:	
		storedglobvar.write('real(wp),intent(inout) :: q(idarray(1):idarray(2),neq),&'+'\n')
		if varstored!= {}:
			storedglobvar.write('                        qst(idarray(1):idarray(2),neqst)'+'\n')	
		else:	
			storedglobvar.write('                        qst(1)'+'\n')	

	if dim == 3:
		storedglobvarF.write('real(wp),intent(inout) :: q(1-hlo:nx+hlo,1-hlo:ny+hlo,1-hlo:nz+hlo,neq),&'+'\n')
		storedglobvarF.write('                        rhs(1-hlo:nx+hlo,1-hlo:ny+hlo,1-hlo:nz+hlo,neq),&'+'\n')
		if varstored!= {}:
			storedglobvarF.write('                        qst(1-hlo:nx+hlo,1-hlo:ny+hlo,1-hlo:nz+hlo,neqst)'+'\n')	
		else:	
			storedglobvarF.write('                        qst(1)'+'\n')	
	elif dim == 2 :
		storedglobvarF.write('real(wp),intent(inout) :: q(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&'+'\n')
		storedglobvarF.write('                        rhs(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&'+'\n')
		if varstored!= {}:
			storedglobvarF.write('                        qst(1-hlo:nx+hlo,1-hlo:ny+hlo,neqst)'+'\n')	
		else:	
			storedglobvarF.write('                        qst(1)'+'\n')	
	else:	
		storedglobvarF.write('real(wp),intent(inout) :: q(1-hlo:nx+hlo,neq),&'+'\n')
		storedglobvarF.write('                        rhs(1-hlo:nx+hlo,neq),&'+'\n')
		if varstored!= {}:
			storedglobvarF.write('                        qst(1-hlo:nx+hlo,neqst)'+'\n')	
		else:	
			storedglobvarF.write('                        qst(1)'+'\n')	


	if dim == 3:
		globalvarRhs.write('real(wp),intent(inout) :: q(idarray(1):idarray(2),idarray(3):idarray(4),idarray(5):idarray(6),neq),&'+'\n')
		globalvarRhs.write('                        rhs(idarray(1):idarray(2),idarray(3):idarray(4),idarray(5):idarray(6),neq),&'+'\n')	
		if varstored!= {}:
			globalvarRhs.write('                        qst(idarray(1):idarray(2),idarray(3):idarray(4),idarray(5):idarray(6),neqst)'+'\n')	
		else:	
			globalvarRhs.write('                        qst(1)'+'\n')	
	elif dim == 2 :
		globalvarRhs.write('real(wp),intent(inout) :: q(idarray(1):idarray(2),idarray(3):idarray(4),neq),&'+'\n')
		globalvarRhs.write('                        rhs(idarray(1):idarray(2),idarray(3):idarray(4),neq),&'+'\n')	
		if varstored!= {}:
			globalvarRhs.write('                        qst(idarray(1):idarray(2),idarray(3):idarray(4),neqst)'+'\n')	
		else:	
			globalvarRhs.write('                        qst(1)'+'\n')	
	else:	
		globalvarRhs.write('real(wp),intent(inout) :: q(idarray(1):idarray(2),neq),&'+'\n')
		globalvarRhs.write('                        rhs(idarray(1):idarray(2),neq),&'+'\n')
		if varstored!= {}:
			globalvarRhs.write('                        qst(idarray(1):idarray(2),neqst)'+'\n')	
		else:	
			globalvarRhs.write('                        qst(1)'+'\n')	

		# vglob.write('!dir$ ASSUME_ALIGNED rhs: 64'+'\n')
		# vglob.write('!dir$ ASSUME_ALIGNED q: 64'+'\n')

	if dim == 3:
		globalvarRk3.write('real(wp),intent(inout) :: q(1-hlo:nx+hlo,1-hlo:ny+hlo,1-hlo:nz+hlo,neq),&'+'\n')
		globalvarRk3.write('                         q1(1-hlo:nx+hlo,1-hlo:ny+hlo,1-hlo:nz+hlo,neq),&'+'\n')
		globalvarRk3.write('                         q2(1-hlo:nx+hlo,1-hlo:ny+hlo,1-hlo:nz+hlo,neq),&'+'\n')
		globalvarRk3.write('                        rhs(1-hlo:nx+hlo,1-hlo:ny+hlo,1-hlo:nz+hlo,neq),&'+'\n')	
		if varstored!= {}:
			globalvarRk3.write('                        qst(1-hlo:nx+hlo,1-hlo:ny+hlo,1-hlo:nz+hlo,neqst)'+'\n')	
		else:	
			globalvarRk3.write('                        qst(1)'+'\n')	
	elif dim == 2 :
		globalvarRk3.write('real(wp),intent(inout) ::  q(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&'+'\n')
		globalvarRk3.write('                          q1(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&'+'\n')
		globalvarRk3.write('                          q2(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&'+'\n')
		globalvarRk3.write('                         rhs(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&'+'\n')			
		if varstored!= {}:
			globalvarRk3.write('                        qst(1-hlo:nx+hlo,1-hlo:ny+hlo,neqst)'+'\n')	
		else:	
			globalvarRk3.write('                        qst(1)'+'\n')	
	else:
		globalvarRk3.write('real(wp),intent(inout) ::  q(1-hlo:nx+hlo,neq),&'+'\n')
		globalvarRk3.write('                          q1(1-hlo:nx+hlo,neq),&'+'\n')
		globalvarRk3.write('                          q2(1-hlo:nx+hlo,neq),&'+'\n')
		globalvarRk3.write('                         rhs(1-hlo:nx+hlo,neq),&'+'\n')	
		if varstored!= {}:
			globalvarRk3.write('                        qst(1-hlo:nx+hlo,neqst)'+'\n')	
		else:	
			globalvarRk3.write('                        qst(1)'+'\n')	
		# vglob.write('!dir$ ASSUME_ALIGNED q: 64'+'\n')
		# vglob.write('!dir$ ASSUME_ALIGNED q1: 64'+'\n')
		# vglob.write('!dir$ ASSUME_ALIGNED q2: 64'+'\n')								
		# vglob.write('!dir$ ASSUME_ALIGNED rhs: 64'+'\n')		
		
	if dim == 3:
		globalvarFlt.write('real(wp),intent(inout) :: q(1-hlo:nx+hlo,1-hlo:ny+hlo,1-hlo:nz+hlo,neq),&'+'\n')
		globalvarFlt.write('                         q2(1-hlo:nx+hlo,1-hlo:ny+hlo,1-hlo:nz+hlo,neq)'+'\n')			
	elif dim == 2 :
		globalvarFlt.write('real(wp),intent(inout) ::  q(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&'+'\n')			
		globalvarFlt.write('                          q2(1-hlo:nx+hlo,1-hlo:ny+hlo,neq)'+'\n')
				
	else:
		globalvarFlt.write('real(wp),intent(inout) ::  q(1-hlo:nx+hlo,neq),&'+'\n')			
		globalvarFlt.write('                          q2(1-hlo:nx+hlo,neq)'+'\n')
			
		# vglob.write('!dir$ ASSUME_ALIGNED q2: 64'+'\n')			
		# vglob.write('!dir$ ASSUME_ALIGNED q: 64'+'\n')			
	
# define derivatives operators :

der1 = {}
der2 = {}

dder1 = {}
dder2 = {}

dirvars = ['x','y','z']
for v in dirvars:
	der1[v] = re.compile('(?<=\[)([^]]+)(?=\]_1'+v+'{1}[\W])')
	for vv in dirvars:
		der2[v+vv] = re.compile('(?<=\[)([^]]+)(?=\]_2['+v+vv+']{2}[\W])')


for v in dirvars:
	dder1[v] = re.compile('(?<=\{)([^}]+)(?=\}_1'+v+'{1}[\W])')
	for vv in dirvars:
		dder2[v+vv] = re.compile('(?<=\[)([^]]+)(?=\]_2['+v+vv+']{2}[\W])')

der1n = {}
der2n = {}

dder1n = {}
dder2n = {}

dirvars = ['x','y','z']
for v in dirvars:
	der1n[v] = re.compile('(\[[^]]+\]_1'+v+'{1}[\W])')
	for vv in dirvars:
		der2n[v+vv] = re.compile('(\[[^]]+\]_2['+v+vv+']{2}[\W])')		

for v in dirvars:
	dder1n[v] = re.compile('(\{[^}]+\}_1'+v+'{1}[\W])')
	for vv in dirvars:
		dder2n[v+vv] = re.compile('(\{[^}]+\}_2['+v+vv+']{2}[\W])')		

import sys

def genSymbDer1(exp,output,locvar,
				order=2,stencil=3,
				indi='i',indj='j',indk='k',
				vname='symb',
				history={'x':{'var':[],'symb':[]},'y':{'var':[],'symb':[]},'z':{'var':[],'symb':[]}},
				dhistory={'x':[[],[],[],[]],'y':[[],[],[],[]],'z':[[],[],[],[]]}):

		from genRhs import dim

		vname  = vname.strip()

		for v in dirvars:

				indiri = indi
				indirj = indj
				indirk = indk

				subdname = []							
				for i,derexp in enumerate(re.finditer(der1[v], exp)):	

					iprv_d = i
					ihist  = None

					for itprv_d,prv_d in enumerate(list(re.finditer(der1[v], exp))[:i]):			
						if prv_d.group(0).replace(" ", "") == derexp.group(0).replace(" ", ""):							
							iprv_d = itprv_d
							break
			
					for ithist,exphist in enumerate(history[v]['symb']):
						if(derexp.group(0).replace(" ", "") == exphist.replace(" ", "")):									
							ihist = ithist
							break

					if iprv_d == i and (ihist == None):	 # new derivative (not already computed...)
						history[v]['symb'].append(derexp.group(0))
						history[v]['var'].append('d1_'+vname+'_d'+v+'_'+str(i)+'_')

						if(dim < 3):
							if (v=='z' and derexp.group(0) != []):
								nsymb = list(re.finditer(der1n[v], exp))
								exception('wrong symbolic equation : ask for z derivatives ('+str(nsymb[i].group(0))+') with dim ='+str(dim),message='error')

						elif(dim < 2):		
							if  (v=='y' and derexp.group(0) != []):
								nsymb = list(re.finditer(der1n[v], exp))
								exception('wrong symbolic equation : ask for y derivatives ('+str(nsymb[i].group(0))+') with dim ='+str(dim),message='error')	

						dsubdname = {'x':[[],[]],'y':[[],[]],'z':[[],[]]}

						# generate and cross derivatives:

						for vv in dirvars:	
							
							for di,dderexp in enumerate(re.finditer(dder1[vv], derexp.group(0))):
							
								diprv_d = di
								dihist  = None
			
								for ditprv_d,dprv_d in enumerate(list(re.finditer(dder1[vv], dderexp.group(0)))[:di]):						
									if dprv_d.group(0).replace(" ", "") == dderexp.group(0).replace(" ", ""):							
										diprv_d = ditprv_d
										break
					
		
								for dithist,dexphist in enumerate(dhistory[v][2]):										
									if (dderexp.group(0).replace(" ", "") == dexphist.replace(" ", "")) and (vv == dhistory[v][1][dithist]):								
										dihist = dithist
										break
								if diprv_d == di and (dihist == None):	 # new derivative (not already computed...)
									dhistory[v][2].append(dderexp.group(0))	
									dhistory[v][1].append(vv)				

									vlocnam = 'd2_'+vname+'_d'+v+'d'+vv+'_'+str(i)+'_'+str(di)+'_'

									if(dim < 3):
										if (vv=='z' and derexp.group(0) != []):
											nsymb = list(re.finditer(dder1n[vv], exp))
											exception('wrong symbolic equation : ask for z derivatives ('+str(nsymb[di].group(0))+') with dim ='+str(dim),message='error')

									elif(dim < 2):		
										if  (vv=='y' and derexp.group(0) != []):
											nsymb = list(re.finditer(dder1n[vv], exp))
											exception('wrong symbolic equation : ask for y derivatives ('+str(nsymb[di].group(0))+') with dim ='+str(dim),message='error')			

									hlo = int((stencil-1)/2)
						

									if (np.mod((stencil-1),2)==1) : 
										raise ValueError('Wrong stencil length for centrered finite differences (should be odd value)  ! ',stencil)
						
									if ( vv == 'x' ) or ( vv == 'y' and (dim >= 2) ) or (vv == 'z' and (dim >= 3)):

										for shift in range(-hlo,0):
					
												if   v == 'x':
													indiri = indi+'{:+d}'.format(shift)
												elif v == 'y' and (dim >= 2):
													indirj = indj+'{:+d}'.format(shift)
												elif v == 'z' and dim >= 3:
													indirk = indk+'{:+d}'.format(shift)	
					
												ddname = genVname(vlocnam,      i=indiri,j=indirj,k=indirk)
												[dexprnbg,dvnamebg] = genNbg(dderexp.group(0), vv , stencil, i=indiri,j=indirj,k=indirk,vname=ddname+'_')
												locvar.append(dvnamebg)
												createNbg(dvnamebg,dexprnbg,output)
												output.write(der(dvnamebg,order,stencil,varname=ddname))
												locvar.append([ddname])
								
												if   vv == 'x':
													output.write(ddname+' = '+ddname.strip()+'*param_float(1)'+'\n\n')
												elif vv == 'y':	
													output.write(ddname+' = '+ddname.strip()+'*param_float(2)'+'\n\n')
												elif vv == 'z':
													output.write(ddname+' = '+ddname.strip()+'*param_float(3)'+'\n\n')
									
										for shift in range(1,hlo+1):
					
												if   v == 'x':
													indiri = indi+'{:+d}'.format(shift)
												elif v == 'y' and (dim >= 2):
													indirj = indj+'{:+d}'.format(shift)
												elif v == 'z' and dim >= 3:
													indirk = indk+'{:+d}'.format(shift)	
					
												ddname = genVname(vlocnam,      i=indiri,j=indirj,k=indirk)
												[dexprnbg,dvnamebg] = genNbg(dderexp.group(0), vv , stencil, i=indiri,j=indirj,k=indirk,vname=ddname+'_')
												locvar.append(dvnamebg)
												createNbg(dvnamebg,dexprnbg,output)
												output.write(der(dvnamebg,order,stencil,varname=ddname))
												locvar.append([ddname])

												if   vv == 'x':
													output.write(ddname+' = '+ddname.strip()+'*param_float(1)'+'\n\n')
												elif vv == 'y':	
													output.write(ddname+' = '+ddname.strip()+'*param_float(2)'+'\n\n')
												elif vv == 'z':
													output.write(ddname+' = '+ddname.strip()+'*param_float(3)'+'\n\n')
									dhistory[v][3].append(vlocnam)		
								
								if(dihist != None): 
									dsubptr = dhistory[v][3][dihist]						
								else:	
									dsubptr = 'd2_'+vname+'_d'+v+'d'+vv+'_'+str(i)+'_'+str(diprv_d)+'_'													
								dsubdname[vv][0].append(dsubptr)
								dsubdname[vv][1].append(dderexp.group(0))	

							if( dsubdname[vv] != [[],[]]):								
								for ri,rderexp in enumerate(re.finditer(der1n[v], exp)):	
									rptr = re.escape(rderexp.group(0))
									rexp = rderexp.group(0)
									for di,dderexp in enumerate(re.finditer(dder1n[vv], rexp)):		
										dptr = re.escape(dderexp.group(0))		
										if ( '{'+dsubdname[vv][1][di]+'}_1'+vv.strip() == dderexp.group(0).strip()):
											rexp = re.sub(dptr,dsubdname[vv][0][di] , rexp)		
									exp = re.sub(rptr,rexp , exp)															

						if( dsubdname == {'x':[[],[]],'y':[[],[]],'z':[[],[]]}):
							dname =  genVname('d1_'+vname+'_d'+v+'_'+str(i)+'_',indi,indj,indk)
							[exprnbg,vnamebg] = genNbg(derexp.group(0),v,stencil,vname='d1_'+vname+'_d'+v+'_'+str(i)+'_')
							locvar.append(vnamebg)
			
							createNbg(vnamebg,exprnbg,output)
							output.write(der(vnamebg,order,stencil,varname=dname))
							locvar.append([dname])
							if   v == 'x':
								output.write(dname + ' = '+ dname +'*param_float(1)'+'\n\n')
							elif v == 'y':	
								output.write(dname + ' = '+ dname +'*param_float(2)'+'\n\n')
							elif v == 'z':
								output.write(dname + ' = '+ dname +'*param_float(3)'+'\n\n')
						
				for i,derexp in enumerate(re.finditer(der1n[v], exp)):	

					iprv_d = i
					ihist  = None

					for itprv_d,prv_d in enumerate(list(re.finditer(der1n[v], exp))[:i]):						
						if prv_d.group(0).replace(" ", "") == derexp.group(0).replace(" ", ""):							
							iprv_d = itprv_d
							break					


					for ithist,exphist in enumerate(history[v]['symb']):
						if(re.findall(der1[v],derexp.group(0))[0].replace(" ", "") == exphist.replace(" ", "")):						
							ihist = ithist
							break

					if iprv_d == i and (ihist == None):	 # new derivative (not already computed...)

						history[v]['symb'].append(derexp.group(0))					
						subptr = re.findall(der1[v],derexp.group(0))
						
						dname =  genVname('d1_'+vname+'_d'+v+'_'+str(i)+'_',indi,indj,indk)
						[exprnbg,vnamebg] = genNbg(subptr[0],v,stencil,vname='d1_'+vname+'_d'+v+'_'+str(i)+'_')
						locvar.append(vnamebg)
		
						createNbg(vnamebg,exprnbg,output)
						output.write(der(vnamebg,order,stencil,varname=dname))
						locvar.append([dname])
						if   v == 'x':
							output.write(dname + ' = '+ dname +'*param_float(1)'+'\n\n')
						elif v == 'y':	
							output.write(dname + ' = '+ dname +'*param_float(2)'+'\n\n')
						elif v == 'z':
							output.write(dname + ' = '+ dname +'*param_float(3)'+'\n\n')
												
						history[v]['var'].append('d1_'+vname+'_d'+v+'_'+str(i)+'_')

					if(ihist != None): 
						subptr = history[v]['var'][ihist]						
					else:	
						subptr = 'd1_'+vname+'_d'+v+'_'+str(iprv_d)+'_'													
					subdname.append(subptr)			
							# 
	
				for i,derexp in enumerate(re.finditer(der1n[v], exp)):		
					ptr = re.escape(derexp.group(0))
					exp = re.sub(ptr,subdname[i] ,exp)

		return [exp,locvar,history]

def genSymbDer2(exp,output,locvar,
				order=2,stencil=3,
				indi='i',indj='j',indk='k',
				vname='symb',
				history={'x':[[],[]],'y':[[],[]],'z':[[],[]]}):

		from genRhs import dim

		vname = vname.strip()
		indiri = indi
		indirj = indj
		indirk = indk

		for v in dirvars:

				subdname = []							
				for i,derexp in enumerate(re.finditer(der2[v+v], exp)):	

					iprv_d = i
					ihist  = None

					for itprv_d,prv_d in enumerate(list(re.finditer(der2[v+v], exp))[:i]):						
						if prv_d.group(0).replace(" ", "") == derexp.group(0).replace(" ", ""):							
							iprv_d = itprv_d
							break

					# if history[dir][v][1] != []:						
					for ithist,exphist in enumerate(history[v][1]):
							if(derexp.group(0).replace(" ", "") == exphist.replace(" ", "")):								
								ihist = ithist
								break

					if iprv_d == i and (ihist == None):	 # new derivative (not already computed...)
						history[v][1].append(derexp.group(0))
						if(dim < 3):
							if (v=='z' and derexp.group(0) != []):
								nsymb = list(re.finditer(der2n[v+v], exp))
								exception('wrong symbolic equation : ask for z derivatives ('+str(nsymb[i].group(0))+') with dim ='+str(dim),message='error')

						elif(dim < 2):		
							if  (v=='y' and derexp.group(0) != []):
								nsymb = list(re.finditer(der2n[v+v], exp))
								exception('wrong symbolic equation : ask for y derivatives ('+str(nsymb[i].group(0))+') with dim ='+str(dim),message='error')

	
						hlo = int((stencil-1)/2)
			
						if (np.mod((stencil-1),2)==1) : 
							raise ValueError('Wrong stencil length for centrered finite differences (should be odd value)  ! ',stencil)
			
	
						dname = genVname('d2_'+vname+'_d'+v+'_'+str(i)+'_',      i=indiri,j=indirj,k=indirk)
						[exprnbg,vnamebg] = genNbg(derexp.group(0), v , stencil, i=indiri,j=indirj,k=indirk,vname=dname+'_')
						locvar.append(vnamebg)
						createNbg(vnamebg,exprnbg,output,'second')
						output.write(dder(vnamebg,order,stencil,varname=dname))
						locvar.append([dname])
	
						if   v == 'x':
							output.write(dname+' = '+dname.strip()+'*param_float(1)*param_float(1)'+'\n\n')
						elif v == 'y':	
							output.write(dname+' = '+dname.strip()+'*param_float(2)*param_float(2)'+'\n\n')
						elif v == 'z':
								output.write(dname+' = '+dname.strip()+'*param_float(3)*param_float(3)'+'\n\n')
					
						
						history[v][0].append('d2_'+vname+'_d'+v+'_'+str(i)+'_')		
													

					if(ihist != None): 
						subptr = history[v][0][ihist]						
					else:	
						subptr = 'd2_'+vname+'_d'+v+'_'+str(iprv_d)+'_'													
					subdname.append(subptr)			
	

				for i,derexp in enumerate(re.finditer(der2n[v+v], exp)):		
					ptr = re.escape(derexp.group(0))					
					exp = re.sub(ptr,subdname[i] ,exp)

				# print(history)	
		return [exp,locvar,history]

def genSymbDer1_bc(bcdic,exp,output,locvar,
				   order=2,stencil=3,
				   indi='i',indj='j',indk='k',
				   vname='symb',
				   history={'x':{'var':[],'symb':[]},'y':{'var':[],'symb':[]},'z':{'var':[],'symb':[]}},
				   dhistory={'x':[[],[],[],[]],'y':[[],[],[],[]],'z':[[],[],[],[]]}):

		from genRhs import dim

		dirBc_v   = '    '
		dirBc_vv  = '    '
		bclayer_v = int((stencil-1)/2)
		bclayer_vv= int((stencil-1)/2)
		indbc_v   = int((stencil-1)/2)
		indbc_vv  = int((stencil-1)/2)

		vname = vname.strip()
		bc  = False

		for v in dirvars:

				indiri = indi
				indirj = indj
				indirk = indk

				subdname = []							
				for i,derexp in enumerate(re.finditer(der1[v], exp)):	

					iprv_d = i
					ihist  = None

					for itprv_d,prv_d in enumerate(list(re.finditer(der1[v], exp))[:i]):			
						if prv_d.group(0).replace(" ", "") == derexp.group(0).replace(" ", ""):							
							iprv_d = itprv_d
							break
			
					for ithist,exphist in enumerate(history[v]['symb']):
						if(derexp.group(0).replace(" ", "") == exphist.replace(" ", "")):									
							ihist = ithist
							break

					if iprv_d == i and (ihist == None):	 # new derivative (not already computed...)
						history[v]['symb'].append(derexp.group(0))
						history[v]['var'].append('d1_'+vname+'_d'+v+'_'+str(i)+'_')

						if(dim < 3):
							if (v=='z' and derexp.group(0) != []):
								nsymb = list(re.finditer(der1n[v], exp))
								exception('wrong symbolic equation : ask for z derivatives ('+str(nsymb[i].group(0))+') with dim ='+str(dim),message='error')

						elif(dim < 2):		
							if  (v=='y' and derexp.group(0) != []):
								nsymb = list(re.finditer(der1n[v], exp))
								exception('wrong symbolic equation : ask for y derivatives ('+str(nsymb[i].group(0))+') with dim ='+str(dim),message='error')

						dsubdname = {'x':[[],[]],'y':[[],[]],'z':[[],[]]}

						# generate second and cross derivatives:

						for vv in dirvars:	
							
							for di,dderexp in enumerate(re.finditer(dder1[vv], derexp.group(0))):
							
								diprv_d = di
								dihist  = None
			
								for ditprv_d,dprv_d in enumerate(list(re.finditer(dder1[vv], dderexp.group(0)))[:di]):						
									if dprv_d.group(0).replace(" ", "") == dderexp.group(0).replace(" ", ""):							
										diprv_d = ditprv_d
										break
					
		
								for dithist,dexphist in enumerate(dhistory[v][2]):										
									if (dderexp.group(0).replace(" ", "") == dexphist.replace(" ", "")) and (vv == dhistory[v][1][dithist]):								
										dihist = dithist
										break
								if diprv_d == di and (dihist == None):	 # new derivative (not already computed...)
									dhistory[v][2].append(dderexp.group(0))	
									dhistory[v][1].append(vv)				

									vlocnam = 'd2_'+vname+'_d'+v+'d'+vv+'_'+str(i)+'_'+str(di)+'_'

									hlo = int((stencil-1)/2)
									bc = False


									if(dim < 3):
										if (vv=='z' and derexp.group(0) != []):
											nsymb = list(re.finditer(dder1n[vv], exp))
											exception('wrong symbolic equation : ask for z derivatives ('+str(nsymb[di].group(0))+') with dim ='+str(dim),message='error')

									elif(dim < 2):		
										if  (vv=='y' and derexp.group(0) != []):
											nsymb = list(re.finditer(dder1n[vv], exp))
											exception('wrong symbolic equation : ask for y derivatives ('+str(nsymb[di].group(0))+') with dim ='+str(dim),message='error')
				

									hlo = int((stencil-1)/2)
						
									if (np.mod((stencil-1),2)==1) : 
										raise ValueError('Wrong stencil length for centrered finite differences (should be odd value)  ! ',stencil)
						
									if ( vv == 'x' ) or ( vv == 'y' and (dim >= 2) ) or (vv == 'z' and (dim >= 3)):


										# print("derivatives outer ",v," inner ",vv,hlo)
										if ( v == 'x') and (bcdic['i']['dir'] != None):	
											dirBc_v   = bcdic['i']['dir']
											indbc_v   = bcdic['i']['indbc']		
											bclayer_v = indbc_v				
											hlo       = indbc_v
										if ( v == 'y') and (bcdic['j']['dir'] != None):	
											dirBc_v   = bcdic['j']['dir']
											indbc_v   = bcdic['j']['indbc']		
											bclayer_v = indbc_v				
											hlo       = indbc_v
										if ( v == 'z') and (bcdic['k']['dir'] != None):	
											dirBc_v   = bcdic['k']['dir']
											indbc_v   = bcdic['k']['indbc']		
											bclayer_v = indbc_v				
											hlo       = indbc_v									
										# print("derivatives outer ",v," inner ",vv,hlo)

										# dirBc_vv   = dirBc_v
										# bclayer_vv = int((stencil-1)/2)

										if(hlo == 0):
											offsetdir = 1
											if( dirBc_v[1:] == 'max'):
												offsetdir = -1
											for shift in range(0,3):	
													if   v == 'x':
														indiri = indi+'{:+d}'.format(offsetdir*shift)
													elif v == 'y' and (dim >= 2):
														indirj = indj+'{:+d}'.format(offsetdir*shift)
													elif dim >= 3:
														indirk = indk+'{:+d}'.format(offsetdir*shift)	

													if ( vv == 'x') and (bcdic['i']['dir'] != None):							
														bclayer_vv = bcdic['i']['indbc']
														dirBc_vv   = bcdic['i']['dir']
														bc      = True
													if ( vv == 'y') and (bcdic['j']['dir'] != None):							
														bclayer_vv = bcdic['j']['indbc']
														dirBc_vv   = bcdic['j']['dir']	
														bc      = True	
													if ( vv == 'z') and (bcdic['k']['dir'] != None):							
														bclayer_vv = bcdic['k']['indbc']
														dirBc_vv   = bcdic['k']['dir']
														bc      = True															

													if ( vv == 'x') and (dirBc_v[0] == 'i'):							
														bclayer_vv = indbc_v + shift
														bc      = True
													if ( vv == 'y') and (dirBc_v[0] == 'j'):
														bclayer_vv = indbc_v + shift	
														bc      = True	
													if ( vv == 'z') and (dirBc_v[0] == 'k'):
														bclayer_vv = indbc_v + shift
														bc      = True	
						
													ddname = genVname(vlocnam,
																      i=indiri,j=indirj,k=indirk)
													[dexprnbg,dvnamebg] = genNbg(dderexp.group(0), vv , stencil, 
																	  i=indiri,j=indirj,k=indirk,vname=ddname+'_',dirBc=dirBc_vv,indbc=bclayer_vv)
													locvar.append(dvnamebg)
													createNbg(dvnamebg,dexprnbg,output,indbc=bclayer_vv,bc=bc)
													output.write(der(dvnamebg,order,stencil,varname=ddname,dirBC=dirBc_vv,indbc=bclayer_vv,bc=bc))
													locvar.append([ddname])
						
													if   vv == 'x':
														output.write(ddname+' = '+ddname.strip()+'*param_float(1)'+'\n\n')
													elif vv == 'y':	
														output.write(ddname+' = '+ddname.strip()+'*param_float(2)'+'\n\n')
													elif vv == 'z':
														output.write(ddname+' = '+ddname.strip()+'*param_float(3)'+'\n\n')


										for shift in range(-hlo,0):
					
												if   v == 'x':
													indiri = indi+'{:+d}'.format(shift)
												elif v == 'y' and (dim >= 2):
													indirj = indj+'{:+d}'.format(shift)
												elif v == 'z' and dim >= 3:
													indirk = indk+'{:+d}'.format(shift)
												
												if ( vv == 'x') and (bcdic['i']['dir'] != None):								
													if(v == 'x'): 
														if dirBc_v[1:] == 'max':
															bclayer_vv = indbc_v + abs(shift)	
														else:	
															bclayer_vv = indbc_v + shift
													else:
														bclayer_vv = bcdic['i']['indbc']
													dirBc_vv   = bcdic['i']['dir']				
													bc      = True

													
												if ( vv == 'y') and (bcdic['j']['dir'] != None):
													if(v == 'y'): 
														if dirBc_v[1:] == 'max':
															bclayer_vv = indbc_v + abs(shift)	
														else:	
															bclayer_vv = indbc_v + shift	
													else:
														bclayer_vv = bcdic['j']['indbc']
													dirBc_vv   = bcdic['j']['dir']			
													bc      = True


												if ( vv == 'z') and (bcdic['k']['dir'] != None):
													if(v == 'z'): 
														if dirBc_v[1:] == 'max':
															bclayer_vv = indbc_v + abs(shift)	
														else:	
															bclayer_vv = indbc_v + shift	
													else:
														bclayer_vv = bcdic['k']['indbc']
													dirBc_vv   = bcdic['k']['dir']		 		
													bc      = True

												ddname = genVname(vlocnam, i=indiri,j=indirj,k=indirk)
												[dexprnbg,dvnamebg] = genNbg(dderexp.group(0), vv , stencil, 
																  i=indiri,j=indirj,k=indirk,vname=ddname+'_',dirBc=dirBc_vv,indbc=bclayer_vv)
												locvar.append(dvnamebg)
												createNbg(dvnamebg,dexprnbg,output,indbc=bclayer_vv,bc=bc)									
												output.write(der(dvnamebg,order,stencil,varname=ddname,dirBC=dirBc_vv,indbc=bclayer_vv,bc=bc))
												locvar.append([ddname])
								
												if   vv == 'x':
													output.write(ddname+' = '+ddname.strip()+'*param_float(1)'+'\n\n')
												elif vv == 'y':	
													output.write(ddname+' = '+ddname.strip()+'*param_float(2)'+'\n\n')
												elif vv == 'z':
													output.write(ddname+' = '+ddname.strip()+'*param_float(3)'+'\n\n')
									
										for shift in range(1,hlo+1):
					
												if   v == 'x':
													indiri = indi+'{:+d}'.format(shift)
												elif v == 'y' and (dim >= 2):
													indirj = indj+'{:+d}'.format(shift)
												elif v == 'z' and dim >= 3:
													indirk = indk+'{:+d}'.format(shift)	

												if ( vv == 'x') and (bcdic['i']['dir'] != None):							
													if(v == 'x'): 
														if dirBc_v[1:] == 'max':
															bclayer_vv = indbc_v - abs(shift)	
														else:	
															bclayer_vv = indbc_v + shift
													else:
														bclayer_vv = bcdic['i']['indbc']
													dirBc_vv   = bcdic['i']['dir']				
													bc      = True
												if ( vv == 'y') and (bcdic['j']['dir'] != None):
													if(v == 'y'): 
														if dirBc_v[1:] == 'max':
															bclayer_vv = indbc_v - abs(shift)	
														else:	
															bclayer_vv = indbc_v + shift
													else:
														bclayer_vv = bcdic['j']['indbc']	
													dirBc_vv   = bcdic['j']['dir']			
													bc      = True	
												if ( vv == 'z') and (bcdic['k']['dir'] != None):
													if(v == 'z'): 
														if dirBc_v[1:] == 'max':
															bclayer_vv = indbc_v - abs(shift)	
														else:	
															bclayer_vv = indbc_v + shift
													else:
														bclayer_vv = bcdic['k']['indbc']
													dirBc_vv   = bcdic['k']['dir']				
													bc      = True

												ddname = genVname(vlocnam, i=indiri,j=indirj,k=indirk)
												[dexprnbg,dvnamebg] = genNbg(dderexp.group(0), vv , stencil, 
																  i=indiri,j=indirj,k=indirk,vname=ddname+'_',dirBc=dirBc_vv,indbc=bclayer_vv)
												locvar.append(dvnamebg)
												createNbg(dvnamebg,dexprnbg,output,indbc=bclayer_vv,bc=bc)
												output.write(der(dvnamebg,order,stencil,varname=ddname,dirBC=dirBc_vv,indbc=bclayer_vv,bc=bc))
												locvar.append([ddname])
								
												if   vv == 'x':
													output.write(ddname+' = '+ddname.strip()+'*param_float(1)'+'\n\n')
												elif vv == 'y':	
													output.write(ddname+' = '+ddname.strip()+'*param_float(2)'+'\n\n')
												elif vv == 'z':
													output.write(ddname+' = '+ddname.strip()+'*param_float(3)'+'\n\n')
									dhistory[v][3].append(vlocnam)		
								
								if(dihist != None): 
									dsubptr = dhistory[v][3][dihist]						
								else:	
									dsubptr = 'd2_'+vname+'_d'+v+'d'+vv+'_'+str(i)+'_'+str(diprv_d)+'_'													
								dsubdname[vv][0].append(dsubptr)
								dsubdname[vv][1].append(dderexp.group(0))	

							if( dsubdname[vv] != [[],[]]):								
								for ri,rderexp in enumerate(re.finditer(der1n[v], exp)):	
									rptr = re.escape(rderexp.group(0))
									rexp = rderexp.group(0)
									for di,dderexp in enumerate(re.finditer(dder1n[vv], rexp)):		
										dptr = re.escape(dderexp.group(0))		
										if ( '{'+dsubdname[vv][1][di]+'}_1'+vv.strip() == dderexp.group(0).strip()):
											rexp = re.sub(dptr,dsubdname[vv][0][di] , rexp)		
									exp = re.sub(rptr,rexp , exp)															

						if( dsubdname == {'x':[[],[]],'y':[[],[]],'z':[[],[]]}):


							bc = False

							if ( v == 'x') and (bcdic['i']['dir'] != None):
								dirBc_v   = bcdic['i']['dir']
								indbc_v   = bcdic['i']['indbc']			
								bc = True
							if ( v == 'y') and (bcdic['j']['dir'] != None):
								dirBc_v   = bcdic['j']['dir']
								indbc_v   = bcdic['j']['indbc']									
								bc = True
							if ( v == 'z') and (bcdic['k']['dir'] != None):
								dirBc_v   = bcdic['k']['dir']
								indbc_v   = bcdic['k']['indbc']	
								bc = True		
		
							dname =  genVname('d1_'+vname+'_d'+v+'_'+str(i)+'_', i=indi,j=indj,k=indk)
							[exprnbg,vnamebg] = genNbg(derexp.group(0),v,stencil,i=indi,j=indj,k=indk,
												vname='d1_'+vname+'_d'+v+'_'+str(i)+'_',dirBc=dirBc_v,indbc=indbc_v)
							locvar.append(vnamebg)
			
							createNbg(vnamebg,exprnbg,output,indbc=indbc_v,bc=bc)
							output.write(der(vnamebg,order,stencil,varname=dname,dirBC=dirBc_v,indbc=indbc_v,bc=bc))
							locvar.append([dname])
							if   v == 'x':
								output.write(dname + ' = '+ dname +'*param_float(1)'+'\n\n')
							elif v == 'y':	
								output.write(dname + ' = '+ dname +'*param_float(2)'+'\n\n')
							elif v == 'z':
								output.write(dname + ' = '+ dname +'*param_float(3)'+'\n\n')
						
				for i,derexp in enumerate(re.finditer(der1n[v], exp)):	

					iprv_d = i
					ihist  = None

					for itprv_d,prv_d in enumerate(list(re.finditer(der1n[v], exp))[:i]):						
						if prv_d.group(0).replace(" ", "") == derexp.group(0).replace(" ", ""):							
							iprv_d = itprv_d
							break					


					for ithist,exphist in enumerate(history[v]['symb']):
						if(re.findall(der1[v],derexp.group(0))[0].replace(" ", "") == exphist.replace(" ", "")):						
							ihist = ithist
							break

					if iprv_d == i and (ihist == None):	 # new derivative (not already computed...)

						history[v]['symb'].append(derexp.group(0))					
						subptr = re.findall(der1[v],derexp.group(0))

						bc = False
						if ( v == 'x') and (bcdic['i']['dir'] != None):
							dirBc_v   = bcdic['i']['dir']
							indbc_v   = bcdic['i']['indbc']			
							bc = True
						if ( v == 'y') and (bcdic['j']['dir'] != None):
							dirBc_v   = bcdic['j']['dir']
							indbc_v   = bcdic['j']['indbc']									
							bc = True
						if ( v == 'z') and (bcdic['k']['dir'] != None):
							dirBc_v   = bcdic['k']['dir']
							indbc_v   = bcdic['k']['indbc']	
							bc = True		


						dname =  genVname('d1_'+vname+'_d'+v+'_'+str(i)+'_', i=indi,j=indj,k=indk)
						[exprnbg,vnamebg] = genNbg(subptr[0],v,stencil,i=indi,j=indj,k=indk,
											vname='d1_'+vname+'_d'+v+'_'+str(i)+'_',dirBc=dirBc_v,indbc=indbc_v)
						locvar.append(vnamebg)
		
						createNbg(vnamebg,exprnbg,output,indbc=indbc_v,bc=bc)
						output.write(der(vnamebg,order,stencil,varname=dname,dirBC=dirBc_v,indbc=indbc_v,bc=bc))
						locvar.append([dname])
						if   v == 'x':
							output.write(dname + ' = '+ dname +'*param_float(1)'+'\n\n')
						elif v == 'y':	
							output.write(dname + ' = '+ dname +'*param_float(2)'+'\n\n')
						elif v == 'z':
							output.write(dname + ' = '+ dname +'*param_float(3)'+'\n\n')


						history[v]['var'].append('d1_'+vname+'_d'+v+'_'+str(i)+'_')

					if(ihist != None): 
						subptr = history[v]['var'][ihist]						
					else:	
						subptr = 'd1_'+vname+'_d'+v+'_'+str(iprv_d)+'_'													
					subdname.append(subptr)			
							# 
	
				for i,derexp in enumerate(re.finditer(der1n[v], exp)):		
					ptr = re.escape(derexp.group(0))
					exp = re.sub(ptr,subdname[i] ,exp)

		return [exp,locvar,history]

def genSymbDer2_bc(bcdic,exp,output,locvar,
				   order=2,stencil=3,
				   indi='i',indj='j',indk='k',
				   vname='symb',
				   history={'x':[[],[]],'y':[[],[]],'z':[[],[]]}):

		from genRhs import dim

		dirBc   = '    '
		bclayer = int((stencil-1)/2)

		vname = vname.strip()
		indiri = indi
		indirj = indj
		indirk = indk

		for v in dirvars:

				subdname = []							
				for i,derexp in enumerate(re.finditer(der2[v+v], exp)):	

					iprv_d = i
					ihist  = None

					for itprv_d,prv_d in enumerate(list(re.finditer(der2[v+v], exp))[:i]):						
						if prv_d.group(0).replace(" ", "") == derexp.group(0).replace(" ", ""):							
							iprv_d = itprv_d
							break

					# if history[dir][v][1] != []:						
					for ithist,exphist in enumerate(history[v][1]):
							if(derexp.group(0).replace(" ", "") == exphist.replace(" ", "")):								
								ihist = ithist
								break

					if iprv_d == i and (ihist == None):	 # new derivative (not already computed...)
						history[v][1].append(derexp.group(0))

						vlocnam = 'd2_'+vname+'_d'+v+'_'+str(i)+'_'

						if(dim < 3):
							if (v=='z' and derexp.group(0) != []):
								nsymb = list(re.finditer(der2n[v+v], exp))
								exception('wrong symbolic equation : ask for z derivatives ('+str(nsymb[i].group(0))+') with dim ='+str(dim),message='error')

						elif(dim < 2):		
							if  (v=='y' and derexp.group(0) != []):
								nsymb = list(re.finditer(der2n[v+v], exp))
								exception('wrong symbolic equation : ask for y derivatives ('+str(nsymb[i].group(0))+') with dim ='+str(dim),message='error')

	
						hlo = int((stencil-1)/2)

						bclayer = hlo
						bc = False
										
						if ( v == 'x') and (bcdic['i']['dir'] != None):
							dirBc = bcdic['i']['dir']							
							bclayer = min(bcdic['i']['indbc'],hlo)
							bc = True									
						if ( v == 'y') and (bcdic['j']['dir'] != None):
							dirBc = bcdic['j']['dir']
							bclayer = min(bcdic['j']['indbc'],hlo)	
							bc = True	
						if ( v == 'z') and (bcdic['k']['dir'] != None):
							dirBc = bcdic['k']['dir']
							bclayer = min(bcdic['k']['indbc'],hlo)	
							bc = True	
			
						if (np.mod((stencil-1),2)==1) : 
							raise ValueError('Wrong stencil length for centrered finite differences (should be odd value)  ! ',stencil)
			
	
						dname = genVname(vlocnam,   i=indiri,j=indirj,k=indirk)
						[exprnbg,vnamebg] = genNbg(derexp.group(0), v , stencil, i=indiri,j=indirj,k=indirk,vname=dname+'_',dirBc=dirBc,indbc=bclayer,der='second')
						locvar.append(vnamebg)
						createNbg(vnamebg,exprnbg,output,'second',indbc=bclayer,bc=bc)
						output.write(dder(vnamebg,order,stencil,varname=dname,indbc=bclayer,bc=bc))
						locvar.append([dname])
	
						if   v == 'x':
							output.write(dname+' = '+dname.strip()+'*param_float(1)*param_float(1)'+'\n\n')
						elif v == 'y':	
							output.write(dname+' = '+dname.strip()+'*param_float(2)*param_float(2)'+'\n\n')
						elif v == 'z':
								output.write(dname+' = '+dname.strip()+'*param_float(3)*param_float(3)'+'\n\n')
					
						
						history[v][0].append(vlocnam)		
													

					if(ihist != None): 
						subptr = history[v][0][ihist]						
					else:	
						subptr = 'd2_'+vname+'_d'+v+'_'+str(iprv_d)+'_'												
					subdname.append(subptr)			
	

				for i,derexp in enumerate(re.finditer(der2n[v+v], exp)):		
					ptr = re.escape(derexp.group(0))					
					exp = re.sub(ptr,subdname[i] ,exp)

				# print(history)	
		return [exp,locvar,history]

def rhsinfo(rhs):

	genBC_calls(rhs)

	wp     	     = rhs.wp     	     
	dim     	 = rhs.dim     	     
	stencil 	 = rhs.stencil 	     
	order        = rhs.order         
	coefficients = rhs.coefficients  
	varname      = rhs.varname   
	varsolved    = rhs.varsolved
	varstored    = rhs.varstored    
	hlo_rhs      = rhs.hlo_rhs  
	bc_info      = rhs.bc_info     
	varbc        = rhs.varbc


	instpath = os.environ['INSTALLPATH']
	rhsinf = open(instpath+'/pymod/rhsinfo.py','w')

	rhsinf.write('wp = \''+str(wp)+'\'\n')
	rhsinf.write('dim = '+str(dim)+'\n')
	rhsinf.write('stencil = '+str(stencil)+'\n')
	rhsinf.write('order = '+str(order)+'\n')
	rhsinf.write('coefficients = '+str(coefficients)+'\n')
	rhsinf.write('varname = '+str(varname)+'\n')
	rhsinf.write('varsolved = '+str(varsolved)+'\n')	
	rhsinf.write('varstored = '+str(varstored)+'\n')
	rhsinf.write('hlo_rhs = '+str(hlo_rhs)+'\n')
	rhsinf.write('bc_info = '+str(bc_info)+'\n')
	rhsinf.write('varbc = '+str(varbc)+'\n')

def gendtype():

	from genRhs import incPATH,wp

	dtype = open(incPATH+'/dtypes.h','w')

	if   wp == 'float64':
		dtype.write('integer,parameter :: wp = kind(0.0D0) ! working precision')
	elif wp == 'float32':
		dtype.write('integer,parameter :: wp = kind(0.0E0) ! working precision')
	else:
		import sys
		exception('gendtype) working precision not supported'+'\n'
		      '                  given value = '+str(wp)+' (set in genRhs.py)',message='error')


def genFilter(stencil,order,nvar,dirBC='',indbc='',fltbeg=2,rhs=None):

	from genRhs import incPATH

	dim     = rhs.dim

	if indbc == '':

		fd = fltDic[stencil][order]
	
		hlo = int((stencil-1)/2)
		
		for dir in ['x','y','z']:
			open(incPATH+'updateFilter_'+dir,'w') # set to empty
			open(incPATH+'Filter_'+dir,'w')       # set to empty
	
		fh = {}
		dirs = ['x']
		if dim == 2:
			dirs.append('y')		
		if dim == 3:
			dirs.append('y')	
			dirs.append('z')	
	
		for dir in dirs:
	
			updateflt = open(incPATH+'updateFilter_'+dir,'w')
			updateflt.write('')
		
			fh[dir] = open(incPATH+'Filter_'+dir,'w')
			fh[dir].write('')
	
			for nv in range(1,nvar+1):
				flt0 = ''
			
				if   dim == 1 and (dir in ['x']):
					updateflt.write('q(i,'+str(nv)+') = q(i,'+str(nv)+') - param_float(5)*q2(i,'+str(nv)+')\n\n')
					flt0 = 'q2(i,'+str(nv)+') ='	
				elif dim == 2 and (dir in ['x','y']):	
					updateflt.write('q(i,j,'+str(nv)+') = q(i,j,'+str(nv)+') - param_float(5)*q2(i,j,'+str(nv)+')\n\n')
					flt0 = 'q2(i,j,'+str(nv)+') ='	
				elif dim == 3 and (dir in ['x','y','z']):
					updateflt.write('q(i,j,k,'+str(nv)+') = q(i,j,k,'+str(nv)+') - param_float(5)*q2(i,j,k,'+str(nv)+')\n\n')
					flt0 = 'q2(i,j,k,'+str(nv)+') ='	
	
		
				if dir == 'x':
					flt = flt0						
					for shift in range(-hlo,hlo+1):
						if flt[-1]== '=':
							if dim == 1:
								flt = flt + str(fd[shift+hlo])+'_wp*'+'q(i'+'{:+d}'.format(shift)+','+str(nv)+') &\n            '
							elif dim == 2:	
								flt = flt + str(fd[shift+hlo])+'_wp*'+'q(i'+'{:+d}'.format(shift)+',j,'+str(nv)+') &\n            '
							else:
								flt = flt + str(fd[shift+hlo])+'_wp*'+'q(i'+'{:+d}'.format(shift)+',j,k,'+str(nv)+') &\n            '
						else:	
							if dim == 1:
								flt = flt + '{:+}'.format(fd[shift+hlo])+'_wp*'+'q(i'+'{:+d}'.format(shift)+','+str(nv)+') &\n            '
							elif dim == 2:
								flt = flt + '{:+}'.format(fd[shift+hlo])+'_wp*'+'q(i'+'{:+d}'.format(shift)+',j,'+str(nv)+') &\n            '
							else:	
								flt = flt + '{:+}'.format(fd[shift+hlo])+'_wp*'+'q(i'+'{:+d}'.format(shift)+',j,k,'+str(nv)+') &\n            '
	
				elif dir == 'y' and (dim >=2):	
					flt = flt0					
					for shift in range(-hlo,hlo+1):
						if flt[-1]== '=':
							if dim == 1:
								flt = ' '
							elif dim == 2:	
								flt = flt + str(fd[shift+hlo])+'_wp*'+'q(i,j'+'{:+d}'.format(shift)+','+str(nv)+') &\n            '
							else:
								flt = flt + str(fd[shift+hlo])+'_wp*'+'q(i,j'+'{:+d}'.format(shift)+',k,'+str(nv)+') &\n            '
						else:	
							if dim == 1:
								flt = ' '
							elif dim == 2:
								flt = flt + '{:+}'.format(fd[shift+hlo])+'_wp*'+'q(i,j'+'{:+d}'.format(shift)+','+str(nv)+') &\n            '
							else:	
								flt = flt + '{:+}'.format(fd[shift+hlo])+'_wp*'+'q(i,j'+'{:+d}'.format(shift)+',k,'+str(nv)+') &\n            '
	
				elif dir == 'z' and (dim >=3):
					flt = flt0						
					for shift in range(-hlo,hlo+1):
						if flt[-1]== '=':
							if dim == 1:
								flt = ' '
							elif dim == 2:	
								flt = ' '
							else:
								flt = flt + str(fd[shift+hlo])+'_wp*'+'q(i,j,k'+'{:+d}'.format(shift)+','+str(nv)+') &\n            '
						else:	
							if dim == 1:
								flt = ' '
							elif dim == 2:
								flt = ' '
							else:	
								flt = flt + '{:+}'.format(fd[shift+hlo])+'_wp*'+'q(i,j,k'+'{:+d}'.format(shift)+','+str(nv)+') &\n            '
								
				rhs.hlo_rhs = max(rhs.hlo_rhs,hlo)
						
				flt = flt.rstrip()
				fh[dir].write(flt[:-1]+'\n\n')
	else:
		
		axes  = {'i':'x','j':'y','k':'z'}

		# fltbeg = 2

		if indbc == fltbeg:		
			up    = open(incPATH+'update_filterbc_'+axes[dirBC[0]]+'.f90','w') # set to empty
			fltbc = open(incPATH+'filterbc_'+axes[dirBC[0]]+'.f90'       ,'w') # set to empty	
			if dim > 2:
				up.write('!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) \n')	
				fltbc.write('!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) \n')
			else:
				up.write('!$OMP DO SCHEDULE(GUIDED,4)  \n')	
				fltbc.write('!$OMP DO SCHEDULE(GUIDED,4)  \n')					

			loop_create('begin',up   , bc=dirBC[0])
			loop_create('begin',fltbc, bc=dirBC[0])

			up.write('\n')
			fltbc.write('\n')

		elif indbc > fltbeg:	
			up    = open(incPATH+'update_filterbc_'+axes[dirBC[0]]+'.f90','a+') # set to empty
			fltbc = open(incPATH+'filterbc_'+axes[dirBC[0]]+'.f90'       ,'a+') # set to empty		

		if rhs.stencil == 1:
			exception('RHS must be defined before BCs',message='error')

		else:	
			hlo_rhs = rhs.hlo_rhs		

		# flt BC:
		
		fd = fltDicBc[indbc]

		indi = 'i'
		indj = 'j'
		indk = 'k'

		indiq = 'i'
		indjq = 'j'
		indkq = 'k'
		dirscan = 1
		if dirBC == 'i1': 
			indi = '1-'+str(hlo_rhs)
			indiq= indi+'shiftmp'
			indi = indi +'layertmp'

		if dirBC == 'imax': 
			indi = 'nx+'+str(hlo_rhs)
			indiq= indi+'shiftmp'
			indi = indi +'layertmp'
			dirscan = -1

		if dirBC == 'j1': 
			indj = '1-'+str(hlo_rhs)
			indjq= indj +'shiftmp'
			indj = indj +'layertmp'
	
		if dirBC == 'jmax': 
			indj = 'ny+'+str(hlo_rhs)
			indjq= indj +'shiftmp'
			indj = indj +'layertmp'
			dirscan = -1

		if dirBC == 'k1': 
			indk = '1-'+str(hlo_rhs)
			indkq= indk+'shiftmp'				
			indk = indk +'layertmp'
		
		if dirBC == 'kmax': 
			indk = 'nz+'+str(hlo_rhs)
			indkq= indk +'shiftmp'			
			indk = indk +'layertmp'
			dirscan = -1

		if dim == 3:	
			index  = indi+','+indj+','+indk+','+'nvtmp'	
			indexq = indiq+','+indjq+','+indkq+','+'nvtmp'	
		elif dim == 2:
			index  = indi+','+indj+','+'nvtmp'	
			indexq = indiq+','+indjq+','+'nvtmp'				
		elif dim == 1:
			index  = indi+','+'nvtmp'		
			indexq = indiq+','+'nvtmp'

		for nv in range(1,nvar+1):
			flt0 = ''
			indexout  = index.replace('nvtmp',str(nv))
			indexqout = indexq.replace('nvtmp',str(nv))
			indexout  = indexout.replace('layertmp','{:+d}'.format(indbc*dirscan))

			up.write('q('+indexout+') = q('+indexout+') - param_float(5)*q2('+indexout+')\n')
			flt0 =  'q2('+indexout+') ='		
			for shift in range(0,len(fd)):
				indexqoutout = indexqout.replace('shiftmp','{:+d}'.format(shift*dirscan))

				if flt0[-1] == '=':
					flt0 = flt0 + str(fd[shift])+'_wp*'+'q('+indexqoutout+') &\n            '
				else:
					flt0 = flt0 + '{:+}'.format(fd[shift])+'_wp*'+'q('+indexqoutout+') &\n            '	
			flt0 = flt0.rstrip()
			fltbc.write(flt0[:-1]+'\n\n')	

		up.write('\n')
		if indbc == hlo_rhs-1:
			loop_create('end',up   , bc=dirBC[0])
			loop_create('end',fltbc, bc=dirBC[0])
			up.write('\n'+'!$OMP END DO \n')
			fltbc.write('\n'+'!$OMP END DO \n')
		up.close()
		fltbc.close()	

def genrk3(nvar,rhs=None,bc=[False,[]],rk3=None):

	from genRhs import incPATH

	dim = rhs.dim
	
	if not rk3:
		rk3 = open(incPATH+'includeRK3.f90','w')

	modvar = range(1,nvar+1)
	if bc[0]:
		modvar = bc[1]	
	
	if dim == 3:
		for nv in modvar:
			rk3.write('q2(i,j,k,'+str(nv)+') = q1(i,j,k,'+str(nv)+')'+'\n')
			rk3.write('q1(i,j,k,'+str(nv)+') = param_float(fadrTIMESTEP)*rk2(nrk)*rhs(i,j,k,'+str(nv)+') + q1(i,j,k,'+str(nv)+')'+'\n')
			#rk3.write('q (i,j,k,'+str(nv)+') = param_float(fadrTIMESTEP)*rk1(nrk)*rhs(i,j,k,'+str(nv)+') + q2(i,j,k,'+str(nv)+')'+'\n\n')
			
		# genP2C(output,'primitive')


	elif dim == 2 :
		for nv in modvar:
			rk3.write('q2(i,j,'+str(nv)+') = q1(i,j,'+str(nv)+')'+'\n')
			rk3.write('q1(i,j,'+str(nv)+') = param_float(fadrTIMESTEP)*rk2(nrk)*rhs(i,j,'+str(nv)+') + q1(i,j,'+str(nv)+')'+'\n')
			#rk3.write('q (i,j,'+str(nv)+') = param_float(fadrTIMESTEP)*rk1(nrk)*rhs(i,j,'+str(nv)+') + q2(i,j,'+str(nv)+')'+'\n\n')

		# genP2C(output,'primitive')	
	else:	
		for nv in modvar:
			rk3.write('q2(i,'+str(nv)+') = q1(i,'+str(nv)+')'+'\n')
			rk3.write('q1(i,'+str(nv)+') = param_float(fadrTIMESTEP)*rk2(nrk)*rhs(i,'+str(nv)+') + q1(i,'+str(nv)+')'+'\n')
			#rk3.write('q (i,'+str(nv)+') = param_float(fadrTIMESTEP)*rk1(nrk)*rhs(i,'+str(nv)+') + q2(i,'+str(nv)+')'+'\n\n')	

		# genP2C(output,'primitive')

def genrk3update(nvar,rhs=None,bc=[False,[]], updaterk3=None):

	from genRhs import incPATH
	
	dim     = rhs.dim
	consvar = rhs.consvar

	if not updaterk3:
		updaterk3   = open(incPATH+'includeRK3update.f90','w')
	
	modvar = range(1,nvar+1)
	if bc[0]:
		modvar = bc[1]

	if dim == 3:
		for nv in modvar:
			updaterk3.write('q (i,j,k,'+str(nv)+') = param_float(fadrTIMESTEP)*rk1(nrk)*rhs(i,j,k,'+str(nv)+') + q2(i,j,k,'+str(nv)+')'+'\n\n')
			
		if consvar != []: genP2C(updaterk3,'primitive',rhs=rhs,bc=bc)


	elif dim == 2 :
		for nv in modvar:
			updaterk3.write('q (i,j,'+str(nv)+') = param_float(fadrTIMESTEP)*rk1(nrk)*rhs(i,j,'+str(nv)+') + q2(i,j,'+str(nv)+')'+'\n\n')

		if consvar != []: genP2C(updaterk3,'primitive',rhs=rhs,bc=bc)
	else:	
		for nv in modvar:
			updaterk3.write('q (i,'+str(nv)+') = param_float(fadrTIMESTEP)*rk1(nrk)*rhs(i,'+str(nv)+') + q2(i,'+str(nv)+')'+'\n\n')	

		if consvar != []: genP2C(updaterk3,'primitive',rhs=rhs,bc=bc)

def genP2C(output,conserv,rhs=None,bc=[False,[]]):

	dim       = rhs.dim
	varsolved = rhs.varsolved
	varname   = rhs.varname
	consvar   = rhs.consvar

	modvar = range(1,len(varsolved)+1)

	if bc[0]:
		modvar = bc[1]

	if conserv == 'conservative':
		if dim == 3:			
			output.write('q1(i,j,k,' +str(varname['rho'])+') = q(i,j,k,'+str(varname['rho'])+')'+'\n')
			for nv in modvar:
				if nv in consvar:
					output.write('q1(i,j,k,' +str(nv)+') = q(i,j,k,'+str(nv)+')*q(i,j,k,'+str(varname['rho'])+')'+'\n')
				elif nv!=varname['rho'] :
					output.write('q1(i,j,k,' +str(nv)+') = q(i,j,k,'+str(nv)+')'+'\n')
		elif dim == 2:
			output.write('q1(i,j,' +str(varname['rho'])+') = q(i,j,'+str(varname['rho'])+')'+'\n')
			for nv in modvar:
				if nv in consvar:
					output.write('q1(i,j,' +str(nv)+') = q(i,j,'+str(nv)+')*q(i,j,'+str(varname['rho'])+')'+'\n')
				elif nv!=varname['rho'] :
					output.write('q1(i,j,' +str(nv)+') = q(i,j,'+str(nv)+')'+'\n')

		elif dim == 1:	
			output.write('q1(i,' +str(varname['rho'])+') = q(i,'+str(varname['rho'])+')'+'\n')
			for nv in modvar:
				if nv in consvar:
					output.write('q1(i,' +str(nv)+') = q(i,'+str(nv)+')*q(i,'+str(varname['rho'])+')'+'\n')
				elif nv!=varname['rho'] :
					output.write('q1(i,' +str(nv)+') = q(i,'+str(nv)+')'+'\n')
	
	elif conserv == 'primitive':
		if dim == 3:
			output.write('one_over_rho = 1.0_wp/q(i,j,k,'+str(varname['rho'])+')\n')
			for nv in modvar:
				if nv in consvar: output.write('q(i,j,k,'+str(nv)+') = q(i,j,k,'+str(nv)+')*one_over_rho'+'\n')	
		elif dim == 2:
			output.write('one_over_rho = 1.0_wp/q(i,j,'+str(varname['rho'])+')\n')
			for nv in modvar:
				if nv in consvar: output.write('q(i,j,'+str(nv)+') = q(i,j,'+str(nv)+')*one_over_rho'+'\n')
		elif dim == 1:
			output.write('one_over_rho = 1.0_wp/q(i,'+str(varname['rho'])+')\n')	
			for nv in modvar:
				if nv in consvar: output.write('q(i,'+str(nv)+') = q(i,'+str(nv)+')*one_over_rho'+'\n')

	elif conserv == 'standard':
		if dim == 3:			
			for nv in varsolved:
				output.write('q1(i,j,k,'+str(varname[nv])+') = q(i,j,k,'+str(varname[nv])+')'+'\n')			
	
		elif dim == 2:
			for nv in varsolved:
				output.write('q1(i,j,'+str(varname[nv])+') = q(i,j,'+str(varname[nv])+')'+'\n')

		elif dim == 1:	
			for nv in varsolved:
				output.write('q1(i,'+str(varname[nv])+') = q(i,'+str(varname[nv])+')'+'\n')

def geninit(output,nvar,rhs=None):

	dim = rhs.dim

	if dim == 3:
		for nv in range(1,nvar+1):
			output.write(  'q(i,j,k,'+str(nv)+') = 0.0_wp'+'\n')
			output.write( 'q1(i,j,k,'+str(nv)+') = 0.0_wp'+'\n')
			output.write( 'q2(i,j,k,'+str(nv)+') = 0.0_wp'+'\n')
			output.write('rhs(i,j,k,'+str(nv)+') = 0.0_wp'+'\n')

	elif dim == 2 :
		for nv in range(1,nvar+1):
			output.write(  'q(i,j,'+str(nv)+') = 0.0_wp'+'\n')
			output.write( 'q1(i,j,'+str(nv)+') = 0.0_wp'+'\n')
			output.write( 'q2(i,j,'+str(nv)+') = 0.0_wp'+'\n')
			output.write('rhs(i,j,'+str(nv)+') = 0.0_wp'+'\n')

	else:	
		for nv in range(1,nvar+1):
			output.write(  'q(i,'+str(nv)+') = 0.0_wp'+'\n')
			output.write( 'q1(i,'+str(nv)+') = 0.0_wp'+'\n')
			output.write( 'q2(i,'+str(nv)+') = 0.0_wp'+'\n')
			output.write('rhs(i,'+str(nv)+') = 0.0_wp'+'\n')

def genbcsrc(nvar,rhs=None):
	
	from genRhs import incPATH

	if rhs == None:
		exception("BC can't be generated before RHS",message='error')

	else:	
		hlo_rhs = rhs.hlo_rhs
		order   = rhs.order
		stencil = rhs.stencil 
		dim     = rhs.dim

	fh = {}
	
	for dir in ['x','y','z']:

		fh[dir] = open(incPATH+'periodic_'+dir,'w')
		fh[dir].write('')

		for nv in range(1,nvar+1):					
			bcexp_m = ''
			bcexp_p = ''

			if dir == 'x':
						
				for shift in range(0,hlo_rhs):					
					if dim == 1:
						bcexp_m = bcexp_m + 'q(0'+'{:+d}'.format(-shift)+','+str(nv)+') = ' + 'q(nx'+'{:+d}'.format(-shift)+','+str(nv)+')\n'
						bcexp_p = bcexp_p + 'q(nx+1'+'{:+d}'.format(shift)+','+str(nv)+') = ' + 'q(1'+'{:+d}'.format(shift)+','+str(nv)+')\n'
					elif dim == 2:
						bcexp_m = bcexp_m + 'q(0'+'{:+d}'.format(-shift)+',j,'+str(nv)+') = ' + 'q(nx'+'{:+d}'.format(-shift)+',j,'+str(nv)+')\n'
						bcexp_p = bcexp_p + 'q(nx+1'+'{:+d}'.format(shift)+',j,'+str(nv)+') = ' + 'q(1'+'{:+d}'.format(shift)+',j,'+str(nv)+')\n'
					else:	
						bcexp_m = bcexp_m + 'q(0'+'{:+d}'.format(-shift)+',j,k,'+str(nv)+') = ' + 'q(nx'+'{:+d}'.format(-shift)+',j,k,'+str(nv)+')\n'
						bcexp_p = bcexp_p + 'q(nx+1'+'{:+d}'.format(shift)+',j,k,'+str(nv)+') = ' + 'q(1'+'{:+d}'.format(shift)+',j,k,'+str(nv)+')\n'

			if dir == 'y':
						
				for shift in range(0,hlo_rhs):										
					if dim == 2:
						bcexp_m = bcexp_m + 'q(i,0'   +'{:+d}'.format(-shift)+','+str(nv)+') = ' + 'q(i,ny'+'{:+d}'.format(-shift)+','+str(nv)+')\n'
						bcexp_p = bcexp_p + 'q(i,ny+1'+'{:+d}'.format(shift)+','+str(nv)+') = ' + 'q(i,1'+'{:+d}'.format(shift)+','+str(nv)+')\n'
					elif dim == 3:	
						bcexp_m = bcexp_m + 'q(i,0'   +'{:+d}'.format(-shift)+',k,'+str(nv)+') = ' + 'q(i,ny'+'{:+d}'.format(-shift)+',k,'+str(nv)+')\n'
						bcexp_p = bcexp_p + 'q(i,ny+1'+'{:+d}'.format(shift)+',k,'+str(nv)+') = ' + 'q(i,1'+'{:+d}'.format(shift)+',k,'+str(nv)+')\n'
			if dir == 'z':
						
				for shift in range(0,hlo_rhs):										
					if dim == 3:
						bcexp_m = bcexp_m + 'q(i,j,0'   +'{:+d}'.format(-shift)+','+str(nv)+') = ' + 'q(i,j,nz'+'{:+d}'.format(-shift)+','+str(nv)+')\n'
						bcexp_p = bcexp_p + 'q(i,j,nz+1'+'{:+d}'.format(shift)+','+str(nv)+') = ' + 'q(i,j,1'+'{:+d}'.format(shift)+','+str(nv)+')\n'

								
			fh[dir].write(bcexp_m[:-1]+'\n\n')
			fh[dir].write(bcexp_p[:-1]+'\n\n')

def create_bcsubroutine(fname,fctname,locvname,loopname,indrange,phy=None):
	
	from genRhs import incPATH
	
	if phy:

		bcScheme_template = open(incPATH[:-13]+'template_genbc_'+phy+'.for','r')
		
		lines = bcScheme_template.readlines()
		
		edit = [[],[]]					

		edit[0].append('BCNAME');edit[1].append(fctname)
		edit[0].append('!PHYBCLOCVAR_'+phy);edit[1].append('#include "'+locvname+'"')
		edit[0].append('!PHYBCLOOPS_'+phy);edit[1].append('#include "'+loopname+'"')		
		edit[0].append('!INDBC_RANGEi');edit[1].append(indrange['i'])
		edit[0].append('!INDBC_RANGEj');edit[1].append(indrange['j'])
		edit[0].append('!INDBC_RANGEk');edit[1].append(indrange['k'])

		fname.write("\n\n")
		for l in lines:		
			lnew = l		
			for i,ed in enumerate(edit[0]):						
				lnew = lnew.replace(ed,edit[1][i])						
			fname.write(lnew)
		
		fname.write("\n\n")			

	else:	

		bcScheme_template = open(incPATH[:-13]+'template_bcscheme.for','r')
		
		lines = bcScheme_template.readlines()
		
		edit = [[],[]]					
		
		edit[0].append('BCNAME');edit[1].append(fctname)
		edit[0].append('!boundarySchemeLOCVAR');edit[1].append('#include "'+locvname+'"')
		edit[0].append('!boundarySchemeLOOPS');edit[1].append('#include "'+loopname+'"')				
		edit[0].append('!INDBC_RANGEi');edit[1].append(indrange['i'])
		edit[0].append('!INDBC_RANGEj');edit[1].append(indrange['j'])
		edit[0].append('!INDBC_RANGEk');edit[1].append(indrange['k'])
		
		fname.write("\n\n")
		for l in lines:		
			lnew = l		
			for i,ed in enumerate(edit[0]):						
				lnew = lnew.replace(ed,edit[1][i])						
			fname.write(lnew)
		
		fname.write("\n\n")		

def create_bccalls(fname,fctname,fctcall):
	
	from genRhs import incPATH
	
	bcScheme_template = open(incPATH[:-13]+'template_bccall.for','r')
	
	lines = bcScheme_template.readlines()
	
	edit = [[],[]]					
	
	edit[0].append('BCNAME');edit[1].append(fctname)
	edit[0].append('!boundarySchemeCALL');edit[1].append(fctcall)

	
	fname.write("\n\n")
	for l in lines:		
		lnew = l		
		for i,ed in enumerate(edit[0]):						
			lnew = lnew.replace(ed,edit[1][i])						
		fname.write(lnew)
	
	fname.write("\n\n")			

def gen_eqns_bc(Eqns,output,localvar,
	            eqname,Order=2,Stencil=3,
	            indi    ='i',indj = 'j',indk = 'k',
	            DirDic  = {'i':{'dir':None,'indbc':None},'j':{'dir':None,'indbc':None},'k':{'dir':None,'indbc':None}},
	            vname   = 'symb',
	            update  = False,
	            updateq = False,
	            updatest= False,
	            updateqbc=False):
					
					from genRhs import dim
					
					indiri = indi
					indirj = indj
					indirk = indk

					locvar = []
		
					output.write(comment('Start building layers for BC : '+str(DirDic['i']['dir'])+' '+str(DirDic['j']['dir'])+' '+str(DirDic['k']['dir'])))

					output.write(comment('BC layer: '+str(DirDic['i']['indbc'])+' '+str(DirDic['j']['indbc'])+' '+str(DirDic['k']['indbc'])))

					bc     = 'all'
					edge   = 'all'
					corner = 'all'

					fbc     = False
					fedge   = False
					fcorner = False

					for d in DirDic:
						if DirDic[d]['dir'] != None:
							if not fbc:
								bc  = DirDic[d]['dir'][0]
								fbc = True
							elif not fedge:
								edge  = DirDic[d]['dir'][0]
								fedge = True	
							elif not fcorner:
								corner  = DirDic[d]['dir'][0]
								fcorner = True																	


					loop_create('begin',output,bc=bc,edge=edge,corner=corner)
					
					history={'x':{'var':[],'symb':[]},'y':{'var':[],'symb':[]},'z':{'var':[],'symb':[]}}
					history2={'x':[[],[]],'y':[[],[]],'z':[[],[]]}
					dhistory={'x':[[],[],[],[]],'y':[[],[],[],[]],'z':[[],[],[],[]]}
	
					if   DirDic['i']['dir'] == 'i1':
						indiri = indi+'{:+d}'.format(DirDic['i']['indbc'])
					if DirDic['j']['dir'] == 'j1' and (dim >= 2):
						indirj = indj+'{:+d}'.format(DirDic['j']['indbc'])
					if DirDic['k']['dir'] == 'k1' and (dim >= 3):
						indirk = indk+'{:+d}'.format(DirDic['k']['indbc'])	
					
					if   DirDic['i']['dir'] == 'imax':
						indiri = indi+'{:+d}'.format(-DirDic['i']['indbc'])
					if DirDic['j']['dir'] == 'jmax' and (dim >= 2):
						indirj = indj+'{:+d}'.format(-DirDic['j']['indbc'])
					if DirDic['k']['dir'] == 'kmax' and (dim >= 3):
						indirk = indk+'{:+d}'.format(-DirDic['k']['indbc'])	

					for eqn in  Eqns.keys():
						if Eqns[eqn].replace(' ',''):
							vnametmp = vname[eqn].strip() # '_' + str(bcnum)
	
							output.write(comment('building source terms in RHS for layer '+str(DirDic['i']['indbc'])+' '+str(DirDic['j']['indbc'])+' '+str(DirDic['k']['indbc'])+' '+eqname[eqn]))
						
							op  = Eqns[eqn].replace(" ", "")+'\n'
							
							output.write('!'.ljust(60,'~')+'\n')
							output.write('!'+'\n')
							output.write('! '+op)
							output.write('!'+'\n')	
							output.write('!'.ljust(60,'~')+'\n\n')	
							
							# generates BC layers :			


							[Out,locvar,history]  = genSymbDer1_bc(DirDic,Eqns[eqn],output,locvar,order=Order,stencil=Stencil,indi=indiri,indj=indirj,indk=indirk,vname=vnametmp,history=history,dhistory=dhistory)			
							[Out,locvar,history2] = genSymbDer2_bc(DirDic,Out,output,locvar,order=4,stencil=5,indi=indiri,indj=indirj,indk=indirk,vname=vnametmp,history=history2)
							
							output.write(comment('Update BC terms for layer '+str(DirDic['i']['indbc'])+' '+str(DirDic['j']['indbc'])+' '+str(DirDic['k']['indbc'])+' '+eqname[eqn]))
	
							if updateq:
								exp_rhs =  op_to_dNami(Out,i=indiri,j=indirj,k=indirk) 
								output.write(updateRHS(eqn,exp_rhs,i=indiri,j=indirj,k=indirk,update=update,updateq=True)+'\n\n')
							elif updatest:
								exp_stored = op_to_dNami(Out,i=indiri,j=indirj,k=indirk) 
								output.write(updateStored(eqn,exp_stored,i=indiri,j=indirj,k=indirk,update=update)+'\n\n')	
							elif updateqbc:
								exp_qbc = op_to_dNami(Out,i=indiri,j=indirj,k=indirk) 
								output.write(updateVarbc(eqn,exp_qbc,i=indiri,j=indirj,k=indirk,update=update)+'\n\n')									
							else:
								exp_rhs = ' - '
								exp_rhs = exp_rhs + ' ( ' + op_to_dNami(Out,i=indiri,j=indirj,k=indirk) + ' ) '
								output.write(updateRHS(eqn,exp_rhs,i=indiri,j=indirj,k=indirk,update=update)+'\n\n')

					loop_create('end',output,bc=bc,edge=edge,corner=corner)
			
					tmpvar = ''
					for var1 in locvar: 
						for var2 in var1:
							if not tmpvar: 
								tmpvar = tmpvar + ' ' +var2
							else:
								tmpvar = tmpvar + ',' + var2
								
						tmpvar = tmpvar + ' &\n            '
					
					tmpvar = tmpvar.rstrip()
					if(tmpvar !=''):
						localvar.write('\n\n real(wp) :: ')						
						localvar.write(tmpvar[:-1])

						



def color(str,message='com'):
	if message == 'com':
		return '\033[1;40;94m'+str+'\033[0m'
	elif message == 'error':
		return '\033[1;40;41m'+str+'\033[0m'

def exception(str,message='com'):
	if message == 'com':
		print(color('[info]'+str,message='com'))
	elif message == 'error':
		print(color('[error] '+str,message='error'))
		import sys
		sys.exit()		
