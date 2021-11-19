# =============================================================================
# 2D NS Curvilinear Equations ig/vdw
# =============================================================================

# number of dimensions
dim = 2


#coefficients ////////////////////////////////////////////////////////////////

coefficients = {'mu'         : 1, 
                'kappa'      : 2,  
                'gamma_m1'   : 3,
                'CvInv'      : 4,
                'Cv'         : 5,     
                'u_wall'     : 6,
                'T_0'        : 7,
                'gamma'      : 8,
                'u_0'        : 9,
                'v_0'        :10, 
                'Rgas'       :11,
                'p_0'        :12,
                'sigma'      :13,
                'L_x'        :14,
                'L_y'        :15,
                'eta_1'      :16,
                'eta_2'      :17,
                'eta_3'      :18  
                }                


# unknowns to march in time ///////////////////////////////////////////////////
varname      = {'rho' : 1,
                'u'   : 2,
                'v'   : 3,
                'et'  : 4, 
                }

varsolved = ['rho','u','v','et']

# of which
consvar      = [2,3,4] # are conservative variables

# derived local variables /////////////////////////////////////////////////////

# IG EOS:

varloc       = {'p': 'gamma_m1*rho*(e)',
                'e': '(et-0.5_wp*(u*u+v*v))',#+w*w))',
                'T': 'CvInv*(e)',
                'Utilde':'(  detady * u - dksidy * v ) ',
                'Vtilde':'( -detadx * u + dksidx * v ) ',
                'x_xi' : '    (J*detady)  ',            
                'x_eta': ' (- J*dksidy) ',              
                'y_xi' : ' (- J*detadx) ',             
                'y_eta': '   (J*dksidx)  '}

# (x,y,z)       --> Computational Domain Coordinate System (arbitrary curvilinear space)
# (xi,eta,zeta) --> Original Cartesian Coordinate System

# -- Same stored var for both

varstored    = {
         'ksi' : {'symb': 'ksi', 
                         'ind':1 ,
                         'static': True},
         'eta' : {'symb': 'eta', 
                         'ind':2 ,
                         'static': True},
         'dksidx' : {'symb': ' [ ksi ]_1x ', 
                         'ind':3 ,
                         'static': True},
         'dksidy' : {'symb': ' [ ksi ]_1y ', 
                         'ind':4 , 
                         'static': True},
         'detadx' : {'symb': ' [ eta ]_1x ', 
                         'ind':5 , 
                         'static': True},
         'detady' : {'symb': ' [ eta ]_1y ', 
                         'ind':6 ,
                         'static': True},
         'J'      : {'symb': '1.0_wp  / ( dksidx * detady - detadx * dksidy  )', 
                         'ind':7,
                         'static': True},
         'omg'    : {'symb': ' J * ( detady * [v]_1x - detadx * [v]_1y + dksidy * [u]_1x - dksidx * [u]_1y  ) ', 
                         'ind':8,
                         'static': False},

         'dJx_xidx' : {'symb': ' [ {eta}_1y ]_1x ', 
                         'ind':9 ,
                         'static': True},
         'dJy_xidy' : {'symb': ' - [ {eta}_1x ]_1y ', 
                         'ind':10 ,
                         'static': True},
         'P' : {'symb': 'p', 
                         'ind':11 ,
                         'static': False}                                                         

         # 'x_xi'   : {'symb': ' detady ', 
         #                 'ind': 9, 
         #                 'static': True}, 
         # 'x_eta'  : {'symb': ' - dksidy ', 
         #                 'ind': 10, 
         #                 'static': True},
         # 'y_xi'   : {'symb': ' - detadx ', 
         #                 'ind': 11, 
         #                 'static': True},
         # 'y_eta'   : {'symb': ' dksidx ', 
         #                 'ind': 12, 
         #                 'static': True},


         # 'Jm1'    : {'symb': '1.0_wp/J ', 
         #                 'ind':9,
         #                 'static': True},
         # 'Utilde'      : {'symb': '(  detady * u - detadx * v ) ', 
         #                 'ind':9 ,
         #                 'static': False},
         # 'Vtilde'      : {'symb': '( -dksidy * u + dksidx * v )  ', 
         #                 'ind':10 ,
         #                 'static': False},                
         }


# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file
rhsname = {'rho' : 'd(rho)/dt',
           'u'   : 'd(rho u)/dt',
           'v'   : 'd(rho v)/dt',
           'et'  : 'd(rho et)/dt'}

# .. name tags to use for intermediate variables created by the constructor

vnamesrc_divFconv = {'rho'  : 'FluRconv',
                     'u'    : 'FluXconv',
                     'v'    : 'FluYconv',
                     'et'   : 'FluEconv'}
  
vnamesrc_divFdif = {'rho'  : 'FluRdif',
                    'u'    : 'FluXdif',
                    'v'    : 'FluYdif',
                    'et'   : 'FluEdif'}

vnamesrc_bc  = {'rho': 'bc_rho',
                'u'  : 'bc_u',
                'v'  : 'bc_v',
                # 'w'  : 'bc_w',
                'et' : 'bc_et  '}                      

# RHS terms ///////////////////////////////////////////////////////////////////

divFx = {
    'rho' : ' J  * ( [ rho * Utilde                   ]_1x )', 
    'u'   : ' J  * ( [ rho * u * Utilde + detady*p    ]_1x )', 
    'v'   : ' J  * ( [ rho * v * Utilde - dksidy*p    ]_1x )', 
    'et'  : ' J  * ( [ (rho * et  + p )  *  Utilde    ]_1x )', 
    }

divFy = {
    'rho' : ' J  * ( [ rho * Vtilde                 ]_1y )', 
    'u'   : ' J  * ( [ rho * u * Vtilde - detadx*p  ]_1y )', 
    'v'   : ' J  * ( [ rho * v * Vtilde + dksidx*p  ]_1y )', 
    'et'  : ' J  * ( [ (rho * et  + p )  *  Vtilde  ]_1y )', 
    }


# -- Assemble both 

divF = {}

for key in divFx.keys():
    divF[key] = divFx[key] + ' + ' + divFy[key]

from genNSBC import LODI, dNami_to_sympy
import sympy as sym
sym.init_printing(use_latex=True,wrap_line=False)

Ftild      = {}

Ftild['x'] = divFx
Ftild['y'] = divFy

# RHS Diffusion //////////////////////////////////////////////////////////////////////////

from NS_curvi2D import dNamiFlx_x,dNamiFlx_y


divViscFx = {
    #'rho' : ' - J  * ( [ ' + dNamiFlx_x[0] + '  ]_1x )', 
    'u'   : '  - J  * ( [ J**(-1)*(' + dNamiFlx_x[1] + ')  ]_1x )', 
    'v'   : '  - J  * ( [ J**(-1)*(' + dNamiFlx_x[2] + ')  ]_1x )', 
    'et'  : '  - J  * ( [ J**(-1)*(' + dNamiFlx_x[3] + ')  ]_1x )', 
    }

divViscFy = {
    #'rho' : ' - J  * ( [ ' + dNamiFlx_y[0] + '  ]_1y )', 
    'u'   : ' - J  * ( [ J**(-1)*( ' + dNamiFlx_y[1] + ')  ]_1y )', 
    'v'   : ' - J  * ( [ J**(-1)*( ' + dNamiFlx_y[2] + ')  ]_1y )', 
    'et'  : ' - J  * ( [ J**(-1)*( ' + dNamiFlx_y[3] + ')  ]_1y )', 
    }


# -- Assemble both 

divViscF = {}

for key in ['u','v','et']:
    divViscF[key] = divViscFx[key] + ' + ' + divViscFy[key]

 
for k in divViscFx:
    Ftild['x'][k] = Ftild['x'][k] + ' + ('+divViscFx[k]+')'
    Ftild['y'][k] = Ftild['y'][k] + ' + ('+divViscFy[k]+')'

print("Done generating bulk equations.")


###########################################################################
#
# Physical Boundary Conditions:
#
###########################################################################


#template for Wall moving at u_0 in x direction

    # etcmp = {'j1'   : ' ( Cv*T_min + 0.5_wp * (u*u + v*v )  ) ' ,
    #          'jmax' : ' ( Cv*T_max + 0.5_wp * (u*u + v*v )  ) '}
    
    # PhyBcFx = {'rho' : 'rho*u_0         ',
    #            'u'   : 'rho*u_0*u_0  + p',
    #            'v'   : 'rho*u_0*v       ',
    #            'et'  : '(rho*et + p)*u_0 '}
    
    # PhyBcFy = {'rho' : 'rho*v         ',
    #            'u'   : 'rho*v*u_0       ', 
    #            'v'   : 'rho*v*v  + p  ', 
    #            'et'  : '(rho*et + p)*v '} 
                                
                                
    # PhyBcFz = {'rho' : 'rho*w         ',
    #            'u'   : 'rho*w*u_0       ',
    #            'v'   : 'rho*w*v       ',
    #            'et'  : '(rho*et + p)*w '}
    
    
    
    # # Src_phybc_rhs = {'rho' : '[ '+PhyBcFx['rho']+' ]_1x' + ' + ' + 'J*[ -detadx * u + dksidx * v ]_1y ' }
    #                   # 'u'   : '[ '+PhyBcFx['u']  +' ]_1x' + ' + ' + '[ '+PhyBcFy['u']  +' ]_1y' + ' + ' + '[ '+PhyBcFz['u']  +' ]_1z ',
    #                   # 'v'   : '[ '+PhyBcFx['v']  +' ]_1x' + ' + ' + '[ '+PhyBcFy['v']  +' ]_1y' + ' + ' + '[ '+PhyBcFz['v']  +' ]_1z ',
    #                   # 'w'   : '[ '+PhyBcFx['w']  +' ]_1x' + ' + ' + '[ '+PhyBcFy['w']  +' ]_1y' + ' + ' + '[ '+PhyBcFz['w']  +' ]_1z ',
    #                   #'et' : ' [ '+PhyBcFx['et'] +' '+ Fx['et'] +' ]_1x' + ' + ' + '[ '+PhyBcFy['et'] +' '+ Fy['et'] +' ]_1y ' }
    
    # Src_phybc_rhs = {'rho' : '  J * ( rho *[  Vtilde                 ]_1y + Vtilde * [ rho ]_1y )' }
    #                   # 'u'   : '[ '+PhyBcFx['u']  +' ]_1x' + ' + ' + '[ '+PhyBcFy['u']  +' ]_1y' + ' + ' + '[ '+PhyBcFz['u']  +' ]_1z ',
    #                   # 'v'   : '[ '+PhyBcFx['v']  +' ]_1x' + ' + ' + '[ '+PhyBcFy['v']  +' ]_1y' + ' + ' + '[ '+PhyBcFz['v']  +' ]_1z ',
    #                   # 'w'   : '[ '+PhyBcFx['w']  +' ]_1x' + ' + ' + '[ '+PhyBcFy['w']  +' ]_1y' + ' + ' + '[ '+PhyBcFz['w']  +' ]_1z ',
    #                   #'et' : ' [ '+PhyBcFx['et'] +' '+ Fx['et'] +' ]_1x' + ' + ' + '[ '+PhyBcFy['et'] +' '+ Fy['et'] +' ]_1y ' }
    
    
    # Src_phybc_qmax   = { 'u'   : 'u_0',
    #                      'v'   : '0.0_wp',
    #                      'et'  : etcmp['jmax']
    #                     }
    
    # Src_phybc_q1   = { 'u'   : 'u_0',
    #                    'v'   : '0.0_wp',
    #                   'et'   : etcmp['j1']}
    

# Characteristic Boundary Conditions:                  

from CharsForConsLaw import characteristics

Char = {} 
Char = characteristics('Euler')


print("Done computing characteristic equations in x and y direction.")


x_xi,y_xi,z_xi          = sym.symbols(['x_xi','y_xi','z_xi'],Real=True)    
x_eta,y_eta,z_eta       = sym.symbols(['x_eta','y_eta','z_eta'],Real=True)  
x_zeta,y_zeta,z_zeta    = sym.symbols(['x_zeta','y_zeta','z_zeta'],Real=True)  

J = sym.Symbol('J',real=True)

betax = sym.Symbol('beta_x',real=True,positive = True)
betay = sym.Symbol('beta_y',real=True,positive = True)

betaswp = {'x':[(x_xi**2+ x_eta**2),betax**2],'y':[(y_xi**2+ y_eta**2),betay**2]}


x      = sym.Symbol('x')
y      = sym.Symbol('y')
z      = sym.Symbol('z')
t      = sym.Symbol('t')

rho   = sym.Function('rho')(x,y,z,t)
u     = sym.Function('u')  (x,y,z,t)
v     = sym.Function('v')  (x,y,z,t)
w     = sym.Function('w')  (x,y,z,t)
et    = sym.Function('et') (x,y,z,t)
p     = sym.Function('p')  (x,y,z,t)

rho   = sym.Symbol('rho')
u     = sym.Symbol('u')  
v     = sym.Symbol('v')  
w     = sym.Symbol('w')  
et    = sym.Symbol('et') 
p     = sym.Symbol('p')  
gamma = sym.Symbol('gamma')

Q_CS = sym.Matrix([[rho],
                   [rho*u],
                   [rho*v],
                   [0.5*rho*(u*u+v*v)+1.0/(gamma-1)*p]])
                   # [rho*w],
                   # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])

Q    = sym.Matrix([[rho],
                   [u],
                   [v],
                   # [w],
                   [p]])


Pmatrix = Q_CS.jacobian(Q)

print("Done computing Prim/Cons transfer Matrix.")

from genNSBC import sympy2dNami

# i1 non-reflective outfow:

Li_BC_i1sym = Char['xi'][0].copy()
Li_BC_i1sym = Li_BC_i1sym.applyfunc(lambda x:sym.simplify(x).subs(betaswp['x'][0],betaswp['x'][1]))

Li_BC_i1 = []
for i in Li_BC_i1sym:
    Li_BC_i1.append(sympy2dNami(i))

  

# Li_BC_i1[3] = 'sigma*c*(1-M_i1*M_i1)/L_x*( p - p_0)'

# i1 inflow:

Li_BC_i1[3] = 'eta_1*rho*c*c*( 1.0_wp - M_i1*M_i1)/L_x*( u - u_0) '
Li_BC_i1[1] = 'eta_2*rho*c*Rgas*(T - T_0)/L_x '            
Li_BC_i1[0] = 'eta_3*c*(v - v_0)/L_x '                        

print('Li i1 done')

#imax non-reflecive outflow:

Li_BC_imaxsym = Char['xi'][0].copy()
Li_BC_imaxsym = Li_BC_imaxsym.applyfunc(lambda x:sym.simplify(x).subs(betaswp['x'][0],betaswp['x'][1]))

Li_BC_imax = []
for i in Li_BC_imaxsym:
    Li_BC_imax.append(sympy2dNami(i))

Li_BC_imax[2] = 'sigma*c*(1-M_imax*M_imax)/L_x*( p - p_0)'


print('Li imax done')

#j1 non-reflective outfow:

Mi_BC_j1sym = Char['eta'][0].copy()
Mi_BC_j1sym = Mi_BC_j1sym.applyfunc(lambda x:sym.simplify(x).subs(betaswp['y'][0],betaswp['y'][1]))

Mi_BC_j1 = []
for i in Mi_BC_j1sym:
    Mi_BC_j1.append(sympy2dNami(i))

# Mi_BC_j1[3] = 'sigma*c*(1-M_j1*M_j1)/L_y*( p - p_0)'

#j1 inflow:

Mi_BC_j1[3] = 'eta_1*rho*c*c*( 1.0_wp - M_j1*M_j1)/L_y*( v - v_0) '
Mi_BC_j1[1] = 'eta_2*rho*c*Rgas*(T - T_0)/L_y '           
Mi_BC_j1[0] = 'eta_3*c*(u - u_0)/L_y '                       

print('Mi j1 done')

#jmax non-reflective outfow:

Mi_BC_jmaxsym = Char['eta'][0].copy()
Mi_BC_jmaxsym = Mi_BC_jmaxsym.applyfunc(lambda x:sym.simplify(x).subs(betaswp['y'][0],betaswp['y'][1]))

Mi_BC_jmax = []
for i in Mi_BC_jmaxsym:
    Mi_BC_jmax.append(sympy2dNami(i))

Mi_BC_jmax[2] = 'sigma*c*(1-M_jmax*M_jmax)/L_y*( p - p_0)'

print('Mi jmax done')

curvi = {'x':'xi','y':'eta'}

Pleft_inv = {}
for d in ['x','y']:
  Pleft_inv[d] = Char[curvi[d]][2].inv()

def setLiNames(bc,dir='x'):  
    
    pre = 'L'

    if dir == 'y': pre = 'M'
    if dir == 'z': pre = 'N'

    L0 = sym.Symbol(pre+'0_'+bc) # u-c
    L1 = sym.Symbol(pre+'1_'+bc) # u+c
    L2 = sym.Symbol(pre+'2_'+bc) # v
    L3 = sym.Symbol(pre+'3_'+bc) # u

    return [L0,L1,L2,L3]

src_rhsBC = {}
tdir = {'x':'y','y':'x'}

lodiRHSdNami = {}
for bc in ['1','max']:
  for dir in ['x','y']:

    idir = 'i'
    if dir == 'y':idir = 'j'
    if dir == 'z':idir = 'j'
    

    [L0,L1,L2,L3] = setLiNames(idir+bc,dir)

    lodiRHS       = Pmatrix*(Pleft_inv[dir].applyfunc(lambda x: sym.simplify(sym.factor(x)))*sym.Matrix([L0,L1,L2,L3]))
    # lodiRHS       = (Pleft_inv[dir].applyfunc(lambda x: sym.simplify(sym.factor(x)))*sym.Matrix([L0,L1,L2,L3]))
     
    lodiRHS       = lodiRHS.applyfunc(lambda x:J*   (sym.expand(sym.simplify(x.subs(betaswp[dir][0],betaswp[dir][1]))))) #+ Ftild[tdir[dir]]

    lodiRHSdNami[idir+bc] = []
    for symb in lodiRHS:
        eqdNami = sympy2dNami(sym.simplify(symb)).replace('c','c_'+idir+bc)
        eqdNami = eqdNami.replace('beta_'+dir,'beta_'+idir+bc)        
        lodiRHSdNami[idir+bc].append(eqdNami)

    src_rhsBC[idir+bc] = {'rho':lodiRHSdNami[idir+bc][0] + '+ ('+Ftild[tdir[dir]]['rho']+' )',
                          'u'  :lodiRHSdNami[idir+bc][1] + '+ ('+Ftild[tdir[dir]]['u'  ]+' )',
                          'v'  :lodiRHSdNami[idir+bc][2] + '+ ('+Ftild[tdir[dir]]['v'  ]+' )',
                          # 'w'  :,.replace('c','c_'+idir+bc)
                          'et' :lodiRHSdNami[idir+bc][3] + '+ ('+Ftild[tdir[dir]]['et']+' )'}


src_rhsBC['i1j1'] = {'rho':lodiRHSdNami['i1'][0] + '+ ('+lodiRHSdNami['j1'][0]+' )',
                     'u'  :lodiRHSdNami['i1'][1] + '+ ('+lodiRHSdNami['j1'][1]+' )',
                     'v'  :lodiRHSdNami['i1'][2] + '+ ('+lodiRHSdNami['j1'][2]+' )',
                     # 'w'  :,.replace('c','c_'+idir+bc)
                     'et' :lodiRHSdNami['i1'][3] + '+ ('+lodiRHSdNami['j1'][3]+' )'}

src_rhsBC['i1jmax'] = {'rho':lodiRHSdNami['i1'][0] + '+ ('+lodiRHSdNami['jmax'][0]+' )',
                     'u'  :lodiRHSdNami['i1'][1] + '+ ('+lodiRHSdNami['jmax'][1]+' )',
                     'v'  :lodiRHSdNami['i1'][2] + '+ ('+lodiRHSdNami['jmax'][2]+' )',
                     # 'w'  :,.replace('c','c_'+idir+bc)
                     'et' :lodiRHSdNami['i1'][3] + '+ ('+lodiRHSdNami['jmax'][3]+' )'}     


src_rhsBC['imaxj1'] = {'rho':lodiRHSdNami['imax'][0] + '+ ('+lodiRHSdNami['j1'][0]+' )',
                     'u'  :lodiRHSdNami['imax'][1] + '+ ('+lodiRHSdNami['j1'][1]+' )',
                     'v'  :lodiRHSdNami['imax'][2] + '+ ('+lodiRHSdNami['j1'][2]+' )',
                     # 'w'  :,.replace('c','c_'+idir+bc)
                     'et' :lodiRHSdNami['imax'][3] + '+ ('+lodiRHSdNami['j1'][3]+' )'}

src_rhsBC['imaxjmax'] = {'rho':lodiRHSdNami['imax'][0] + '+ ('+lodiRHSdNami['jmax'][0]+' )',
                     'u'  :lodiRHSdNami['imax'][1] + '+ ('+lodiRHSdNami['jmax'][1]+' )',
                     'v'  :lodiRHSdNami['imax'][2] + '+ ('+lodiRHSdNami['jmax'][2]+' )',
                     # 'w'  :,.replace('c','c_'+idir+bc)
                     'et' :lodiRHSdNami['imax'][3] + '+ ('+lodiRHSdNami['jmax'][3]+' )'}


varbc = {'M_i1'   : {'symb'  : ' (J*(Utilde))/sqrt(gamma*Rgas*T) ', 
                    'ind'   : 1 ,
                    'static': False,
                    'face'  : 'i1'},                         
         'M_imax' : {'symb'  : ' (J*(Utilde))/sqrt(gamma*Rgas*T) ', 
                    'ind'   : 2 ,
                    'static': False,
                    'face'  : 'imax'},
         'M_j1'   : {'symb'  : ' (J*(Vtilde))/sqrt(gamma*Rgas*T) ', 
                    'ind'   : 1 ,
                    'static': False,
                    'face'  : 'j1'},                    
         'M_jmax' : {'symb'  : ' (J*(Vtilde))/sqrt(gamma*Rgas*T) ', 
                    'ind'   : 2 ,
                    'static': False,
                    'face'  : 'jmax'},
         'c_i1'   : {'symb'   : 'sqrt(gamma*p/rho) ', 
                    'ind'   : 3 ,
                    'static': False,
                    'face'  : 'i1'},             
         'c_imax'   : {'symb'   : 'sqrt(gamma*p/rho) ', 
                    'ind'   : 4 ,
                    'static': False,
                    'face'  : 'imax'},
         'c_j1'   : {'symb'   : 'sqrt(gamma*p/rho) ', 
                    'ind'   : 3 ,
                    'static': False,
                    'face'  : 'j1'},                   
         'c_jmax'   : {'symb'   : 'sqrt(gamma*p/rho) ', 
                    'ind'   : 4 ,
                    'static': False,
                    'face'  : 'jmax'},
         'beta_i1'  : {'symb': ' sqrt(  (x_xi)**2 + (x_eta)**2 )  ', 
                    'ind'   : 5 ,
                    'static': True,
                    'face'  : 'i1'},
         'beta_imax': {'symb': ' sqrt(  (x_xi)**2 + (x_eta)**2 )  ', 
                    'ind'   : 6 ,
                    'static': True,
                    'face'  : 'imax'},                    
         'beta_j1'  : {'symb': ' sqrt(  (y_xi)**2 + (y_eta)**2 )  ', 
                    'ind':    5 ,
                    'static': True,
                    'face'  : 'j1'},                              
         'beta_jmax': {'symb': ' sqrt(  (y_xi)**2 + (y_eta)**2 )  ', 
                    'ind':    6 ,
                    'static': True,
                    'face'  : 'jmax'},

          # Li's definition i1: 
         'L0_i1' : {'symb'   : Li_BC_i1[0].replace('c','c_i1').replace('beta_x','beta_i1 ')+' ', 
                     'ind'   : 7,
                     'static': False,
                     'face'  : 'i1'}, 
            
         'L1_i1' : {'symb'   : Li_BC_i1[1].replace('c','c_i1').replace('beta_x','beta_i1 ')+' ', 
                     'ind'   : 8 ,
                     'static': False,
                     'face'  : 'i1'},
         'L2_i1' : {'symb'   : Li_BC_i1[2].replace('c','c_i1').replace('beta_x','beta_i1 ')+' ', 
                     'ind'   : 9 ,
                     'static': False,
                     'face'  : 'i1'},            
         'L3_i1' : {'symb'   : Li_BC_i1[3].replace('c','c_i1').replace('beta_x','beta_i1 ')+' ', 
                     'ind'   : 10 ,
                     'static': False,
                     'face'  : 'i1'},            
          #Li's definition imax: 
         'L0_imax' : {'symb'   : Li_BC_imax[0].replace('c','c_imax').replace('beta_x','beta_imax ')+' ', 
                     'ind'   : 11 ,
                     'static': False,
                     'face'  : 'imax'},            
         'L1_imax' : {'symb'   : Li_BC_imax[1].replace('c','c_imax').replace('beta_x','beta_imax ')+' ', 
                     'ind'   : 12 ,
                     'static': False,
                     'face'  : 'imax'},
         'L2_imax' : {'symb'   : Li_BC_imax[2].replace('c','c_imax').replace('beta_x','beta_imax ')+' ', 
                     'ind'   : 13  ,
                     'static': False,
                     'face'  : 'imax'},            
         'L3_imax' : {'symb'   : Li_BC_imax[3].replace('c','c_imax').replace('beta_x','beta_imax ')+' ', 
                     'ind'   : 14 ,
                     'static': False,
                     'face'  : 'imax'},   

          # Mi's definition jmax: 
         'M0_j1' : {'symb'   : Mi_BC_j1[0].replace('c','c_j1').replace('beta_y','beta_j1 ')+' ', 
                     'ind'   : 7 ,
                     'static': False,
                     'face'  : 'j1'},            
         'M1_j1' : {'symb'   : Mi_BC_j1[1].replace('c','c_j1').replace('beta_y','beta_j1 ')+' ', 
                     'ind'   : 8 ,
                     'static': False,
                     'face'  : 'j1'},
         'M2_j1' : {'symb'   : Mi_BC_j1[2].replace('c','c_j1').replace('beta_y','beta_j1 ')+' ', 
                     'ind'   : 9 ,
                     'static': False,
                     'face'  : 'j1'},            
         'M3_j1' : {'symb'   : Mi_BC_j1[3].replace('c','c_j1').replace('beta_y','beta_j1 ')+' ', 
                     'ind'   : 10 ,
                     'static': False,
                     'face'  : 'j1'},        
          # Mi's definition jmax: 
         'M0_jmax' : {'symb'   : Mi_BC_jmax[0].replace('c','c_jmax').replace('beta_y','beta_jmax ')+' ', 
                     'ind'   : 11 ,
                     'static': False,
                     'face'  : 'jmax'},            
         'M1_jmax' : {'symb'   : Mi_BC_jmax[1].replace('c','c_jmax').replace('beta_y','beta_jmax ')+' ', 
                     'ind'   : 12 ,
                     'static': False,
                     'face'  : 'jmax'},
         'M2_jmax' : {'symb'   : Mi_BC_jmax[2].replace('c','c_jmax').replace('beta_y','beta_jmax ')+' ', 
                     'ind'   : 13 ,
                     'static': False,
                     'face'  : 'jmax'},            
         'M3_jmax' : {'symb'   : Mi_BC_jmax[3].replace('c','c_jmax').replace('beta_y','beta_jmax ')+' ', 
                     'ind'   : 14 ,
                     'static': False,
                     'face'  : 'jmax'}


          }           
