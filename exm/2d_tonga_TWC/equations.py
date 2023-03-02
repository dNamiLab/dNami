# -----------------------------------------------------------------------------
# 2D Two-way coupled long-wave equations on a sphere for the Tonga explosion case  
# 
#                                                          -sw 15-JUN-20
# -----------------------------------------------------------------------------
from functions import *

# NB: theta is the polar angle, phi is the azimuthal angle
#     x -> theta, y -> phi

# number of dimensions
dim = 2

# air and water or just water
iAtmos = True 
 
# coefficients ////////////////////////////////////////////////////////////////

coefficients = {}

coefficients['ft'] = 1 
coefficients['one_over_Rad'] = 2 
coefficients['rhow']         = 3 
coefficients['one_over_Rew'] = 4 

if iAtmos:
    coefficients['gamma'] = 5 
    coefficients['eta2']  = 6 

# unknowns to march in time ///////////////////////////////////////////////////

varname      = {}

# -- Add variables for water layer  

varname['eta1'] = 1
varname['vt'  ] = 2
varname['vp'  ] = 3
                
varsolved = ['eta1','vt','vp']

# -- Add variables for atmosphere

if iAtmos:
    varname['rhoa'] = 4
    varname['ut']   = 5
    varname['up']   = 6
    varname['pa']   = 7

    varsolved += ['rhoa', 'ut', 'up', 'pa']

# of which
consvar      = []  

# derived local variables /////////////////////////////////////////////////////

varloc = {} 

varloc['hw']    = '(eta1 - eta0)'                 #water depth 

if iAtmos:
    varloc['eta1a'] = '(eta1 + max(eta0,0.0_wp))' #eta1 gradient taking into account ground
    varloc['ha']    = '(eta2 - eta1a)'            #atmospheric thickness    
    varloc['ca2']   = 'gamma * pa/ rhoa'           #local speed of sound

# -- Same stored var for both

varstored    = {
         'eta0' : {'symb': ' eta0 ',       # elevation map 
                         'ind':1,
                         'static': True},
         'sint' : {'symb': '  sint  ',     # stored at the start 
                         'ind':2,
                         'static': True},
         'cost' : {'symb': '  cost ',      # stored at the start 
                         'ind':3,
                         'static': True},
         'sintinv' : {'symb': ' sintinv ', # stored at the start 
                         'ind':4,
                         'static': True},
         'hmax' : {'symb': ' hmax  ',      # historical max height difference 
                         'ind':5,
                         'static': True},
         'tmax' : {'symb': ' tmax  ',      # historical time 
                         'ind':6,
                         'static': True},
         }

varsidx = len(varstored.keys()) + 1
for j,var in enumerate(varsolved):
    varstored['S' + var] = { 'symb' : 'S'+ var , 'ind' : varsidx  + j , 'static' : True}

# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file
rhsname = {}

for var in list(varname.keys()):
    rhsname[var]  = 'd('+var+')/dt'

# .. name tags to use for intermediate variables created by the constructor

vnamesrc_Adv  = {}

for var in varsolved:
    vnamesrc_Adv[var] = 'Flux_AdvF_' + var

vnamesrc_S = {}

for var in varsolved:
    vnamesrc_S[var] = 'Flux_S_' + var

vnamesrc_C = {}

for var in varsolved:
    vnamesrc_C[var] = 'Flux_C_' + var

vnamesrc_Fss = {}

for var in varsolved:
    vnamesrc_Fss[var] = 'Flux_Fss_' + var

# RHS terms ///////////////////////////////////////////////////////////////////

# -- Define some functions to help with the syntax

# --------- Start filling ... 

Adv  = {}
S    = {}
C    = {}
Fss  = {}

# ----- Filling Adv 

# -- water
Adv['eta1'] = 'vt *' + grad1d('hw',0)  + ' + vp *'  + grad1d('hw',1) 
Adv['vt']   = 'vt *' + grad1d('vt',0)  + ' + vp *'  + grad1d('vt',1) 
Adv['vp']   = 'vt *' + grad1d('vp',0)  + ' + vp *'  + grad1d('vp',1) 

if iAtmos:
    # -- air
    Adv['rhoa']  = 'ut *' + grad1d('rhoa',0)  + ' + up *'  + grad1d('rhoa',1) 
    Adv['ut']    = 'ut *' + grad1d('ut'  ,0)  + ' + up *'  + grad1d('ut'  ,1) 
    Adv['up']    = 'ut *' + grad1d('up'  ,0)  + ' + up *'  + grad1d('up'  ,1) 
    Adv['pa']    = 'ut *' + grad1d('pa'  ,0)  + ' + up *'  + grad1d('pa'  ,1) 


# ----- Filling source terms (Written on LHS!) 

# -- water
S['eta1']  = 'hw * ( ' + div1d('vt', 'vp')  + ')' 
S['vt']    = grad1d('eta1', 0)  
S['vp']    = grad1d('eta1', 1)  

if iAtmos:
    #  -- psi
    psi = '1.0_wp/ha * (' + div1d('(ha*ut + hw*vt)', '(ha*up + hw*vp)') + ')'
    ## -- Coupling term 
    S['vt']   += '+ 1.0_wp/rhow *' + grad1d('rhoa*ha',0)  
    S['vp']   += '+ 1.0_wp/rhow *' + grad1d('rhoa*ha',1)  
    ## -- air 
    S['rhoa'] = 'rhoa *' + psi 
    S['ut']   = grad1d('eta1a', 0) + ' + 1.0_wp/(rhoa*ha) *' + grad1d('pa*ha',0)  
    S['up']   = grad1d('eta1a', 1) + ' + 1.0_wp/(rhoa*ha) *' + grad1d('pa*ha',1) 
    S['pa']   = 'rhoa * ca2 * ' + psi  

#
for var in varsolved:
    C[var] = ' - ( ft * S' + var + ')'

# -- Create stored forcing variables
varsidx = len(varstored.keys()) + 1
for j,var in enumerate(varsolved):

    if var in Adv.keys():
        a = Adv[var] 
    else:
        a = ' 0.0_wp' 

    if var in S.keys():
        b = S[var]
    else:
        b = ' 0.0_wp' 

    if var in C.keys():
        c = C[var]
    else:
        c = ' 0.0_wp' 

    varstored['rhs_' + var] = { 'symb' : ' - (' + a + ' + ' + b + '+' + c + ')' , 'ind' : varsidx  + j , 'static' :False}

# -- Add forcing to components 

for var in varsolved:
    Fss[var] = ' rhs_' + var

