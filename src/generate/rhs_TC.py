# =============================================================================
# 2D euler equations ig/vdw
# =============================================================================

# number of dimensions
dim = 2

# switches
iVDW = False 
#iVDW = True 

iChar = False 
#iChar = False 

#Critical compressibility
if iVDW:
    Zc = '0.375_wp*'
    invZc = '1.0_wp/0.375_wp*'
else:
    Zc = ''
    invZc = ''

# coefficients ////////////////////////////////////////////////////////////////

coefficients = {
        'delta' : 1,
        'mub'   : 2,
        'mu'    : 3,
        'u0'    : 4,
        'u1'    : 5,
        'Tw'    : 6,
        }

# unknowns to march in time ///////////////////////////////////////////////////
varname      = {'rho' : 1,
                'ur'  : 2,
                'ut'  : 3,
                'et'  : 4, 
                }

varsolved = ['rho','ur','ut','et']

# of which
consvar      = [2,3,4] # are conservative variables

# derived local variables /////////////////////////////////////////////////////

if iVDW:
    varloc = { 'e' : ' (et - 0.5_wp*(ur*ur + ut*ut) ) ',
               'p' : ' (8.0_wp * rho / (3.0_wp-rho) * delta * ( (e) + 9.0_wp/8.0_wp  * rho ) - 3.0_wp * rho*rho) ',
               'T' : '( ( (p) + 3.0_wp * rho**2 ) /8.0_wp / rho * (3.0_wp-rho) ) ', 
               'c' : ' ((T)*(1.0_wp + delta)*( 1.0_wp / ( 1.0_wp - rho / 3.0_wp) )**2 - 2.0_wp*rho*9.0_wp/8.0_wp)**0.5_wp ',
               'dedp' : '( (3.0_wp - rho)/(8.0_wp*rho*delta)  )', #AT FIXED RHO
               'cp_o_av' : invZc+' rho * (c)*(c)*(dedp)',
               } 
else:
    varloc = { 'e' : ' (et - 0.5_wp*(ur*ur + ut*ut) ) ',
               'p' : 'delta*rho* ( e )',
               'c' : ' ( ( 1.0_wp + delta ) * p / rho  )**0.5_wp ',
               'dedp' : '1.0_wp/delta/rho', #AT FIXED RHO
               'cp_o_av' : invZc+' rho * (c)*(c)*(dedp)',
               'dTdp'    : '( 1.0_wp / rho )', 
               'dTdrho'    : '( -p/rho/rho  )', 
               }

# -- Same stored var for both


varstored    = {
         'r'   : {'symb': ' r ',                   # FILL IN COMPUTE 
                         'ind':1,
                         'static': True},
         'ksi' : {'symb': ' ksi ',                 # FILL IN COMPUTE 
                         'ind':2,
                         'static': True},
         'eta' : {'symb': ' eta  ',                # FILL IN COMPUTE
                         'ind':3,
                         'static': True},
         'divV' : {'symb': '  1./r * [r * ur]_1x + 1./r * [ut]_1y  ',  #COMPUTE NEW EVERYTIME STEP 
                         'ind':4,
                         'static': False },
         'durdr' : {'symb': ' [ur]_1x ',  #COMPUTE NEW EVERYTIME STEP 
                         'ind':5,
                         'static': False },
         'durdt' : {'symb': ' [ur]_1y ',  #COMPUTE NEW EVERYTIME STEP 
                         'ind':6,
                         'static': False },
         'dutdr' : {'symb': ' [ut]_1x ',  #COMPUTE NEW EVERYTIME STEP 
                         'ind':7,
                         'static': False },
         'dutdt' : {'symb': ' [ut]_1y ',  #COMPUTE NEW EVERYTIME STEP 
                         'ind':8,
                         'static': False },
         }

# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file
rhsname = {'rho' : 'd(rho)/dt',
           'ur'   : 'd(rho ur)/dt',
           'ut'   : 'd(rho ut)/dt',
           'et'   : 'd(rho et)/dt'}

# .. name tags to use for intermediate variables created by the constructor

vnamesrc_divF = {'rho'  : 'FluRx',
                  'ur'    : 'FluXx',
                  'ut'    : 'FluYx',
                  'et'   : 'FluEx'}

vnamesrc_visc  = {'rho'  : 'FluRv',
                  'ur'    : 'FluXv',
                  'ut'    : 'FluYv',
                  'et'   : 'FluEv'}
  
# RHS terms ///////////////////////////////////////////////////////////////////

Fr    = {  
        'rho' : ' (r*rho*ur             ) ', 
        'ur'  : ' (r*(rho*ur*ur + '+Zc+'p)) ', 
        'ut'  : ' (r*rho*ur*ut          ) ', 
        'et'  : ' (r*ur*(rho*et + '+Zc+'p)) '
        }
 
Ft    = {  
        'rho' : ' (rho*ut               ) ', 
        'ur'  : ' (rho*ur*ut            ) ', 
        'ut'  : ' (rho*ut*ut  + '+Zc+'p ) ', 
        'et'  : ' (ut*(rho*et + '+Zc+'p)) '
        }

divF = { 
        'rho' : ' 1.0_wp/r * ( [ '+Fr['rho']+'  ]_1x + [ '+Ft['rho']+' ]_1y ) + 0.0_wp                ', 
        'ur'  : ' 1.0_wp/r * ( [ '+Fr['ur'] +'  ]_1x + [ '+Ft['ur']+ ' ]_1y ) - ( ut*ut + '+Zc+'p )/r ', 
        'ut'  : ' 1.0_wp/r * ( [ '+Fr['ut'] +'  ]_1x + [ '+Ft['ut']+ ' ]_1y ) + ( ur*ut           )/r ', 
        'et'  : ' 1.0_wp/r * ( [ '+Fr['et'] +'  ]_1x + [ '+Ft['et']+ ' ]_1y ) +0.0_wp                 '
        }

Drr = ' ( durdr - 1.0_wp/3.0_wp * divV ) '
Drt = ' ( 1.0_wp / (2.0_wp * r) * ( durdt - ut  ) + 1.0_wp /2.0_wp * dutdr ) '
Dtt = ' ( 1.0_wp / r * ( dutdt + ur ) - 1.0_wp/3.0_wp * divV ) '

visc_shear= { 
        'rho' : ' 0.0_wp  ', 
        'ur'  : ' - 2.0_wp * mu *( 1.0_wp / r * [r*'+Drr+']_1x + 1.0_wp / r * ['+Drt+']_1y - '+Dtt+'/r               ) ', 
        'ut'  : ' - 2.0_wp * mu *( 1.0_wp / r * [r*'+Drt+']_1x + 1.0_wp / r * ['+Dtt+']_1y + '+Drt+'/r               ) ', 
        'et'  : ' - 2.0_wp * mu *( 1.0_wp / r * [r*(ur*'+Drr+'+ut*'+Drt+')]_1x + 1.0_wp / r * [ur*'+Drt+'+ut*'+Dtt+']_1y ) ', 
        }

visc_bulk = { 
        'rho' : ' 0.0_wp  ', 
        'ur'  : ' - mub * ( 1.0_wp /r * [ r * divV  ]_1x - divV/r  ) ', 
        'ut'  : ' - mub * ( 1.0_wp /r * [     divV  ]_1y ) ', 
        'et'  : ' - mub * ( 1.0_wp /r * [ r * ur * divV  ]_1x + 1.0_wp /r * [ ut * divV ]_1y ) ', 
        }

visc = {
        }

# -- Add both terms
for key in ['rho', 'ur', 'ut', 'et']:
    visc[key] = visc_shear[key] + ' + ' + visc_bulk[key]


# BOUNDARY CONDITIONS ////////////////////////////////////////////////////////

locname_bc  = {'rho' : 'bc_rho',
               'ur'  : 'bc_ur',
               'ut'  : 'bc_ut',
               'et'  : 'bc_et  '}                

# ====================================================== I / r ======================================================================

# ------------- i1

src_phybc_i1 = {
        'rho' :' rho * [ur]_1x + 1.0/r * ( rho * [ut]_1y + u0 * [rho]_1y ) ',
        }

q_i1 = {
        'ur' : '0.0_wp ', 
        'ut' : ' u0 ',
        'et' : ' Tw/delta + 0.5_wp * ( ur*ur + ut*ut  )', 
        }

# ------------- imax

src_phybc_imax = {
        'rho' :' rho * [ur]_1x + 1.0/r * ( rho * [ut]_1y + u1 * [rho]_1y ) ',
        }

q_imax = {
        'ur' : '0.0_wp ', 
        'ut' : ' u1 ',
        'et' : ' Tw/delta + 0.5_wp * ( ur*ur + ut*ut  )', 
        }
