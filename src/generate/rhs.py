# =============================================================================
# 2D euler equations ig/vdw
# =============================================================================

# number of dimensions
dim = 2

# switches
#iVDW = False 
iVDW = True 

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
        }

# unknowns to march in time ///////////////////////////////////////////////////
varname      = {'rho' : 1,
                'u'  : 2,
                'v'  : 3,
                'et'  : 4, 
                }

varsolved = ['rho','u','v','et']

# of which
consvar      = [2,3,4] # are conservative variables

# derived local variables /////////////////////////////////////////////////////

if iVDW:
    varloc = { 'e' : ' (et - 0.5_wp*(u*u +v*v ) ) ',
               'p' : ' (8.0_wp * rho / (3.0_wp-rho) * delta * ( (e) + 9.0_wp/8.0_wp  * rho ) - 3.0_wp * rho*rho) ',
               } 
else:
    varloc = { 'e' : ' (et - 0.5_wp*(u*u + v*v) ) ',
               'p' : 'delta*rho* ( e )',
               }

# -- Same stored var for both

varstored    = {
         'divV' : {'symb': '   [u]_1x  + [v]_1y ',  #COMPUTE NEW EVERYTIME STEP 
                         'ind':1,
                         'static': False },
         }

# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file
rhsname = {'rho' : 'd(rho)/dt',
           'u'   : 'd(rho u)/dt',
           'v'   : 'd(rho v)/dt',
           'et'   : 'd(rho et)/dt'}

# .. name tags to use for intermediate variables created by the constructor

vnamesrc_divF = {'rho'  : 'FluRx',
                  'u'    : 'FluXx',
                  'v'    : 'FluYx',
                  'et'   : 'FluEx'}

vnamesrc_visc  = {'rho'  : 'FluRv',
                  'u'    : 'FluXv',
                  'v'    : 'FluYv',
                  'et'   : 'FluEv'}
  
# RHS terms ///////////////////////////////////////////////////////////////////

divFx = {  'rho' : ' [rho*u                ]_1x ', 
           'u'   : ' [rho*u*u + '+Zc+'p    ]_1x ', 
           'v'   : ' [rho*u*v              ]_1x ', 
           'et'  : ' [u *(rho*et + '+Zc+'p)]_1x ' }
 
divFy = {  'rho' : ' [rho*v                ]_1y ', 
           'u'   : ' [rho*v*u              ]_1y ', 
           'v'   : ' [rho*v*v + '+Zc+'p    ]_1y ', 
           'et'  : ' [v *(rho*et + '+Zc+'p)]_1y ' }

divF = { 
        }

for key in divFx.keys():
    divF[key] = divFx[key] + ' + ' + divFy[key]

viscx = { 
        'rho' : ' 0.0_wp  ', 
        'u'  : ' - mub * ( [ divV  ]_1x  ) ', 
        'v'  : ' 0.0_wp  ', 
        'et'  : ' - mub * ( [ u * divV ]_1x  ) ', 
        }

viscy = { 
        'rho' : ' 0.0_wp  ', 
        'u'  : ' 0.0_wp  ', 
        'v'  : ' - mub * ( [ divV  ]_1y ) ', 
        'et'  : ' - mub * ( [ v * divV ]_1y ) ', 
        }

visc = { 
        'rho' : ' 0.0_wp  ', 
        'u'  : ' - mub * ( [ divV  ]_1x  ) ', 
        'v'  : ' - mub * ( [ divV  ]_1y ) ', 
        'et'  : ' - mub * ( [ u * divV ]_1x + [ v * divV ]_1y ) ', 
        }

