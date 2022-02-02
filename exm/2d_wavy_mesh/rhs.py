# =============================================================================
# 2D euler equations ig/vdw
# =============================================================================

# number of dimensions
dim = 2

# coefficients ////////////////////////////////////////////////////////////////

coefficients = {
        'delta' : 1
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

varloc = { 'e' : '(et - 0.5_wp*(u*u + v*v) )',
           'p' : 'delta*rho* ( e )',
           'c' : ' ( ( 1.0_wp + delta ) * p / rho  )**0.5_wp ',
           }

# -- Same stored var for both

varstored    = {
         'ksi' : {'symb': "ksi", 
                         'ind':1 ,
                         'static': True}, # metric terms
         'eta' : {'symb': "eta", 
                         'ind':2 ,
                         'static': True}, # metric terms
         'dksidx' : {'symb': ' [ksi]_1x ', 
                         'ind':3 ,
                         'static': True}, # metric derivative
         'dksidy' : {'symb': ' [ksi]_1y ', 
                         'ind':4 , 
                         'static': True}, # metric derivative
         'detadx' : {'symb': ' [eta]_1x ', 
                         'ind':5 , 
                         'static': True}, # metric derivative
         'detady' : {'symb': ' [eta]_1y ', 
                         'ind':6 ,
                         'static': True}, # metric derivative
         'Jm1'    : {'symb': '1.0_wp  / ( dksidx * detady - detadx * dksidy  )', 
                         'ind':7,
                         'static': True}, # coordinate transformation jacobian 
         'omg'    : {'symb': ' Jm1 * ( detady * [v]_1x - detadx * [v]_1y + dksidy * [u]_1x - dksidx * [u]_1y  ) ', 
                         'ind':8,
                         'static': False}, # vorticity in curvilinear coordinates
         }


# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file
rhsname = {'rho' : 'd(rho)/dt',
           'u'   : 'd(rho u)/dt',
           'v'   : 'd(rho v)/dt',
           'et'   : 'd(rho et)/dt'}

# .. name tags to use for intermediate variables created by the constructor

vnamesrc_divFx = {'rho'  : 'FluRx',
                  'u'    : 'FluXx',
                  'v'    : 'FluYx',
                  'et'   : 'FluEx'}
  
vnamesrc_divFy = {'rho'  : 'FluRy',
                  'u'    : 'FluXy',
                  'v'    : 'FluYy',
                  'et'   : 'FluEy'}

# RHS terms ///////////////////////////////////////////////////////////////////

Vel = {
        'U' : ' (detady * u - dksidy * v  ) ', 
        'V' : ' (-detadx *u + dksidx * v  ) ', 
        }

divFx = {
    'rho' : ' Jm1  * ( [ rho * '+Vel['U']+'                 ]_1x )', 
    'u'   : ' Jm1  * ( [ rho * u * '+Vel['U']+' + detady*p  ]_1x )', 
    'v'   : ' Jm1  * ( [ rho * v * '+Vel['U']+' - dksidy*p  ]_1x )', 
    'et'  : ' Jm1  * ( [ (rho * et  + p )  *  '+Vel['U']+'  ]_1x )', 
    }

divFy = {
    'rho' : ' Jm1  * ( [ rho * '+Vel['V']+'                 ]_1y )', 
    'u'   : ' Jm1  * ( [ rho * u * '+Vel['V']+' - detadx*p  ]_1y )', 
    'v'   : ' Jm1  * ( [ rho * v * '+Vel['V']+' + dksidx*p  ]_1y )', 
    'et'  : ' Jm1  * ( [ (rho * et  + p )  *  '+Vel['V']+'  ]_1y )', 
    }

# -- Assemble both directions 
divF = {}

for key in divFx.keys():
    divF[key] = divFx[key] + ' + ' + divFy[key]

