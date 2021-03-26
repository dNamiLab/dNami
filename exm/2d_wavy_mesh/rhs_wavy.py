# =============================================================================
# 2D euler equations ig/vdw
# =============================================================================

# number of dimensions
dim = 2

#Strong or weak conservative form
iStrongConserv = False # False 

#Critical compressibility
Zc = ''
invZc = ''

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
                         'static': True},
         'eta' : {'symb': "eta", 
                         'ind':2 ,
                         'static': True},
         'dksidx' : {'symb': ' [ksi]_1x ', 
                         'ind':3 ,
                         'static': True},
         'dksidy' : {'symb': ' [ksi]_1y ', 
                         'ind':4 , 
                         'static': True},
         'detadx' : {'symb': ' [eta]_1x ', 
                         'ind':5 , 
                         'static': True},
         'detady' : {'symb': ' [eta]_1y ', 
                         'ind':6 ,
                         'static': True},
         'Jm1'    : {'symb': '1.0_wp  / ( dksidx * detady - detadx * dksidy  )', 
                         'ind':7,
                         'static': True},
         'U'      : {'symb': '(  detady * u - detadx * v ) ', 
                         'ind':8 ,
                         'static': False},
         'V'      : {'symb': '( -dksidy * u + dksidx * v )  ', 
                         'ind':9 ,
                         'static': False},
         'omg'    : {'symb': ' Jm1 * ( detady * [v]_1x - detadx * [v]_1y + dksidy * [u]_1x - dksidx * [u]_1y  ) ', 
                         'ind':10,
                         'static': False},
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

if iStrongConserv:

    Vel = {
            'U' : ' U ', 
            'V' : ' V ', 
            }

    divFx = {
        'rho' : ' Jm1  * ( [ rho * '+Vel['U']+'                 ]_1x )', 
        'u'   : ' Jm1  * ( [ rho * u * '+Vel['U']+' + detady*p  ]_1x )', 
        'v'   : ' Jm1  * ( [ rho * v * '+Vel['U']+' - detadx*p  ]_1x )', 
        'et'  : ' Jm1  * ( [ (rho * et  + p )  *  '+Vel['U']+'  ]_1x )', 
        }

    divFy = {
        'rho' : ' Jm1  * ( [ rho * '+Vel['V']+'                 ]_1y )', 
        'u'   : ' Jm1  * ( [ rho * u * '+Vel['V']+' - dksidx*p  ]_1y )', 
        'v'   : ' Jm1  * ( [ rho * v * '+Vel['V']+' + dksidy*p  ]_1y )', 
        'et'  : ' Jm1  * ( [ (rho * et  + p )  *  '+Vel['V']+'  ]_1y )', 
        }

else:

    Fx    = {  'rho' : ' (rho*u                ) ', 
               'u'   : ' (rho*u*u + '+Zc+'p    ) ', 
               'v'   : ' (rho*u*v              ) ', 
               'et'  : ' (u *(rho*et + '+Zc+'p)) ' }
     
    Fy    = {  'rho' : ' (rho*v                ) ', 
               'u'   : ' (rho*v*u              ) ', 
               'v'   : ' (rho*v*v + '+Zc+'p    ) ', 
               'et'  : ' (v *(rho*et + '+Zc+'p)) ' }

    divFx = {  'rho' : ' Jm1 * (   detady *  ['+Fx['rho']+']_1x   - dksidy * [ '+Fy['rho']+' ]_1x ) ', 
               'u'   : ' Jm1 * (   detady *  ['+Fx['u']+'  ]_1x   - dksidy * [ '+Fy['u']+'   ]_1x ) ', 
               'v'   : ' Jm1 * (   detady *  ['+Fx['v']+'  ]_1x   - dksidy * [ '+Fy['v']+'   ]_1x ) ', 
               'et'  : ' Jm1 * (   detady *  ['+Fx['et']+' ]_1x   - dksidy * [ '+Fy['et']+'  ]_1x ) ' }
     
    divFy = {  'rho' : ' Jm1 * ( - detadx * [ '+Fx['rho']+' ]_1y  + dksidx * [ '+Fy['rho']+' ]_1y ) ', 
               'u'   : ' Jm1 * ( - detadx * [ '+Fx['u']+'   ]_1y  + dksidx * [ '+Fy['u']+'   ]_1y ) ', 
               'v'   : ' Jm1 * ( - detadx * [ '+Fx['v']+'   ]_1y  + dksidx * [ '+Fy['v']+'   ]_1y ) ', 
               'et'  : ' Jm1 * ( - detadx * [ '+Fx['et']+'  ]_1y  + dksidx * [ '+Fy['et']+'  ]_1y ) ' }

# -- Assemble both 

divF = {}

for key in divFx.keys():
    divF[key] = divFx[key] + ' + ' + divFy[key]

