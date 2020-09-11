# =============================================================================
# 1D euler equations in cartesian coordiantes ig/vdw
#                                                          -     sw 01-SEP-20
# =============================================================================

# number of dimensions
dim = 1

# coefficients ////////////////////////////////////////////////////////////////

coefficients = {
        'delta' : 1,
        }

# unknowns to march in time ///////////////////////////////////////////////////
varname      = {'rho' : 1,
                'u'   : 2,
                'et'  : 3, 
                }

varsolved = ['rho','u','et']

# of which
consvar      = [2,3] # are conservative variables

# derived local variables /////////////////////////////////////////////////////

varloc = { 'e' : ' (et - 0.5_wp*u*u) ',
           'p' : 'delta*rho* ( e )',
           }

# -- Same stored var for both

varbc = {
        }

varstored    = {
         }

# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file
rhsname = {'rho'  : 'd(rho)/dt',
           'u'    : 'd(rho u)/dt',
           'et'   : 'd(rho et)/dt'}

# .. name tags to use for intermediate variables created by the constructor

vnamesrc_divF = {'rho'  : 'FluRx',
                  'u'    : 'FluXx',
                  'et'   : 'FluEx'}

# RHS terms ///////////////////////////////////////////////////////////////////

F    = {  
        'rho' : ' (rho*u         ) ', 
        'u'   : ' (rho*u*u + p   ) ', 
        'et'  : ' (u*(rho*et + p)) '
        }
 
divF = { 
        'rho' : ' [ '+F['rho']+'  ]_1x  ', 
        'u'   : ' [ '+F['u'] +'   ]_1x  ', 
        'et'  : ' [ '+F['et'] +'  ]_1x  '
        }

