# =============================================================================
#
# 1D kdv equations 
#
# =============================================================================

# number of dimensions
dim = 1

# coefficients ////////////////////////////////////////////////////////////////

coefficients = {
        'epsilon' : 1,
        'mu' : 2,
        }

# unknowns to march in time ///////////////////////////////////////////////////
varname      = {'rho' : 1,
                }

varsolved = ['rho']

# of which
consvar      = [] # are conservative variables

# derived local variables /////////////////////////////////////////////////////

varloc = { 
           }

# -- Same stored var for both

varbc = {
        }

varstored    = {
         'rho_2x'   : {'symb': ' [rho]_2xx ',                   # FILL IN COMPUTE 
                         'ind':1,
                         'static': False },
         }

# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file
rhsname = {'rho'  : 'd(rho)/dt',
           }

# .. name tags to use for intermediate variables created by the constructor

vnamesrc_divF = {'rho'  : 'FluRx',
                  }

# RHS terms ///////////////////////////////////////////////////////////////////

 
divF = { 
        'rho' : ' epsilon * rho * [ rho ]_1x + mu * [ rho_2x ]_1x ', 
        }

