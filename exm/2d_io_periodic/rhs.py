# =============================================================================
# 3D navier stokes equations  
# =============================================================================

# number of dimensions
dim = 2 

# coefficients ////////////////////////////////////////////////////////////////
coefficients = {
                }

# unknowns to march in time ///////////////////////////////////////////////////
varname      = {'a' : 1,
                'b' : 2,}
                  
varsolved = ['a', 'b']

# of which
consvar      = [] # are conservative variables         

# derived local variables /////////////////////////////////////////////////////

varloc       = {
                }


# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file
rhsname = { 'a': 'd(a)/dt',  
            'b': 'd(b)/dt'
            }
       

# .. name tags to use for intermediate variables created by the constructor
locname_f =  {'a' : 'f_a',
              'b' : 'f_b  '}

locname_bc  = {'a' : 'bc_a',
               'b' : 'bc_b',
               }                

# RHS terms ///////////////////////////////////////////////////////////////////

F = {'a' : ' 0.0_wp  ',
     'b' : ' 0.0_wp  '}

