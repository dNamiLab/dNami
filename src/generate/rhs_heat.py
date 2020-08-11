# =============================================================================
# 2D curvilinear grid generation 
# =============================================================================

# number of dimensions
dim = 2

# coefficients ////////////////////////////////////////////////////////////////

coefficients = {'alpha':1}

# unknowns to march in time ///////////////////////////////////////////////////
varname      = {'T' : 1}

varsolved = ['T']

# of which
consvar = [] # are conservative variables

# derived local variables /////////////////////////////////////////////////////

varbc = {'T_i1'  : {'symb'  : '', 
                    'ind'   :1 ,
                    'static': True,
                    'face'  :'i'},
         'T_imax': {'symb'  : '', 
                    'ind'   :2 ,
                    'static': True,
                    'face'  :'i'},
         'T_j1'  : {'symb'  : '', 
                    'ind'   :1 ,
                    'static': True,
                    'face'  :'j'},
         'T_jmax': {'symb'  : '', 
                    'ind'   :2 ,
                    'static': True,
                    'face'  :'j'}}


# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file
rhsname = {'T' : 'd(T)/dt'}

# .. name tags to use for intermediate variables created by the constructor

vnamesrc = {'T'  : 'F_T'}

# RHS terms ///////////////////////////////////////////////////////////////////

Lap = {  'T' : ' - alpha * ( [T]_2xx + [T]_2yy ) ' }

#Lap = {  'T' : ' - alpha * ( [ {T}_1x ]_1x + [ {T}_1y ]_1y ) ' }

# BOUNDARY CONDITIONS ////////////////////////////////////////////////////////

locname_bc  = {'T' : 'bc_T',
               }                

# -- Set each of the values 

phybc_i1   = { 'T': 'T_i1'}
phybc_imax = { 'T': 'T_imax'}
phybc_j1   = { 'T': 'T_j1'}
phybc_jmax = { 'T': 'T_jmax'}

