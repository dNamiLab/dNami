# =============================================================================
# 2D Euler with bulk viscosity // **vdW** EoS
# =============================================================================

# number of dimensions
dim = 2

# coefficients
a = '3.0_wp'
b = '0.3333333333333333_wp'
coefficients = {'visc'       : 1, 
                'kappa'      : 2,  
                'gamma_m1'   : 3,
                'CvInv'      : 4}

# unknowns to march in time
varname      = {'rho' : 1,
                  'u' : 2, 
                  'v' : 3, 
                  'et': 4}
# of which
consvar      = [2,3,4] # are conservative variables

# local variables
varloc       = {'p': 'rho*gamma_m1/(1.0_wp-'+b+'*rho)*(e +'+a+'*rho)-'+a+'*rho*rho',
                'e': 'et-0.5_wp*(u*u+v*v)',
                'T': 'CvInv*(e+'+a+'*rho)'}

# for comments in the src
rhsname = {'rho' : 'd(rho)/dt',  
             'u' : 'd(rho u)/dt',
             'v' : 'd(rho v)/dt', 
             'et': 'd(rho et)/dt'} 

# name tags to use for variable constructor
vnamesrc_div = {'rho': 'FluR',
                'u'  : 'FluX',
                'v'  : 'FluY',
                'et' : 'FluE'}

vnamesrc_dif = {'u'  : 'difU',
                'v'  : 'difV',
                'et' : 'difE'}            

# specify the PDEs
# IMPORTANT SYNTAX RULE: put a space after derivative symbol []_1* 

# Euler terms in -div(F)

divF = {'rho' : ' [rho*u]_1x + [rho*v]_1y ',
        'u'   : ' [rho*u*u + p]_1x  + [rho*u*v]_1y ',
        'v'   : ' [rho*u*v]_1x  + [rho*v*v + p]_1y ',
        'et'  : ' [u*(rho*et + p)]_1x + [v*(rho*et + p)]_1y '}

# Bulk viscosity & thermal conductivity terms

# diff = {'u'    : '- visc*( [u]_2xx + [u]_2yy + [ {v}_1x ]_1y )',
#         'et'   : '- visc*( [u]_1x * [u]_1x  +  [v]_1y * [v]_1y  +  2.0_wp * [u]_1x * [v]_1y + u*( [u]_2xx + [ {v}_1x + rho  ]_1y )  +  v*( [v]_2yy + [ {u}_1x ]_1y ) )'}

# diff = {'u'    : '- visc*( [ {v}_1x ]_1y )',
#         'et'   : '- visc*( [u]_1y + [ {v}_1x ]_1y )'}

# diff = {'et'   : '- visc*(  u*( [u]_2xx + [ {v}_1x ]_1y )   +  [v]_1y * [v]_1y  )'}

# diff = {'et'   : '  [ {v}_1x ]_1y + [ {u}_1x ]_1y + [ u ]_2xx '}

# diff = {'et'   : ' [ {u}_1x ]_1y '}

# diff = {'u'    : '- visc*( [u]_2xx + [u]_2yy )',
#         'v'    : '- visc*( [v]_2xx + [v]_2yy )',
#         'et'   : '- visc*( u*( [u]_2xx + [ {v}_1x ]_1y )  +  v*( [v]_2yy + [ {u}_1y ]_1x ) + [u]_1x * [u]_1x  +  [v]_1y * [v]_1y  +  2.0_wp * [u]_1x * [v]_1y )'}

diff = {'u'    : '- visc*( [u]_2xx + [u]_2yy )',
        'v'    : '- visc*( [v]_2xx + [v]_2yy )',
        'et'   : '- visc*( u*( [u]_2xx + [ {v}_1x ]_1y )  +  v*( [v]_2yy + [ {u}_1y ]_1x + [u]_1x * [u]_1x  +  [v]_1y * [v]_1y  +  2.0_wp * [u]_1x * [v]_1y ) )'}

# diff['et'] = diff['et'] + '- kappa* ( [T]_2xx + [T]_2yy )'