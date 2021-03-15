# =============================================================================
# 3D Navier-Stokes equations : 3D counter-flow simulation in skew-symmetric form with a uniform grid
# =============================================================================

# number of dimensions
dim = 3
iVDW = False 

#Critical compressibility
if iVDW:
    Zc = '0.375_wp*'
    invZc = '1.0_wp/0.375_wp*'
else:
    Zc = ''
    invZc = ''

# coefficients ////////////////////////////////////////////////////////////////

coefficients = {
        'ReInv' : 1,
        'M0_sq' : 2,
        'delta' : 3,
        'kappa' : 4,
        'mu'  : 5,
        'mub' : 6,
        'gamma' : 7,
        'Tw'   : 8,
        }

# unknowns to march in time ///////////////////////////////////////////////////
varname      = {
        'rho' : 1,
        'u'   : 2,
        'v'   : 3,
        'w'   : 4,
        'et'  : 5, 
        }

varsolved = ['rho','u','v','w','et']

# of which
consvar      = [2,3,4,5] # are conservative variables
                       # i.e. advancing  
                       #    | rho 
                       # d  | rho*u 
                       # dt | rho*v
                       # dt | rho*w
                       #    | rho*et

# derived local variables /////////////////////////////////////////////////////
if iVDW:
    varloc = { 'e' : ' (et - 0.5_wp*(u*u + v*v + w*w) ) ',
               'p' : ' (8.0_wp * rho / (3.0_wp-rho) * delta * ( (e) + 9.0_wp/8.0_wp  * rho ) - 3.0_wp * rho*rho) ',
               } 
else:
    varloc = { 'e' : ' (et - 0.5_wp*(u*u + v*v + w*w ) ) ',
               'p' : 'delta*rho* ( e )',
               'T' : 'delta*M0_sq*e*gamma'
              }

varstored    = {
         'divV' : {'symb': '   [u]_1x  + [v]_1y + [w]_1z ',  #Computed every timestep
                         'ind':1,
                         'static': False },
        'c' : {'symb' : ' ( gamma*p / rho )**0.5 ',
              'ind':2, 'static' : False},
        'M' : {'symb' : ' ( u*u + v*v + w*w )**0.5 /  c ',
              'ind':3, 'static' : False},
        'phi' : {'symb' : ' ( phi ) ', 'ind':4, 'static': True},
        }
last = 4
# Statistics
variables = ['rho', 'u', 'v', 'w', 'M', 'T', 'et']
avgs = ['%s_mean' % var for var in variables]
rhs = [' %s' % var1 + ' + %s ' % var2 for (var1,var2) in zip(avgs, variables)]
stats = [{'symb' : ' %s ' % rh, 'ind':n+1+last, 'static': False} for n, rh in enumerate(rhs)]
stats = dict(zip(avgs, stats))
varstored.update(stats)

# names to give to the constructor ////////////////////////////////////////////

# .. for comments in the Fortran file

rhsname = {'rho' : 'd(rho)/dt',
           'u'   : 'd(rho u)/dt',
           'v'   : 'd(rho v)/dt',
           'w'   : 'd(rho w)/dt',
           'et'   : 'd(rho et)/dt'}

# .. name tags to use for intermediate variables created by the constructor

vname_Euler = {'rho'   : 'Euler_r',
                'u'    : 'Euler_u',
                'v'    : 'Euler_v',
                'w'    : 'Euler_w',
                'et'   : 'Euler_et'}

vname_Diffusive = {'rho'  : 'Diff_r',
                'u'    : 'Diff_u',
                'v'    : 'Diff_v',
                'w'    : 'Diff_w',
                'et'   : 'Diff_et'}
  
# RHS terms ///////////////////////////////////////////////////////////////////

# -- Inviscid contribution:
dilatation = '[u]_1x + [v]_1y + [w]_1z'

Euler= {
       'rho' : ' 0.5*( rho*( '+ dilatation +' ) + u*[rho]_1x + v*[rho]_1y + w*[rho]_1z + [ rho*u ]_1x + [ rho*v ]_1y + [ rho*w ]_1z ) ', 
       'u'   : ' 0.5*( rho*u*( '+ dilatation +' ) + u*[ rho*u ]_1x + v*[ rho*u ]_1y + w*[ rho*u ]_1z + 2.0*[p]_1x + [ rho*u*u ]_1x + [ rho*u*v ]_1y + [ rho*u*w ]_1z ) ', 
       'v'   : ' 0.5*( rho*v*( '+ dilatation +' ) + u*[ rho*v ]_1x + v*[ rho*v ]_1y + w*[ rho*v ]_1z + 2.0*[p]_1y + [ rho*v*u ]_1x + [ rho*v*v ]_1y + [ rho*v*w ]_1z ) ', 
       'w'   : ' 0.5*( rho*w*( '+ dilatation +' ) + u*[ rho*w ]_1x + v*[ rho*w ]_1y + w*[ rho*w ]_1z + 2.0*[p]_1z + [ rho*w*u ]_1x + [ rho*w*v ]_1y + [ rho*w*w ]_1z ) ',
       'et'  : ' 0.5*( rho*et*( '+ dilatation +' ) + u*[ rho*et ]_1x + v*[ rho*et ]_1y + w*[ rho*et ]_1z + 2.0*( [ p*u ]_1x + [ p*v ]_1y + [ p*w ]_1z ) + [ rho*et*u ]_1x + [ rho*et*v ]_1y + [ rho*et*w ]_1z ) ' 
        }

# Add body forcing term in the streamwise direction
Euler['u'] += ' - phi '
Euler['et'] += ' - u*phi '

Diffusive = {
          'u' : ' ',
          'v' : ' ',
          'w' : ' ',
          'et' : ' ',
}

Shear_visc = {
            'u' : ' ( -mu*ReInv / 3.0 )*( 4.0*[u]_2xx + 3.0*[u]_2yy + 3.0*[u]_2zz + 3.0*[ {v}_1y ]_1x - 2.0*[ {v}_1x ]_1y + 3.0*[ {w}_1z ]_1x - 2.0*[ {w}_1x ]_1z ) ',
            'v' : ' ( -mu*ReInv / 3.0 )*( 3.0*[v]_2xx + 4.0*[v]_2yy + 3.0*[v]_2zz + 3.0*[ {u}_1y ]_1x - 2.0*[ {u}_1x ]_1y + 3.0*[ {w}_1z ]_1y - 2.0*[ {w}_1y ]_1z ) ',
            'w' : ' ( -mu*ReInv / 3.0 )*( 3.0*[w]_2xx + 3.0*[w]_2yy + 4.0*[w]_2zz + 3.0*[ {u}_1z ]_1x - 2.0*[ {u}_1x ]_1z + 3.0*[ {v}_1y ]_1z - 2.0*[ {v}_1z ]_1y ) ',
            'et': ' ( -mu*ReInv / 3.0 )*( u*( 4.0*[u]_2xx + 3.0*[u]_2yy + 3.0*[u]_2zz + 3.0*[ {v}_1y ]_1x - 2.0*[ {v}_1x ]_1y + 3.0*[ {w}_1z ]_1x - 2.0*[ {w}_1x ]_1z ) +\
                                          v*( 3.0*[v]_2xx + 4.0*[v]_2yy + 3.0*[v]_2zz + 3.0*[ {u}_1x ]_1y - 2.0*[ {u}_1y ]_1x + 3.0*[ {w}_1z ]_1y - 2.0*[ {w}_1y ]_1z ) +\
                                          w*( 3.0*[w]_2xx + 3.0*[w]_2yy + 4.0*[w]_2zz + 3.0*[ {u}_1x ]_1z - 2.0*[ {u}_1z ]_1x + 3.0*[ {v}_1y ]_1z - 2.0*[ {v}_1z ]_1y ) +\
                                          4.0*( [u]_1x * [u]_1x - [u]_1x * [v]_1y - [u]_1x * [w]_1z + [v]_1y * [v]_1y - [v]_1y * [w]_1z + [w]_1z * [w]_1z ) + \
                                          6.0*( [u]_1y * [v]_1x + [u]_1z * [w]_1x + [v]_1z * [w]_1y ) + \
                                          3.0*( [u]_1y * [u]_1y + [u]_1z * [u]_1z + [v]_1x * [v]_1x + [v]_1z * [v]_1z + [w]_1x * [w]_1x + [w]_1y * [w]_1y ) ) '
            }

Bulk_visc = {
            'u' : ' - mub*ReInv * ( [divV]_1x  ) ',
            'v' : ' - mub*ReInv * ( [divV]_1y  ) ',
            'w' : ' - mub*ReInv * ( [divV]_1z  ) ',
            'et': ' - mub*ReInv * ( u*[divV]_1x +  v*[divV]_1y + w*[divV]_1z + divV*divV ) '
}

# Combine the terms
Diffusive['u'] = Diffusive['u'] + Shear_visc['u'] + Bulk_visc['u']
Diffusive['v'] = Diffusive['v'] + Shear_visc['v'] + Bulk_visc['v']
Diffusive['w'] = Diffusive['w'] + Shear_visc['w'] + Bulk_visc['w']
Diffusive['et'] = Diffusive['et'] + Shear_visc['et'] + Bulk_visc['et']
# Add heat flux dq_i/dx_i
Diffusive['et'] = Diffusive['et'] + '- ( mu*kappa )*( [T]_2xx + [T]_2yy + [T]_2zz )'

# Boundary conditions for Couette flow (3D)
locname_bc  = {'rho' : 'bc_rho',
               'u'  : 'bc_u',
               'v'  : 'bc_v',
               'w'  : 'bc_w',
               'et'  : 'bc_et  '}  

## Bottom wall j=1
# Set u_i = 0, constant temperature no-slip wall
q_j1 = {
        'u' : ' 0.0 ', 
        'v' : ' 0.0 ',
        'w' : ' 0.0 ',
        'et' : ' Tw/(M0_sq * delta * gamma)', 
        }

rhs_j1 = {'rho' : ' [ rho*u ]_1x + [ rho*v ]_1y + [ rho*w ]_1z ' }

## Top wall j=jmax, u_i = 0, constant temperature no-slip wall
q_jmax = {
        'u' : ' 0.0 ', 
        'v' : ' 0.0 ',
        'w' : ' 0.0 ',
        'et' : ' Tw/(M0_sq * delta * gamma) ', 
        }

rhs_jmax = {'rho' : ' [ rho*u ]_1x + [ rho*v ]_1y + [ rho*w ]_1z ' }
