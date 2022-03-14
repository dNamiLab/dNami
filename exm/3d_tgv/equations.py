# =============================================================================
# 3D navier stokes equations  
# =============================================================================

# number of dimensions
dim = 3 

# coefficients ////////////////////////////////////////////////////////////////
coefficients = {'visc'       : 1, 
                'kappa'      : 2,  
                'gamma_m1'   : 3,
                'CvInv'      : 4,
                'u_0'        : 5}

# unknowns to march in time ///////////////////////////////////////////////////
varname      = {'rho' : 1,
                  'u' : 2,  
                  'v' : 3, 
                  'w' : 4,
                  'et': 5}
                  
varsolved = ['rho','u','v','w','et']

# of which
consvar      = [2,3,4,5] # are conservative variables         

# derived local variables /////////////////////////////////////////////////////

# -- Two divergence forms for inside and outside first derivative
divops  = ' [u]_1x + [v]_1y + [w]_1z '
ddivops = ' {u}_1x + {v}_1y + {w}_1z '

varloc       = {'p': 'gamma_m1*rho*(e)',
                'e': '(et-0.5_wp*(u*u+v*v+w*w))',
                'T': 'CvInv*(e)'}


# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file
rhsname = {'rho' : 'd(rho)/dt'  
          ,  'u' : 'd(rho u)/dt',
             'v' : 'd(rho v)/dt', 
             'w' : 'd(rho w)/dt', 
             'et': 'd(rho et)/dt'}
       

# .. name tags to use for intermediate variables created by the constructor
locname_dif = {'rho': 'dif_rho',
               'u'  : 'dif_rhou',
               'v'  : 'dif_rhov',
               'w'  : 'dif_rhow',
               'et' : 'dif_et  '}

locname_conv = {'rho': 'conv_rho',
                'u'  : 'conv_rhou',
                'v'  : 'conv_rhov',
                'w'  : 'conv_rhow',
                'et' : 'conv_et  '}

locname_bc  = {'rho': 'bc_rho',
               'u'  : 'bc_u',
               'v'  : 'bc_v',
               'w'  : 'bc_w',
               'et' : 'bc_et  '}                

# RHS terms ///////////////////////////////////////////////////////////////////

# Euler 

Fx = {'rho' : 'rho*u         ',
      'u'   : 'rho*u*u  + p  ',
      'v'   : 'rho*u*v       ',
      'w'   : 'rho*u*w       ',
      'et'  : '(rho*et + p)*u '}

Fy = {'rho' : 'rho*v         ',
      'u'   : 'rho*v*u       ', 
      'v'   : 'rho*v*v  + p  ', 
      'w'   : 'rho*v*w       ', 
      'et'  : '(rho*et + p)*v '} 

Fz = {'rho' : 'rho*w         ',
      'u'   : 'rho*w*u       ',
      'v'   : 'rho*w*v       ',
      'w'   : 'rho*w*w  + p  ',
      'et'  : '(rho*et + p)*w '}

# -- Conservative formulation

Src_conv = {'rho' : '[ '+Fx['rho']+' ]_1x' + ' + ' + '[ '+Fy['rho']+' ]_1y' + ' + ' + '[ '+Fz['rho']+' ]_1z ',
            'u'   : '[ '+Fx['u']  +' ]_1x' + ' + ' + '[ '+Fy['u']  +' ]_1y' + ' + ' + '[ '+Fz['u']  +' ]_1z ',
            'v'   : '[ '+Fx['v']  +' ]_1x' + ' + ' + '[ '+Fy['v']  +' ]_1y' + ' + ' + '[ '+Fz['v']  +' ]_1z ',
            'w'   : '[ '+Fx['w']  +' ]_1x' + ' + ' + '[ '+Fy['w']  +' ]_1y' + ' + ' + '[ '+Fz['w']  +' ]_1z ',
            'et' : ' [ '+Fx['et'] +' ]_1x' + ' + ' + '[ '+Fy['et'] +' ]_1y' + ' + ' + '[ '+Fz['et'] +' ]_1z '}

# -- Skew symmetric formulation

Src_skew = {'rho' : '0.5_wp*( [rho*u]_1x + [rho*v]_1y + [rho*w]_1z '
                   '+   u*[rho]_1x + v*[rho]_1y + w*[rho]_1z '
                   '+ rho*( '+divops+' ) )',
        'u'   : '0.5_wp*( [rho*u*u]_1x + [rho*u*v]_1y + [rho*u*w]_1z '
                   '+   u*[rho*u]_1x + v*[rho*u]_1y + w*[rho*u]_1z '
                   '+rho*u*( '+divops+' ) ) + [p]_1x ',
        'v'   : '0.5_wp*( [rho*u*v]_1x + [rho*v*v]_1y + [rho*v*w]_1z '
                   '+   u*[rho*v]_1x + v*[rho*v]_1y + w*[rho*v]_1z '
                   '+rho*v*( '+divops+' ) ) + [p]_1y ',
        'w'   : '0.5_wp*( [rho*w*u]_1x + [rho*w*v]_1y + [rho*w*w]_1z '
                   '+   u*[rho*w]_1x + v*[rho*w]_1y + w*[rho*w]_1z '
                   '+rho*w*( '+divops+' ) ) + [p]_1z ',
        'et'  : '0.5_wp*( [rho*et*u]_1x + [rho*et*v]_1y + [rho*et*w]_1z '
                   '+   u*[rho*et]_1x + v*[rho*et]_1y + w*[rho*et]_1z '
                   '+rho*et*( '+divops+') ) + [p*u]_1x + [p*v]_1y + [p*w]_1z '}

#Navier-Stokes Diffusive terms only 

Fx = {'u'   : ' - visc *( 2.0_wp * {u}_1x - 2.0_wp/3.0_wp * ( '+ ddivops +'  ) )',
      'v'   : ' - visc *( {u}_1y + {v}_1x )', 
      'w'   : ' - visc *( {u}_1z + {w}_1x )',
      
      'et'  : ' - kappa*( {T}_1x ) '
              ' - u*(visc *( 2.0_wp * {u}_1x - 2.0_wp/3.0_wp * ( '+ ddivops +'  )))'
              ' - v*(visc *( {u}_1y + {v}_1x ))'
              ' - w*(visc *( {u}_1z + {w}_1x ))'}

Fy = {'u'   : ' - visc *( {u}_1y + {v}_1x )  ',
      'v'   : ' - visc *( 2.0_wp * {v}_1y - 2.0_wp/3.0_wp * ( '+ ddivops +'  ) )', 
      'w'   : ' - visc *( {v}_1z + {w}_1y )  ',
      'et'  : ' - kappa*( {T}_1y )'
              ' - u*(visc *( {u}_1y + {v}_1x ))'
              ' - v*(visc *( 2.0_wp * {v}_1y - 2.0_wp/3.0_wp * ( '+ ddivops +'  )))'
              ' - w*(visc *( {v}_1z + {w}_1y ))'}
       
Fz = {'u'   : ' - visc *( {u}_1z + {w}_1x )',
      'v'   : ' - visc *( {v}_1z + {w}_1y )', 
      'w'   : ' - visc *( 2.0_wp * {w}_1z - 2.0_wp/3.0_wp * ( '+ ddivops +'  ) )',
      'et'  : ' - kappa*( {T}_1z )'
              ' - u*(visc *( {u}_1z + {w}_1x ))'
              ' - v*(visc *( {v}_1z + {w}_1y ))'
              ' - w*(visc *( 2.0_wp * {w}_1z - 2.0_wp/3.0_wp * ( '+ ddivops +'  )))'}

# -- Divergence formulation

Src_dif  = {'u'   : '[ '+Fx['u']  +' ]_1x' + ' + ' + '[ '+Fy['u']  +' ]_1y' + ' + ' + '[ '+Fz['u']  +' ]_1z ',
            'v'   : '[ '+Fx['v']  +' ]_1x' + ' + ' + '[ '+Fy['v']  +' ]_1y' + ' + ' + '[ '+Fz['v']  +' ]_1z ',
            'w'   : '[ '+Fx['w']  +' ]_1x' + ' + ' + '[ '+Fy['w']  +' ]_1y' + ' + ' + '[ '+Fz['w']  +' ]_1z ',
            'et' : ' [ '+Fx['et'] +' ]_1x' + ' + ' + '[ '+Fy['et'] +' ]_1y' + ' + ' + '[ '+Fz['et'] +' ]_1z '}

# -- Laplacian formulation

Src_Laplace = {'u'   : '- visc* ( [u]_2xx + [u]_2yy + [u]_2zz ) - visc*0.3333333333333333_wp*( [ '+ ddivops+' ]_1x )  ',
     'v'   : '- visc* ( [v]_2xx + [v]_2yy + [v]_2zz ) - visc*0.3333333333333333_wp*( [ '+ ddivops+' ]_1y )  ',
     'w'   : '- visc* ( [w]_2xx + [w]_2yy + [w]_2zz ) - visc*0.3333333333333333_wp*( [ '+ ddivops+' ]_1z )  ',
     'et'  : '- visc* ( ( [u]_2xx + [u]_2yy + [u]_2zz )*u '
             ' +        ( [v]_2xx + [v]_2yy + [v]_2zz )*v '
             ' +        ( [w]_2xx + [w]_2yy + [w]_2zz )*w )'
             '- visc*0.3333333333333333_wp*('
             '  [ '+ ddivops+' ]_1x * u '
             '+ [ '+ ddivops+' ]_1y * v '
             '+ [ '+ ddivops+' ]_1z * w )'
             '- visc* ('
             ' ( 2.0_wp * [u]_1x * [u]_1x ) '
             '+( [u]_1y + [v]_1x ) * [v]_1x '
             '+( [u]_1z + [w]_1x ) * [w]_1x '
             '+( 2.0_wp * [v]_1y * [v]_1y ) '
             '+( [v]_1x + [u]_1y ) * [u]_1y '
             '+( [v]_1z + [w]_1y ) * [w]_1y '
             '+( 2.0_wp * [w]_1z * [w]_1z ) '
             '+( [w]_1x + [u]_1z ) * [u]_1z '
             '+( [w]_1y + [v]_1z ) * [v]_1z '
             ')'
             '+ 0.6666666666666666_wp*visc* ( ' + divops+' )*( ' + divops+' ) '} 

Src_Laplace['et'] = Src_Laplace['et'] + '- kappa* ( [T]_2xx + [T]_2yy + [T]_2zz )'                          
