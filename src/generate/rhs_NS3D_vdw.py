
dim = 3 # 1 :1D, 2 : 2D, 3 : 3D

# Coefficitents needed in the pdes

a = '3.0_wp'
b = '0.3333333333333333_wp'

coefficients = {'visc'       : 1, 
                'kappa'      : 2,  
                'gamma_m1'   : 3,
                'CvInv'      : 4}

varname      = {'rho' : 1,
                  'u' : 2, 
                  'v' : 3, 
                  'w' : 4,
                  'et': 5}

consvar      = [2,3,4,5]            

varloc       = {'p': 'rho*gamma_m1/(1.0_wp-'+b+'*rho)*(e +'+a+'*rho)-'+a+'*rho*rho',
                'e': 'et-0.5_wp*(u*u+v*v+w*w)',
                'T': 'CvInv*(e+'+a+'*rho)'}

# These names are for the comments in the src

rhsname = {'rho' : 'd(rho)/dt'  
          ,  'u' : 'd(rho u)/dt',
             'v' : 'd(rho v)/dt', 
             'w' : 'd(rho w)/dt', 
             'et': 'd(rho et)/dt'}
       

vname = {'rho': {'x':'Flx_rho ','y':'Fly_rho ','z':'Flz_rho '},
         'u'  : {'x':'Flx_rhou','y':'Fly_rhou','z':'Flz_rhou'},
         'v'  : {'x':'Flx_rhov','y':'Fly_rhov','z':'Flz_rhov'},
         'w'  : {'x':'Flx_rhow','y':'Fly_rhow','z':'Flz_rhow'},
         'et' : {'x':'Flx_et  ','y':'Fly_et  ','z':'Flz_et  '}}


vnamesrc = {'rho': 'src_rho',
            'u'  : 'src_rhou',
            'v'  : 'src_rhov',
            'w'  : 'src_rhow',
            'et' : 'src_et  '}

# Start specifying the PDEs

# IMPORTANT SYNTAX RULE :
# put a space after derivative symbol []_1* 


divops = '[u]_1x + [v]_1y + [w]_1z'

# Euler RHS        

Fx = {'rho' : 'rho*u         ',
      'u'   : 'rho*u*u  + p  ',
      'v'   : 'rho*u*v       ',
      'w'   : 'rho*u*w       ',
      'et'  : '(rho*et + p)*u '}

Fy = {'rho' : 'rho*v         ',
      'u'   : 'rho*v*u       ', 
      'v'   : 'rho*v*v  + p  ', 
      'w'   : 'rho*v*w       ', 
      'et'  : '(rho*et + p)*v'} 
                             # 
       # 
Fz = {'rho' : 'rho*w         ',
      'u'   : 'rho*w*u       ',
      'v'   : 'rho*w*v       ',
      'w'   : 'rho*w*w  + p  ',
      'et'  : '(rho*et + p)*w'}

Flxconv = { 'x': Fx , 'y': Fy, 'z': Fz }

# Navier-Stokes Diffusive termes only

Fx = {'u'   : ' - visc *( 2.0_wp * [u]_1x - 2.0_wp/3.0_wp * ( '+ divops +'  ) )',
      'v'   : ' - visc *( [u]_1y + [v]_1x )', 
      'w'   : ' - visc *( [u]_1z + [w]_1x )',
      
      'et'  : ' - kappa*( [T]_1x ) '
              ' - u*(visc *( 2.0_wp * [u]_1x - 2.0_wp/3.0_wp * ( '+ divops +'  )))'
              ' - v*(visc *( [u]_1y + [v]_1x ))'
              ' - w*(visc *( [u]_1z + [w]_1x ))'}

Fy = {'u'   : ' - visc *( [u]_1y + [v]_1x )  ',
      'v'   : ' - visc *( 2.0_wp * [v]_1y - 2.0_wp/3.0_wp * ( '+ divops +'  ) )', 
      'w'   : ' - visc *( [v]_1z + [w]_1y )  ',
      'et'  : ' - kappa*( [T]_1y )'
              ' - u*(visc *( [u]_1y + [v]_1x ))'
              ' - v*(visc *( 2.0_wp * [v]_1y - 2.0_wp/3.0_wp * ( '+ divops +'  )))'
              ' - w*(visc *( [v]_1z + [w]_1y ))'}
       
Fz = {'u'   : ' - visc *( [u]_1z + [w]_1x )',
      'v'   : ' - visc *( [v]_1z + [w]_1y )', 
      'w'   : ' - visc *( 2.0_wp * [w]_1z - 2.0_wp/3.0_wp * ( '+ divops +'  ) )',
      'et'  : ' - kappa*( [T]_1z )'
              ' - u*(visc *( [u]_1z + [w]_1x ))'
              ' - v*(visc *( [v]_1z + [w]_1y ))'
              ' - w*(visc *( 2.0_wp * [w]_1z - 2.0_wp/3.0_wp * ( '+ divops +'  )))'}

Flxdif = { 'x': Fx , 'y': Fy, 'z': Fz }

ddivops = '{u}_1x + {v}_1y + {w}_1z'

# Navier-Stokes Diffusuve termes only (Laplacian forme)

SLaplace = {'u'   : '- visc* ( [u]_2xx + [u]_2yy + [u]_2zz ) - visc*0.3333333333333333_wp*( [ '+ ddivops+' ]_1x )  ',
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

SLaplace['et'] = SLaplace['et'] + '- kappa* ( [T]_2xx + [T]_2yy + [T]_2zz )'


# Navier-Stokes (Full set of eqs)

# Fx = {'rho' : 'rho*u         ',
#       'u'   : 'rho*u*u  + p    - visc *( 2.0_wp * [u]_1x - 2.0_wp/3.0_wp * ( '+ divops +'  ) )',
#       'v'   : 'rho*u*v         - visc *( [u]_1y + [v]_1x )', 
#       'w'   : 'rho*u*w         - visc *( [u]_1z + [w]_1x )',
#       'et'  : '(rho*et + p)*u  - kappa*( [T]_1x ) '
#                               '- u*(visc *( 2.0_wp * [u]_1x - 2.0_wp/3.0_wp * ( '+ divops +'  )))'
#                               '- v*(visc *( [u]_1y + [v]_1x ))'
#                               '- w*(visc *( [u]_1z + [w]_1x ))'}

# Fy = {'rho' : 'rho*v         ',
#       'u'   : 'rho*v*u         - visc *( [u]_1y + [v]_1x )  ',
#       'v'   : 'rho*v*v  + p    - visc *( 2.0_wp * [v]_1y - 2.0_wp/3.0_wp * ( '+ divops +'  ) )', 
#       'w'   : 'rho*v*w         - visc *( [v]_1z + [w]_1y )  ',
#       'et'  : '(rho*et + p)*v  - kappa*( [T]_1y )'
#                               '- u*(visc *( [u]_1y + [v]_1x ))'
#                               '- v*(visc *( 2.0_wp * [v]_1y - 2.0_wp/3.0_wp * ( '+ divops +'  )))'
#                               '- w*(visc *( [v]_1z + [w]_1y ))'}
       
# Fz = {'rho' : 'rho*w         ',
#       'u'   : 'rho*w*u        - visc *( [u]_1z + [w]_1x )',
#       'v'   : 'rho*w*v        - visc *( [v]_1z + [w]_1y )', 
#       'w'   : 'rho*w*w  + p   - visc *( 2.0_wp * [w]_1z - 2.0_wp/3.0_wp * ( '+ divops +'  ) )',
#       'et'  : '(rho*et + p)*w - kappa*( [T]_1z )'
#                              '- u*(visc *( [u]_1z + [w]_1x ))'
#                              '- v*(visc *( [v]_1z + [w]_1y ))'
#                              '- w*(visc *( 2.0_wp * [w]_1z - 2.0_wp/3.0_wp * ( '+ divops +'  )))'}

# Flx = { 'x': Fx , 'y': Fy, 'z': Fz }                             