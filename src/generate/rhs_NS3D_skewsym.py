
dim = 3 # 1 :1D, 2 : 2D, 3 : 3D

# Coefficitents needed in the pdes

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

varloc       = {'p': 'gamma_m1*rho*(e)',
                'e': 'et-0.5_wp*(u*u+v*v+w*w)',
                'T': 'CvInv*(e)'}
varstored    = {}

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

# Start specifying the PDEs

# IMPORTANT SYNTAX RULE :
# put a space after derivative symbol []_1* 


divops = ' [u]_1x + [v]_1y + [w]_1z '

# Euler Skew-symetric forme RHS  

Src_skew = {'rho' : '0.5_wp*( [rho*u]_1x + [rho*v]_1y + [rho*w]_1z '
                   '+   u*[rho]_1x + v*[rho]_1y + w*[rho]_1z '
                   '+ rho*( '+divops+' ) )',
        'u'   : '0.5_wp*( [rho*u*u]_1x + [rho*u*v]_1y + [rho*u*w]_1z '
                   '+   u*[rho*u]_1x + v*[rho*u]_1y + w*[rho*u]_1z '
                   '+rho*u*( '+divops+' ) ) + [p]_1x ',
        'v'   : '0.5_wp*( [rho*u*v]_1x + [rho*v*v]_1y + [rho*v*w]_1z '
                   '+   u*[rho*v]_1x + v*[rho*v]_1y + w*[rho*v]_1z '
                   '+rho*v*( '+divops+') ) + [p]_1y ',
        'w'   : '0.5_wp*( [rho*w*u]_1x + [rho*w*v]_1y + [rho*w*w]_1z '
                   '+   u*[rho*w]_1x + v*[rho*w]_1y + w*[rho*w]_1z '
                   '+rho*w*( '+divops+') ) + [p]_1z ',
        'et'  : '0.5_wp*( [rho*et*u]_1x + [rho*et*v]_1y + [rho*et*w]_1z '
                   '+   u*[rho*et]_1x + v*[rho*et]_1y + w*[rho*et]_1z '
                   '+rho*et*( '+divops+') ) + [p*u]_1x + [p*v]_1y + [p*w]_1z '}

# Skew = {'rho' : '0.5_wp*( [rho*u]_1x + [rho*v]_1y + [rho*w]_1z '
#                    '+   u*[rho]_1x + v*[rho]_1y + w*[rho]_1z '
#                    '+ rho*( '+divops+' ) )',
#         'u'   : '0.5_wp*( [rho*u*u]_1x + [rho*u*v]_1y + [rho*u*w]_1z '
#                    '+   u*[rho*u]_1x + v*[rho*u]_1y + w*[rho*u]_1z '
#                    '+rho*u*( '+divops+' ) ) + [p]_1x ',
#         'v'   : '0.5_wp*( [rho*u*v]_1x + [rho*v*v]_1y + [rho*v*w]_1z '
#                    '+   u*[rho*v]_1x + v*[rho*v]_1y + w*[rho*v]_1z '
#                    '+rho*v*( '+divops+') ) + [p]_1y ',
#         'w'   : '0.5_wp*( [rho*w*u]_1x + [rho*w*v]_1y + [rho*w*w]_1z '
#                    '+   u*[rho*w]_1x + v*[rho*w]_1y + w*[rho*w]_1z '
#                    '+rho*w*( '+divops+') ) + [p]_1z ',
#         'et'  : '0.5_wp*( [(rho*et+p)*u]_1x + [(rho*et+p)*v]_1y + [(rho*et+p)*w]_1z '
#                    '+   u*[rho*et+p]_1x + v*[rho*et+p]_1y + w*[rho*et+p]_1z '
#                    '+(rho*et+p)*( '+divops+') ) '}

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