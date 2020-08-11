
dim = 2 # 1 :1D, 2 : 2D, 3 : 3D

# Coefficitents needed in the pdes

coefficients = {'visc'       : 1, 
                'kappa'      : 2,  
                'gamma_m1'   : 3,
                'CvInv'      : 4}

varname      = {'rho' : 1,
                  'u' : 2, 
                  'v' : 3, 
                  'et': 4}

consvar      = [2,3,4]            

varloc       = {'p': 'gamma_m1*rho*(e)',
                'e': 'et-0.5_wp*(u*u+v*v)',
                'T': 'CvInv*(e)'}

# These names are for the comments in the src

rhsname = {'rho' : 'd(rho)/dt'  
          ,  'u' : 'd(rho u)/dt',
             'v' : 'd(rho v)/dt', 
             'et': 'd(rho et)/dt'}
       

vname = {'rho': {'x':'Flx_rho ','y':'Fly_rho ','z':'Flz_rho '},
         'u'  : {'x':'Flx_rhou','y':'Fly_rhou','z':'Flz_rhou'},
         'v'  : {'x':'Flx_rhov','y':'Fly_rhov','z':'Flz_rhov'},
         'et' : {'x':'Flx_et  ','y':'Fly_et  ','z':'Flz_et  '}}


vnamesrc = {'rho': 'src_rho',
            'u'  : 'src_rhou',
            'v'  : 'src_rhov',
            'et' : 'src_et  '}

# Start specifying the PDEs

# IMPORTANT SYNTAX RULE :
# put a space after derivative symbol []_1* 

# Euler RHS        

Fx = {'rho' : 'rho*u         ',
      'u'   : 'rho*u*u  + p  ',
      'v'   : 'rho*u*v       ',
      'et'  : '(rho*et + p)*u '}
                              # 
                              # 
                              # 

Fy = {'rho' : 'rho*v         ',
      'u'   : 'rho*v*u       ', 
      'v'   : 'rho*v*v  + p  ', 
      'et'  : '(rho*et + p)*v'} 


Flxconv = { 'x': Fx , 'y': Fy}


divops = '[u]_1x + [v]_1y'

# Navier-Stokes Diffusive termes only

Fx = {'u'   : ' - visc *( 2.0_wp * [u]_1x - 2.0_wp/3.0_wp * ( '+ divops +'  ) )',
      'v'   : ' - visc *( [u]_1y + [v]_1x )',       
      'et'  : ' - kappa*( [T]_1x ) '
              ' - u*(visc *( 2.0_wp * [u]_1x - 2.0_wp/3.0_wp * ( '+ divops +'  )))'
              ' - v*(visc *( [u]_1y + [v]_1x ))'
              }

Fy = {'u'   : ' - visc *( [u]_1y + [v]_1x )  ',
      'v'   : ' - visc *( 2.0_wp * [v]_1y - 2.0_wp/3.0_wp * ( '+ divops +'  ) )',       
      'et'  : ' - kappa*( [T]_1y )'
              ' - u*(visc *( [u]_1y + [v]_1x ))'
              ' - v*(visc *( 2.0_wp * [v]_1y - 2.0_wp/3.0_wp * ( '+ divops +'  )))'
              }


Flxdif = { 'x': Fx , 'y': Fy}

ddivops = '{u}_1x + {v}_1y'

# Navier-Stokes Diffusuve termes only (Laplacian forme)

SLaplace = {'u'   : '- visc* ( [u]_2xx + [u]_2yy ) - visc*0.3333333333333333_wp*( [ '+ ddivops+' ]_1x )  ',
       'v'   : '- visc* ( [v]_2xx + [v]_2yy ) - visc*0.3333333333333333_wp*( [ '+ ddivops+' ]_1y )  ',         
       'et'  : '- visc* ( ( [u]_2xx + [u]_2yy  )*u '
               ' +        ( [v]_2xx + [v]_2yy  )*v '
               ' )'
               '- visc*0.3333333333333333_wp*('
               '  [ '+ ddivops+' ]_1x * u '
               '+ [ '+ ddivops+' ]_1y * v '
               '  )'
               '- visc* ('
               ' ( 2.0_wp * [u]_1x * [u]_1x ) '
               '+( [u]_1y + [v]_1x ) * [v]_1x '
               '  '
               '+( 2.0_wp * [v]_1y * [v]_1y ) '
               '+( [v]_1x + [u]_1y ) * [u]_1y '
               '   '
               '   '
               '   '
               '   '
               ')'
               '+ 0.6666666666666666_wp*visc* ( ' + divops+' )*( ' + divops+' ) '} 

SLaplace['et'] = SLaplace['et'] + '- kappa* ( [T]_2xx + [T]_2yy )'
