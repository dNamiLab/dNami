
dim = 2 # 1 :1D, 2 : 2D, 3 : 3D

# Coefficitents needed in the pdes

coefficients = {'visc'       : 1, 
                'kappa'      : 2,  
                'gamma_m1'   : 3,
                'CvInv'      : 4,
                'u_0'        : 5}

varname      = {'rho' : 1,
                  'u' : 2,  
                  'v' : 3, 
                  'et': 4}
                  
varsolved = ['rho','u','v','et']

varstored    = {'divU' : {'symb': ' [u]_1x + [v]_1y ', 
                          'ind':1 }}
                 # 'e'   : {'symb': ' ( et-0.5_wp*(u*u+v*v+w*w) ) ', 
                 #         'ind':7 },
                 # 'p' : {'symb': '( gamma_m1*rho*(e) ) ', 
                 #        'ind': 8 }}         

divops  = ' [u]_1x + [v]_1y  '
ddivops = ' {u}_1x + {v}_1y  '

# divops  = 'divU'
# ddivops = 'divU'

consvar      = [2,3,4]            

# VdW EOS:
a = '3.0_wp'
b = '0.3333333333333333_wp'

# varloc       = {'e': 'et-0.5_wp*(u*u+v*v+w*w)',
#                 'p': 'rho*gamma_m1/(1.0_wp-'+b+'*rho)*(e +'+a+'*rho)-'+a+'*rho*rho',                
#                 'T': 'CvInv*(e+'+a+'*rho)'}
# IG EOS:
varloc       = {'p': 'gamma_m1*rho*(e)',
                'e': '(et-0.5_wp*(u*u+v*v))',
                'T': 'CvInv*(e)'}
# WARNING: the order of definition matters:
# example:
# varstored = {'p' : {'symb': 'gamma_m1*rho*(e)',
#                     'ind': 6 },
#                     'e' : {'symb': 'et-0.5_wp*(u*u+v*v+w*w)',
#                      'ind':7 }}
# IS NOT VALID because 'p' uses 'e' which appears after in varstored dictionary...


# IMPORTANT limitation: stored variables cannot appear in derivatives (for now...), which makes them mostly useless...

varstored    = {} # This dictionary is compulsory (declare it as empty if necessary.)

# varstored    = {'e' : {'symb': ' ( et-0.5_wp*(u*u+v*v+w*w) ) ', 
#                          'ind':6 },
#                 'p' : {'symb': '( gamma_m1*rho*(e) ) ', 
#                         'ind': 7 }}

# varstored    = {#'e' : {'symb': ' ( et-0.5_wp*(u*u+v*v+w*w) ) ', 
#                          'ind':6 },
#                 'p' : {'symb': '( rho*gamma_m1/(1.0_wp-'+b+'*rho)*(e +'+a+'*rho)-'+a+'*rho*rho) ', 
#                        'ind': 6 }}#,
#                 'T' : {'symb' : 'CvInv*(e+'+a+'*rho)',
#                         'ind': 8}}

# varstored = { #'e' : {'symb': ' ( et-0.5_wp*(u*u+v*v+w*w) ) ', 
#                          # 'ind':6 },
#               # 'p' : {'symb': '( gamma_m1*rho*(e) ) ', 
#               #           'ind': 7 },
#                # 'T':  {'symb': 'CvInv*(e)', 'ind' : 7},           
#               'visc' : {'symb' : '1_Re * (2.0_wp*sqrt(T)**(3)/(T+0.5_wp) )', 'ind' : 6}}

# These names are for the comments in the src

rhsname = {'rho' : 'd(rho)/dt'  
          ,  'u' : 'd(rho u)/dt',
             'v' : 'd(rho v)/dt', 
             'et': 'd(rho et)/dt'}
       

# These names are for the variables in the fortran src

locname_dif = {'rho': 'dif_rho',
               'u'  : 'dif_rhou',
               'v'  : 'dif_rhov',
               'et' : 'dif_et  '}

locname_conv = {'rho': 'conv_rho',
                'u'  : 'conv_rhou',
                'v'  : 'conv_rhov',
                'et' : 'conv_et  '}

locname_bc  = {'rho': 'bc_rho',
               'u'  : 'bc_u',
               'v'  : 'bc_v',
               'et' : 'bc_et  '}                

# Start specifying the PDEs

# IMPORTANT SYNTAX RULE :
# put a space after derivative symbol []_1* 

# Euler RHS        

Fx = {'rho' : 'rho*u         ',
      'u'   : 'rho*u*u  + p  ',
      'v'   : 'rho*u*v       ',
      'et'  : '(rho*et + p)*u '}

Fy = {'rho' : 'rho*v         ',
      'u'   : 'rho*v*u       ', 
      'v'   : 'rho*v*v  + p  ',  
      'et'  : '(rho*et + p)*v '} 
                             # 
       # 
Fz = {'rho' : 'rho*w         ',
      'u'   : 'rho*w*u       ',
      'v'   : 'rho*w*v       ',
      'et'  : '(rho*et + p)*w '}

Src_conv = {'rho' : '[ '+Fx['rho']+' ]_1x' + ' + ' + '[ '+Fy['rho']+' ]_1y ',
            'u'   : '[ '+Fx['u']  +' ]_1x' + ' + ' + '[ '+Fy['u']  +' ]_1y ',
            'v'   : '[ '+Fx['v']  +' ]_1x' + ' + ' + '[ '+Fy['v']  +' ]_1y ',
            'et' : ' [ '+Fx['et'] +' ]_1x' + ' + ' + '[ '+Fy['et'] +' ]_1y '}

# Navier-Stokes Diffusive termes only


Fx = {'u'   : ' [ {u}_1y ]_1y '}

Fy = {'u'   : ' - visc *( {u}_1y + {v}_1x )  ',
      'v'   : ' - visc *( 2.0_wp * {v}_1y - 2.0_wp/3.0_wp * ( '+ ddivops +'  ) )', 
      'et'  : ' - kappa*( {T}_1y )'
              ' - u*(visc *( {u}_1y + {v}_1x ))'
              ' - v*(visc *( 2.0_wp * {v}_1y - 2.0_wp/3.0_wp * ( '+ ddivops +'  )))'}
       
Fz = {'u'   : ' - visc *( {u}_1z + {w}_1x )',
      'v'   : ' - visc *( {v}_1z + {w}_1y )', 
      'et'  : ' - kappa*( {T}_1z )'
              ' - u*(visc *( {u}_1z + {w}_1x ))'
              ' - v*(visc *( {v}_1z + {w}_1y ))'
              ' - w*(visc *( 2.0_wp * {w}_1z - 2.0_wp/3.0_wp * ( '+ ddivops +'  )))'}

# Src_dif  = {'u'   : '[ '+Fx['u']  +' ]_1x' + ' + ' + '[ '+Fy['u']  +' ]_1y ' ,
#             'v'   : '[ '+Fx['v']  +' ]_1x' + ' + ' + '[ '+Fy['v']  +' ]_1y ' ,
#             'et' : ' [ '+Fx['et'] +' ]_1x' + ' + ' + '[ '+Fy['et'] +' ]_1y ' }


Src_skew = {'rho' : '0.5_wp*( [rho*u]_1x + [rho*v]_1y + [rho*w]_1z '
                   '+   u*[rho]_1x + v*[rho]_1y + w*[rho]_1z '
                   '+ rho*( '+divops+' ) )',
        'u'   : '0.5_wp*( [rho*u*u]_1x + [rho*u*v]_1y + [rho*u*w]_1z '
                   '+   u*[rho*u]_1x + v*[rho*u]_1y + w*[rho*u]_1z '
                   '+rho*u*( '+divops+' ) ) + [p]_1x ',
        'v'   : '0.5_wp*( [rho*u*v]_1x + [rho*v*v]_1y + [rho*v*w]_1z '
                   '+   u*[rho*v]_1x + v*[rho*v]_1y + w*[rho*v]_1z '
                   '+rho*v*( '+divops+' ) ) + [p]_1y ',
        'et'  : '0.5_wp*( [rho*et*u]_1x + [rho*et*v]_1y + [rho*et*w]_1z '
                   '+   u*[rho*et]_1x + v*[rho*et]_1y + w*[rho*et]_1z '
                   '+rho*et*( '+divops+') ) + [p*u]_1x + [p*v]_1y + [p*w]_1z '}

# Navier-Stokes Diffusive terms only (Laplacian form)

Src_Laplace = {'u'   : '- visc* ( [u]_2xx + [u]_2yy + [u]_2zz ) - visc*0.3333333333333333_wp*( [ '+ ddivops+' ]_1x )  ',
     'v'   : '- visc* ( [v]_2xx + [v]_2yy + [v]_2zz ) - visc*0.3333333333333333_wp*( [ '+ ddivops+' ]_1y )  ',
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


# Physical Boundary Conditions (Exemple Wall moving at u_0 in x direction):

PhyBcFx = {'rho' : 'rho*u_0         ',
           'u'   : 'rho*u_0*u_0  + p',
           'v'   : 'rho*u_0*v       ',
           'et'  : '(rho*et + p)*u_0 '}

PhyBcFy = {'rho' : 'rho*v         ',
           'u'   : 'rho*v*u_0       ', 
           'v'   : 'rho*v*v  + p  ', 
           'et'  : '(rho*et + p)*v '} 
                            # 
       # 
PhyBcFz = {'rho' : 'rho*w         ',
           'u'   : 'rho*w*u_0       ',
           'v'   : 'rho*w*v       ',
           'et'  : '(rho*et + p)*w '}



# Src_phybc_rhs = {'rho' : '[ '+PhyBcFx['rho']+' ]_1x' + ' + ' + '[ '+PhyBcFy['rho']+' ]_1y ' ,
#                   # 'u'   : '[ '+PhyBcFx['u']  +' ]_1x' + ' + ' + '[ '+PhyBcFy['u']  +' ]_1y' + ' + ' + '[ '+PhyBcFz['u']  +' ]_1z ',
#                   # 'v'   : '[ '+PhyBcFx['v']  +' ]_1x' + ' + ' + '[ '+PhyBcFy['v']  +' ]_1y' + ' + ' + '[ '+PhyBcFz['v']  +' ]_1z ',
#                   # 'w'   : '[ '+PhyBcFx['w']  +' ]_1x' + ' + ' + '[ '+PhyBcFy['w']  +' ]_1y' + ' + ' + '[ '+PhyBcFz['w']  +' ]_1z ',
#                   'et' : ' [ '+PhyBcFx['et'] +' '+ Fx['et'] +' ]_1x' + ' + ' + '[ '+PhyBcFy['et'] +' '+ Fy['et'] +' ]_1y ' }


Src_phybc_q   = { 'u'   : 'u_0',
                  'v'   : '0.0_wp'
                }
