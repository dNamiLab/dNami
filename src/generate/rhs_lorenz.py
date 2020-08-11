
dim = 3 # 1 :1D, 2 : 2D, 3 : 3D

# Coefficitents needed in the pdes

coefficients = {'sigma'       : 1, 
                'beta'        : 2,  
                'rho'         : 3}

varname      = {'p'   : 1}
                  
varsolved = ['p']

varstored    = {'x' : {'symb': ' ', 'ind':1 },
                'y' : {'symb': ' ', 'ind':2 },
                'z' : {'symb': ' ', 'ind':3 }
                }
         
# varstored    = {}        
consvar      = []            
varloc       = {}


# These names are for the comments in the src

rhsname = {'p' : 'd(p)/dt'}

# These names are for the variables in the fortran src

locname_dif = {'p': 'conv_p'}

F_x = {'x':  'sigma*(x-y)',
       'y':  'rho*x-y-x*z',
       'z':'-beta*z + x*y'}


Src_conv = {'p' : ' [ ( '+ F_x['x'] + ' )* p ]_1x + ' + ' [ ( '+ F_x['y'] +' )* p  ]_1y + '+ ' [ ( '+ F_x['z'] +')* p  ]_1z '}

