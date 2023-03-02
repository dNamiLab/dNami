# -----------------------------------------------------------------------------
# Divergence and gradient on spherical shell 
#                                                          -sw 15-JUN-20
# -----------------------------------------------------------------------------

def div1d(v1, v2):
    v1mid = v1.replace('{', '[').replace('}',']')
    res = 'one_over_Rad*sintinv * ( [ '+v1+'  ]_1x * sint  + cost *( '+v1mid+' ) + [ '+v2+' ]_1y )' 
    return res

def grad1d(v1, rturn=None):
    res1 = 'one_over_Rad      * ( ['+v1+']_1x )' 
    res2 = 'one_over_Rad*sintinv * ( ['+v1+']_1y )' 
    if rturn == None:
        return [res1,res2]
    elif rturn == 0:
        return res1
    elif rturn == 1:
        return res2

