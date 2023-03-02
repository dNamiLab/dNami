# -----------------------------------------------------------------------------
# 2D Two-way coupled long-wave equations on a sphere for the Tonga explosion case  
#
#                                                          -sw 15-JUN-20
# -----------------------------------------------------------------------------
import sys
import re
import numpy as np
from genKer import genrk3, genrk3update, genFilter, genBC, append_Rhs
import os 

wp = 'float64'
from equations import *

# Set install path (needed by genKer)
instpath = os.environ['INSTALLPATH']
incPATH  = instpath+'/src_for/includes/gen/'

hlo_glob = 6

def main():

    from genKer import rhs_info    
    rhs = rhs_info(dim,wp,hlo_glob,incPATH,varsolved,varname,  
                 consvar=consvar,varloc=varloc,varstored=varstored,
                 coefficients=coefficients,premult='phi')
# Generate LHS:
    genrk3(len(varsolved)      ,rhs=rhs)
    genrk3update(len(varsolved),rhs=rhs)

# Generate RHS:
    append_Rhs(Adv , 5,4, rhsname, vnamesrc_Adv, update=False,rhs=rhs,stored=True)                           
    append_Rhs(S   , 5,4, rhsname, vnamesrc_S   , update=True, rhs=rhs,stored=False)                           
    append_Rhs(C   , 5,4, rhsname, vnamesrc_C   , update=True, rhs=rhs,stored=False)                           
    append_Rhs(Fss , 5,4, rhsname, vnamesrc_Fss , update=True, rhs=rhs,stored=False)                           

# Generate Filters (if required):      
    genFilter(13,8, len(varsolved),rhs=rhs)

    # Extract RHS info:
    rhs.export()

if __name__ == '__main__':
    main()
