import sys
import re
import numpy as np
from genKer import rhsinfo, genrk3, genrk3update, genFilter, genBC, append_Rhs, genbcsrc
import os 

wp = 'float64'
iVisc = True  
  
from  rhs import *

# Set install path (needed by genKer)
instpath = os.environ['INSTALLPATH']
incPATH  = instpath+'/src_for/includes/gen/'

def main():
      
    from genKer import rhs_info    
    rhs = rhs_info()

# Generate LHS:
    genrk3(len(varsolved)      ,rhs=rhs) 
    genrk3update(len(varsolved),rhs=rhs)

# Generate RHS:
    Save_eqns = {'divF':divF.copy()}
    if iVisc:
        Save_eqns = {'divF':divF.copy(), 'visc':visc.copy()}

    append_Rhs(divF, 5,4, rhsname, vnamesrc_divF, update=False,rhs=rhs,stored=True)                           
    if iVisc:
        append_Rhs(visc , 5,4, rhsname, vnamesrc_visc , update=True ,rhs=rhs,stored=False)                           

# Generate Filters (if required):      
    genFilter(11,10, len(varsolved),rhs=rhs)

    # Extract RHS info:
    rhsinfo(rhs)

if __name__ == '__main__':
    main()
