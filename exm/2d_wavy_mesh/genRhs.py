import sys
import re
import numpy as np
from genKer import rhsinfo, genrk3, genrk3update, genFilter, genBC, append_Rhs, genbcsrc
import os 

wp = 'float64'
  
from rhs import *

# Set install path (needed by genKer)
instpath = os.environ['INSTALLPATH']
incPATH  = instpath+'/src_for/includes/gen/'

hlo_glob = 5
#hlo_glob = 4
#hlo_glob = 3
#hlo_glob = 2
#hlo_glob = 1

def main():
      
    from genKer import rhs_info    
    rhs = rhs_info()

# Generate LHS:
    genrk3(len(varsolved)      ,rhs=rhs) 
    genrk3update(len(varsolved),rhs=rhs)

# Generate RHS:
    append_Rhs(divF, 11,10, rhsname, vnamesrc_divFx, update=False,rhs=rhs,stored=True)                           
    #append_Rhs(divF, 9,8, rhsname, vnamesrc_divFx, update=False,rhs=rhs,stored=True)                           
    #append_Rhs(divF, 7,6, rhsname, vnamesrc_divFx, update=False,rhs=rhs,stored=True)                           
    #append_Rhs(divF, 5,4, rhsname, vnamesrc_divFx, update=False,rhs=rhs,stored=True)                           
    #append_Rhs(divF, 3,2, rhsname, vnamesrc_divFx, update=False,rhs=rhs,stored=True)                           

# Generate Filters (if required):      
    genFilter(11,10, len(varsolved),rhs=rhs)

    # Extract RHS info:
    rhsinfo(rhs)

if __name__ == '__main__':
    main()
