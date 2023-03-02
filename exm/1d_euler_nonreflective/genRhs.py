import sys
import re
import numpy as np
from genKer import genrk3, genrk3update, genFilter, genBC, append_Rhs
import os 

# Set working precision
wp = 'float64'
  
from equations import *

# Set install path (needed by genKer)
instpath = os.environ['INSTALLPATH']
incPATH  = instpath+'/src_for/includes/gen/'

# Specify global number of halos
hlo_glob = 5

def main():
      
    from genKer import rhs_info    
    # rhs = rhs_info(dim,wp,hlo_glob,incPATH,varsolved,varname, 
    #                consvar=consvar,varstored=varstored,varloc=varloc,varbc=varbc,
    #                coefficients=coefficients) 
    
    rhs = rhs_info(dim,wp,hlo_glob,incPATH,varsolved,varname,  
                   consvar=consvar,varloc=varloc,
                   coefficients=coefficients)

# Generate LHS:
    genrk3(len(varsolved)      ,rhs=rhs) 
    genrk3update(len(varsolved),rhs=rhs)

# Generate RHS:
    append_Rhs(divF, 5,4, rhsname, vnamesrc_divF, update=False,rhs=rhs,stored=True)                           
# Generate Filters (if required):      
    genFilter(11,10, len(varsolved),rhs=rhs)

# Progressive stencil/order adjustement from domain to boundary 
    genBC(divF,3,2, rhsname , vnamesrc_divF, update=False,rhs=rhs)

# Boundary conditions on d(q)/dt 
    #i1
    genBC(src_phybc_wave_i1,3,2,rhsname , vnamesrc_divF, setbc=[True,{'char':{'i1':['rhs']}}]  , update=False,rhs=rhs)
    #imax
    genBC(src_phybc_wave_imax,3,2,rhsname ,vnamesrc_divF, setbc=[True,{'char':{'imax':['rhs']}}]  , update=False,rhs=rhs)

    # Extract RHS info:
    rhs.export()

if __name__ == '__main__':
    main()
