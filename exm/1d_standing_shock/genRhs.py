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
hlo_glob = 5

def main():
      
    from genKer import rhs_info    
    rhs = rhs_info(dim,wp,hlo_glob,incPATH,varsolved,varname, 
                   consvar=consvar,varloc=varloc,
                   coefficients=coefficients)

# Generate LHS:
    genrk3(len(varsolved)      ,rhs=rhs) 
    genrk3update(len(varsolved),rhs=rhs)

# Generate RHS:
    append_Rhs(divF, 11,10, rhsname, vnamesrc_divF, update=False,rhs=rhs,stored=True)                           
    append_Rhs(visc,  5, 4, rhsname, vnamesrc_visc, update=True,rhs=rhs,stored=False)

# Generate Filters (if required):      
    genFilter(11,10, len(varsolved),rhs=rhs)

# RHS where full stencil doesnt fit:
    genBC(divF  ,11,10, rhsname , vnamesrc_divF, update=False,rhs=rhs)
    genBC(visc  , 5, 4, rhsname , vnamesrc_visc, update=True,rhs=rhs)

# Boundary conditions
    #i1
    genBC(src_phybc_wave_i1,11,10,rhsname , vnamesrc_divF, setbc=[True,{'char':{'i1':['rhs']}}]  , update=False,rhs=rhs)
    #imax
    genBC(src_phybc_wave_imax,11,10,rhsname ,vnamesrc_divF, setbc=[True,{'char':{'imax':['rhs']}}]  , update=False,rhs=rhs)

    # Extract RHS info:
    rhs.export()

if __name__ == '__main__':
    main()
