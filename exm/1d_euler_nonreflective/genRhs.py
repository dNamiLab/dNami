import sys
import re
import numpy as np
from genKer import rhsinfo, genrk3, genrk3update, genFilter, genBC, append_Rhs, genbcsrc
import os 

# Set working precision
wp = 'float64'
  
from rhs import *

# Set install path (needed by genKer)
instpath = os.environ['INSTALLPATH']
incPATH  = instpath+'/src_for/includes/gen/'

# Specify global number of halos
hlo_glob = 5

def main():
      
    from genKer import rhs_info    
    rhs = rhs_info()

# Generate LHS:
    genrk3(len(varsolved)      ,rhs=rhs) 
    genrk3update(len(varsolved),rhs=rhs)

# Generate RHS:
    Save_eqns = {'divF':divF.copy()}
    append_Rhs(divF, 5,4, rhsname, vnamesrc_divF, update=False,rhs=rhs,stored=True)                           

# Generate Filters (if required):      
    genFilter(11,10, len(varsolved),rhs=rhs)

# Progressive stencil/order adjustement from domain to boundary 
    genBC(Save_eqns['divF']  ,3,2, rhsname , vnamesrc_divF, update=False,rhs=rhs)

# Boundary conditions on d(q)/dt 
    #i1
    genBC(src_phybc_wave_i1,3,2,rhsname , vnamesrc_divF, setbc=[True,{'char':{'i1':['rhs']}}]  , update=False,rhs=rhs)
    #imax
    genBC(src_phybc_wave_imax,3,2,rhsname ,vnamesrc_divF, setbc=[True,{'char':{'imax':['rhs']}}]  , update=False,rhs=rhs)

    # Extract RHS info:
    rhsinfo(rhs)

if __name__ == '__main__':
    main()
