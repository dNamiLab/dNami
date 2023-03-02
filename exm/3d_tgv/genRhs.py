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

# - Number of halo points 
hlo_glob = 6 

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
# Adding convective and diffusive terms to the RHS with different schemes
      append_Rhs(Src_conv, 13, 4, rhsname,locname_conv,update=False,rhs=rhs)                           
      append_Rhs(Src_dif , 5 , 4 , rhsname,locname_dif ,update=True ,rhs=rhs,stored=True)

# Generate Filters (if required):      
      genFilter(13,8, len(varsolved),rhs=rhs)

# Extract RHS info:
      rhs.export()

if __name__ == '__main__':
    main()
