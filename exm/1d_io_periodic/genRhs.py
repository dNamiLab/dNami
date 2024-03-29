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
      rhs = rhs_info()

# Generate LHS:
      genrk3(len(varsolved)      ,rhs=rhs)
      genrk3update(len(varsolved),rhs=rhs)

# Generate RHS:
      append_Rhs(F, 13, 4, rhsname,locname_f,update=False,rhs=rhs)

# Extract RHS info:
      rhs.export()

if __name__ == '__main__':
    main()
