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

# - Number of halo points 
hlo_glob = 5 

def main(): 
      
      from genKer import rhs_info    
      rhs = rhs_info()

# Generate LHS:
      genrk3(len(varsolved)      ,rhs=rhs) 
      genrk3update(len(varsolved),rhs=rhs)

# Generate RHS:
      Fbc = F.copy()
      append_Rhs(F,5,4, rhsname,locname_f,update=False,rhs=rhs)                           

# Generate Bcs:
      genBC(Fbc  ,3,2, rhsname , locname_f, update=False,rhs=rhs)

      genBC(qbc,3,2,rhsname , locname_bc, setbc=[True,{'fakebc':{'i1'      :['q']}}]  , update=False,rhs=rhs)
      genBC(qbc,3,2,rhsname , locname_bc, setbc=[True,{'fakebc':{'imax'    :['q']}}]  , update=False,rhs=rhs)
      genBC(qbc,3,2,rhsname , locname_bc, setbc=[True,{'fakebc':{'j1'      :['q']}}]  , update=False,rhs=rhs)
      genBC(qbc,3,2,rhsname , locname_bc, setbc=[True,{'fakebc':{'jmax'    :['q']}}]  , update=False,rhs=rhs)
      genBC(qbc,3,2,rhsname , locname_bc, setbc=[True,{'fakebc':{'k1'      :['q']}}]  , update=False,rhs=rhs)
      genBC(qbc,3,2,rhsname , locname_bc, setbc=[True,{'fakebc':{'kmax'    :['q']}}]  , update=False,rhs=rhs)

# Extract RHS info:
      rhsinfo(rhs)

if __name__ == '__main__':
    main()
