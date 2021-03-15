import sys
import re
import numpy as np
from genKer import rhsinfo, genrk3, genrk3update, genFilter, genBC, append_Rhs, genbcsrc
import os 

wp = 'float64'
  
# from  rhs_NS2D_test import *
from  rhs_NS3D import *

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
      Save_eqns = {}
      Save_eqns = {'Src_conv':Src_conv.copy(),'Src_dif': Src_dif.copy()}
      # append_Rhs(Src_conv, 11, 10, rhsname,locname_conv,update=False,rhs=rhs)                           
      # append_Rhs(Src_dif , 5 , 4 , rhsname,locname_dif ,update=True ,rhs=rhs,stored=True)

      append_Rhs(Src_conv, 13, 4, rhsname,locname_conv,update=False,rhs=rhs)                           
      append_Rhs(Src_dif , 5 , 4 , rhsname,locname_dif ,update=True ,rhs=rhs,stored=True)

# Generate Filters (if required):      
      genFilter(13,8, len(varsolved),rhs=rhs)

# Generate BCs:
      # genBC(Save_eqns['Src_conv'] ,11,10,rhsname , locname_conv, update=False,rhs=rhs)
      # genBC(Save_eqns['Src_dif']  ,5 ,4 ,rhsname , locname_dif , update=True ,rhs=rhs,stored=True)
      # # j1
      # genBC(Src_phybc_rhs ,5,4,rhsname , locname_bc, setbc=[True,{'wall':{'j1':['rhs']}}]  , update=False,rhs=rhs)
      # genBC(Src_phybc_q1  ,5,4,rhsname , locname_bc, setbc=[True,{'wall':{'j1':['q'  ]}}]  , update=False,rhs=rhs)
      # # jmax
      # genBC(Src_phybc_rhs ,5,4,rhsname , locname_bc, setbc=[True,{'wall':{'jmax':['rhs']}}], update=False,rhs=rhs)
      # genBC(Src_phybc_qmax,5,4,rhsname , locname_bc, setbc=[True,{'wall':{'jmax':['q'  ]}}], update=False,rhs=rhs)
      # # i1
      # genBC(Src_phybc_rhs ,5,4,rhsname , locname_bc, setbc=[True,{'wall':{'i1':['rhs']}}]  , update=False,rhs=rhs)
      # genBC(Src_phybc_q1  ,5,4,rhsname , locname_bc, setbc=[True,{'wall':{'i1':['q'  ]}}]  , update=False,rhs=rhs)
      # # imax
      # genBC(Src_phybc_rhs ,5,4,rhsname , locname_bc, setbc=[True,{'wall':{'imax':['rhs']}}], update=False,rhs=rhs)
      # genBC(Src_phybc_qmax,5,4,rhsname , locname_bc, setbc=[True,{'wall':{'imax':['q'  ]}}], update=False,rhs=rhs)      

      # # genBC(bc_test ,5,4,rhsname , locname_conv, update=False,rhs=rhs)
# Extract RHS info:
      rhsinfo(rhs)

if __name__ == '__main__':
    main()
