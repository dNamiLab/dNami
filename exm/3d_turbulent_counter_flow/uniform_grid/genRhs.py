import sys
import re
import numpy as np
from genKer import rhsinfo, genrk3, genrk3update, genFilter, genBC, append_Rhs, genbcsrc, genrk_Williamson
import os 

wp = 'float64'
  
from  rhs import *

# Set install path (needed by genKer)
instpath = os.environ['INSTALLPATH']
incPATH  = instpath+'/src_for/includes/gen/'

def main():
   
    from genKer import rhs_info    
    rhs = rhs_info()
    RK_type = 'Williamson'
# Generate LHS:
    if RK_type == 'standard':
        genrk3(len(varsolved)      ,rhs=rhs) 
        genrk3update(len(varsolved),rhs=rhs)
    else:
        genrk_Williamson(len(varsolved),rhs=rhs, order=3, SSP=False)

# Generate RHS:

    Save_eqns = {'Euler':Euler.copy(), 'Diffusive':Diffusive.copy()}


    # Add terms to the right hand side with the following options:
    #
    #  Symbolic equations (see rhs.py)   Nb. of points in stencil   Order of scheme
    #             |     _________________|                                 |
    #             |    |   ________________________________________________|
    #             |    |  |
    #             v    v  v 
    append_Rhs(Euler, 13,4, rhsname, vname_Euler, update=False,rhs=rhs,stored=True)                           
    append_Rhs(Diffusive, 5,4, rhsname, vname_Diffusive, update=True ,rhs=rhs)                         

    # Generate Filters (if required):      
    genFilter(13,8, len(varsolved),rhs=rhs)

    # Save_eqns = {'Euler':Euler.copy(), 'Diffusive':Diffusive.copy()}


    genBC(Save_eqns['Euler'] ,13,4, rhsname , vname_Euler, update=False,rhs=rhs,stored=True)
    genBC(Save_eqns['Diffusive']  ,5,4, rhsname , vname_Diffusive , update=True ,rhs=rhs)

    # Boundary conditions
    # (Equations, Stencil, order, rhsname, vname, setbc, update, rhs, stored)
    # Lower wall j=1
    genBC(q_j1,5,4,rhsname , locname_bc, setbc=[True,{'bottom_wall':{'j1':['q']}}]  , update=False,rhs=rhs, stored=True)
    genBC(rhs_j1,5,4,rhsname , locname_bc, setbc=[True,{'bottom_wall':{'j1':['rhs']}}]  , update=False,rhs=rhs, stored=True)

    # # Top wall j=jmax
    genBC(q_jmax,5,4,rhsname ,locname_bc , setbc=[True,{'top_wall':{'jmax':['q']}}]  , update=False,rhs=rhs, stored=True)
    genBC(rhs_jmax,5,4,rhsname ,locname_bc , setbc=[True,{'top_wall':{'jmax':['rhs']}}]  , update=False,rhs=rhs, stored=True)
    # Extract RHS info:
    rhsinfo(rhs)

if __name__ == '__main__':
    main()
