import sys
import re
import numpy as np
from genKer import rhsinfo, genrk3, genrk3update, genFilter, genBC, append_Rhs, genbcsrc
import os 

wp = 'float64'
from  rhs import *

# Set install path (needed by genKer)
instpath = os.environ['INSTALLPATH']
incPATH  = instpath+'/src_for/includes/gen/'

hlo_glob = 5 

def main():
      
    from genKer import rhs_info    
    rhs = rhs_info()

# Generate LHS:
    genrk3(len(varsolved)      ,rhs=rhs) 
    genrk3update(len(varsolved),rhs=rhs)

# Generate RHS:
    append_Rhs(divF, 9,8, rhsname, vnamesrc_divF, update=False,rhs=rhs)                           

# Generate Filters (if required):      
    genFilter(11,10, len(varsolved),rhs=rhs)

# Boundary conditions
    genBC(divF  ,9,8, rhsname , vnamesrc_divF, update=False,rhs=rhs)

    # ------------ I 
    #i1
    genBC(src_phybc_wave_i1,5,4,rhsname , locname_bc, setbc=[True,{'char':{'i1':['rhs']}}]  , update=False,rhs=rhs)
    #imax
    genBC(src_phybc_wave_imax,5,4,rhsname ,locname_bc , setbc=[True,{'char':{'imax':['rhs']}}]  , update=False,rhs=rhs)

    # ------------ J
    #j1
    genBC(src_phybc_wave_j1,5,4,rhsname , locname_bc, setbc=[True,{'char':{'j1':['rhs']}}]  , update=False,rhs=rhs)
    #jmax
    genBC(src_phybc_wave_jmax,5,4,rhsname ,locname_bc , setbc=[True,{'char':{'jmax':['rhs']}}]  , update=False,rhs=rhs)

    # ------------ Corners
    # - i1/j1 
    genBC(src_phybc_wave_i1j1,5,4,rhsname , locname_bc, setbc=[True,{'cornerchar':{'i1j1':['rhs']}}]  , update=False,rhs=rhs)
    # - i1/jmax 
    genBC(src_phybc_wave_i1jmax,5,4,rhsname , locname_bc, setbc=[True,{'cornerchar':{'i1jmax':['rhs']}}]  , update=False,rhs=rhs)
    # - imax/j1 
    genBC(src_phybc_wave_imaxj1,5,4,rhsname , locname_bc, setbc=[True,{'cornerchar':{'imaxj1':['rhs']}}]  , update=False,rhs=rhs)
    # - imax/jmax 
    genBC(src_phybc_wave_imaxjmax,5,4,rhsname , locname_bc, setbc=[True,{'cornerchar':{'imaxjmax':['rhs']}}]  , update=False,rhs=rhs)

    # Extract RHS info:
    rhsinfo(rhs)

if __name__ == '__main__':
    main()
