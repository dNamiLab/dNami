

!***********************************************************
!                                                           
! Start building layers for BC : i1 jmax None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 4 4 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 4 4 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_1m5p4m2nyp5m4k = q(1-5+4-2,ny+5-4,indvars(2))

d1_divV_dx_0_1m5p4m1nyp5m4k = q(1-5+4-1,ny+5-4,indvars(2))

d1_divV_dx_0_1m5p4p1nyp5m4k = q(1-5+4+1,ny+5-4,indvars(2))

d1_divV_dx_0_1m5p4p2nyp5m4k = q(1-5+4+2,ny+5-4,indvars(2))

d1_divV_dx_0_1m5p4nyp5m4k = 0.08333333333333333_wp*d1_divV_dx_0_1m5p4m2nyp5m4k-&
          0.666666666666667_wp*d1_divV_dx_0_1m5p4m1nyp5m4k+&
          0.666666666666667_wp*d1_divV_dx_0_1m5p4p1nyp5m4k-&
          0.08333333333333333_wp*d1_divV_dx_0_1m5p4p2nyp5m4k

d1_divV_dx_0_1m5p4nyp5m4k = d1_divV_dx_0_1m5p4nyp5m4k*param_float(1)

d1_divV_dy_0_1m5p4nyp5m4m2k = q(1-5+4,ny+5-4-2,indvars(3))

d1_divV_dy_0_1m5p4nyp5m4m1k = q(1-5+4,ny+5-4-1,indvars(3))

d1_divV_dy_0_1m5p4nyp5m4p1k = q(1-5+4,ny+5-4+1,indvars(3))

d1_divV_dy_0_1m5p4nyp5m4p2k = q(1-5+4,ny+5-4+2,indvars(3))

d1_divV_dy_0_1m5p4nyp5m4k = 0.08333333333333333_wp*d1_divV_dy_0_1m5p4nyp5m4m2k-&
          0.666666666666667_wp*d1_divV_dy_0_1m5p4nyp5m4m1k+&
          0.666666666666667_wp*d1_divV_dy_0_1m5p4nyp5m4p1k-&
          0.08333333333333333_wp*d1_divV_dy_0_1m5p4nyp5m4p2k

d1_divV_dy_0_1m5p4nyp5m4k = d1_divV_dy_0_1m5p4nyp5m4k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 4 4 None divV ******************
!                                                           
!***********************************************************


qst(1-5+4,ny+5-4,indvarsst(1)) =  d1_divV_dx_0_1m5p4nyp5m4k+d1_divV_dy_0_1m5p4nyp5m4k

