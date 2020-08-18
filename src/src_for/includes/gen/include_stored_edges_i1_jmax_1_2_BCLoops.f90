

!***********************************************************
!                                                           
! Start building layers for BC : i1 jmax None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 1 2 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 1 2 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_1m5p1m1nyp5m2k = q(1-5+1-1,ny+5-2,indvars(2))

d1_divV_dx_0_1m5p1p1nyp5m2k = q(1-5+1+1,ny+5-2,indvars(2))

d1_divV_dx_0_1m5p1nyp5m2k = -&
          0.5_wp*d1_divV_dx_0_1m5p1m1nyp5m2k+&
          0.5_wp*d1_divV_dx_0_1m5p1p1nyp5m2k

d1_divV_dx_0_1m5p1nyp5m2k = d1_divV_dx_0_1m5p1nyp5m2k*param_float(1)

d1_divV_dy_0_1m5p1nyp5m2m2k = q(1-5+1,ny+5-2-2,indvars(3))

d1_divV_dy_0_1m5p1nyp5m2m1k = q(1-5+1,ny+5-2-1,indvars(3))

d1_divV_dy_0_1m5p1nyp5m2p1k = q(1-5+1,ny+5-2+1,indvars(3))

d1_divV_dy_0_1m5p1nyp5m2p2k = q(1-5+1,ny+5-2+2,indvars(3))

d1_divV_dy_0_1m5p1nyp5m2k = 0.08333333333333333_wp*d1_divV_dy_0_1m5p1nyp5m2m2k-&
          0.666666666666667_wp*d1_divV_dy_0_1m5p1nyp5m2m1k+&
          0.666666666666667_wp*d1_divV_dy_0_1m5p1nyp5m2p1k-&
          0.08333333333333333_wp*d1_divV_dy_0_1m5p1nyp5m2p2k

d1_divV_dy_0_1m5p1nyp5m2k = d1_divV_dy_0_1m5p1nyp5m2k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 2 None divV ******************
!                                                           
!***********************************************************


qst(1-5+1,ny+5-2,indvarsst(1)) =  d1_divV_dx_0_1m5p1nyp5m2k+d1_divV_dy_0_1m5p1nyp5m2k

