

!***********************************************************
!                                                           
! Start building layers for BC : i1 jmax None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 3 3 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 3 3 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_1m5p3m2nyp5m3k = q(1-5+3-2,ny+5-3,indvars(2))

d1_divV_dx_0_1m5p3m1nyp5m3k = q(1-5+3-1,ny+5-3,indvars(2))

d1_divV_dx_0_1m5p3p1nyp5m3k = q(1-5+3+1,ny+5-3,indvars(2))

d1_divV_dx_0_1m5p3p2nyp5m3k = q(1-5+3+2,ny+5-3,indvars(2))

d1_divV_dx_0_1m5p3nyp5m3k = 0.08333333333333333_wp*d1_divV_dx_0_1m5p3m2nyp5m3k-&
          0.666666666666667_wp*d1_divV_dx_0_1m5p3m1nyp5m3k+&
          0.666666666666667_wp*d1_divV_dx_0_1m5p3p1nyp5m3k-&
          0.08333333333333333_wp*d1_divV_dx_0_1m5p3p2nyp5m3k

d1_divV_dx_0_1m5p3nyp5m3k = d1_divV_dx_0_1m5p3nyp5m3k*param_float(1)

d1_divV_dy_0_1m5p3nyp5m3m2k = q(1-5+3,ny+5-3-2,indvars(3))

d1_divV_dy_0_1m5p3nyp5m3m1k = q(1-5+3,ny+5-3-1,indvars(3))

d1_divV_dy_0_1m5p3nyp5m3p1k = q(1-5+3,ny+5-3+1,indvars(3))

d1_divV_dy_0_1m5p3nyp5m3p2k = q(1-5+3,ny+5-3+2,indvars(3))

d1_divV_dy_0_1m5p3nyp5m3k = 0.08333333333333333_wp*d1_divV_dy_0_1m5p3nyp5m3m2k-&
          0.666666666666667_wp*d1_divV_dy_0_1m5p3nyp5m3m1k+&
          0.666666666666667_wp*d1_divV_dy_0_1m5p3nyp5m3p1k-&
          0.08333333333333333_wp*d1_divV_dy_0_1m5p3nyp5m3p2k

d1_divV_dy_0_1m5p3nyp5m3k = d1_divV_dy_0_1m5p3nyp5m3k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 3 3 None divV ******************
!                                                           
!***********************************************************


qst(1-5+3,ny+5-3,indvarsst(1)) =  d1_divV_dx_0_1m5p3nyp5m3k+d1_divV_dy_0_1m5p3nyp5m3k

