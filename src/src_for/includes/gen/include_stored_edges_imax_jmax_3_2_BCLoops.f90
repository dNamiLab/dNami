

!***********************************************************
!                                                           
! Start building layers for BC : imax jmax None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 3 2 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 3 2 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_nxp5m3m2nyp5m2k = q(nx+5-3-2,ny+5-2,indvars(2))

d1_divV_dx_0_nxp5m3m1nyp5m2k = q(nx+5-3-1,ny+5-2,indvars(2))

d1_divV_dx_0_nxp5m3p1nyp5m2k = q(nx+5-3+1,ny+5-2,indvars(2))

d1_divV_dx_0_nxp5m3p2nyp5m2k = q(nx+5-3+2,ny+5-2,indvars(2))

d1_divV_dx_0_nxp5m3nyp5m2k = 0.08333333333333333_wp*d1_divV_dx_0_nxp5m3m2nyp5m2k-&
          0.666666666666667_wp*d1_divV_dx_0_nxp5m3m1nyp5m2k+&
          0.666666666666667_wp*d1_divV_dx_0_nxp5m3p1nyp5m2k-&
          0.08333333333333333_wp*d1_divV_dx_0_nxp5m3p2nyp5m2k

d1_divV_dx_0_nxp5m3nyp5m2k = d1_divV_dx_0_nxp5m3nyp5m2k*param_float(1)

d1_divV_dy_0_nxp5m3nyp5m2m2k = q(nx+5-3,ny+5-2-2,indvars(3))

d1_divV_dy_0_nxp5m3nyp5m2m1k = q(nx+5-3,ny+5-2-1,indvars(3))

d1_divV_dy_0_nxp5m3nyp5m2p1k = q(nx+5-3,ny+5-2+1,indvars(3))

d1_divV_dy_0_nxp5m3nyp5m2p2k = q(nx+5-3,ny+5-2+2,indvars(3))

d1_divV_dy_0_nxp5m3nyp5m2k = 0.08333333333333333_wp*d1_divV_dy_0_nxp5m3nyp5m2m2k-&
          0.666666666666667_wp*d1_divV_dy_0_nxp5m3nyp5m2m1k+&
          0.666666666666667_wp*d1_divV_dy_0_nxp5m3nyp5m2p1k-&
          0.08333333333333333_wp*d1_divV_dy_0_nxp5m3nyp5m2p2k

d1_divV_dy_0_nxp5m3nyp5m2k = d1_divV_dy_0_nxp5m3nyp5m2k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 3 2 None divV ******************
!                                                           
!***********************************************************


qst(nx+5-3,ny+5-2,indvarsst(1)) =  d1_divV_dx_0_nxp5m3nyp5m2k+d1_divV_dy_0_nxp5m3nyp5m2k

