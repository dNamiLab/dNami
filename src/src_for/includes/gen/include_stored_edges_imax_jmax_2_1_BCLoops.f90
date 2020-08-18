

!***********************************************************
!                                                           
! Start building layers for BC : imax jmax None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 2 1 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 2 1 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_nxp5m2m2nyp5m1k = q(nx+5-2-2,ny+5-1,indvars(2))

d1_divV_dx_0_nxp5m2m1nyp5m1k = q(nx+5-2-1,ny+5-1,indvars(2))

d1_divV_dx_0_nxp5m2p1nyp5m1k = q(nx+5-2+1,ny+5-1,indvars(2))

d1_divV_dx_0_nxp5m2p2nyp5m1k = q(nx+5-2+2,ny+5-1,indvars(2))

d1_divV_dx_0_nxp5m2nyp5m1k = 0.08333333333333333_wp*d1_divV_dx_0_nxp5m2m2nyp5m1k-&
          0.666666666666667_wp*d1_divV_dx_0_nxp5m2m1nyp5m1k+&
          0.666666666666667_wp*d1_divV_dx_0_nxp5m2p1nyp5m1k-&
          0.08333333333333333_wp*d1_divV_dx_0_nxp5m2p2nyp5m1k

d1_divV_dx_0_nxp5m2nyp5m1k = d1_divV_dx_0_nxp5m2nyp5m1k*param_float(1)

d1_divV_dy_0_nxp5m2nyp5m1m1k = q(nx+5-2,ny+5-1-1,indvars(3))

d1_divV_dy_0_nxp5m2nyp5m1p1k = q(nx+5-2,ny+5-1+1,indvars(3))

d1_divV_dy_0_nxp5m2nyp5m1k = -&
          0.5_wp*d1_divV_dy_0_nxp5m2nyp5m1m1k+&
          0.5_wp*d1_divV_dy_0_nxp5m2nyp5m1p1k

d1_divV_dy_0_nxp5m2nyp5m1k = d1_divV_dy_0_nxp5m2nyp5m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 2 1 None divV ******************
!                                                           
!***********************************************************


qst(nx+5-2,ny+5-1,indvarsst(1)) =  d1_divV_dx_0_nxp5m2nyp5m1k+d1_divV_dy_0_nxp5m2nyp5m1k

