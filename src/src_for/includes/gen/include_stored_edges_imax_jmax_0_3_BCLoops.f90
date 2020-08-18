

!***********************************************************
!                                                           
! Start building layers for BC : imax jmax None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 3 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_nxp5p0p0nyp5m3k = q(nx+5+0+0,ny+5-3,indvars(2))

d1_divV_dx_0_nxp5p0m1nyp5m3k = q(nx+5+0-1,ny+5-3,indvars(2))

d1_divV_dx_0_nxp5p0m2nyp5m3k = q(nx+5+0-2,ny+5-3,indvars(2))

d1_divV_dx_0_nxp5p0nyp5m3k = 1.5_wp*d1_divV_dx_0_nxp5p0p0nyp5m3k-&
          2.0_wp*d1_divV_dx_0_nxp5p0m1nyp5m3k+&
          0.5_wp*d1_divV_dx_0_nxp5p0m2nyp5m3k

d1_divV_dx_0_nxp5p0nyp5m3k = d1_divV_dx_0_nxp5p0nyp5m3k*param_float(1)

d1_divV_dy_0_nxp5p0nyp5m3m2k = q(nx+5+0,ny+5-3-2,indvars(3))

d1_divV_dy_0_nxp5p0nyp5m3m1k = q(nx+5+0,ny+5-3-1,indvars(3))

d1_divV_dy_0_nxp5p0nyp5m3p1k = q(nx+5+0,ny+5-3+1,indvars(3))

d1_divV_dy_0_nxp5p0nyp5m3p2k = q(nx+5+0,ny+5-3+2,indvars(3))

d1_divV_dy_0_nxp5p0nyp5m3k = 0.08333333333333333_wp*d1_divV_dy_0_nxp5p0nyp5m3m2k-&
          0.666666666666667_wp*d1_divV_dy_0_nxp5p0nyp5m3m1k+&
          0.666666666666667_wp*d1_divV_dy_0_nxp5p0nyp5m3p1k-&
          0.08333333333333333_wp*d1_divV_dy_0_nxp5p0nyp5m3p2k

d1_divV_dy_0_nxp5p0nyp5m3k = d1_divV_dy_0_nxp5p0nyp5m3k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None divV ******************
!                                                           
!***********************************************************


qst(nx+5+0,ny+5-3,indvarsst(1)) =  d1_divV_dx_0_nxp5p0nyp5m3k+d1_divV_dy_0_nxp5p0nyp5m3k

