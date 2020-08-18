

!***********************************************************
!                                                           
! Start building layers for BC : imax jmax None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 4 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 4 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_nxp5p0p0nyp5m4k = q(nx+5+0+0,ny+5-4,indvars(2))

d1_divV_dx_0_nxp5p0m1nyp5m4k = q(nx+5+0-1,ny+5-4,indvars(2))

d1_divV_dx_0_nxp5p0m2nyp5m4k = q(nx+5+0-2,ny+5-4,indvars(2))

d1_divV_dx_0_nxp5p0nyp5m4k = 1.5_wp*d1_divV_dx_0_nxp5p0p0nyp5m4k-&
          2.0_wp*d1_divV_dx_0_nxp5p0m1nyp5m4k+&
          0.5_wp*d1_divV_dx_0_nxp5p0m2nyp5m4k

d1_divV_dx_0_nxp5p0nyp5m4k = d1_divV_dx_0_nxp5p0nyp5m4k*param_float(1)

d1_divV_dy_0_nxp5p0nyp5m4m2k = q(nx+5+0,ny+5-4-2,indvars(3))

d1_divV_dy_0_nxp5p0nyp5m4m1k = q(nx+5+0,ny+5-4-1,indvars(3))

d1_divV_dy_0_nxp5p0nyp5m4p1k = q(nx+5+0,ny+5-4+1,indvars(3))

d1_divV_dy_0_nxp5p0nyp5m4p2k = q(nx+5+0,ny+5-4+2,indvars(3))

d1_divV_dy_0_nxp5p0nyp5m4k = 0.08333333333333333_wp*d1_divV_dy_0_nxp5p0nyp5m4m2k-&
          0.666666666666667_wp*d1_divV_dy_0_nxp5p0nyp5m4m1k+&
          0.666666666666667_wp*d1_divV_dy_0_nxp5p0nyp5m4p1k-&
          0.08333333333333333_wp*d1_divV_dy_0_nxp5p0nyp5m4p2k

d1_divV_dy_0_nxp5p0nyp5m4k = d1_divV_dy_0_nxp5p0nyp5m4k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 4 None divV ******************
!                                                           
!***********************************************************


qst(nx+5+0,ny+5-4,indvarsst(1)) =  d1_divV_dx_0_nxp5p0nyp5m4k+d1_divV_dy_0_nxp5p0nyp5m4k

