

!***********************************************************
!                                                           
! Start building layers for BC : imax jmax None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 4 0 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 4 0 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_nxp5m4m2nyp5p0k = q(nx+5-4-2,ny+5+0,indvars(2))

d1_divV_dx_0_nxp5m4m1nyp5p0k = q(nx+5-4-1,ny+5+0,indvars(2))

d1_divV_dx_0_nxp5m4p1nyp5p0k = q(nx+5-4+1,ny+5+0,indvars(2))

d1_divV_dx_0_nxp5m4p2nyp5p0k = q(nx+5-4+2,ny+5+0,indvars(2))

d1_divV_dx_0_nxp5m4nyp5p0k = 0.08333333333333333_wp*d1_divV_dx_0_nxp5m4m2nyp5p0k-&
          0.666666666666667_wp*d1_divV_dx_0_nxp5m4m1nyp5p0k+&
          0.666666666666667_wp*d1_divV_dx_0_nxp5m4p1nyp5p0k-&
          0.08333333333333333_wp*d1_divV_dx_0_nxp5m4p2nyp5p0k

d1_divV_dx_0_nxp5m4nyp5p0k = d1_divV_dx_0_nxp5m4nyp5p0k*param_float(1)

d1_divV_dy_0_nxp5m4nyp5p0p0k = q(nx+5-4,ny+5+0+0,indvars(3))

d1_divV_dy_0_nxp5m4nyp5p0m1k = q(nx+5-4,ny+5+0-1,indvars(3))

d1_divV_dy_0_nxp5m4nyp5p0m2k = q(nx+5-4,ny+5+0-2,indvars(3))

d1_divV_dy_0_nxp5m4nyp5p0k = 1.5_wp*d1_divV_dy_0_nxp5m4nyp5p0p0k-&
          2.0_wp*d1_divV_dy_0_nxp5m4nyp5p0m1k+&
          0.5_wp*d1_divV_dy_0_nxp5m4nyp5p0m2k

d1_divV_dy_0_nxp5m4nyp5p0k = d1_divV_dy_0_nxp5m4nyp5p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 4 0 None divV ******************
!                                                           
!***********************************************************


qst(nx+5-4,ny+5+0,indvarsst(1)) =  d1_divV_dx_0_nxp5m4nyp5p0k+d1_divV_dy_0_nxp5m4nyp5p0k

