

!***********************************************************
!                                                           
! Start building layers for BC : imax j1 None **************
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

d1_divV_dx_0_nxp5m3m21m5p3k = q(nx+5-3-2,1-5+3,indvars(2))

d1_divV_dx_0_nxp5m3m11m5p3k = q(nx+5-3-1,1-5+3,indvars(2))

d1_divV_dx_0_nxp5m3p11m5p3k = q(nx+5-3+1,1-5+3,indvars(2))

d1_divV_dx_0_nxp5m3p21m5p3k = q(nx+5-3+2,1-5+3,indvars(2))

d1_divV_dx_0_nxp5m31m5p3k = 0.08333333333333333_wp*d1_divV_dx_0_nxp5m3m21m5p3k-&
          0.666666666666667_wp*d1_divV_dx_0_nxp5m3m11m5p3k+&
          0.666666666666667_wp*d1_divV_dx_0_nxp5m3p11m5p3k-&
          0.08333333333333333_wp*d1_divV_dx_0_nxp5m3p21m5p3k

d1_divV_dx_0_nxp5m31m5p3k = d1_divV_dx_0_nxp5m31m5p3k*param_float(1)

d1_divV_dy_0_nxp5m31m5p3m2k = q(nx+5-3,1-5+3-2,indvars(3))

d1_divV_dy_0_nxp5m31m5p3m1k = q(nx+5-3,1-5+3-1,indvars(3))

d1_divV_dy_0_nxp5m31m5p3p1k = q(nx+5-3,1-5+3+1,indvars(3))

d1_divV_dy_0_nxp5m31m5p3p2k = q(nx+5-3,1-5+3+2,indvars(3))

d1_divV_dy_0_nxp5m31m5p3k = 0.08333333333333333_wp*d1_divV_dy_0_nxp5m31m5p3m2k-&
          0.666666666666667_wp*d1_divV_dy_0_nxp5m31m5p3m1k+&
          0.666666666666667_wp*d1_divV_dy_0_nxp5m31m5p3p1k-&
          0.08333333333333333_wp*d1_divV_dy_0_nxp5m31m5p3p2k

d1_divV_dy_0_nxp5m31m5p3k = d1_divV_dy_0_nxp5m31m5p3k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 3 3 None divV ******************
!                                                           
!***********************************************************


qst(nx+5-3,1-5+3,indvarsst(1)) =  d1_divV_dx_0_nxp5m31m5p3k+d1_divV_dy_0_nxp5m31m5p3k

