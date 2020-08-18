

!***********************************************************
!                                                           
! Start building layers for BC : imax j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 1 0 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 1 0 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_nxp5m1m11m5p0k = q(nx+5-1-1,1-5+0,indvars(2))

d1_divV_dx_0_nxp5m1p11m5p0k = q(nx+5-1+1,1-5+0,indvars(2))

d1_divV_dx_0_nxp5m11m5p0k = -&
          0.5_wp*d1_divV_dx_0_nxp5m1m11m5p0k+&
          0.5_wp*d1_divV_dx_0_nxp5m1p11m5p0k

d1_divV_dx_0_nxp5m11m5p0k = d1_divV_dx_0_nxp5m11m5p0k*param_float(1)

d1_divV_dy_0_nxp5m11m5p0p0k = q(nx+5-1,1-5+0+0,indvars(3))

d1_divV_dy_0_nxp5m11m5p0p1k = q(nx+5-1,1-5+0+1,indvars(3))

d1_divV_dy_0_nxp5m11m5p0p2k = q(nx+5-1,1-5+0+2,indvars(3))

d1_divV_dy_0_nxp5m11m5p0k = -&
          1.5_wp*d1_divV_dy_0_nxp5m11m5p0p0k+&
          2.0_wp*d1_divV_dy_0_nxp5m11m5p0p1k-&
          0.5_wp*d1_divV_dy_0_nxp5m11m5p0p2k

d1_divV_dy_0_nxp5m11m5p0k = d1_divV_dy_0_nxp5m11m5p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 0 None divV ******************
!                                                           
!***********************************************************


qst(nx+5-1,1-5+0,indvarsst(1)) =  d1_divV_dx_0_nxp5m11m5p0k+d1_divV_dy_0_nxp5m11m5p0k

