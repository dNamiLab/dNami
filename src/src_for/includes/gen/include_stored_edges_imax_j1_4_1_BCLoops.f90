

!***********************************************************
!                                                           
! Start building layers for BC : imax j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 4 1 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 4 1 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_nxp5m4m21m5p1k = q(nx+5-4-2,1-5+1,indvars(2))

d1_divV_dx_0_nxp5m4m11m5p1k = q(nx+5-4-1,1-5+1,indvars(2))

d1_divV_dx_0_nxp5m4p11m5p1k = q(nx+5-4+1,1-5+1,indvars(2))

d1_divV_dx_0_nxp5m4p21m5p1k = q(nx+5-4+2,1-5+1,indvars(2))

d1_divV_dx_0_nxp5m41m5p1k = 0.08333333333333333_wp*d1_divV_dx_0_nxp5m4m21m5p1k-&
          0.666666666666667_wp*d1_divV_dx_0_nxp5m4m11m5p1k+&
          0.666666666666667_wp*d1_divV_dx_0_nxp5m4p11m5p1k-&
          0.08333333333333333_wp*d1_divV_dx_0_nxp5m4p21m5p1k

d1_divV_dx_0_nxp5m41m5p1k = d1_divV_dx_0_nxp5m41m5p1k*param_float(1)

d1_divV_dy_0_nxp5m41m5p1m1k = q(nx+5-4,1-5+1-1,indvars(3))

d1_divV_dy_0_nxp5m41m5p1p1k = q(nx+5-4,1-5+1+1,indvars(3))

d1_divV_dy_0_nxp5m41m5p1k = -&
          0.5_wp*d1_divV_dy_0_nxp5m41m5p1m1k+&
          0.5_wp*d1_divV_dy_0_nxp5m41m5p1p1k

d1_divV_dy_0_nxp5m41m5p1k = d1_divV_dy_0_nxp5m41m5p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 4 1 None divV ******************
!                                                           
!***********************************************************


qst(nx+5-4,1-5+1,indvarsst(1)) =  d1_divV_dx_0_nxp5m41m5p1k+d1_divV_dy_0_nxp5m41m5p1k

