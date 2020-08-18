

!***********************************************************
!                                                           
! Start building layers for BC : imax j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 1 4 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 1 4 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_nxp5m1m11m5p4k = q(nx+5-1-1,1-5+4,indvars(2))

d1_divV_dx_0_nxp5m1p11m5p4k = q(nx+5-1+1,1-5+4,indvars(2))

d1_divV_dx_0_nxp5m11m5p4k = -&
          0.5_wp*d1_divV_dx_0_nxp5m1m11m5p4k+&
          0.5_wp*d1_divV_dx_0_nxp5m1p11m5p4k

d1_divV_dx_0_nxp5m11m5p4k = d1_divV_dx_0_nxp5m11m5p4k*param_float(1)

d1_divV_dy_0_nxp5m11m5p4m2k = q(nx+5-1,1-5+4-2,indvars(3))

d1_divV_dy_0_nxp5m11m5p4m1k = q(nx+5-1,1-5+4-1,indvars(3))

d1_divV_dy_0_nxp5m11m5p4p1k = q(nx+5-1,1-5+4+1,indvars(3))

d1_divV_dy_0_nxp5m11m5p4p2k = q(nx+5-1,1-5+4+2,indvars(3))

d1_divV_dy_0_nxp5m11m5p4k = 0.08333333333333333_wp*d1_divV_dy_0_nxp5m11m5p4m2k-&
          0.666666666666667_wp*d1_divV_dy_0_nxp5m11m5p4m1k+&
          0.666666666666667_wp*d1_divV_dy_0_nxp5m11m5p4p1k-&
          0.08333333333333333_wp*d1_divV_dy_0_nxp5m11m5p4p2k

d1_divV_dy_0_nxp5m11m5p4k = d1_divV_dy_0_nxp5m11m5p4k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 4 None divV ******************
!                                                           
!***********************************************************


qst(nx+5-1,1-5+4,indvarsst(1)) =  d1_divV_dx_0_nxp5m11m5p4k+d1_divV_dy_0_nxp5m11m5p4k

