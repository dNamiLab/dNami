

!***********************************************************
!                                                           
! Start building layers for BC : imax j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 2 2 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 2 2 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_nxp5m2m21m5p2k = q(nx+5-2-2,1-5+2,indvars(2))

d1_divV_dx_0_nxp5m2m11m5p2k = q(nx+5-2-1,1-5+2,indvars(2))

d1_divV_dx_0_nxp5m2p11m5p2k = q(nx+5-2+1,1-5+2,indvars(2))

d1_divV_dx_0_nxp5m2p21m5p2k = q(nx+5-2+2,1-5+2,indvars(2))

d1_divV_dx_0_nxp5m21m5p2k = 0.08333333333333333_wp*d1_divV_dx_0_nxp5m2m21m5p2k-&
          0.666666666666667_wp*d1_divV_dx_0_nxp5m2m11m5p2k+&
          0.666666666666667_wp*d1_divV_dx_0_nxp5m2p11m5p2k-&
          0.08333333333333333_wp*d1_divV_dx_0_nxp5m2p21m5p2k

d1_divV_dx_0_nxp5m21m5p2k = d1_divV_dx_0_nxp5m21m5p2k*param_float(1)

d1_divV_dy_0_nxp5m21m5p2m2k = q(nx+5-2,1-5+2-2,indvars(3))

d1_divV_dy_0_nxp5m21m5p2m1k = q(nx+5-2,1-5+2-1,indvars(3))

d1_divV_dy_0_nxp5m21m5p2p1k = q(nx+5-2,1-5+2+1,indvars(3))

d1_divV_dy_0_nxp5m21m5p2p2k = q(nx+5-2,1-5+2+2,indvars(3))

d1_divV_dy_0_nxp5m21m5p2k = 0.08333333333333333_wp*d1_divV_dy_0_nxp5m21m5p2m2k-&
          0.666666666666667_wp*d1_divV_dy_0_nxp5m21m5p2m1k+&
          0.666666666666667_wp*d1_divV_dy_0_nxp5m21m5p2p1k-&
          0.08333333333333333_wp*d1_divV_dy_0_nxp5m21m5p2p2k

d1_divV_dy_0_nxp5m21m5p2k = d1_divV_dy_0_nxp5m21m5p2k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 2 2 None divV ******************
!                                                           
!***********************************************************


qst(nx+5-2,1-5+2,indvarsst(1)) =  d1_divV_dx_0_nxp5m21m5p2k+d1_divV_dy_0_nxp5m21m5p2k

