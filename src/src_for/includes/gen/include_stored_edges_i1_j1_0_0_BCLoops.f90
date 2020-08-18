

!***********************************************************
!                                                           
! Start building layers for BC : i1 j1 None ****************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 0 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_1m5p0p01m5p0k = q(1-5+0+0,1-5+0,indvars(2))

d1_divV_dx_0_1m5p0p11m5p0k = q(1-5+0+1,1-5+0,indvars(2))

d1_divV_dx_0_1m5p0p21m5p0k = q(1-5+0+2,1-5+0,indvars(2))

d1_divV_dx_0_1m5p01m5p0k = -&
          1.5_wp*d1_divV_dx_0_1m5p0p01m5p0k+&
          2.0_wp*d1_divV_dx_0_1m5p0p11m5p0k-&
          0.5_wp*d1_divV_dx_0_1m5p0p21m5p0k

d1_divV_dx_0_1m5p01m5p0k = d1_divV_dx_0_1m5p01m5p0k*param_float(1)

d1_divV_dy_0_1m5p01m5p0p0k = q(1-5+0,1-5+0+0,indvars(3))

d1_divV_dy_0_1m5p01m5p0p1k = q(1-5+0,1-5+0+1,indvars(3))

d1_divV_dy_0_1m5p01m5p0p2k = q(1-5+0,1-5+0+2,indvars(3))

d1_divV_dy_0_1m5p01m5p0k = -&
          1.5_wp*d1_divV_dy_0_1m5p01m5p0p0k+&
          2.0_wp*d1_divV_dy_0_1m5p01m5p0p1k-&
          0.5_wp*d1_divV_dy_0_1m5p01m5p0p2k

d1_divV_dy_0_1m5p01m5p0k = d1_divV_dy_0_1m5p01m5p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None divV ******************
!                                                           
!***********************************************************


qst(1-5+0,1-5+0,indvarsst(1)) =  d1_divV_dx_0_1m5p01m5p0k+d1_divV_dy_0_1m5p01m5p0k

