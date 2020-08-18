

!***********************************************************
!                                                           
! Start building layers for BC : i1 j1 None ****************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 1 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_1m5p0p01m5p1k = q(1-5+0+0,1-5+1,indvars(2))

d1_divV_dx_0_1m5p0p11m5p1k = q(1-5+0+1,1-5+1,indvars(2))

d1_divV_dx_0_1m5p0p21m5p1k = q(1-5+0+2,1-5+1,indvars(2))

d1_divV_dx_0_1m5p01m5p1k = -&
          1.5_wp*d1_divV_dx_0_1m5p0p01m5p1k+&
          2.0_wp*d1_divV_dx_0_1m5p0p11m5p1k-&
          0.5_wp*d1_divV_dx_0_1m5p0p21m5p1k

d1_divV_dx_0_1m5p01m5p1k = d1_divV_dx_0_1m5p01m5p1k*param_float(1)

d1_divV_dy_0_1m5p01m5p1m1k = q(1-5+0,1-5+1-1,indvars(3))

d1_divV_dy_0_1m5p01m5p1p1k = q(1-5+0,1-5+1+1,indvars(3))

d1_divV_dy_0_1m5p01m5p1k = -&
          0.5_wp*d1_divV_dy_0_1m5p01m5p1m1k+&
          0.5_wp*d1_divV_dy_0_1m5p01m5p1p1k

d1_divV_dy_0_1m5p01m5p1k = d1_divV_dy_0_1m5p01m5p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None divV ******************
!                                                           
!***********************************************************


qst(1-5+0,1-5+1,indvarsst(1)) =  d1_divV_dx_0_1m5p01m5p1k+d1_divV_dy_0_1m5p01m5p1k

