

!***********************************************************
!                                                           
! Start building layers for BC : i1 j1 None ****************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 4 4 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 4 4 None divV *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_1m5p4m21m5p4k = q(1-5+4-2,1-5+4,indvars(2))

d1_divV_dx_0_1m5p4m11m5p4k = q(1-5+4-1,1-5+4,indvars(2))

d1_divV_dx_0_1m5p4p11m5p4k = q(1-5+4+1,1-5+4,indvars(2))

d1_divV_dx_0_1m5p4p21m5p4k = q(1-5+4+2,1-5+4,indvars(2))

d1_divV_dx_0_1m5p41m5p4k = 0.08333333333333333_wp*d1_divV_dx_0_1m5p4m21m5p4k-&
          0.666666666666667_wp*d1_divV_dx_0_1m5p4m11m5p4k+&
          0.666666666666667_wp*d1_divV_dx_0_1m5p4p11m5p4k-&
          0.08333333333333333_wp*d1_divV_dx_0_1m5p4p21m5p4k

d1_divV_dx_0_1m5p41m5p4k = d1_divV_dx_0_1m5p41m5p4k*param_float(1)

d1_divV_dy_0_1m5p41m5p4m2k = q(1-5+4,1-5+4-2,indvars(3))

d1_divV_dy_0_1m5p41m5p4m1k = q(1-5+4,1-5+4-1,indvars(3))

d1_divV_dy_0_1m5p41m5p4p1k = q(1-5+4,1-5+4+1,indvars(3))

d1_divV_dy_0_1m5p41m5p4p2k = q(1-5+4,1-5+4+2,indvars(3))

d1_divV_dy_0_1m5p41m5p4k = 0.08333333333333333_wp*d1_divV_dy_0_1m5p41m5p4m2k-&
          0.666666666666667_wp*d1_divV_dy_0_1m5p41m5p4m1k+&
          0.666666666666667_wp*d1_divV_dy_0_1m5p41m5p4p1k-&
          0.08333333333333333_wp*d1_divV_dy_0_1m5p41m5p4p2k

d1_divV_dy_0_1m5p41m5p4k = d1_divV_dy_0_1m5p41m5p4k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 4 4 None divV ******************
!                                                           
!***********************************************************


qst(1-5+4,1-5+4,indvarsst(1)) =  d1_divV_dx_0_1m5p41m5p4k+d1_divV_dy_0_1m5p41m5p4k

