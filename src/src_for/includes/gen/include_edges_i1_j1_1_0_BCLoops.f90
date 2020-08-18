

!***********************************************************
!                                                           
! Start building layers for BC : i1 j1 None ****************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 1 0 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 1 0 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u]_1x+[rho*v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluRx_dx_0_1m5p1m11m5p0k = q(1-5+1-1,1-5+0,indvars(1))*q(1-5+1-1,1-5+0,indvars(2))

d1_FluRx_dx_0_1m5p1p11m5p0k = q(1-5+1+1,1-5+0,indvars(1))*q(1-5+1+1,1-5+0,indvars(2))

d1_FluRx_dx_0_1m5p11m5p0k = -&
          0.5_wp*d1_FluRx_dx_0_1m5p1m11m5p0k+&
          0.5_wp*d1_FluRx_dx_0_1m5p1p11m5p0k

d1_FluRx_dx_0_1m5p11m5p0k = d1_FluRx_dx_0_1m5p11m5p0k*param_float(1)

d1_FluRx_dy_0_1m5p11m5p0p0k = q(1-5+1,1-5+0+0,indvars(1))*q(1-5+1,1-5+0+0,indvars(3))

d1_FluRx_dy_0_1m5p11m5p0p1k = q(1-5+1,1-5+0+1,indvars(1))*q(1-5+1,1-5+0+1,indvars(3))

d1_FluRx_dy_0_1m5p11m5p0p2k = q(1-5+1,1-5+0+2,indvars(1))*q(1-5+1,1-5+0+2,indvars(3))

d1_FluRx_dy_0_1m5p11m5p0k = -&
          1.5_wp*d1_FluRx_dy_0_1m5p11m5p0p0k+&
          2.0_wp*d1_FluRx_dy_0_1m5p11m5p0p1k-&
          0.5_wp*d1_FluRx_dy_0_1m5p11m5p0p2k

d1_FluRx_dy_0_1m5p11m5p0k = d1_FluRx_dy_0_1m5p11m5p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 0 None d(rho)/dt *************
!                                                           
!***********************************************************


rhs(1-5+1,1-5+0,indvars(1)) =   -  ( d1_FluRx_dx_0_1m5p11m5p0k+d1_FluRx_dy_0_1m5p11m5p0k ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 0 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u*u+0.375_wp*p]_1x+[rho*v*u]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluXx_dx_0_1m5p1m11m5p0k = q(1-5+1-1,1-5+0,indvars(1))*q(1-5+1-1,1-5+0,indvars(2))*q(1-5+1-1,1-5+0,indvars(2))+&
                    0.375_wp*(8.0_wp*q(1-5+1-1,1-5+0,indvars(1))/(3.0_wp-&
                    q(1-5+1-1,1-5+0,indvars(1)))*param_float(1 + 5)*(((q(1-5+1-1,1-5+0,indvars(4))-&
                    0.5_wp*(q(1-5+1-1,1-5+0,indvars(2))*q(1-5+1-1,1-5+0,indvars(2))+&
                    q(1-5+1-1,1-5+0,indvars(3))*q(1-5+1-1,1-5+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+1-1,1-5+0,indvars(1)))-&
                    3.0_wp*q(1-5+1-1,1-5+0,indvars(1))*q(1-5+1-1,1-5+0,indvars(1)))

d1_FluXx_dx_0_1m5p1p11m5p0k = q(1-5+1+1,1-5+0,indvars(1))*q(1-5+1+1,1-5+0,indvars(2))*q(1-5+1+1,1-5+0,indvars(2))+&
                    0.375_wp*(8.0_wp*q(1-5+1+1,1-5+0,indvars(1))/(3.0_wp-&
                    q(1-5+1+1,1-5+0,indvars(1)))*param_float(1 + 5)*(((q(1-5+1+1,1-5+0,indvars(4))-&
                    0.5_wp*(q(1-5+1+1,1-5+0,indvars(2))*q(1-5+1+1,1-5+0,indvars(2))+&
                    q(1-5+1+1,1-5+0,indvars(3))*q(1-5+1+1,1-5+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+1+1,1-5+0,indvars(1)))-&
                    3.0_wp*q(1-5+1+1,1-5+0,indvars(1))*q(1-5+1+1,1-5+0,indvars(1)))

d1_FluXx_dx_0_1m5p11m5p0k = -&
          0.5_wp*d1_FluXx_dx_0_1m5p1m11m5p0k+&
          0.5_wp*d1_FluXx_dx_0_1m5p1p11m5p0k

d1_FluXx_dx_0_1m5p11m5p0k = d1_FluXx_dx_0_1m5p11m5p0k*param_float(1)

d1_FluXx_dy_0_1m5p11m5p0p0k = q(1-5+1,1-5+0+0,indvars(1))*q(1-5+1,1-5+0+0,indvars(3))*q(1-5+1,1-5+0+0,indvars(2))

d1_FluXx_dy_0_1m5p11m5p0p1k = q(1-5+1,1-5+0+1,indvars(1))*q(1-5+1,1-5+0+1,indvars(3))*q(1-5+1,1-5+0+1,indvars(2))

d1_FluXx_dy_0_1m5p11m5p0p2k = q(1-5+1,1-5+0+2,indvars(1))*q(1-5+1,1-5+0+2,indvars(3))*q(1-5+1,1-5+0+2,indvars(2))

d1_FluXx_dy_0_1m5p11m5p0k = -&
          1.5_wp*d1_FluXx_dy_0_1m5p11m5p0p0k+&
          2.0_wp*d1_FluXx_dy_0_1m5p11m5p0p1k-&
          0.5_wp*d1_FluXx_dy_0_1m5p11m5p0p2k

d1_FluXx_dy_0_1m5p11m5p0k = d1_FluXx_dy_0_1m5p11m5p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 0 None d(rho u)/dt ***********
!                                                           
!***********************************************************


rhs(1-5+1,1-5+0,indvars(2)) =   -  ( d1_FluXx_dx_0_1m5p11m5p0k+d1_FluXx_dy_0_1m5p11m5p0k ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 0 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u*v]_1x+[rho*v*v+0.375_wp*p]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluYx_dx_0_1m5p1m11m5p0k = q(1-5+1-1,1-5+0,indvars(1))*q(1-5+1-1,1-5+0,indvars(2))*q(1-5+1-1,1-5+0,indvars(3))

d1_FluYx_dx_0_1m5p1p11m5p0k = q(1-5+1+1,1-5+0,indvars(1))*q(1-5+1+1,1-5+0,indvars(2))*q(1-5+1+1,1-5+0,indvars(3))

d1_FluYx_dx_0_1m5p11m5p0k = -&
          0.5_wp*d1_FluYx_dx_0_1m5p1m11m5p0k+&
          0.5_wp*d1_FluYx_dx_0_1m5p1p11m5p0k

d1_FluYx_dx_0_1m5p11m5p0k = d1_FluYx_dx_0_1m5p11m5p0k*param_float(1)

d1_FluYx_dy_0_1m5p11m5p0p0k = q(1-5+1,1-5+0+0,indvars(1))*q(1-5+1,1-5+0+0,indvars(3))*q(1-5+1,1-5+0+0,indvars(3))+&
                    0.375_wp*(8.0_wp*q(1-5+1,1-5+0+0,indvars(1))/(3.0_wp-&
                    q(1-5+1,1-5+0+0,indvars(1)))*param_float(1 + 5)*(((q(1-5+1,1-5+0+0,indvars(4))-&
                    0.5_wp*(q(1-5+1,1-5+0+0,indvars(2))*q(1-5+1,1-5+0+0,indvars(2))+&
                    q(1-5+1,1-5+0+0,indvars(3))*q(1-5+1,1-5+0+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+1,1-5+0+0,indvars(1)))-&
                    3.0_wp*q(1-5+1,1-5+0+0,indvars(1))*q(1-5+1,1-5+0+0,indvars(1)))

d1_FluYx_dy_0_1m5p11m5p0p1k = q(1-5+1,1-5+0+1,indvars(1))*q(1-5+1,1-5+0+1,indvars(3))*q(1-5+1,1-5+0+1,indvars(3))+&
                    0.375_wp*(8.0_wp*q(1-5+1,1-5+0+1,indvars(1))/(3.0_wp-&
                    q(1-5+1,1-5+0+1,indvars(1)))*param_float(1 + 5)*(((q(1-5+1,1-5+0+1,indvars(4))-&
                    0.5_wp*(q(1-5+1,1-5+0+1,indvars(2))*q(1-5+1,1-5+0+1,indvars(2))+&
                    q(1-5+1,1-5+0+1,indvars(3))*q(1-5+1,1-5+0+1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+1,1-5+0+1,indvars(1)))-&
                    3.0_wp*q(1-5+1,1-5+0+1,indvars(1))*q(1-5+1,1-5+0+1,indvars(1)))

d1_FluYx_dy_0_1m5p11m5p0p2k = q(1-5+1,1-5+0+2,indvars(1))*q(1-5+1,1-5+0+2,indvars(3))*q(1-5+1,1-5+0+2,indvars(3))+&
                    0.375_wp*(8.0_wp*q(1-5+1,1-5+0+2,indvars(1))/(3.0_wp-&
                    q(1-5+1,1-5+0+2,indvars(1)))*param_float(1 + 5)*(((q(1-5+1,1-5+0+2,indvars(4))-&
                    0.5_wp*(q(1-5+1,1-5+0+2,indvars(2))*q(1-5+1,1-5+0+2,indvars(2))+&
                    q(1-5+1,1-5+0+2,indvars(3))*q(1-5+1,1-5+0+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+1,1-5+0+2,indvars(1)))-&
                    3.0_wp*q(1-5+1,1-5+0+2,indvars(1))*q(1-5+1,1-5+0+2,indvars(1)))

d1_FluYx_dy_0_1m5p11m5p0k = -&
          1.5_wp*d1_FluYx_dy_0_1m5p11m5p0p0k+&
          2.0_wp*d1_FluYx_dy_0_1m5p11m5p0p1k-&
          0.5_wp*d1_FluYx_dy_0_1m5p11m5p0p2k

d1_FluYx_dy_0_1m5p11m5p0k = d1_FluYx_dy_0_1m5p11m5p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 0 None d(rho v)/dt ***********
!                                                           
!***********************************************************


rhs(1-5+1,1-5+0,indvars(3)) =   -  ( d1_FluYx_dx_0_1m5p11m5p0k+d1_FluYx_dy_0_1m5p11m5p0k ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 0 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u*(rho*et+0.375_wp*p)]_1x+[v*(rho*et+0.375_wp*p)]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluEx_dx_0_1m5p1m11m5p0k = q(1-5+1-1,1-5+0,indvars(2))*(q(1-5+1-1,1-5+0,indvars(1))*q(1-5+1-1,1-5+0,indvars(4))+&
                    0.375_wp*(8.0_wp*q(1-5+1-1,1-5+0,indvars(1))/(3.0_wp-&
                    q(1-5+1-1,1-5+0,indvars(1)))*param_float(1 + 5)*(((q(1-5+1-1,1-5+0,indvars(4))-&
                    0.5_wp*(q(1-5+1-1,1-5+0,indvars(2))*q(1-5+1-1,1-5+0,indvars(2))+&
                    q(1-5+1-1,1-5+0,indvars(3))*q(1-5+1-1,1-5+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+1-1,1-5+0,indvars(1)))-&
                    3.0_wp*q(1-5+1-1,1-5+0,indvars(1))*q(1-5+1-1,1-5+0,indvars(1))))

d1_FluEx_dx_0_1m5p1p11m5p0k = q(1-5+1+1,1-5+0,indvars(2))*(q(1-5+1+1,1-5+0,indvars(1))*q(1-5+1+1,1-5+0,indvars(4))+&
                    0.375_wp*(8.0_wp*q(1-5+1+1,1-5+0,indvars(1))/(3.0_wp-&
                    q(1-5+1+1,1-5+0,indvars(1)))*param_float(1 + 5)*(((q(1-5+1+1,1-5+0,indvars(4))-&
                    0.5_wp*(q(1-5+1+1,1-5+0,indvars(2))*q(1-5+1+1,1-5+0,indvars(2))+&
                    q(1-5+1+1,1-5+0,indvars(3))*q(1-5+1+1,1-5+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+1+1,1-5+0,indvars(1)))-&
                    3.0_wp*q(1-5+1+1,1-5+0,indvars(1))*q(1-5+1+1,1-5+0,indvars(1))))

d1_FluEx_dx_0_1m5p11m5p0k = -&
          0.5_wp*d1_FluEx_dx_0_1m5p1m11m5p0k+&
          0.5_wp*d1_FluEx_dx_0_1m5p1p11m5p0k

d1_FluEx_dx_0_1m5p11m5p0k = d1_FluEx_dx_0_1m5p11m5p0k*param_float(1)

d1_FluEx_dy_0_1m5p11m5p0p0k = q(1-5+1,1-5+0+0,indvars(3))*(q(1-5+1,1-5+0+0,indvars(1))*q(1-5+1,1-5+0+0,indvars(4))+&
                    0.375_wp*(8.0_wp*q(1-5+1,1-5+0+0,indvars(1))/(3.0_wp-&
                    q(1-5+1,1-5+0+0,indvars(1)))*param_float(1 + 5)*(((q(1-5+1,1-5+0+0,indvars(4))-&
                    0.5_wp*(q(1-5+1,1-5+0+0,indvars(2))*q(1-5+1,1-5+0+0,indvars(2))+&
                    q(1-5+1,1-5+0+0,indvars(3))*q(1-5+1,1-5+0+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+1,1-5+0+0,indvars(1)))-&
                    3.0_wp*q(1-5+1,1-5+0+0,indvars(1))*q(1-5+1,1-5+0+0,indvars(1))))

d1_FluEx_dy_0_1m5p11m5p0p1k = q(1-5+1,1-5+0+1,indvars(3))*(q(1-5+1,1-5+0+1,indvars(1))*q(1-5+1,1-5+0+1,indvars(4))+&
                    0.375_wp*(8.0_wp*q(1-5+1,1-5+0+1,indvars(1))/(3.0_wp-&
                    q(1-5+1,1-5+0+1,indvars(1)))*param_float(1 + 5)*(((q(1-5+1,1-5+0+1,indvars(4))-&
                    0.5_wp*(q(1-5+1,1-5+0+1,indvars(2))*q(1-5+1,1-5+0+1,indvars(2))+&
                    q(1-5+1,1-5+0+1,indvars(3))*q(1-5+1,1-5+0+1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+1,1-5+0+1,indvars(1)))-&
                    3.0_wp*q(1-5+1,1-5+0+1,indvars(1))*q(1-5+1,1-5+0+1,indvars(1))))

d1_FluEx_dy_0_1m5p11m5p0p2k = q(1-5+1,1-5+0+2,indvars(3))*(q(1-5+1,1-5+0+2,indvars(1))*q(1-5+1,1-5+0+2,indvars(4))+&
                    0.375_wp*(8.0_wp*q(1-5+1,1-5+0+2,indvars(1))/(3.0_wp-&
                    q(1-5+1,1-5+0+2,indvars(1)))*param_float(1 + 5)*(((q(1-5+1,1-5+0+2,indvars(4))-&
                    0.5_wp*(q(1-5+1,1-5+0+2,indvars(2))*q(1-5+1,1-5+0+2,indvars(2))+&
                    q(1-5+1,1-5+0+2,indvars(3))*q(1-5+1,1-5+0+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+1,1-5+0+2,indvars(1)))-&
                    3.0_wp*q(1-5+1,1-5+0+2,indvars(1))*q(1-5+1,1-5+0+2,indvars(1))))

d1_FluEx_dy_0_1m5p11m5p0k = -&
          1.5_wp*d1_FluEx_dy_0_1m5p11m5p0p0k+&
          2.0_wp*d1_FluEx_dy_0_1m5p11m5p0p1k-&
          0.5_wp*d1_FluEx_dy_0_1m5p11m5p0p2k

d1_FluEx_dy_0_1m5p11m5p0k = d1_FluEx_dy_0_1m5p11m5p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 0 None d(rho et)/dt **********
!                                                           
!***********************************************************


rhs(1-5+1,1-5+0,indvars(4)) =   -  ( d1_FluEx_dx_0_1m5p11m5p0k+d1_FluEx_dy_0_1m5p11m5p0k ) 

