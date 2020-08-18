

!***********************************************************
!                                                           
! Start building layers for BC : imax jmax None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 3 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u]_1x+[rho*v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluRx_dx_0_nxp5p0p0nyp5m3k = q(nx+5+0+0,ny+5-3,indvars(1))*q(nx+5+0+0,ny+5-3,indvars(2))

d1_FluRx_dx_0_nxp5p0m1nyp5m3k = q(nx+5+0-1,ny+5-3,indvars(1))*q(nx+5+0-1,ny+5-3,indvars(2))

d1_FluRx_dx_0_nxp5p0m2nyp5m3k = q(nx+5+0-2,ny+5-3,indvars(1))*q(nx+5+0-2,ny+5-3,indvars(2))

d1_FluRx_dx_0_nxp5p0nyp5m3k = 1.5_wp*d1_FluRx_dx_0_nxp5p0p0nyp5m3k-&
          2.0_wp*d1_FluRx_dx_0_nxp5p0m1nyp5m3k+&
          0.5_wp*d1_FluRx_dx_0_nxp5p0m2nyp5m3k

d1_FluRx_dx_0_nxp5p0nyp5m3k = d1_FluRx_dx_0_nxp5p0nyp5m3k*param_float(1)

d1_FluRx_dy_0_nxp5p0nyp5m3m2k = q(nx+5+0,ny+5-3-2,indvars(1))*q(nx+5+0,ny+5-3-2,indvars(3))

d1_FluRx_dy_0_nxp5p0nyp5m3m1k = q(nx+5+0,ny+5-3-1,indvars(1))*q(nx+5+0,ny+5-3-1,indvars(3))

d1_FluRx_dy_0_nxp5p0nyp5m3p1k = q(nx+5+0,ny+5-3+1,indvars(1))*q(nx+5+0,ny+5-3+1,indvars(3))

d1_FluRx_dy_0_nxp5p0nyp5m3p2k = q(nx+5+0,ny+5-3+2,indvars(1))*q(nx+5+0,ny+5-3+2,indvars(3))

d1_FluRx_dy_0_nxp5p0nyp5m3k = 0.08333333333333333_wp*d1_FluRx_dy_0_nxp5p0nyp5m3m2k-&
          0.666666666666667_wp*d1_FluRx_dy_0_nxp5p0nyp5m3m1k+&
          0.666666666666667_wp*d1_FluRx_dy_0_nxp5p0nyp5m3p1k-&
          0.08333333333333333_wp*d1_FluRx_dy_0_nxp5p0nyp5m3p2k

d1_FluRx_dy_0_nxp5p0nyp5m3k = d1_FluRx_dy_0_nxp5p0nyp5m3k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None d(rho)/dt *************
!                                                           
!***********************************************************


rhs(nx+5+0,ny+5-3,indvars(1)) =   -  ( d1_FluRx_dx_0_nxp5p0nyp5m3k+d1_FluRx_dy_0_nxp5p0nyp5m3k ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u*u+0.375_wp*p]_1x+[rho*v*u]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluXx_dx_0_nxp5p0p0nyp5m3k = q(nx+5+0+0,ny+5-3,indvars(1))*q(nx+5+0+0,ny+5-3,indvars(2))*q(nx+5+0+0,ny+5-3,indvars(2))+&
                    0.375_wp*(8.0_wp*q(nx+5+0+0,ny+5-3,indvars(1))/(3.0_wp-&
                    q(nx+5+0+0,ny+5-3,indvars(1)))*param_float(1 + 5)*(((q(nx+5+0+0,ny+5-3,indvars(4))-&
                    0.5_wp*(q(nx+5+0+0,ny+5-3,indvars(2))*q(nx+5+0+0,ny+5-3,indvars(2))+&
                    q(nx+5+0+0,ny+5-3,indvars(3))*q(nx+5+0+0,ny+5-3,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5+0+0,ny+5-3,indvars(1)))-&
                    3.0_wp*q(nx+5+0+0,ny+5-3,indvars(1))*q(nx+5+0+0,ny+5-3,indvars(1)))

d1_FluXx_dx_0_nxp5p0m1nyp5m3k = q(nx+5+0-1,ny+5-3,indvars(1))*q(nx+5+0-1,ny+5-3,indvars(2))*q(nx+5+0-1,ny+5-3,indvars(2))+&
                    0.375_wp*(8.0_wp*q(nx+5+0-1,ny+5-3,indvars(1))/(3.0_wp-&
                    q(nx+5+0-1,ny+5-3,indvars(1)))*param_float(1 + 5)*(((q(nx+5+0-1,ny+5-3,indvars(4))-&
                    0.5_wp*(q(nx+5+0-1,ny+5-3,indvars(2))*q(nx+5+0-1,ny+5-3,indvars(2))+&
                    q(nx+5+0-1,ny+5-3,indvars(3))*q(nx+5+0-1,ny+5-3,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5+0-1,ny+5-3,indvars(1)))-&
                    3.0_wp*q(nx+5+0-1,ny+5-3,indvars(1))*q(nx+5+0-1,ny+5-3,indvars(1)))

d1_FluXx_dx_0_nxp5p0m2nyp5m3k = q(nx+5+0-2,ny+5-3,indvars(1))*q(nx+5+0-2,ny+5-3,indvars(2))*q(nx+5+0-2,ny+5-3,indvars(2))+&
                    0.375_wp*(8.0_wp*q(nx+5+0-2,ny+5-3,indvars(1))/(3.0_wp-&
                    q(nx+5+0-2,ny+5-3,indvars(1)))*param_float(1 + 5)*(((q(nx+5+0-2,ny+5-3,indvars(4))-&
                    0.5_wp*(q(nx+5+0-2,ny+5-3,indvars(2))*q(nx+5+0-2,ny+5-3,indvars(2))+&
                    q(nx+5+0-2,ny+5-3,indvars(3))*q(nx+5+0-2,ny+5-3,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5+0-2,ny+5-3,indvars(1)))-&
                    3.0_wp*q(nx+5+0-2,ny+5-3,indvars(1))*q(nx+5+0-2,ny+5-3,indvars(1)))

d1_FluXx_dx_0_nxp5p0nyp5m3k = 1.5_wp*d1_FluXx_dx_0_nxp5p0p0nyp5m3k-&
          2.0_wp*d1_FluXx_dx_0_nxp5p0m1nyp5m3k+&
          0.5_wp*d1_FluXx_dx_0_nxp5p0m2nyp5m3k

d1_FluXx_dx_0_nxp5p0nyp5m3k = d1_FluXx_dx_0_nxp5p0nyp5m3k*param_float(1)

d1_FluXx_dy_0_nxp5p0nyp5m3m2k = q(nx+5+0,ny+5-3-2,indvars(1))*q(nx+5+0,ny+5-3-2,indvars(3))*q(nx+5+0,ny+5-3-2,indvars(2))

d1_FluXx_dy_0_nxp5p0nyp5m3m1k = q(nx+5+0,ny+5-3-1,indvars(1))*q(nx+5+0,ny+5-3-1,indvars(3))*q(nx+5+0,ny+5-3-1,indvars(2))

d1_FluXx_dy_0_nxp5p0nyp5m3p1k = q(nx+5+0,ny+5-3+1,indvars(1))*q(nx+5+0,ny+5-3+1,indvars(3))*q(nx+5+0,ny+5-3+1,indvars(2))

d1_FluXx_dy_0_nxp5p0nyp5m3p2k = q(nx+5+0,ny+5-3+2,indvars(1))*q(nx+5+0,ny+5-3+2,indvars(3))*q(nx+5+0,ny+5-3+2,indvars(2))

d1_FluXx_dy_0_nxp5p0nyp5m3k = 0.08333333333333333_wp*d1_FluXx_dy_0_nxp5p0nyp5m3m2k-&
          0.666666666666667_wp*d1_FluXx_dy_0_nxp5p0nyp5m3m1k+&
          0.666666666666667_wp*d1_FluXx_dy_0_nxp5p0nyp5m3p1k-&
          0.08333333333333333_wp*d1_FluXx_dy_0_nxp5p0nyp5m3p2k

d1_FluXx_dy_0_nxp5p0nyp5m3k = d1_FluXx_dy_0_nxp5p0nyp5m3k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None d(rho u)/dt ***********
!                                                           
!***********************************************************


rhs(nx+5+0,ny+5-3,indvars(2)) =   -  ( d1_FluXx_dx_0_nxp5p0nyp5m3k+d1_FluXx_dy_0_nxp5p0nyp5m3k ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u*v]_1x+[rho*v*v+0.375_wp*p]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluYx_dx_0_nxp5p0p0nyp5m3k = q(nx+5+0+0,ny+5-3,indvars(1))*q(nx+5+0+0,ny+5-3,indvars(2))*q(nx+5+0+0,ny+5-3,indvars(3))

d1_FluYx_dx_0_nxp5p0m1nyp5m3k = q(nx+5+0-1,ny+5-3,indvars(1))*q(nx+5+0-1,ny+5-3,indvars(2))*q(nx+5+0-1,ny+5-3,indvars(3))

d1_FluYx_dx_0_nxp5p0m2nyp5m3k = q(nx+5+0-2,ny+5-3,indvars(1))*q(nx+5+0-2,ny+5-3,indvars(2))*q(nx+5+0-2,ny+5-3,indvars(3))

d1_FluYx_dx_0_nxp5p0nyp5m3k = 1.5_wp*d1_FluYx_dx_0_nxp5p0p0nyp5m3k-&
          2.0_wp*d1_FluYx_dx_0_nxp5p0m1nyp5m3k+&
          0.5_wp*d1_FluYx_dx_0_nxp5p0m2nyp5m3k

d1_FluYx_dx_0_nxp5p0nyp5m3k = d1_FluYx_dx_0_nxp5p0nyp5m3k*param_float(1)

d1_FluYx_dy_0_nxp5p0nyp5m3m2k = q(nx+5+0,ny+5-3-2,indvars(1))*q(nx+5+0,ny+5-3-2,indvars(3))*q(nx+5+0,ny+5-3-2,indvars(3))+&
                    0.375_wp*(8.0_wp*q(nx+5+0,ny+5-3-2,indvars(1))/(3.0_wp-&
                    q(nx+5+0,ny+5-3-2,indvars(1)))*param_float(1 + 5)*(((q(nx+5+0,ny+5-3-2,indvars(4))-&
                    0.5_wp*(q(nx+5+0,ny+5-3-2,indvars(2))*q(nx+5+0,ny+5-3-2,indvars(2))+&
                    q(nx+5+0,ny+5-3-2,indvars(3))*q(nx+5+0,ny+5-3-2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5+0,ny+5-3-2,indvars(1)))-&
                    3.0_wp*q(nx+5+0,ny+5-3-2,indvars(1))*q(nx+5+0,ny+5-3-2,indvars(1)))

d1_FluYx_dy_0_nxp5p0nyp5m3m1k = q(nx+5+0,ny+5-3-1,indvars(1))*q(nx+5+0,ny+5-3-1,indvars(3))*q(nx+5+0,ny+5-3-1,indvars(3))+&
                    0.375_wp*(8.0_wp*q(nx+5+0,ny+5-3-1,indvars(1))/(3.0_wp-&
                    q(nx+5+0,ny+5-3-1,indvars(1)))*param_float(1 + 5)*(((q(nx+5+0,ny+5-3-1,indvars(4))-&
                    0.5_wp*(q(nx+5+0,ny+5-3-1,indvars(2))*q(nx+5+0,ny+5-3-1,indvars(2))+&
                    q(nx+5+0,ny+5-3-1,indvars(3))*q(nx+5+0,ny+5-3-1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5+0,ny+5-3-1,indvars(1)))-&
                    3.0_wp*q(nx+5+0,ny+5-3-1,indvars(1))*q(nx+5+0,ny+5-3-1,indvars(1)))

d1_FluYx_dy_0_nxp5p0nyp5m3p1k = q(nx+5+0,ny+5-3+1,indvars(1))*q(nx+5+0,ny+5-3+1,indvars(3))*q(nx+5+0,ny+5-3+1,indvars(3))+&
                    0.375_wp*(8.0_wp*q(nx+5+0,ny+5-3+1,indvars(1))/(3.0_wp-&
                    q(nx+5+0,ny+5-3+1,indvars(1)))*param_float(1 + 5)*(((q(nx+5+0,ny+5-3+1,indvars(4))-&
                    0.5_wp*(q(nx+5+0,ny+5-3+1,indvars(2))*q(nx+5+0,ny+5-3+1,indvars(2))+&
                    q(nx+5+0,ny+5-3+1,indvars(3))*q(nx+5+0,ny+5-3+1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5+0,ny+5-3+1,indvars(1)))-&
                    3.0_wp*q(nx+5+0,ny+5-3+1,indvars(1))*q(nx+5+0,ny+5-3+1,indvars(1)))

d1_FluYx_dy_0_nxp5p0nyp5m3p2k = q(nx+5+0,ny+5-3+2,indvars(1))*q(nx+5+0,ny+5-3+2,indvars(3))*q(nx+5+0,ny+5-3+2,indvars(3))+&
                    0.375_wp*(8.0_wp*q(nx+5+0,ny+5-3+2,indvars(1))/(3.0_wp-&
                    q(nx+5+0,ny+5-3+2,indvars(1)))*param_float(1 + 5)*(((q(nx+5+0,ny+5-3+2,indvars(4))-&
                    0.5_wp*(q(nx+5+0,ny+5-3+2,indvars(2))*q(nx+5+0,ny+5-3+2,indvars(2))+&
                    q(nx+5+0,ny+5-3+2,indvars(3))*q(nx+5+0,ny+5-3+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5+0,ny+5-3+2,indvars(1)))-&
                    3.0_wp*q(nx+5+0,ny+5-3+2,indvars(1))*q(nx+5+0,ny+5-3+2,indvars(1)))

d1_FluYx_dy_0_nxp5p0nyp5m3k = 0.08333333333333333_wp*d1_FluYx_dy_0_nxp5p0nyp5m3m2k-&
          0.666666666666667_wp*d1_FluYx_dy_0_nxp5p0nyp5m3m1k+&
          0.666666666666667_wp*d1_FluYx_dy_0_nxp5p0nyp5m3p1k-&
          0.08333333333333333_wp*d1_FluYx_dy_0_nxp5p0nyp5m3p2k

d1_FluYx_dy_0_nxp5p0nyp5m3k = d1_FluYx_dy_0_nxp5p0nyp5m3k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None d(rho v)/dt ***********
!                                                           
!***********************************************************


rhs(nx+5+0,ny+5-3,indvars(3)) =   -  ( d1_FluYx_dx_0_nxp5p0nyp5m3k+d1_FluYx_dy_0_nxp5p0nyp5m3k ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u*(rho*et+0.375_wp*p)]_1x+[v*(rho*et+0.375_wp*p)]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluEx_dx_0_nxp5p0p0nyp5m3k = q(nx+5+0+0,ny+5-3,indvars(2))*(q(nx+5+0+0,ny+5-3,indvars(1))*q(nx+5+0+0,ny+5-3,indvars(4))+&
                    0.375_wp*(8.0_wp*q(nx+5+0+0,ny+5-3,indvars(1))/(3.0_wp-&
                    q(nx+5+0+0,ny+5-3,indvars(1)))*param_float(1 + 5)*(((q(nx+5+0+0,ny+5-3,indvars(4))-&
                    0.5_wp*(q(nx+5+0+0,ny+5-3,indvars(2))*q(nx+5+0+0,ny+5-3,indvars(2))+&
                    q(nx+5+0+0,ny+5-3,indvars(3))*q(nx+5+0+0,ny+5-3,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5+0+0,ny+5-3,indvars(1)))-&
                    3.0_wp*q(nx+5+0+0,ny+5-3,indvars(1))*q(nx+5+0+0,ny+5-3,indvars(1))))

d1_FluEx_dx_0_nxp5p0m1nyp5m3k = q(nx+5+0-1,ny+5-3,indvars(2))*(q(nx+5+0-1,ny+5-3,indvars(1))*q(nx+5+0-1,ny+5-3,indvars(4))+&
                    0.375_wp*(8.0_wp*q(nx+5+0-1,ny+5-3,indvars(1))/(3.0_wp-&
                    q(nx+5+0-1,ny+5-3,indvars(1)))*param_float(1 + 5)*(((q(nx+5+0-1,ny+5-3,indvars(4))-&
                    0.5_wp*(q(nx+5+0-1,ny+5-3,indvars(2))*q(nx+5+0-1,ny+5-3,indvars(2))+&
                    q(nx+5+0-1,ny+5-3,indvars(3))*q(nx+5+0-1,ny+5-3,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5+0-1,ny+5-3,indvars(1)))-&
                    3.0_wp*q(nx+5+0-1,ny+5-3,indvars(1))*q(nx+5+0-1,ny+5-3,indvars(1))))

d1_FluEx_dx_0_nxp5p0m2nyp5m3k = q(nx+5+0-2,ny+5-3,indvars(2))*(q(nx+5+0-2,ny+5-3,indvars(1))*q(nx+5+0-2,ny+5-3,indvars(4))+&
                    0.375_wp*(8.0_wp*q(nx+5+0-2,ny+5-3,indvars(1))/(3.0_wp-&
                    q(nx+5+0-2,ny+5-3,indvars(1)))*param_float(1 + 5)*(((q(nx+5+0-2,ny+5-3,indvars(4))-&
                    0.5_wp*(q(nx+5+0-2,ny+5-3,indvars(2))*q(nx+5+0-2,ny+5-3,indvars(2))+&
                    q(nx+5+0-2,ny+5-3,indvars(3))*q(nx+5+0-2,ny+5-3,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5+0-2,ny+5-3,indvars(1)))-&
                    3.0_wp*q(nx+5+0-2,ny+5-3,indvars(1))*q(nx+5+0-2,ny+5-3,indvars(1))))

d1_FluEx_dx_0_nxp5p0nyp5m3k = 1.5_wp*d1_FluEx_dx_0_nxp5p0p0nyp5m3k-&
          2.0_wp*d1_FluEx_dx_0_nxp5p0m1nyp5m3k+&
          0.5_wp*d1_FluEx_dx_0_nxp5p0m2nyp5m3k

d1_FluEx_dx_0_nxp5p0nyp5m3k = d1_FluEx_dx_0_nxp5p0nyp5m3k*param_float(1)

d1_FluEx_dy_0_nxp5p0nyp5m3m2k = q(nx+5+0,ny+5-3-2,indvars(3))*(q(nx+5+0,ny+5-3-2,indvars(1))*q(nx+5+0,ny+5-3-2,indvars(4))+&
                    0.375_wp*(8.0_wp*q(nx+5+0,ny+5-3-2,indvars(1))/(3.0_wp-&
                    q(nx+5+0,ny+5-3-2,indvars(1)))*param_float(1 + 5)*(((q(nx+5+0,ny+5-3-2,indvars(4))-&
                    0.5_wp*(q(nx+5+0,ny+5-3-2,indvars(2))*q(nx+5+0,ny+5-3-2,indvars(2))+&
                    q(nx+5+0,ny+5-3-2,indvars(3))*q(nx+5+0,ny+5-3-2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5+0,ny+5-3-2,indvars(1)))-&
                    3.0_wp*q(nx+5+0,ny+5-3-2,indvars(1))*q(nx+5+0,ny+5-3-2,indvars(1))))

d1_FluEx_dy_0_nxp5p0nyp5m3m1k = q(nx+5+0,ny+5-3-1,indvars(3))*(q(nx+5+0,ny+5-3-1,indvars(1))*q(nx+5+0,ny+5-3-1,indvars(4))+&
                    0.375_wp*(8.0_wp*q(nx+5+0,ny+5-3-1,indvars(1))/(3.0_wp-&
                    q(nx+5+0,ny+5-3-1,indvars(1)))*param_float(1 + 5)*(((q(nx+5+0,ny+5-3-1,indvars(4))-&
                    0.5_wp*(q(nx+5+0,ny+5-3-1,indvars(2))*q(nx+5+0,ny+5-3-1,indvars(2))+&
                    q(nx+5+0,ny+5-3-1,indvars(3))*q(nx+5+0,ny+5-3-1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5+0,ny+5-3-1,indvars(1)))-&
                    3.0_wp*q(nx+5+0,ny+5-3-1,indvars(1))*q(nx+5+0,ny+5-3-1,indvars(1))))

d1_FluEx_dy_0_nxp5p0nyp5m3p1k = q(nx+5+0,ny+5-3+1,indvars(3))*(q(nx+5+0,ny+5-3+1,indvars(1))*q(nx+5+0,ny+5-3+1,indvars(4))+&
                    0.375_wp*(8.0_wp*q(nx+5+0,ny+5-3+1,indvars(1))/(3.0_wp-&
                    q(nx+5+0,ny+5-3+1,indvars(1)))*param_float(1 + 5)*(((q(nx+5+0,ny+5-3+1,indvars(4))-&
                    0.5_wp*(q(nx+5+0,ny+5-3+1,indvars(2))*q(nx+5+0,ny+5-3+1,indvars(2))+&
                    q(nx+5+0,ny+5-3+1,indvars(3))*q(nx+5+0,ny+5-3+1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5+0,ny+5-3+1,indvars(1)))-&
                    3.0_wp*q(nx+5+0,ny+5-3+1,indvars(1))*q(nx+5+0,ny+5-3+1,indvars(1))))

d1_FluEx_dy_0_nxp5p0nyp5m3p2k = q(nx+5+0,ny+5-3+2,indvars(3))*(q(nx+5+0,ny+5-3+2,indvars(1))*q(nx+5+0,ny+5-3+2,indvars(4))+&
                    0.375_wp*(8.0_wp*q(nx+5+0,ny+5-3+2,indvars(1))/(3.0_wp-&
                    q(nx+5+0,ny+5-3+2,indvars(1)))*param_float(1 + 5)*(((q(nx+5+0,ny+5-3+2,indvars(4))-&
                    0.5_wp*(q(nx+5+0,ny+5-3+2,indvars(2))*q(nx+5+0,ny+5-3+2,indvars(2))+&
                    q(nx+5+0,ny+5-3+2,indvars(3))*q(nx+5+0,ny+5-3+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5+0,ny+5-3+2,indvars(1)))-&
                    3.0_wp*q(nx+5+0,ny+5-3+2,indvars(1))*q(nx+5+0,ny+5-3+2,indvars(1))))

d1_FluEx_dy_0_nxp5p0nyp5m3k = 0.08333333333333333_wp*d1_FluEx_dy_0_nxp5p0nyp5m3m2k-&
          0.666666666666667_wp*d1_FluEx_dy_0_nxp5p0nyp5m3m1k+&
          0.666666666666667_wp*d1_FluEx_dy_0_nxp5p0nyp5m3p1k-&
          0.08333333333333333_wp*d1_FluEx_dy_0_nxp5p0nyp5m3p2k

d1_FluEx_dy_0_nxp5p0nyp5m3k = d1_FluEx_dy_0_nxp5p0nyp5m3k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None d(rho et)/dt **********
!                                                           
!***********************************************************


rhs(nx+5+0,ny+5-3,indvars(4)) =   -  ( d1_FluEx_dx_0_nxp5p0nyp5m3k+d1_FluEx_dy_0_nxp5p0nyp5m3k ) 

