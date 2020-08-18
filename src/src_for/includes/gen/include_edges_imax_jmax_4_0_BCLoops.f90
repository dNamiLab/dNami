

!***********************************************************
!                                                           
! Start building layers for BC : imax jmax None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 4 0 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 4 0 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u]_1x+[rho*v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluRx_dx_0_nxp5m4m2nyp5p0k = q(nx+5-4-2,ny+5+0,indvars(1))*q(nx+5-4-2,ny+5+0,indvars(2))

d1_FluRx_dx_0_nxp5m4m1nyp5p0k = q(nx+5-4-1,ny+5+0,indvars(1))*q(nx+5-4-1,ny+5+0,indvars(2))

d1_FluRx_dx_0_nxp5m4p1nyp5p0k = q(nx+5-4+1,ny+5+0,indvars(1))*q(nx+5-4+1,ny+5+0,indvars(2))

d1_FluRx_dx_0_nxp5m4p2nyp5p0k = q(nx+5-4+2,ny+5+0,indvars(1))*q(nx+5-4+2,ny+5+0,indvars(2))

d1_FluRx_dx_0_nxp5m4nyp5p0k = 0.08333333333333333_wp*d1_FluRx_dx_0_nxp5m4m2nyp5p0k-&
          0.666666666666667_wp*d1_FluRx_dx_0_nxp5m4m1nyp5p0k+&
          0.666666666666667_wp*d1_FluRx_dx_0_nxp5m4p1nyp5p0k-&
          0.08333333333333333_wp*d1_FluRx_dx_0_nxp5m4p2nyp5p0k

d1_FluRx_dx_0_nxp5m4nyp5p0k = d1_FluRx_dx_0_nxp5m4nyp5p0k*param_float(1)

d1_FluRx_dy_0_nxp5m4nyp5p0p0k = q(nx+5-4,ny+5+0+0,indvars(1))*q(nx+5-4,ny+5+0+0,indvars(3))

d1_FluRx_dy_0_nxp5m4nyp5p0m1k = q(nx+5-4,ny+5+0-1,indvars(1))*q(nx+5-4,ny+5+0-1,indvars(3))

d1_FluRx_dy_0_nxp5m4nyp5p0m2k = q(nx+5-4,ny+5+0-2,indvars(1))*q(nx+5-4,ny+5+0-2,indvars(3))

d1_FluRx_dy_0_nxp5m4nyp5p0k = 1.5_wp*d1_FluRx_dy_0_nxp5m4nyp5p0p0k-&
          2.0_wp*d1_FluRx_dy_0_nxp5m4nyp5p0m1k+&
          0.5_wp*d1_FluRx_dy_0_nxp5m4nyp5p0m2k

d1_FluRx_dy_0_nxp5m4nyp5p0k = d1_FluRx_dy_0_nxp5m4nyp5p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 4 0 None d(rho)/dt *************
!                                                           
!***********************************************************


rhs(nx+5-4,ny+5+0,indvars(1)) =   -  ( d1_FluRx_dx_0_nxp5m4nyp5p0k+d1_FluRx_dy_0_nxp5m4nyp5p0k ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 4 0 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u*u+0.375_wp*p]_1x+[rho*v*u]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluXx_dx_0_nxp5m4m2nyp5p0k = q(nx+5-4-2,ny+5+0,indvars(1))*q(nx+5-4-2,ny+5+0,indvars(2))*q(nx+5-4-2,ny+5+0,indvars(2))+&
                    0.375_wp*(8.0_wp*q(nx+5-4-2,ny+5+0,indvars(1))/(3.0_wp-&
                    q(nx+5-4-2,ny+5+0,indvars(1)))*param_float(1 + 5)*(((q(nx+5-4-2,ny+5+0,indvars(4))-&
                    0.5_wp*(q(nx+5-4-2,ny+5+0,indvars(2))*q(nx+5-4-2,ny+5+0,indvars(2))+&
                    q(nx+5-4-2,ny+5+0,indvars(3))*q(nx+5-4-2,ny+5+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5-4-2,ny+5+0,indvars(1)))-&
                    3.0_wp*q(nx+5-4-2,ny+5+0,indvars(1))*q(nx+5-4-2,ny+5+0,indvars(1)))

d1_FluXx_dx_0_nxp5m4m1nyp5p0k = q(nx+5-4-1,ny+5+0,indvars(1))*q(nx+5-4-1,ny+5+0,indvars(2))*q(nx+5-4-1,ny+5+0,indvars(2))+&
                    0.375_wp*(8.0_wp*q(nx+5-4-1,ny+5+0,indvars(1))/(3.0_wp-&
                    q(nx+5-4-1,ny+5+0,indvars(1)))*param_float(1 + 5)*(((q(nx+5-4-1,ny+5+0,indvars(4))-&
                    0.5_wp*(q(nx+5-4-1,ny+5+0,indvars(2))*q(nx+5-4-1,ny+5+0,indvars(2))+&
                    q(nx+5-4-1,ny+5+0,indvars(3))*q(nx+5-4-1,ny+5+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5-4-1,ny+5+0,indvars(1)))-&
                    3.0_wp*q(nx+5-4-1,ny+5+0,indvars(1))*q(nx+5-4-1,ny+5+0,indvars(1)))

d1_FluXx_dx_0_nxp5m4p1nyp5p0k = q(nx+5-4+1,ny+5+0,indvars(1))*q(nx+5-4+1,ny+5+0,indvars(2))*q(nx+5-4+1,ny+5+0,indvars(2))+&
                    0.375_wp*(8.0_wp*q(nx+5-4+1,ny+5+0,indvars(1))/(3.0_wp-&
                    q(nx+5-4+1,ny+5+0,indvars(1)))*param_float(1 + 5)*(((q(nx+5-4+1,ny+5+0,indvars(4))-&
                    0.5_wp*(q(nx+5-4+1,ny+5+0,indvars(2))*q(nx+5-4+1,ny+5+0,indvars(2))+&
                    q(nx+5-4+1,ny+5+0,indvars(3))*q(nx+5-4+1,ny+5+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5-4+1,ny+5+0,indvars(1)))-&
                    3.0_wp*q(nx+5-4+1,ny+5+0,indvars(1))*q(nx+5-4+1,ny+5+0,indvars(1)))

d1_FluXx_dx_0_nxp5m4p2nyp5p0k = q(nx+5-4+2,ny+5+0,indvars(1))*q(nx+5-4+2,ny+5+0,indvars(2))*q(nx+5-4+2,ny+5+0,indvars(2))+&
                    0.375_wp*(8.0_wp*q(nx+5-4+2,ny+5+0,indvars(1))/(3.0_wp-&
                    q(nx+5-4+2,ny+5+0,indvars(1)))*param_float(1 + 5)*(((q(nx+5-4+2,ny+5+0,indvars(4))-&
                    0.5_wp*(q(nx+5-4+2,ny+5+0,indvars(2))*q(nx+5-4+2,ny+5+0,indvars(2))+&
                    q(nx+5-4+2,ny+5+0,indvars(3))*q(nx+5-4+2,ny+5+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5-4+2,ny+5+0,indvars(1)))-&
                    3.0_wp*q(nx+5-4+2,ny+5+0,indvars(1))*q(nx+5-4+2,ny+5+0,indvars(1)))

d1_FluXx_dx_0_nxp5m4nyp5p0k = 0.08333333333333333_wp*d1_FluXx_dx_0_nxp5m4m2nyp5p0k-&
          0.666666666666667_wp*d1_FluXx_dx_0_nxp5m4m1nyp5p0k+&
          0.666666666666667_wp*d1_FluXx_dx_0_nxp5m4p1nyp5p0k-&
          0.08333333333333333_wp*d1_FluXx_dx_0_nxp5m4p2nyp5p0k

d1_FluXx_dx_0_nxp5m4nyp5p0k = d1_FluXx_dx_0_nxp5m4nyp5p0k*param_float(1)

d1_FluXx_dy_0_nxp5m4nyp5p0p0k = q(nx+5-4,ny+5+0+0,indvars(1))*q(nx+5-4,ny+5+0+0,indvars(3))*q(nx+5-4,ny+5+0+0,indvars(2))

d1_FluXx_dy_0_nxp5m4nyp5p0m1k = q(nx+5-4,ny+5+0-1,indvars(1))*q(nx+5-4,ny+5+0-1,indvars(3))*q(nx+5-4,ny+5+0-1,indvars(2))

d1_FluXx_dy_0_nxp5m4nyp5p0m2k = q(nx+5-4,ny+5+0-2,indvars(1))*q(nx+5-4,ny+5+0-2,indvars(3))*q(nx+5-4,ny+5+0-2,indvars(2))

d1_FluXx_dy_0_nxp5m4nyp5p0k = 1.5_wp*d1_FluXx_dy_0_nxp5m4nyp5p0p0k-&
          2.0_wp*d1_FluXx_dy_0_nxp5m4nyp5p0m1k+&
          0.5_wp*d1_FluXx_dy_0_nxp5m4nyp5p0m2k

d1_FluXx_dy_0_nxp5m4nyp5p0k = d1_FluXx_dy_0_nxp5m4nyp5p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 4 0 None d(rho u)/dt ***********
!                                                           
!***********************************************************


rhs(nx+5-4,ny+5+0,indvars(2)) =   -  ( d1_FluXx_dx_0_nxp5m4nyp5p0k+d1_FluXx_dy_0_nxp5m4nyp5p0k ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 4 0 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u*v]_1x+[rho*v*v+0.375_wp*p]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluYx_dx_0_nxp5m4m2nyp5p0k = q(nx+5-4-2,ny+5+0,indvars(1))*q(nx+5-4-2,ny+5+0,indvars(2))*q(nx+5-4-2,ny+5+0,indvars(3))

d1_FluYx_dx_0_nxp5m4m1nyp5p0k = q(nx+5-4-1,ny+5+0,indvars(1))*q(nx+5-4-1,ny+5+0,indvars(2))*q(nx+5-4-1,ny+5+0,indvars(3))

d1_FluYx_dx_0_nxp5m4p1nyp5p0k = q(nx+5-4+1,ny+5+0,indvars(1))*q(nx+5-4+1,ny+5+0,indvars(2))*q(nx+5-4+1,ny+5+0,indvars(3))

d1_FluYx_dx_0_nxp5m4p2nyp5p0k = q(nx+5-4+2,ny+5+0,indvars(1))*q(nx+5-4+2,ny+5+0,indvars(2))*q(nx+5-4+2,ny+5+0,indvars(3))

d1_FluYx_dx_0_nxp5m4nyp5p0k = 0.08333333333333333_wp*d1_FluYx_dx_0_nxp5m4m2nyp5p0k-&
          0.666666666666667_wp*d1_FluYx_dx_0_nxp5m4m1nyp5p0k+&
          0.666666666666667_wp*d1_FluYx_dx_0_nxp5m4p1nyp5p0k-&
          0.08333333333333333_wp*d1_FluYx_dx_0_nxp5m4p2nyp5p0k

d1_FluYx_dx_0_nxp5m4nyp5p0k = d1_FluYx_dx_0_nxp5m4nyp5p0k*param_float(1)

d1_FluYx_dy_0_nxp5m4nyp5p0p0k = q(nx+5-4,ny+5+0+0,indvars(1))*q(nx+5-4,ny+5+0+0,indvars(3))*q(nx+5-4,ny+5+0+0,indvars(3))+&
                    0.375_wp*(8.0_wp*q(nx+5-4,ny+5+0+0,indvars(1))/(3.0_wp-&
                    q(nx+5-4,ny+5+0+0,indvars(1)))*param_float(1 + 5)*(((q(nx+5-4,ny+5+0+0,indvars(4))-&
                    0.5_wp*(q(nx+5-4,ny+5+0+0,indvars(2))*q(nx+5-4,ny+5+0+0,indvars(2))+&
                    q(nx+5-4,ny+5+0+0,indvars(3))*q(nx+5-4,ny+5+0+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5-4,ny+5+0+0,indvars(1)))-&
                    3.0_wp*q(nx+5-4,ny+5+0+0,indvars(1))*q(nx+5-4,ny+5+0+0,indvars(1)))

d1_FluYx_dy_0_nxp5m4nyp5p0m1k = q(nx+5-4,ny+5+0-1,indvars(1))*q(nx+5-4,ny+5+0-1,indvars(3))*q(nx+5-4,ny+5+0-1,indvars(3))+&
                    0.375_wp*(8.0_wp*q(nx+5-4,ny+5+0-1,indvars(1))/(3.0_wp-&
                    q(nx+5-4,ny+5+0-1,indvars(1)))*param_float(1 + 5)*(((q(nx+5-4,ny+5+0-1,indvars(4))-&
                    0.5_wp*(q(nx+5-4,ny+5+0-1,indvars(2))*q(nx+5-4,ny+5+0-1,indvars(2))+&
                    q(nx+5-4,ny+5+0-1,indvars(3))*q(nx+5-4,ny+5+0-1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5-4,ny+5+0-1,indvars(1)))-&
                    3.0_wp*q(nx+5-4,ny+5+0-1,indvars(1))*q(nx+5-4,ny+5+0-1,indvars(1)))

d1_FluYx_dy_0_nxp5m4nyp5p0m2k = q(nx+5-4,ny+5+0-2,indvars(1))*q(nx+5-4,ny+5+0-2,indvars(3))*q(nx+5-4,ny+5+0-2,indvars(3))+&
                    0.375_wp*(8.0_wp*q(nx+5-4,ny+5+0-2,indvars(1))/(3.0_wp-&
                    q(nx+5-4,ny+5+0-2,indvars(1)))*param_float(1 + 5)*(((q(nx+5-4,ny+5+0-2,indvars(4))-&
                    0.5_wp*(q(nx+5-4,ny+5+0-2,indvars(2))*q(nx+5-4,ny+5+0-2,indvars(2))+&
                    q(nx+5-4,ny+5+0-2,indvars(3))*q(nx+5-4,ny+5+0-2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5-4,ny+5+0-2,indvars(1)))-&
                    3.0_wp*q(nx+5-4,ny+5+0-2,indvars(1))*q(nx+5-4,ny+5+0-2,indvars(1)))

d1_FluYx_dy_0_nxp5m4nyp5p0k = 1.5_wp*d1_FluYx_dy_0_nxp5m4nyp5p0p0k-&
          2.0_wp*d1_FluYx_dy_0_nxp5m4nyp5p0m1k+&
          0.5_wp*d1_FluYx_dy_0_nxp5m4nyp5p0m2k

d1_FluYx_dy_0_nxp5m4nyp5p0k = d1_FluYx_dy_0_nxp5m4nyp5p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 4 0 None d(rho v)/dt ***********
!                                                           
!***********************************************************


rhs(nx+5-4,ny+5+0,indvars(3)) =   -  ( d1_FluYx_dx_0_nxp5m4nyp5p0k+d1_FluYx_dy_0_nxp5m4nyp5p0k ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 4 0 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u*(rho*et+0.375_wp*p)]_1x+[v*(rho*et+0.375_wp*p)]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluEx_dx_0_nxp5m4m2nyp5p0k = q(nx+5-4-2,ny+5+0,indvars(2))*(q(nx+5-4-2,ny+5+0,indvars(1))*q(nx+5-4-2,ny+5+0,indvars(4))+&
                    0.375_wp*(8.0_wp*q(nx+5-4-2,ny+5+0,indvars(1))/(3.0_wp-&
                    q(nx+5-4-2,ny+5+0,indvars(1)))*param_float(1 + 5)*(((q(nx+5-4-2,ny+5+0,indvars(4))-&
                    0.5_wp*(q(nx+5-4-2,ny+5+0,indvars(2))*q(nx+5-4-2,ny+5+0,indvars(2))+&
                    q(nx+5-4-2,ny+5+0,indvars(3))*q(nx+5-4-2,ny+5+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5-4-2,ny+5+0,indvars(1)))-&
                    3.0_wp*q(nx+5-4-2,ny+5+0,indvars(1))*q(nx+5-4-2,ny+5+0,indvars(1))))

d1_FluEx_dx_0_nxp5m4m1nyp5p0k = q(nx+5-4-1,ny+5+0,indvars(2))*(q(nx+5-4-1,ny+5+0,indvars(1))*q(nx+5-4-1,ny+5+0,indvars(4))+&
                    0.375_wp*(8.0_wp*q(nx+5-4-1,ny+5+0,indvars(1))/(3.0_wp-&
                    q(nx+5-4-1,ny+5+0,indvars(1)))*param_float(1 + 5)*(((q(nx+5-4-1,ny+5+0,indvars(4))-&
                    0.5_wp*(q(nx+5-4-1,ny+5+0,indvars(2))*q(nx+5-4-1,ny+5+0,indvars(2))+&
                    q(nx+5-4-1,ny+5+0,indvars(3))*q(nx+5-4-1,ny+5+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5-4-1,ny+5+0,indvars(1)))-&
                    3.0_wp*q(nx+5-4-1,ny+5+0,indvars(1))*q(nx+5-4-1,ny+5+0,indvars(1))))

d1_FluEx_dx_0_nxp5m4p1nyp5p0k = q(nx+5-4+1,ny+5+0,indvars(2))*(q(nx+5-4+1,ny+5+0,indvars(1))*q(nx+5-4+1,ny+5+0,indvars(4))+&
                    0.375_wp*(8.0_wp*q(nx+5-4+1,ny+5+0,indvars(1))/(3.0_wp-&
                    q(nx+5-4+1,ny+5+0,indvars(1)))*param_float(1 + 5)*(((q(nx+5-4+1,ny+5+0,indvars(4))-&
                    0.5_wp*(q(nx+5-4+1,ny+5+0,indvars(2))*q(nx+5-4+1,ny+5+0,indvars(2))+&
                    q(nx+5-4+1,ny+5+0,indvars(3))*q(nx+5-4+1,ny+5+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5-4+1,ny+5+0,indvars(1)))-&
                    3.0_wp*q(nx+5-4+1,ny+5+0,indvars(1))*q(nx+5-4+1,ny+5+0,indvars(1))))

d1_FluEx_dx_0_nxp5m4p2nyp5p0k = q(nx+5-4+2,ny+5+0,indvars(2))*(q(nx+5-4+2,ny+5+0,indvars(1))*q(nx+5-4+2,ny+5+0,indvars(4))+&
                    0.375_wp*(8.0_wp*q(nx+5-4+2,ny+5+0,indvars(1))/(3.0_wp-&
                    q(nx+5-4+2,ny+5+0,indvars(1)))*param_float(1 + 5)*(((q(nx+5-4+2,ny+5+0,indvars(4))-&
                    0.5_wp*(q(nx+5-4+2,ny+5+0,indvars(2))*q(nx+5-4+2,ny+5+0,indvars(2))+&
                    q(nx+5-4+2,ny+5+0,indvars(3))*q(nx+5-4+2,ny+5+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5-4+2,ny+5+0,indvars(1)))-&
                    3.0_wp*q(nx+5-4+2,ny+5+0,indvars(1))*q(nx+5-4+2,ny+5+0,indvars(1))))

d1_FluEx_dx_0_nxp5m4nyp5p0k = 0.08333333333333333_wp*d1_FluEx_dx_0_nxp5m4m2nyp5p0k-&
          0.666666666666667_wp*d1_FluEx_dx_0_nxp5m4m1nyp5p0k+&
          0.666666666666667_wp*d1_FluEx_dx_0_nxp5m4p1nyp5p0k-&
          0.08333333333333333_wp*d1_FluEx_dx_0_nxp5m4p2nyp5p0k

d1_FluEx_dx_0_nxp5m4nyp5p0k = d1_FluEx_dx_0_nxp5m4nyp5p0k*param_float(1)

d1_FluEx_dy_0_nxp5m4nyp5p0p0k = q(nx+5-4,ny+5+0+0,indvars(3))*(q(nx+5-4,ny+5+0+0,indvars(1))*q(nx+5-4,ny+5+0+0,indvars(4))+&
                    0.375_wp*(8.0_wp*q(nx+5-4,ny+5+0+0,indvars(1))/(3.0_wp-&
                    q(nx+5-4,ny+5+0+0,indvars(1)))*param_float(1 + 5)*(((q(nx+5-4,ny+5+0+0,indvars(4))-&
                    0.5_wp*(q(nx+5-4,ny+5+0+0,indvars(2))*q(nx+5-4,ny+5+0+0,indvars(2))+&
                    q(nx+5-4,ny+5+0+0,indvars(3))*q(nx+5-4,ny+5+0+0,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5-4,ny+5+0+0,indvars(1)))-&
                    3.0_wp*q(nx+5-4,ny+5+0+0,indvars(1))*q(nx+5-4,ny+5+0+0,indvars(1))))

d1_FluEx_dy_0_nxp5m4nyp5p0m1k = q(nx+5-4,ny+5+0-1,indvars(3))*(q(nx+5-4,ny+5+0-1,indvars(1))*q(nx+5-4,ny+5+0-1,indvars(4))+&
                    0.375_wp*(8.0_wp*q(nx+5-4,ny+5+0-1,indvars(1))/(3.0_wp-&
                    q(nx+5-4,ny+5+0-1,indvars(1)))*param_float(1 + 5)*(((q(nx+5-4,ny+5+0-1,indvars(4))-&
                    0.5_wp*(q(nx+5-4,ny+5+0-1,indvars(2))*q(nx+5-4,ny+5+0-1,indvars(2))+&
                    q(nx+5-4,ny+5+0-1,indvars(3))*q(nx+5-4,ny+5+0-1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5-4,ny+5+0-1,indvars(1)))-&
                    3.0_wp*q(nx+5-4,ny+5+0-1,indvars(1))*q(nx+5-4,ny+5+0-1,indvars(1))))

d1_FluEx_dy_0_nxp5m4nyp5p0m2k = q(nx+5-4,ny+5+0-2,indvars(3))*(q(nx+5-4,ny+5+0-2,indvars(1))*q(nx+5-4,ny+5+0-2,indvars(4))+&
                    0.375_wp*(8.0_wp*q(nx+5-4,ny+5+0-2,indvars(1))/(3.0_wp-&
                    q(nx+5-4,ny+5+0-2,indvars(1)))*param_float(1 + 5)*(((q(nx+5-4,ny+5+0-2,indvars(4))-&
                    0.5_wp*(q(nx+5-4,ny+5+0-2,indvars(2))*q(nx+5-4,ny+5+0-2,indvars(2))+&
                    q(nx+5-4,ny+5+0-2,indvars(3))*q(nx+5-4,ny+5+0-2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(nx+5-4,ny+5+0-2,indvars(1)))-&
                    3.0_wp*q(nx+5-4,ny+5+0-2,indvars(1))*q(nx+5-4,ny+5+0-2,indvars(1))))

d1_FluEx_dy_0_nxp5m4nyp5p0k = 1.5_wp*d1_FluEx_dy_0_nxp5m4nyp5p0p0k-&
          2.0_wp*d1_FluEx_dy_0_nxp5m4nyp5p0m1k+&
          0.5_wp*d1_FluEx_dy_0_nxp5m4nyp5p0m2k

d1_FluEx_dy_0_nxp5m4nyp5p0k = d1_FluEx_dy_0_nxp5m4nyp5p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 4 0 None d(rho et)/dt **********
!                                                           
!***********************************************************


rhs(nx+5-4,ny+5+0,indvars(4)) =   -  ( d1_FluEx_dx_0_nxp5m4nyp5p0k+d1_FluEx_dy_0_nxp5m4nyp5p0k ) 

