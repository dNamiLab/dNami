

!***********************************************************
!                                                           
! Start building layers for BC : i1 None None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 4 None None ************************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 


!***********************************************************
!                                                           
! building source terms in RHS for layer 4 None None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u]_1x+[rho*v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluRx_dx_0_1m5p4m2jk = q(1-5+4-2,j,indvars(1))*q(1-5+4-2,j,indvars(2))

d1_FluRx_dx_0_1m5p4m1jk = q(1-5+4-1,j,indvars(1))*q(1-5+4-1,j,indvars(2))

d1_FluRx_dx_0_1m5p4p1jk = q(1-5+4+1,j,indvars(1))*q(1-5+4+1,j,indvars(2))

d1_FluRx_dx_0_1m5p4p2jk = q(1-5+4+2,j,indvars(1))*q(1-5+4+2,j,indvars(2))

d1_FluRx_dx_0_1m5p4jk = 0.08333333333333333_wp*d1_FluRx_dx_0_1m5p4m2jk-&
          0.666666666666667_wp*d1_FluRx_dx_0_1m5p4m1jk+&
          0.666666666666667_wp*d1_FluRx_dx_0_1m5p4p1jk-&
          0.08333333333333333_wp*d1_FluRx_dx_0_1m5p4p2jk

d1_FluRx_dx_0_1m5p4jk = d1_FluRx_dx_0_1m5p4jk*param_float(1)

d1_FluRx_dy_0_1m5p4jm2k = q(1-5+4,j-2,indvars(1))*q(1-5+4,j-2,indvars(3))

d1_FluRx_dy_0_1m5p4jm1k = q(1-5+4,j-1,indvars(1))*q(1-5+4,j-1,indvars(3))

d1_FluRx_dy_0_1m5p4jp1k = q(1-5+4,j+1,indvars(1))*q(1-5+4,j+1,indvars(3))

d1_FluRx_dy_0_1m5p4jp2k = q(1-5+4,j+2,indvars(1))*q(1-5+4,j+2,indvars(3))

d1_FluRx_dy_0_1m5p4jk = 0.08333333333333333_wp*d1_FluRx_dy_0_1m5p4jm2k-&
          0.666666666666667_wp*d1_FluRx_dy_0_1m5p4jm1k+&
          0.666666666666667_wp*d1_FluRx_dy_0_1m5p4jp1k-&
          0.08333333333333333_wp*d1_FluRx_dy_0_1m5p4jp2k

d1_FluRx_dy_0_1m5p4jk = d1_FluRx_dy_0_1m5p4jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 4 None None d(rho)/dt **********
!                                                           
!***********************************************************


rhs(1-5+4,j,indvars(1)) =   -  ( d1_FluRx_dx_0_1m5p4jk+d1_FluRx_dy_0_1m5p4jk ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 4 None None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u*u+0.375_wp*p]_1x+[rho*v*u]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluXx_dx_0_1m5p4m2jk = q(1-5+4-2,j,indvars(1))*q(1-5+4-2,j,indvars(2))*q(1-5+4-2,j,indvars(2))+&
                    0.375_wp*(8.0_wp*q(1-5+4-2,j,indvars(1))/(3.0_wp-&
                    q(1-5+4-2,j,indvars(1)))*param_float(1 + 5)*(((q(1-5+4-2,j,indvars(4))-&
                    0.5_wp*(q(1-5+4-2,j,indvars(2))*q(1-5+4-2,j,indvars(2))+&
                    q(1-5+4-2,j,indvars(3))*q(1-5+4-2,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4-2,j,indvars(1)))-&
                    3.0_wp*q(1-5+4-2,j,indvars(1))*q(1-5+4-2,j,indvars(1)))

d1_FluXx_dx_0_1m5p4m1jk = q(1-5+4-1,j,indvars(1))*q(1-5+4-1,j,indvars(2))*q(1-5+4-1,j,indvars(2))+&
                    0.375_wp*(8.0_wp*q(1-5+4-1,j,indvars(1))/(3.0_wp-&
                    q(1-5+4-1,j,indvars(1)))*param_float(1 + 5)*(((q(1-5+4-1,j,indvars(4))-&
                    0.5_wp*(q(1-5+4-1,j,indvars(2))*q(1-5+4-1,j,indvars(2))+&
                    q(1-5+4-1,j,indvars(3))*q(1-5+4-1,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4-1,j,indvars(1)))-&
                    3.0_wp*q(1-5+4-1,j,indvars(1))*q(1-5+4-1,j,indvars(1)))

d1_FluXx_dx_0_1m5p4p1jk = q(1-5+4+1,j,indvars(1))*q(1-5+4+1,j,indvars(2))*q(1-5+4+1,j,indvars(2))+&
                    0.375_wp*(8.0_wp*q(1-5+4+1,j,indvars(1))/(3.0_wp-&
                    q(1-5+4+1,j,indvars(1)))*param_float(1 + 5)*(((q(1-5+4+1,j,indvars(4))-&
                    0.5_wp*(q(1-5+4+1,j,indvars(2))*q(1-5+4+1,j,indvars(2))+&
                    q(1-5+4+1,j,indvars(3))*q(1-5+4+1,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4+1,j,indvars(1)))-&
                    3.0_wp*q(1-5+4+1,j,indvars(1))*q(1-5+4+1,j,indvars(1)))

d1_FluXx_dx_0_1m5p4p2jk = q(1-5+4+2,j,indvars(1))*q(1-5+4+2,j,indvars(2))*q(1-5+4+2,j,indvars(2))+&
                    0.375_wp*(8.0_wp*q(1-5+4+2,j,indvars(1))/(3.0_wp-&
                    q(1-5+4+2,j,indvars(1)))*param_float(1 + 5)*(((q(1-5+4+2,j,indvars(4))-&
                    0.5_wp*(q(1-5+4+2,j,indvars(2))*q(1-5+4+2,j,indvars(2))+&
                    q(1-5+4+2,j,indvars(3))*q(1-5+4+2,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4+2,j,indvars(1)))-&
                    3.0_wp*q(1-5+4+2,j,indvars(1))*q(1-5+4+2,j,indvars(1)))

d1_FluXx_dx_0_1m5p4jk = 0.08333333333333333_wp*d1_FluXx_dx_0_1m5p4m2jk-&
          0.666666666666667_wp*d1_FluXx_dx_0_1m5p4m1jk+&
          0.666666666666667_wp*d1_FluXx_dx_0_1m5p4p1jk-&
          0.08333333333333333_wp*d1_FluXx_dx_0_1m5p4p2jk

d1_FluXx_dx_0_1m5p4jk = d1_FluXx_dx_0_1m5p4jk*param_float(1)

d1_FluXx_dy_0_1m5p4jm2k = q(1-5+4,j-2,indvars(1))*q(1-5+4,j-2,indvars(3))*q(1-5+4,j-2,indvars(2))

d1_FluXx_dy_0_1m5p4jm1k = q(1-5+4,j-1,indvars(1))*q(1-5+4,j-1,indvars(3))*q(1-5+4,j-1,indvars(2))

d1_FluXx_dy_0_1m5p4jp1k = q(1-5+4,j+1,indvars(1))*q(1-5+4,j+1,indvars(3))*q(1-5+4,j+1,indvars(2))

d1_FluXx_dy_0_1m5p4jp2k = q(1-5+4,j+2,indvars(1))*q(1-5+4,j+2,indvars(3))*q(1-5+4,j+2,indvars(2))

d1_FluXx_dy_0_1m5p4jk = 0.08333333333333333_wp*d1_FluXx_dy_0_1m5p4jm2k-&
          0.666666666666667_wp*d1_FluXx_dy_0_1m5p4jm1k+&
          0.666666666666667_wp*d1_FluXx_dy_0_1m5p4jp1k-&
          0.08333333333333333_wp*d1_FluXx_dy_0_1m5p4jp2k

d1_FluXx_dy_0_1m5p4jk = d1_FluXx_dy_0_1m5p4jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 4 None None d(rho u)/dt ********
!                                                           
!***********************************************************


rhs(1-5+4,j,indvars(2)) =   -  ( d1_FluXx_dx_0_1m5p4jk+d1_FluXx_dy_0_1m5p4jk ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 4 None None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u*v]_1x+[rho*v*v+0.375_wp*p]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluYx_dx_0_1m5p4m2jk = q(1-5+4-2,j,indvars(1))*q(1-5+4-2,j,indvars(2))*q(1-5+4-2,j,indvars(3))

d1_FluYx_dx_0_1m5p4m1jk = q(1-5+4-1,j,indvars(1))*q(1-5+4-1,j,indvars(2))*q(1-5+4-1,j,indvars(3))

d1_FluYx_dx_0_1m5p4p1jk = q(1-5+4+1,j,indvars(1))*q(1-5+4+1,j,indvars(2))*q(1-5+4+1,j,indvars(3))

d1_FluYx_dx_0_1m5p4p2jk = q(1-5+4+2,j,indvars(1))*q(1-5+4+2,j,indvars(2))*q(1-5+4+2,j,indvars(3))

d1_FluYx_dx_0_1m5p4jk = 0.08333333333333333_wp*d1_FluYx_dx_0_1m5p4m2jk-&
          0.666666666666667_wp*d1_FluYx_dx_0_1m5p4m1jk+&
          0.666666666666667_wp*d1_FluYx_dx_0_1m5p4p1jk-&
          0.08333333333333333_wp*d1_FluYx_dx_0_1m5p4p2jk

d1_FluYx_dx_0_1m5p4jk = d1_FluYx_dx_0_1m5p4jk*param_float(1)

d1_FluYx_dy_0_1m5p4jm2k = q(1-5+4,j-2,indvars(1))*q(1-5+4,j-2,indvars(3))*q(1-5+4,j-2,indvars(3))+&
                    0.375_wp*(8.0_wp*q(1-5+4,j-2,indvars(1))/(3.0_wp-&
                    q(1-5+4,j-2,indvars(1)))*param_float(1 + 5)*(((q(1-5+4,j-2,indvars(4))-&
                    0.5_wp*(q(1-5+4,j-2,indvars(2))*q(1-5+4,j-2,indvars(2))+&
                    q(1-5+4,j-2,indvars(3))*q(1-5+4,j-2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4,j-2,indvars(1)))-&
                    3.0_wp*q(1-5+4,j-2,indvars(1))*q(1-5+4,j-2,indvars(1)))

d1_FluYx_dy_0_1m5p4jm1k = q(1-5+4,j-1,indvars(1))*q(1-5+4,j-1,indvars(3))*q(1-5+4,j-1,indvars(3))+&
                    0.375_wp*(8.0_wp*q(1-5+4,j-1,indvars(1))/(3.0_wp-&
                    q(1-5+4,j-1,indvars(1)))*param_float(1 + 5)*(((q(1-5+4,j-1,indvars(4))-&
                    0.5_wp*(q(1-5+4,j-1,indvars(2))*q(1-5+4,j-1,indvars(2))+&
                    q(1-5+4,j-1,indvars(3))*q(1-5+4,j-1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4,j-1,indvars(1)))-&
                    3.0_wp*q(1-5+4,j-1,indvars(1))*q(1-5+4,j-1,indvars(1)))

d1_FluYx_dy_0_1m5p4jp1k = q(1-5+4,j+1,indvars(1))*q(1-5+4,j+1,indvars(3))*q(1-5+4,j+1,indvars(3))+&
                    0.375_wp*(8.0_wp*q(1-5+4,j+1,indvars(1))/(3.0_wp-&
                    q(1-5+4,j+1,indvars(1)))*param_float(1 + 5)*(((q(1-5+4,j+1,indvars(4))-&
                    0.5_wp*(q(1-5+4,j+1,indvars(2))*q(1-5+4,j+1,indvars(2))+&
                    q(1-5+4,j+1,indvars(3))*q(1-5+4,j+1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4,j+1,indvars(1)))-&
                    3.0_wp*q(1-5+4,j+1,indvars(1))*q(1-5+4,j+1,indvars(1)))

d1_FluYx_dy_0_1m5p4jp2k = q(1-5+4,j+2,indvars(1))*q(1-5+4,j+2,indvars(3))*q(1-5+4,j+2,indvars(3))+&
                    0.375_wp*(8.0_wp*q(1-5+4,j+2,indvars(1))/(3.0_wp-&
                    q(1-5+4,j+2,indvars(1)))*param_float(1 + 5)*(((q(1-5+4,j+2,indvars(4))-&
                    0.5_wp*(q(1-5+4,j+2,indvars(2))*q(1-5+4,j+2,indvars(2))+&
                    q(1-5+4,j+2,indvars(3))*q(1-5+4,j+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4,j+2,indvars(1)))-&
                    3.0_wp*q(1-5+4,j+2,indvars(1))*q(1-5+4,j+2,indvars(1)))

d1_FluYx_dy_0_1m5p4jk = 0.08333333333333333_wp*d1_FluYx_dy_0_1m5p4jm2k-&
          0.666666666666667_wp*d1_FluYx_dy_0_1m5p4jm1k+&
          0.666666666666667_wp*d1_FluYx_dy_0_1m5p4jp1k-&
          0.08333333333333333_wp*d1_FluYx_dy_0_1m5p4jp2k

d1_FluYx_dy_0_1m5p4jk = d1_FluYx_dy_0_1m5p4jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 4 None None d(rho v)/dt ********
!                                                           
!***********************************************************


rhs(1-5+4,j,indvars(3)) =   -  ( d1_FluYx_dx_0_1m5p4jk+d1_FluYx_dy_0_1m5p4jk ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 4 None None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u*(rho*et+0.375_wp*p)]_1x+[v*(rho*et+0.375_wp*p)]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluEx_dx_0_1m5p4m2jk = q(1-5+4-2,j,indvars(2))*(q(1-5+4-2,j,indvars(1))*q(1-5+4-2,j,indvars(4))+&
                    0.375_wp*(8.0_wp*q(1-5+4-2,j,indvars(1))/(3.0_wp-&
                    q(1-5+4-2,j,indvars(1)))*param_float(1 + 5)*(((q(1-5+4-2,j,indvars(4))-&
                    0.5_wp*(q(1-5+4-2,j,indvars(2))*q(1-5+4-2,j,indvars(2))+&
                    q(1-5+4-2,j,indvars(3))*q(1-5+4-2,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4-2,j,indvars(1)))-&
                    3.0_wp*q(1-5+4-2,j,indvars(1))*q(1-5+4-2,j,indvars(1))))

d1_FluEx_dx_0_1m5p4m1jk = q(1-5+4-1,j,indvars(2))*(q(1-5+4-1,j,indvars(1))*q(1-5+4-1,j,indvars(4))+&
                    0.375_wp*(8.0_wp*q(1-5+4-1,j,indvars(1))/(3.0_wp-&
                    q(1-5+4-1,j,indvars(1)))*param_float(1 + 5)*(((q(1-5+4-1,j,indvars(4))-&
                    0.5_wp*(q(1-5+4-1,j,indvars(2))*q(1-5+4-1,j,indvars(2))+&
                    q(1-5+4-1,j,indvars(3))*q(1-5+4-1,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4-1,j,indvars(1)))-&
                    3.0_wp*q(1-5+4-1,j,indvars(1))*q(1-5+4-1,j,indvars(1))))

d1_FluEx_dx_0_1m5p4p1jk = q(1-5+4+1,j,indvars(2))*(q(1-5+4+1,j,indvars(1))*q(1-5+4+1,j,indvars(4))+&
                    0.375_wp*(8.0_wp*q(1-5+4+1,j,indvars(1))/(3.0_wp-&
                    q(1-5+4+1,j,indvars(1)))*param_float(1 + 5)*(((q(1-5+4+1,j,indvars(4))-&
                    0.5_wp*(q(1-5+4+1,j,indvars(2))*q(1-5+4+1,j,indvars(2))+&
                    q(1-5+4+1,j,indvars(3))*q(1-5+4+1,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4+1,j,indvars(1)))-&
                    3.0_wp*q(1-5+4+1,j,indvars(1))*q(1-5+4+1,j,indvars(1))))

d1_FluEx_dx_0_1m5p4p2jk = q(1-5+4+2,j,indvars(2))*(q(1-5+4+2,j,indvars(1))*q(1-5+4+2,j,indvars(4))+&
                    0.375_wp*(8.0_wp*q(1-5+4+2,j,indvars(1))/(3.0_wp-&
                    q(1-5+4+2,j,indvars(1)))*param_float(1 + 5)*(((q(1-5+4+2,j,indvars(4))-&
                    0.5_wp*(q(1-5+4+2,j,indvars(2))*q(1-5+4+2,j,indvars(2))+&
                    q(1-5+4+2,j,indvars(3))*q(1-5+4+2,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4+2,j,indvars(1)))-&
                    3.0_wp*q(1-5+4+2,j,indvars(1))*q(1-5+4+2,j,indvars(1))))

d1_FluEx_dx_0_1m5p4jk = 0.08333333333333333_wp*d1_FluEx_dx_0_1m5p4m2jk-&
          0.666666666666667_wp*d1_FluEx_dx_0_1m5p4m1jk+&
          0.666666666666667_wp*d1_FluEx_dx_0_1m5p4p1jk-&
          0.08333333333333333_wp*d1_FluEx_dx_0_1m5p4p2jk

d1_FluEx_dx_0_1m5p4jk = d1_FluEx_dx_0_1m5p4jk*param_float(1)

d1_FluEx_dy_0_1m5p4jm2k = q(1-5+4,j-2,indvars(3))*(q(1-5+4,j-2,indvars(1))*q(1-5+4,j-2,indvars(4))+&
                    0.375_wp*(8.0_wp*q(1-5+4,j-2,indvars(1))/(3.0_wp-&
                    q(1-5+4,j-2,indvars(1)))*param_float(1 + 5)*(((q(1-5+4,j-2,indvars(4))-&
                    0.5_wp*(q(1-5+4,j-2,indvars(2))*q(1-5+4,j-2,indvars(2))+&
                    q(1-5+4,j-2,indvars(3))*q(1-5+4,j-2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4,j-2,indvars(1)))-&
                    3.0_wp*q(1-5+4,j-2,indvars(1))*q(1-5+4,j-2,indvars(1))))

d1_FluEx_dy_0_1m5p4jm1k = q(1-5+4,j-1,indvars(3))*(q(1-5+4,j-1,indvars(1))*q(1-5+4,j-1,indvars(4))+&
                    0.375_wp*(8.0_wp*q(1-5+4,j-1,indvars(1))/(3.0_wp-&
                    q(1-5+4,j-1,indvars(1)))*param_float(1 + 5)*(((q(1-5+4,j-1,indvars(4))-&
                    0.5_wp*(q(1-5+4,j-1,indvars(2))*q(1-5+4,j-1,indvars(2))+&
                    q(1-5+4,j-1,indvars(3))*q(1-5+4,j-1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4,j-1,indvars(1)))-&
                    3.0_wp*q(1-5+4,j-1,indvars(1))*q(1-5+4,j-1,indvars(1))))

d1_FluEx_dy_0_1m5p4jp1k = q(1-5+4,j+1,indvars(3))*(q(1-5+4,j+1,indvars(1))*q(1-5+4,j+1,indvars(4))+&
                    0.375_wp*(8.0_wp*q(1-5+4,j+1,indvars(1))/(3.0_wp-&
                    q(1-5+4,j+1,indvars(1)))*param_float(1 + 5)*(((q(1-5+4,j+1,indvars(4))-&
                    0.5_wp*(q(1-5+4,j+1,indvars(2))*q(1-5+4,j+1,indvars(2))+&
                    q(1-5+4,j+1,indvars(3))*q(1-5+4,j+1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4,j+1,indvars(1)))-&
                    3.0_wp*q(1-5+4,j+1,indvars(1))*q(1-5+4,j+1,indvars(1))))

d1_FluEx_dy_0_1m5p4jp2k = q(1-5+4,j+2,indvars(3))*(q(1-5+4,j+2,indvars(1))*q(1-5+4,j+2,indvars(4))+&
                    0.375_wp*(8.0_wp*q(1-5+4,j+2,indvars(1))/(3.0_wp-&
                    q(1-5+4,j+2,indvars(1)))*param_float(1 + 5)*(((q(1-5+4,j+2,indvars(4))-&
                    0.5_wp*(q(1-5+4,j+2,indvars(2))*q(1-5+4,j+2,indvars(2))+&
                    q(1-5+4,j+2,indvars(3))*q(1-5+4,j+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(1-5+4,j+2,indvars(1)))-&
                    3.0_wp*q(1-5+4,j+2,indvars(1))*q(1-5+4,j+2,indvars(1))))

d1_FluEx_dy_0_1m5p4jk = 0.08333333333333333_wp*d1_FluEx_dy_0_1m5p4jm2k-&
          0.666666666666667_wp*d1_FluEx_dy_0_1m5p4jm1k+&
          0.666666666666667_wp*d1_FluEx_dy_0_1m5p4jp1k-&
          0.08333333333333333_wp*d1_FluEx_dy_0_1m5p4jp2k

d1_FluEx_dy_0_1m5p4jk = d1_FluEx_dy_0_1m5p4jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 4 None None d(rho et)/dt *******
!                                                           
!***********************************************************


rhs(1-5+4,j,indvars(4)) =   -  ( d1_FluEx_dx_0_1m5p4jk+d1_FluEx_dy_0_1m5p4jk ) 

     enddo


!***********************************************************
!                                                           
! Start building layers for BC : i1 None None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 4 None None ************************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 


!***********************************************************
!                                                           
! building source terms in RHS for layer 4 None None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 0.0_wp
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 4 None None d(rho)/dt **********
!                                                           
!***********************************************************


rhs(1-5+4,j,indvars(1)) = rhs(1-5+4,j,indvars(1))  -  ( 0.0_wp ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 4 None None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! -mub*([divV]_1x)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluXv_dx_0_1m5p4m2jk = qst(1-5+4-2,j,indvarsst(1))

d1_FluXv_dx_0_1m5p4m1jk = qst(1-5+4-1,j,indvarsst(1))

d1_FluXv_dx_0_1m5p4p1jk = qst(1-5+4+1,j,indvarsst(1))

d1_FluXv_dx_0_1m5p4p2jk = qst(1-5+4+2,j,indvarsst(1))

d1_FluXv_dx_0_1m5p4jk = 0.08333333333333333_wp*d1_FluXv_dx_0_1m5p4m2jk-&
          0.666666666666667_wp*d1_FluXv_dx_0_1m5p4m1jk+&
          0.666666666666667_wp*d1_FluXv_dx_0_1m5p4p1jk-&
          0.08333333333333333_wp*d1_FluXv_dx_0_1m5p4p2jk

d1_FluXv_dx_0_1m5p4jk = d1_FluXv_dx_0_1m5p4jk*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 4 None None d(rho u)/dt ********
!                                                           
!***********************************************************


rhs(1-5+4,j,indvars(2)) = rhs(1-5+4,j,indvars(2))  -  ( -param_float(2 + 5)*(d1_FluXv_dx_0_1m5p4jk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 4 None None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! -mub*([divV]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluYv_dy_0_1m5p4jm2k = qst(1-5+4,j-2,indvarsst(1))

d1_FluYv_dy_0_1m5p4jm1k = qst(1-5+4,j-1,indvarsst(1))

d1_FluYv_dy_0_1m5p4jp1k = qst(1-5+4,j+1,indvarsst(1))

d1_FluYv_dy_0_1m5p4jp2k = qst(1-5+4,j+2,indvarsst(1))

d1_FluYv_dy_0_1m5p4jk = 0.08333333333333333_wp*d1_FluYv_dy_0_1m5p4jm2k-&
          0.666666666666667_wp*d1_FluYv_dy_0_1m5p4jm1k+&
          0.666666666666667_wp*d1_FluYv_dy_0_1m5p4jp1k-&
          0.08333333333333333_wp*d1_FluYv_dy_0_1m5p4jp2k

d1_FluYv_dy_0_1m5p4jk = d1_FluYv_dy_0_1m5p4jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 4 None None d(rho v)/dt ********
!                                                           
!***********************************************************


rhs(1-5+4,j,indvars(3)) = rhs(1-5+4,j,indvars(3))  -  ( -param_float(2 + 5)*(d1_FluYv_dy_0_1m5p4jk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 4 None None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! -mub*([u*divV]_1x+[v*divV]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluEv_dx_0_1m5p4m2jk = q(1-5+4-2,j,indvars(2))*qst(1-5+4-2,j,indvarsst(1))

d1_FluEv_dx_0_1m5p4m1jk = q(1-5+4-1,j,indvars(2))*qst(1-5+4-1,j,indvarsst(1))

d1_FluEv_dx_0_1m5p4p1jk = q(1-5+4+1,j,indvars(2))*qst(1-5+4+1,j,indvarsst(1))

d1_FluEv_dx_0_1m5p4p2jk = q(1-5+4+2,j,indvars(2))*qst(1-5+4+2,j,indvarsst(1))

d1_FluEv_dx_0_1m5p4jk = 0.08333333333333333_wp*d1_FluEv_dx_0_1m5p4m2jk-&
          0.666666666666667_wp*d1_FluEv_dx_0_1m5p4m1jk+&
          0.666666666666667_wp*d1_FluEv_dx_0_1m5p4p1jk-&
          0.08333333333333333_wp*d1_FluEv_dx_0_1m5p4p2jk

d1_FluEv_dx_0_1m5p4jk = d1_FluEv_dx_0_1m5p4jk*param_float(1)

d1_FluEv_dy_0_1m5p4jm2k = q(1-5+4,j-2,indvars(3))*qst(1-5+4,j-2,indvarsst(1))

d1_FluEv_dy_0_1m5p4jm1k = q(1-5+4,j-1,indvars(3))*qst(1-5+4,j-1,indvarsst(1))

d1_FluEv_dy_0_1m5p4jp1k = q(1-5+4,j+1,indvars(3))*qst(1-5+4,j+1,indvarsst(1))

d1_FluEv_dy_0_1m5p4jp2k = q(1-5+4,j+2,indvars(3))*qst(1-5+4,j+2,indvarsst(1))

d1_FluEv_dy_0_1m5p4jk = 0.08333333333333333_wp*d1_FluEv_dy_0_1m5p4jm2k-&
          0.666666666666667_wp*d1_FluEv_dy_0_1m5p4jm1k+&
          0.666666666666667_wp*d1_FluEv_dy_0_1m5p4jp1k-&
          0.08333333333333333_wp*d1_FluEv_dy_0_1m5p4jp2k

d1_FluEv_dy_0_1m5p4jk = d1_FluEv_dy_0_1m5p4jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 4 None None d(rho et)/dt *******
!                                                           
!***********************************************************


rhs(1-5+4,j,indvars(4)) = rhs(1-5+4,j,indvars(4))  -  ( -param_float(2 + 5)*(d1_FluEv_dx_0_1m5p4jk+&
                    d1_FluEv_dy_0_1m5p4jk) ) 

     enddo
