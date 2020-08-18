

!***********************************************************
!                                                           
! Start building layers for BC : None j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: None 2 None ************************************
!                                                           
!***********************************************************


 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for layer None 2 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u]_1x+[rho*v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluRx_dx_0_im21m5p2k = q(i-2,1-5+2,indvars(1))*q(i-2,1-5+2,indvars(2))

d1_FluRx_dx_0_im11m5p2k = q(i-1,1-5+2,indvars(1))*q(i-1,1-5+2,indvars(2))

d1_FluRx_dx_0_ip11m5p2k = q(i+1,1-5+2,indvars(1))*q(i+1,1-5+2,indvars(2))

d1_FluRx_dx_0_ip21m5p2k = q(i+2,1-5+2,indvars(1))*q(i+2,1-5+2,indvars(2))

d1_FluRx_dx_0_i1m5p2k = 0.08333333333333333_wp*d1_FluRx_dx_0_im21m5p2k-&
          0.666666666666667_wp*d1_FluRx_dx_0_im11m5p2k+&
          0.666666666666667_wp*d1_FluRx_dx_0_ip11m5p2k-&
          0.08333333333333333_wp*d1_FluRx_dx_0_ip21m5p2k

d1_FluRx_dx_0_i1m5p2k = d1_FluRx_dx_0_i1m5p2k*param_float(1)

d1_FluRx_dy_0_i1m5p2m2k = q(i,1-5+2-2,indvars(1))*q(i,1-5+2-2,indvars(3))

d1_FluRx_dy_0_i1m5p2m1k = q(i,1-5+2-1,indvars(1))*q(i,1-5+2-1,indvars(3))

d1_FluRx_dy_0_i1m5p2p1k = q(i,1-5+2+1,indvars(1))*q(i,1-5+2+1,indvars(3))

d1_FluRx_dy_0_i1m5p2p2k = q(i,1-5+2+2,indvars(1))*q(i,1-5+2+2,indvars(3))

d1_FluRx_dy_0_i1m5p2k = 0.08333333333333333_wp*d1_FluRx_dy_0_i1m5p2m2k-&
          0.666666666666667_wp*d1_FluRx_dy_0_i1m5p2m1k+&
          0.666666666666667_wp*d1_FluRx_dy_0_i1m5p2p1k-&
          0.08333333333333333_wp*d1_FluRx_dy_0_i1m5p2p2k

d1_FluRx_dy_0_i1m5p2k = d1_FluRx_dy_0_i1m5p2k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 2 None d(rho)/dt **********
!                                                           
!***********************************************************


rhs(i,1-5+2,indvars(1)) =   -  ( d1_FluRx_dx_0_i1m5p2k+d1_FluRx_dy_0_i1m5p2k ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 2 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u*u+0.375_wp*p]_1x+[rho*v*u]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluXx_dx_0_im21m5p2k = q(i-2,1-5+2,indvars(1))*q(i-2,1-5+2,indvars(2))*q(i-2,1-5+2,indvars(2))+&
                    0.375_wp*(8.0_wp*q(i-2,1-5+2,indvars(1))/(3.0_wp-&
                    q(i-2,1-5+2,indvars(1)))*param_float(1 + 5)*(((q(i-2,1-5+2,indvars(4))-&
                    0.5_wp*(q(i-2,1-5+2,indvars(2))*q(i-2,1-5+2,indvars(2))+&
                    q(i-2,1-5+2,indvars(3))*q(i-2,1-5+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i-2,1-5+2,indvars(1)))-&
                    3.0_wp*q(i-2,1-5+2,indvars(1))*q(i-2,1-5+2,indvars(1)))

d1_FluXx_dx_0_im11m5p2k = q(i-1,1-5+2,indvars(1))*q(i-1,1-5+2,indvars(2))*q(i-1,1-5+2,indvars(2))+&
                    0.375_wp*(8.0_wp*q(i-1,1-5+2,indvars(1))/(3.0_wp-&
                    q(i-1,1-5+2,indvars(1)))*param_float(1 + 5)*(((q(i-1,1-5+2,indvars(4))-&
                    0.5_wp*(q(i-1,1-5+2,indvars(2))*q(i-1,1-5+2,indvars(2))+&
                    q(i-1,1-5+2,indvars(3))*q(i-1,1-5+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i-1,1-5+2,indvars(1)))-&
                    3.0_wp*q(i-1,1-5+2,indvars(1))*q(i-1,1-5+2,indvars(1)))

d1_FluXx_dx_0_ip11m5p2k = q(i+1,1-5+2,indvars(1))*q(i+1,1-5+2,indvars(2))*q(i+1,1-5+2,indvars(2))+&
                    0.375_wp*(8.0_wp*q(i+1,1-5+2,indvars(1))/(3.0_wp-&
                    q(i+1,1-5+2,indvars(1)))*param_float(1 + 5)*(((q(i+1,1-5+2,indvars(4))-&
                    0.5_wp*(q(i+1,1-5+2,indvars(2))*q(i+1,1-5+2,indvars(2))+&
                    q(i+1,1-5+2,indvars(3))*q(i+1,1-5+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i+1,1-5+2,indvars(1)))-&
                    3.0_wp*q(i+1,1-5+2,indvars(1))*q(i+1,1-5+2,indvars(1)))

d1_FluXx_dx_0_ip21m5p2k = q(i+2,1-5+2,indvars(1))*q(i+2,1-5+2,indvars(2))*q(i+2,1-5+2,indvars(2))+&
                    0.375_wp*(8.0_wp*q(i+2,1-5+2,indvars(1))/(3.0_wp-&
                    q(i+2,1-5+2,indvars(1)))*param_float(1 + 5)*(((q(i+2,1-5+2,indvars(4))-&
                    0.5_wp*(q(i+2,1-5+2,indvars(2))*q(i+2,1-5+2,indvars(2))+&
                    q(i+2,1-5+2,indvars(3))*q(i+2,1-5+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i+2,1-5+2,indvars(1)))-&
                    3.0_wp*q(i+2,1-5+2,indvars(1))*q(i+2,1-5+2,indvars(1)))

d1_FluXx_dx_0_i1m5p2k = 0.08333333333333333_wp*d1_FluXx_dx_0_im21m5p2k-&
          0.666666666666667_wp*d1_FluXx_dx_0_im11m5p2k+&
          0.666666666666667_wp*d1_FluXx_dx_0_ip11m5p2k-&
          0.08333333333333333_wp*d1_FluXx_dx_0_ip21m5p2k

d1_FluXx_dx_0_i1m5p2k = d1_FluXx_dx_0_i1m5p2k*param_float(1)

d1_FluXx_dy_0_i1m5p2m2k = q(i,1-5+2-2,indvars(1))*q(i,1-5+2-2,indvars(3))*q(i,1-5+2-2,indvars(2))

d1_FluXx_dy_0_i1m5p2m1k = q(i,1-5+2-1,indvars(1))*q(i,1-5+2-1,indvars(3))*q(i,1-5+2-1,indvars(2))

d1_FluXx_dy_0_i1m5p2p1k = q(i,1-5+2+1,indvars(1))*q(i,1-5+2+1,indvars(3))*q(i,1-5+2+1,indvars(2))

d1_FluXx_dy_0_i1m5p2p2k = q(i,1-5+2+2,indvars(1))*q(i,1-5+2+2,indvars(3))*q(i,1-5+2+2,indvars(2))

d1_FluXx_dy_0_i1m5p2k = 0.08333333333333333_wp*d1_FluXx_dy_0_i1m5p2m2k-&
          0.666666666666667_wp*d1_FluXx_dy_0_i1m5p2m1k+&
          0.666666666666667_wp*d1_FluXx_dy_0_i1m5p2p1k-&
          0.08333333333333333_wp*d1_FluXx_dy_0_i1m5p2p2k

d1_FluXx_dy_0_i1m5p2k = d1_FluXx_dy_0_i1m5p2k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 2 None d(rho u)/dt ********
!                                                           
!***********************************************************


rhs(i,1-5+2,indvars(2)) =   -  ( d1_FluXx_dx_0_i1m5p2k+d1_FluXx_dy_0_i1m5p2k ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 2 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u*v]_1x+[rho*v*v+0.375_wp*p]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluYx_dx_0_im21m5p2k = q(i-2,1-5+2,indvars(1))*q(i-2,1-5+2,indvars(2))*q(i-2,1-5+2,indvars(3))

d1_FluYx_dx_0_im11m5p2k = q(i-1,1-5+2,indvars(1))*q(i-1,1-5+2,indvars(2))*q(i-1,1-5+2,indvars(3))

d1_FluYx_dx_0_ip11m5p2k = q(i+1,1-5+2,indvars(1))*q(i+1,1-5+2,indvars(2))*q(i+1,1-5+2,indvars(3))

d1_FluYx_dx_0_ip21m5p2k = q(i+2,1-5+2,indvars(1))*q(i+2,1-5+2,indvars(2))*q(i+2,1-5+2,indvars(3))

d1_FluYx_dx_0_i1m5p2k = 0.08333333333333333_wp*d1_FluYx_dx_0_im21m5p2k-&
          0.666666666666667_wp*d1_FluYx_dx_0_im11m5p2k+&
          0.666666666666667_wp*d1_FluYx_dx_0_ip11m5p2k-&
          0.08333333333333333_wp*d1_FluYx_dx_0_ip21m5p2k

d1_FluYx_dx_0_i1m5p2k = d1_FluYx_dx_0_i1m5p2k*param_float(1)

d1_FluYx_dy_0_i1m5p2m2k = q(i,1-5+2-2,indvars(1))*q(i,1-5+2-2,indvars(3))*q(i,1-5+2-2,indvars(3))+&
                    0.375_wp*(8.0_wp*q(i,1-5+2-2,indvars(1))/(3.0_wp-&
                    q(i,1-5+2-2,indvars(1)))*param_float(1 + 5)*(((q(i,1-5+2-2,indvars(4))-&
                    0.5_wp*(q(i,1-5+2-2,indvars(2))*q(i,1-5+2-2,indvars(2))+&
                    q(i,1-5+2-2,indvars(3))*q(i,1-5+2-2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,1-5+2-2,indvars(1)))-&
                    3.0_wp*q(i,1-5+2-2,indvars(1))*q(i,1-5+2-2,indvars(1)))

d1_FluYx_dy_0_i1m5p2m1k = q(i,1-5+2-1,indvars(1))*q(i,1-5+2-1,indvars(3))*q(i,1-5+2-1,indvars(3))+&
                    0.375_wp*(8.0_wp*q(i,1-5+2-1,indvars(1))/(3.0_wp-&
                    q(i,1-5+2-1,indvars(1)))*param_float(1 + 5)*(((q(i,1-5+2-1,indvars(4))-&
                    0.5_wp*(q(i,1-5+2-1,indvars(2))*q(i,1-5+2-1,indvars(2))+&
                    q(i,1-5+2-1,indvars(3))*q(i,1-5+2-1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,1-5+2-1,indvars(1)))-&
                    3.0_wp*q(i,1-5+2-1,indvars(1))*q(i,1-5+2-1,indvars(1)))

d1_FluYx_dy_0_i1m5p2p1k = q(i,1-5+2+1,indvars(1))*q(i,1-5+2+1,indvars(3))*q(i,1-5+2+1,indvars(3))+&
                    0.375_wp*(8.0_wp*q(i,1-5+2+1,indvars(1))/(3.0_wp-&
                    q(i,1-5+2+1,indvars(1)))*param_float(1 + 5)*(((q(i,1-5+2+1,indvars(4))-&
                    0.5_wp*(q(i,1-5+2+1,indvars(2))*q(i,1-5+2+1,indvars(2))+&
                    q(i,1-5+2+1,indvars(3))*q(i,1-5+2+1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,1-5+2+1,indvars(1)))-&
                    3.0_wp*q(i,1-5+2+1,indvars(1))*q(i,1-5+2+1,indvars(1)))

d1_FluYx_dy_0_i1m5p2p2k = q(i,1-5+2+2,indvars(1))*q(i,1-5+2+2,indvars(3))*q(i,1-5+2+2,indvars(3))+&
                    0.375_wp*(8.0_wp*q(i,1-5+2+2,indvars(1))/(3.0_wp-&
                    q(i,1-5+2+2,indvars(1)))*param_float(1 + 5)*(((q(i,1-5+2+2,indvars(4))-&
                    0.5_wp*(q(i,1-5+2+2,indvars(2))*q(i,1-5+2+2,indvars(2))+&
                    q(i,1-5+2+2,indvars(3))*q(i,1-5+2+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,1-5+2+2,indvars(1)))-&
                    3.0_wp*q(i,1-5+2+2,indvars(1))*q(i,1-5+2+2,indvars(1)))

d1_FluYx_dy_0_i1m5p2k = 0.08333333333333333_wp*d1_FluYx_dy_0_i1m5p2m2k-&
          0.666666666666667_wp*d1_FluYx_dy_0_i1m5p2m1k+&
          0.666666666666667_wp*d1_FluYx_dy_0_i1m5p2p1k-&
          0.08333333333333333_wp*d1_FluYx_dy_0_i1m5p2p2k

d1_FluYx_dy_0_i1m5p2k = d1_FluYx_dy_0_i1m5p2k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 2 None d(rho v)/dt ********
!                                                           
!***********************************************************


rhs(i,1-5+2,indvars(3)) =   -  ( d1_FluYx_dx_0_i1m5p2k+d1_FluYx_dy_0_i1m5p2k ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 2 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u*(rho*et+0.375_wp*p)]_1x+[v*(rho*et+0.375_wp*p)]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluEx_dx_0_im21m5p2k = q(i-2,1-5+2,indvars(2))*(q(i-2,1-5+2,indvars(1))*q(i-2,1-5+2,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i-2,1-5+2,indvars(1))/(3.0_wp-&
                    q(i-2,1-5+2,indvars(1)))*param_float(1 + 5)*(((q(i-2,1-5+2,indvars(4))-&
                    0.5_wp*(q(i-2,1-5+2,indvars(2))*q(i-2,1-5+2,indvars(2))+&
                    q(i-2,1-5+2,indvars(3))*q(i-2,1-5+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i-2,1-5+2,indvars(1)))-&
                    3.0_wp*q(i-2,1-5+2,indvars(1))*q(i-2,1-5+2,indvars(1))))

d1_FluEx_dx_0_im11m5p2k = q(i-1,1-5+2,indvars(2))*(q(i-1,1-5+2,indvars(1))*q(i-1,1-5+2,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i-1,1-5+2,indvars(1))/(3.0_wp-&
                    q(i-1,1-5+2,indvars(1)))*param_float(1 + 5)*(((q(i-1,1-5+2,indvars(4))-&
                    0.5_wp*(q(i-1,1-5+2,indvars(2))*q(i-1,1-5+2,indvars(2))+&
                    q(i-1,1-5+2,indvars(3))*q(i-1,1-5+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i-1,1-5+2,indvars(1)))-&
                    3.0_wp*q(i-1,1-5+2,indvars(1))*q(i-1,1-5+2,indvars(1))))

d1_FluEx_dx_0_ip11m5p2k = q(i+1,1-5+2,indvars(2))*(q(i+1,1-5+2,indvars(1))*q(i+1,1-5+2,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i+1,1-5+2,indvars(1))/(3.0_wp-&
                    q(i+1,1-5+2,indvars(1)))*param_float(1 + 5)*(((q(i+1,1-5+2,indvars(4))-&
                    0.5_wp*(q(i+1,1-5+2,indvars(2))*q(i+1,1-5+2,indvars(2))+&
                    q(i+1,1-5+2,indvars(3))*q(i+1,1-5+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i+1,1-5+2,indvars(1)))-&
                    3.0_wp*q(i+1,1-5+2,indvars(1))*q(i+1,1-5+2,indvars(1))))

d1_FluEx_dx_0_ip21m5p2k = q(i+2,1-5+2,indvars(2))*(q(i+2,1-5+2,indvars(1))*q(i+2,1-5+2,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i+2,1-5+2,indvars(1))/(3.0_wp-&
                    q(i+2,1-5+2,indvars(1)))*param_float(1 + 5)*(((q(i+2,1-5+2,indvars(4))-&
                    0.5_wp*(q(i+2,1-5+2,indvars(2))*q(i+2,1-5+2,indvars(2))+&
                    q(i+2,1-5+2,indvars(3))*q(i+2,1-5+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i+2,1-5+2,indvars(1)))-&
                    3.0_wp*q(i+2,1-5+2,indvars(1))*q(i+2,1-5+2,indvars(1))))

d1_FluEx_dx_0_i1m5p2k = 0.08333333333333333_wp*d1_FluEx_dx_0_im21m5p2k-&
          0.666666666666667_wp*d1_FluEx_dx_0_im11m5p2k+&
          0.666666666666667_wp*d1_FluEx_dx_0_ip11m5p2k-&
          0.08333333333333333_wp*d1_FluEx_dx_0_ip21m5p2k

d1_FluEx_dx_0_i1m5p2k = d1_FluEx_dx_0_i1m5p2k*param_float(1)

d1_FluEx_dy_0_i1m5p2m2k = q(i,1-5+2-2,indvars(3))*(q(i,1-5+2-2,indvars(1))*q(i,1-5+2-2,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i,1-5+2-2,indvars(1))/(3.0_wp-&
                    q(i,1-5+2-2,indvars(1)))*param_float(1 + 5)*(((q(i,1-5+2-2,indvars(4))-&
                    0.5_wp*(q(i,1-5+2-2,indvars(2))*q(i,1-5+2-2,indvars(2))+&
                    q(i,1-5+2-2,indvars(3))*q(i,1-5+2-2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,1-5+2-2,indvars(1)))-&
                    3.0_wp*q(i,1-5+2-2,indvars(1))*q(i,1-5+2-2,indvars(1))))

d1_FluEx_dy_0_i1m5p2m1k = q(i,1-5+2-1,indvars(3))*(q(i,1-5+2-1,indvars(1))*q(i,1-5+2-1,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i,1-5+2-1,indvars(1))/(3.0_wp-&
                    q(i,1-5+2-1,indvars(1)))*param_float(1 + 5)*(((q(i,1-5+2-1,indvars(4))-&
                    0.5_wp*(q(i,1-5+2-1,indvars(2))*q(i,1-5+2-1,indvars(2))+&
                    q(i,1-5+2-1,indvars(3))*q(i,1-5+2-1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,1-5+2-1,indvars(1)))-&
                    3.0_wp*q(i,1-5+2-1,indvars(1))*q(i,1-5+2-1,indvars(1))))

d1_FluEx_dy_0_i1m5p2p1k = q(i,1-5+2+1,indvars(3))*(q(i,1-5+2+1,indvars(1))*q(i,1-5+2+1,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i,1-5+2+1,indvars(1))/(3.0_wp-&
                    q(i,1-5+2+1,indvars(1)))*param_float(1 + 5)*(((q(i,1-5+2+1,indvars(4))-&
                    0.5_wp*(q(i,1-5+2+1,indvars(2))*q(i,1-5+2+1,indvars(2))+&
                    q(i,1-5+2+1,indvars(3))*q(i,1-5+2+1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,1-5+2+1,indvars(1)))-&
                    3.0_wp*q(i,1-5+2+1,indvars(1))*q(i,1-5+2+1,indvars(1))))

d1_FluEx_dy_0_i1m5p2p2k = q(i,1-5+2+2,indvars(3))*(q(i,1-5+2+2,indvars(1))*q(i,1-5+2+2,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i,1-5+2+2,indvars(1))/(3.0_wp-&
                    q(i,1-5+2+2,indvars(1)))*param_float(1 + 5)*(((q(i,1-5+2+2,indvars(4))-&
                    0.5_wp*(q(i,1-5+2+2,indvars(2))*q(i,1-5+2+2,indvars(2))+&
                    q(i,1-5+2+2,indvars(3))*q(i,1-5+2+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,1-5+2+2,indvars(1)))-&
                    3.0_wp*q(i,1-5+2+2,indvars(1))*q(i,1-5+2+2,indvars(1))))

d1_FluEx_dy_0_i1m5p2k = 0.08333333333333333_wp*d1_FluEx_dy_0_i1m5p2m2k-&
          0.666666666666667_wp*d1_FluEx_dy_0_i1m5p2m1k+&
          0.666666666666667_wp*d1_FluEx_dy_0_i1m5p2p1k-&
          0.08333333333333333_wp*d1_FluEx_dy_0_i1m5p2p2k

d1_FluEx_dy_0_i1m5p2k = d1_FluEx_dy_0_i1m5p2k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 2 None d(rho et)/dt *******
!                                                           
!***********************************************************


rhs(i,1-5+2,indvars(4)) =   -  ( d1_FluEx_dx_0_i1m5p2k+d1_FluEx_dy_0_i1m5p2k ) 

   enddo


!***********************************************************
!                                                           
! Start building layers for BC : None j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: None 2 None ************************************
!                                                           
!***********************************************************


 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for layer None 2 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 0.0_wp
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 2 None d(rho)/dt **********
!                                                           
!***********************************************************


rhs(i,1-5+2,indvars(1)) = rhs(i,1-5+2,indvars(1))  -  ( 0.0_wp ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 2 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! -mub*([divV]_1x)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluXv_dx_0_im21m5p2k = qst(i-2,1-5+2,indvarsst(1))

d1_FluXv_dx_0_im11m5p2k = qst(i-1,1-5+2,indvarsst(1))

d1_FluXv_dx_0_ip11m5p2k = qst(i+1,1-5+2,indvarsst(1))

d1_FluXv_dx_0_ip21m5p2k = qst(i+2,1-5+2,indvarsst(1))

d1_FluXv_dx_0_i1m5p2k = 0.08333333333333333_wp*d1_FluXv_dx_0_im21m5p2k-&
          0.666666666666667_wp*d1_FluXv_dx_0_im11m5p2k+&
          0.666666666666667_wp*d1_FluXv_dx_0_ip11m5p2k-&
          0.08333333333333333_wp*d1_FluXv_dx_0_ip21m5p2k

d1_FluXv_dx_0_i1m5p2k = d1_FluXv_dx_0_i1m5p2k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer None 2 None d(rho u)/dt ********
!                                                           
!***********************************************************


rhs(i,1-5+2,indvars(2)) = rhs(i,1-5+2,indvars(2))  -  ( -param_float(2 + 5)*(d1_FluXv_dx_0_i1m5p2k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 2 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! -mub*([divV]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluYv_dy_0_i1m5p2m2k = qst(i,1-5+2-2,indvarsst(1))

d1_FluYv_dy_0_i1m5p2m1k = qst(i,1-5+2-1,indvarsst(1))

d1_FluYv_dy_0_i1m5p2p1k = qst(i,1-5+2+1,indvarsst(1))

d1_FluYv_dy_0_i1m5p2p2k = qst(i,1-5+2+2,indvarsst(1))

d1_FluYv_dy_0_i1m5p2k = 0.08333333333333333_wp*d1_FluYv_dy_0_i1m5p2m2k-&
          0.666666666666667_wp*d1_FluYv_dy_0_i1m5p2m1k+&
          0.666666666666667_wp*d1_FluYv_dy_0_i1m5p2p1k-&
          0.08333333333333333_wp*d1_FluYv_dy_0_i1m5p2p2k

d1_FluYv_dy_0_i1m5p2k = d1_FluYv_dy_0_i1m5p2k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 2 None d(rho v)/dt ********
!                                                           
!***********************************************************


rhs(i,1-5+2,indvars(3)) = rhs(i,1-5+2,indvars(3))  -  ( -param_float(2 + 5)*(d1_FluYv_dy_0_i1m5p2k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 2 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! -mub*([u*divV]_1x+[v*divV]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluEv_dx_0_im21m5p2k = q(i-2,1-5+2,indvars(2))*qst(i-2,1-5+2,indvarsst(1))

d1_FluEv_dx_0_im11m5p2k = q(i-1,1-5+2,indvars(2))*qst(i-1,1-5+2,indvarsst(1))

d1_FluEv_dx_0_ip11m5p2k = q(i+1,1-5+2,indvars(2))*qst(i+1,1-5+2,indvarsst(1))

d1_FluEv_dx_0_ip21m5p2k = q(i+2,1-5+2,indvars(2))*qst(i+2,1-5+2,indvarsst(1))

d1_FluEv_dx_0_i1m5p2k = 0.08333333333333333_wp*d1_FluEv_dx_0_im21m5p2k-&
          0.666666666666667_wp*d1_FluEv_dx_0_im11m5p2k+&
          0.666666666666667_wp*d1_FluEv_dx_0_ip11m5p2k-&
          0.08333333333333333_wp*d1_FluEv_dx_0_ip21m5p2k

d1_FluEv_dx_0_i1m5p2k = d1_FluEv_dx_0_i1m5p2k*param_float(1)

d1_FluEv_dy_0_i1m5p2m2k = q(i,1-5+2-2,indvars(3))*qst(i,1-5+2-2,indvarsst(1))

d1_FluEv_dy_0_i1m5p2m1k = q(i,1-5+2-1,indvars(3))*qst(i,1-5+2-1,indvarsst(1))

d1_FluEv_dy_0_i1m5p2p1k = q(i,1-5+2+1,indvars(3))*qst(i,1-5+2+1,indvarsst(1))

d1_FluEv_dy_0_i1m5p2p2k = q(i,1-5+2+2,indvars(3))*qst(i,1-5+2+2,indvarsst(1))

d1_FluEv_dy_0_i1m5p2k = 0.08333333333333333_wp*d1_FluEv_dy_0_i1m5p2m2k-&
          0.666666666666667_wp*d1_FluEv_dy_0_i1m5p2m1k+&
          0.666666666666667_wp*d1_FluEv_dy_0_i1m5p2p1k-&
          0.08333333333333333_wp*d1_FluEv_dy_0_i1m5p2p2k

d1_FluEv_dy_0_i1m5p2k = d1_FluEv_dy_0_i1m5p2k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 2 None d(rho et)/dt *******
!                                                           
!***********************************************************


rhs(i,1-5+2,indvars(4)) = rhs(i,1-5+2,indvars(4))  -  ( -param_float(2 + 5)*(d1_FluEv_dx_0_i1m5p2k+&
                    d1_FluEv_dy_0_i1m5p2k) ) 

   enddo
