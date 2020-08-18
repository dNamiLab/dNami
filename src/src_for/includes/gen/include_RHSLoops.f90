

!***********************************************************
!                                                           
! Start building RHS with source terms (2D) ****************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for d(rho)/dt ***************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u]_1x+[rho*v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluRx_dx_0_im2jk = q(i-2,j,indvars(1))*q(i-2,j,indvars(2))

d1_FluRx_dx_0_im1jk = q(i-1,j,indvars(1))*q(i-1,j,indvars(2))

d1_FluRx_dx_0_ip1jk = q(i+1,j,indvars(1))*q(i+1,j,indvars(2))

d1_FluRx_dx_0_ip2jk = q(i+2,j,indvars(1))*q(i+2,j,indvars(2))

d1_FluRx_dx_0_ijk = 0.08333333333333333_wp*d1_FluRx_dx_0_im2jk-&
          0.666666666666667_wp*d1_FluRx_dx_0_im1jk+&
          0.666666666666667_wp*d1_FluRx_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_FluRx_dx_0_ip2jk

d1_FluRx_dx_0_ijk = d1_FluRx_dx_0_ijk*param_float(1)

d1_FluRx_dy_0_ijm2k = q(i,j-2,indvars(1))*q(i,j-2,indvars(3))

d1_FluRx_dy_0_ijm1k = q(i,j-1,indvars(1))*q(i,j-1,indvars(3))

d1_FluRx_dy_0_ijp1k = q(i,j+1,indvars(1))*q(i,j+1,indvars(3))

d1_FluRx_dy_0_ijp2k = q(i,j+2,indvars(1))*q(i,j+2,indvars(3))

d1_FluRx_dy_0_ijk = 0.08333333333333333_wp*d1_FluRx_dy_0_ijm2k-&
          0.666666666666667_wp*d1_FluRx_dy_0_ijm1k+&
          0.666666666666667_wp*d1_FluRx_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_FluRx_dy_0_ijp2k

d1_FluRx_dy_0_ijk = d1_FluRx_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho)/dt ***************************
!                                                           
!***********************************************************


rhs(i,j,indvars(1)) =   -  ( d1_FluRx_dx_0_ijk+d1_FluRx_dy_0_ijk ) 



!***********************************************************
!                                                           
! building source terms in RHS for d(rho u)/dt *************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u*u+0.375_wp*p]_1x+[rho*v*u]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluXx_dx_0_im2jk = q(i-2,j,indvars(1))*q(i-2,j,indvars(2))*q(i-2,j,indvars(2))+&
                    0.375_wp*(8.0_wp*q(i-2,j,indvars(1))/(3.0_wp-&
                    q(i-2,j,indvars(1)))*param_float(1 + 5)*(((q(i-2,j,indvars(4))-&
                    0.5_wp*(q(i-2,j,indvars(2))*q(i-2,j,indvars(2))+&
                    q(i-2,j,indvars(3))*q(i-2,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i-2,j,indvars(1)))-&
                    3.0_wp*q(i-2,j,indvars(1))*q(i-2,j,indvars(1)))

d1_FluXx_dx_0_im1jk = q(i-1,j,indvars(1))*q(i-1,j,indvars(2))*q(i-1,j,indvars(2))+&
                    0.375_wp*(8.0_wp*q(i-1,j,indvars(1))/(3.0_wp-&
                    q(i-1,j,indvars(1)))*param_float(1 + 5)*(((q(i-1,j,indvars(4))-&
                    0.5_wp*(q(i-1,j,indvars(2))*q(i-1,j,indvars(2))+&
                    q(i-1,j,indvars(3))*q(i-1,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i-1,j,indvars(1)))-&
                    3.0_wp*q(i-1,j,indvars(1))*q(i-1,j,indvars(1)))

d1_FluXx_dx_0_ip1jk = q(i+1,j,indvars(1))*q(i+1,j,indvars(2))*q(i+1,j,indvars(2))+&
                    0.375_wp*(8.0_wp*q(i+1,j,indvars(1))/(3.0_wp-&
                    q(i+1,j,indvars(1)))*param_float(1 + 5)*(((q(i+1,j,indvars(4))-&
                    0.5_wp*(q(i+1,j,indvars(2))*q(i+1,j,indvars(2))+&
                    q(i+1,j,indvars(3))*q(i+1,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i+1,j,indvars(1)))-&
                    3.0_wp*q(i+1,j,indvars(1))*q(i+1,j,indvars(1)))

d1_FluXx_dx_0_ip2jk = q(i+2,j,indvars(1))*q(i+2,j,indvars(2))*q(i+2,j,indvars(2))+&
                    0.375_wp*(8.0_wp*q(i+2,j,indvars(1))/(3.0_wp-&
                    q(i+2,j,indvars(1)))*param_float(1 + 5)*(((q(i+2,j,indvars(4))-&
                    0.5_wp*(q(i+2,j,indvars(2))*q(i+2,j,indvars(2))+&
                    q(i+2,j,indvars(3))*q(i+2,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i+2,j,indvars(1)))-&
                    3.0_wp*q(i+2,j,indvars(1))*q(i+2,j,indvars(1)))

d1_FluXx_dx_0_ijk = 0.08333333333333333_wp*d1_FluXx_dx_0_im2jk-&
          0.666666666666667_wp*d1_FluXx_dx_0_im1jk+&
          0.666666666666667_wp*d1_FluXx_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_FluXx_dx_0_ip2jk

d1_FluXx_dx_0_ijk = d1_FluXx_dx_0_ijk*param_float(1)

d1_FluXx_dy_0_ijm2k = q(i,j-2,indvars(1))*q(i,j-2,indvars(3))*q(i,j-2,indvars(2))

d1_FluXx_dy_0_ijm1k = q(i,j-1,indvars(1))*q(i,j-1,indvars(3))*q(i,j-1,indvars(2))

d1_FluXx_dy_0_ijp1k = q(i,j+1,indvars(1))*q(i,j+1,indvars(3))*q(i,j+1,indvars(2))

d1_FluXx_dy_0_ijp2k = q(i,j+2,indvars(1))*q(i,j+2,indvars(3))*q(i,j+2,indvars(2))

d1_FluXx_dy_0_ijk = 0.08333333333333333_wp*d1_FluXx_dy_0_ijm2k-&
          0.666666666666667_wp*d1_FluXx_dy_0_ijm1k+&
          0.666666666666667_wp*d1_FluXx_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_FluXx_dy_0_ijp2k

d1_FluXx_dy_0_ijk = d1_FluXx_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho u)/dt *************************
!                                                           
!***********************************************************


rhs(i,j,indvars(2)) =   -  ( d1_FluXx_dx_0_ijk+d1_FluXx_dy_0_ijk ) 



!***********************************************************
!                                                           
! building source terms in RHS for d(rho v)/dt *************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [rho*u*v]_1x+[rho*v*v+0.375_wp*p]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluYx_dx_0_im2jk = q(i-2,j,indvars(1))*q(i-2,j,indvars(2))*q(i-2,j,indvars(3))

d1_FluYx_dx_0_im1jk = q(i-1,j,indvars(1))*q(i-1,j,indvars(2))*q(i-1,j,indvars(3))

d1_FluYx_dx_0_ip1jk = q(i+1,j,indvars(1))*q(i+1,j,indvars(2))*q(i+1,j,indvars(3))

d1_FluYx_dx_0_ip2jk = q(i+2,j,indvars(1))*q(i+2,j,indvars(2))*q(i+2,j,indvars(3))

d1_FluYx_dx_0_ijk = 0.08333333333333333_wp*d1_FluYx_dx_0_im2jk-&
          0.666666666666667_wp*d1_FluYx_dx_0_im1jk+&
          0.666666666666667_wp*d1_FluYx_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_FluYx_dx_0_ip2jk

d1_FluYx_dx_0_ijk = d1_FluYx_dx_0_ijk*param_float(1)

d1_FluYx_dy_0_ijm2k = q(i,j-2,indvars(1))*q(i,j-2,indvars(3))*q(i,j-2,indvars(3))+&
                    0.375_wp*(8.0_wp*q(i,j-2,indvars(1))/(3.0_wp-&
                    q(i,j-2,indvars(1)))*param_float(1 + 5)*(((q(i,j-2,indvars(4))-&
                    0.5_wp*(q(i,j-2,indvars(2))*q(i,j-2,indvars(2))+&
                    q(i,j-2,indvars(3))*q(i,j-2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,j-2,indvars(1)))-&
                    3.0_wp*q(i,j-2,indvars(1))*q(i,j-2,indvars(1)))

d1_FluYx_dy_0_ijm1k = q(i,j-1,indvars(1))*q(i,j-1,indvars(3))*q(i,j-1,indvars(3))+&
                    0.375_wp*(8.0_wp*q(i,j-1,indvars(1))/(3.0_wp-&
                    q(i,j-1,indvars(1)))*param_float(1 + 5)*(((q(i,j-1,indvars(4))-&
                    0.5_wp*(q(i,j-1,indvars(2))*q(i,j-1,indvars(2))+&
                    q(i,j-1,indvars(3))*q(i,j-1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,j-1,indvars(1)))-&
                    3.0_wp*q(i,j-1,indvars(1))*q(i,j-1,indvars(1)))

d1_FluYx_dy_0_ijp1k = q(i,j+1,indvars(1))*q(i,j+1,indvars(3))*q(i,j+1,indvars(3))+&
                    0.375_wp*(8.0_wp*q(i,j+1,indvars(1))/(3.0_wp-&
                    q(i,j+1,indvars(1)))*param_float(1 + 5)*(((q(i,j+1,indvars(4))-&
                    0.5_wp*(q(i,j+1,indvars(2))*q(i,j+1,indvars(2))+&
                    q(i,j+1,indvars(3))*q(i,j+1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,j+1,indvars(1)))-&
                    3.0_wp*q(i,j+1,indvars(1))*q(i,j+1,indvars(1)))

d1_FluYx_dy_0_ijp2k = q(i,j+2,indvars(1))*q(i,j+2,indvars(3))*q(i,j+2,indvars(3))+&
                    0.375_wp*(8.0_wp*q(i,j+2,indvars(1))/(3.0_wp-&
                    q(i,j+2,indvars(1)))*param_float(1 + 5)*(((q(i,j+2,indvars(4))-&
                    0.5_wp*(q(i,j+2,indvars(2))*q(i,j+2,indvars(2))+&
                    q(i,j+2,indvars(3))*q(i,j+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,j+2,indvars(1)))-&
                    3.0_wp*q(i,j+2,indvars(1))*q(i,j+2,indvars(1)))

d1_FluYx_dy_0_ijk = 0.08333333333333333_wp*d1_FluYx_dy_0_ijm2k-&
          0.666666666666667_wp*d1_FluYx_dy_0_ijm1k+&
          0.666666666666667_wp*d1_FluYx_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_FluYx_dy_0_ijp2k

d1_FluYx_dy_0_ijk = d1_FluYx_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho v)/dt *************************
!                                                           
!***********************************************************


rhs(i,j,indvars(3)) =   -  ( d1_FluYx_dx_0_ijk+d1_FluYx_dy_0_ijk ) 



!***********************************************************
!                                                           
! building source terms in RHS for d(rho et)/dt ************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u*(rho*et+0.375_wp*p)]_1x+[v*(rho*et+0.375_wp*p)]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluEx_dx_0_im2jk = q(i-2,j,indvars(2))*(q(i-2,j,indvars(1))*q(i-2,j,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i-2,j,indvars(1))/(3.0_wp-&
                    q(i-2,j,indvars(1)))*param_float(1 + 5)*(((q(i-2,j,indvars(4))-&
                    0.5_wp*(q(i-2,j,indvars(2))*q(i-2,j,indvars(2))+&
                    q(i-2,j,indvars(3))*q(i-2,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i-2,j,indvars(1)))-&
                    3.0_wp*q(i-2,j,indvars(1))*q(i-2,j,indvars(1))))

d1_FluEx_dx_0_im1jk = q(i-1,j,indvars(2))*(q(i-1,j,indvars(1))*q(i-1,j,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i-1,j,indvars(1))/(3.0_wp-&
                    q(i-1,j,indvars(1)))*param_float(1 + 5)*(((q(i-1,j,indvars(4))-&
                    0.5_wp*(q(i-1,j,indvars(2))*q(i-1,j,indvars(2))+&
                    q(i-1,j,indvars(3))*q(i-1,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i-1,j,indvars(1)))-&
                    3.0_wp*q(i-1,j,indvars(1))*q(i-1,j,indvars(1))))

d1_FluEx_dx_0_ip1jk = q(i+1,j,indvars(2))*(q(i+1,j,indvars(1))*q(i+1,j,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i+1,j,indvars(1))/(3.0_wp-&
                    q(i+1,j,indvars(1)))*param_float(1 + 5)*(((q(i+1,j,indvars(4))-&
                    0.5_wp*(q(i+1,j,indvars(2))*q(i+1,j,indvars(2))+&
                    q(i+1,j,indvars(3))*q(i+1,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i+1,j,indvars(1)))-&
                    3.0_wp*q(i+1,j,indvars(1))*q(i+1,j,indvars(1))))

d1_FluEx_dx_0_ip2jk = q(i+2,j,indvars(2))*(q(i+2,j,indvars(1))*q(i+2,j,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i+2,j,indvars(1))/(3.0_wp-&
                    q(i+2,j,indvars(1)))*param_float(1 + 5)*(((q(i+2,j,indvars(4))-&
                    0.5_wp*(q(i+2,j,indvars(2))*q(i+2,j,indvars(2))+&
                    q(i+2,j,indvars(3))*q(i+2,j,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i+2,j,indvars(1)))-&
                    3.0_wp*q(i+2,j,indvars(1))*q(i+2,j,indvars(1))))

d1_FluEx_dx_0_ijk = 0.08333333333333333_wp*d1_FluEx_dx_0_im2jk-&
          0.666666666666667_wp*d1_FluEx_dx_0_im1jk+&
          0.666666666666667_wp*d1_FluEx_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_FluEx_dx_0_ip2jk

d1_FluEx_dx_0_ijk = d1_FluEx_dx_0_ijk*param_float(1)

d1_FluEx_dy_0_ijm2k = q(i,j-2,indvars(3))*(q(i,j-2,indvars(1))*q(i,j-2,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i,j-2,indvars(1))/(3.0_wp-&
                    q(i,j-2,indvars(1)))*param_float(1 + 5)*(((q(i,j-2,indvars(4))-&
                    0.5_wp*(q(i,j-2,indvars(2))*q(i,j-2,indvars(2))+&
                    q(i,j-2,indvars(3))*q(i,j-2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,j-2,indvars(1)))-&
                    3.0_wp*q(i,j-2,indvars(1))*q(i,j-2,indvars(1))))

d1_FluEx_dy_0_ijm1k = q(i,j-1,indvars(3))*(q(i,j-1,indvars(1))*q(i,j-1,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i,j-1,indvars(1))/(3.0_wp-&
                    q(i,j-1,indvars(1)))*param_float(1 + 5)*(((q(i,j-1,indvars(4))-&
                    0.5_wp*(q(i,j-1,indvars(2))*q(i,j-1,indvars(2))+&
                    q(i,j-1,indvars(3))*q(i,j-1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,j-1,indvars(1)))-&
                    3.0_wp*q(i,j-1,indvars(1))*q(i,j-1,indvars(1))))

d1_FluEx_dy_0_ijp1k = q(i,j+1,indvars(3))*(q(i,j+1,indvars(1))*q(i,j+1,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i,j+1,indvars(1))/(3.0_wp-&
                    q(i,j+1,indvars(1)))*param_float(1 + 5)*(((q(i,j+1,indvars(4))-&
                    0.5_wp*(q(i,j+1,indvars(2))*q(i,j+1,indvars(2))+&
                    q(i,j+1,indvars(3))*q(i,j+1,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,j+1,indvars(1)))-&
                    3.0_wp*q(i,j+1,indvars(1))*q(i,j+1,indvars(1))))

d1_FluEx_dy_0_ijp2k = q(i,j+2,indvars(3))*(q(i,j+2,indvars(1))*q(i,j+2,indvars(4))+&
                    0.375_wp*(8.0_wp*q(i,j+2,indvars(1))/(3.0_wp-&
                    q(i,j+2,indvars(1)))*param_float(1 + 5)*(((q(i,j+2,indvars(4))-&
                    0.5_wp*(q(i,j+2,indvars(2))*q(i,j+2,indvars(2))+&
                    q(i,j+2,indvars(3))*q(i,j+2,indvars(3)))))+&
                    9.0_wp/8.0_wp*q(i,j+2,indvars(1)))-&
                    3.0_wp*q(i,j+2,indvars(1))*q(i,j+2,indvars(1))))

d1_FluEx_dy_0_ijk = 0.08333333333333333_wp*d1_FluEx_dy_0_ijm2k-&
          0.666666666666667_wp*d1_FluEx_dy_0_ijm1k+&
          0.666666666666667_wp*d1_FluEx_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_FluEx_dy_0_ijp2k

d1_FluEx_dy_0_ijk = d1_FluEx_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho et)/dt ************************
!                                                           
!***********************************************************


rhs(i,j,indvars(4)) =   -  ( d1_FluEx_dx_0_ijk+d1_FluEx_dy_0_ijk ) 

     enddo
   enddo


!***********************************************************
!                                                           
! Start building RHS with source terms (2D) ****************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for d(rho)/dt ***************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 0.0_wp
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update RHS terms for d(rho)/dt ***************************
!                                                           
!***********************************************************


rhs(i,j,indvars(1)) = rhs(i,j,indvars(1))  -  ( 0.0_wp ) 



!***********************************************************
!                                                           
! building source terms in RHS for d(rho u)/dt *************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! -mub*([divV]_1x)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluXv_dx_0_im2jk = qst(i-2,j,indvarsst(1))

d1_FluXv_dx_0_im1jk = qst(i-1,j,indvarsst(1))

d1_FluXv_dx_0_ip1jk = qst(i+1,j,indvarsst(1))

d1_FluXv_dx_0_ip2jk = qst(i+2,j,indvarsst(1))

d1_FluXv_dx_0_ijk = 0.08333333333333333_wp*d1_FluXv_dx_0_im2jk-&
          0.666666666666667_wp*d1_FluXv_dx_0_im1jk+&
          0.666666666666667_wp*d1_FluXv_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_FluXv_dx_0_ip2jk

d1_FluXv_dx_0_ijk = d1_FluXv_dx_0_ijk*param_float(1)



!***********************************************************
!                                                           
! Update RHS terms for d(rho u)/dt *************************
!                                                           
!***********************************************************


rhs(i,j,indvars(2)) = rhs(i,j,indvars(2))  -  ( -param_float(2 + 5)*(d1_FluXv_dx_0_ijk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for d(rho v)/dt *************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! -mub*([divV]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluYv_dy_0_ijm2k = qst(i,j-2,indvarsst(1))

d1_FluYv_dy_0_ijm1k = qst(i,j-1,indvarsst(1))

d1_FluYv_dy_0_ijp1k = qst(i,j+1,indvarsst(1))

d1_FluYv_dy_0_ijp2k = qst(i,j+2,indvarsst(1))

d1_FluYv_dy_0_ijk = 0.08333333333333333_wp*d1_FluYv_dy_0_ijm2k-&
          0.666666666666667_wp*d1_FluYv_dy_0_ijm1k+&
          0.666666666666667_wp*d1_FluYv_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_FluYv_dy_0_ijp2k

d1_FluYv_dy_0_ijk = d1_FluYv_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho v)/dt *************************
!                                                           
!***********************************************************


rhs(i,j,indvars(3)) = rhs(i,j,indvars(3))  -  ( -param_float(2 + 5)*(d1_FluYv_dy_0_ijk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for d(rho et)/dt ************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! -mub*([u*divV]_1x+[v*divV]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_FluEv_dx_0_im2jk = q(i-2,j,indvars(2))*qst(i-2,j,indvarsst(1))

d1_FluEv_dx_0_im1jk = q(i-1,j,indvars(2))*qst(i-1,j,indvarsst(1))

d1_FluEv_dx_0_ip1jk = q(i+1,j,indvars(2))*qst(i+1,j,indvarsst(1))

d1_FluEv_dx_0_ip2jk = q(i+2,j,indvars(2))*qst(i+2,j,indvarsst(1))

d1_FluEv_dx_0_ijk = 0.08333333333333333_wp*d1_FluEv_dx_0_im2jk-&
          0.666666666666667_wp*d1_FluEv_dx_0_im1jk+&
          0.666666666666667_wp*d1_FluEv_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_FluEv_dx_0_ip2jk

d1_FluEv_dx_0_ijk = d1_FluEv_dx_0_ijk*param_float(1)

d1_FluEv_dy_0_ijm2k = q(i,j-2,indvars(3))*qst(i,j-2,indvarsst(1))

d1_FluEv_dy_0_ijm1k = q(i,j-1,indvars(3))*qst(i,j-1,indvarsst(1))

d1_FluEv_dy_0_ijp1k = q(i,j+1,indvars(3))*qst(i,j+1,indvarsst(1))

d1_FluEv_dy_0_ijp2k = q(i,j+2,indvars(3))*qst(i,j+2,indvarsst(1))

d1_FluEv_dy_0_ijk = 0.08333333333333333_wp*d1_FluEv_dy_0_ijm2k-&
          0.666666666666667_wp*d1_FluEv_dy_0_ijm1k+&
          0.666666666666667_wp*d1_FluEv_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_FluEv_dy_0_ijp2k

d1_FluEv_dy_0_ijk = d1_FluEv_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho et)/dt ************************
!                                                           
!***********************************************************


rhs(i,j,indvars(4)) = rhs(i,j,indvars(4))  -  ( -param_float(2 + 5)*(d1_FluEv_dx_0_ijk+&
                    d1_FluEv_dy_0_ijk) ) 

     enddo
   enddo
