

!***********************************************************
!                                                           
! Start computing stored variables *************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building stored variable divV ****************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [u]_1x+[v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_divV_dx_0_im2jk = q(i-2,j,indvars(2))

d1_divV_dx_0_im1jk = q(i-1,j,indvars(2))

d1_divV_dx_0_ip1jk = q(i+1,j,indvars(2))

d1_divV_dx_0_ip2jk = q(i+2,j,indvars(2))

d1_divV_dx_0_ijk = 0.08333333333333333_wp*d1_divV_dx_0_im2jk-&
          0.666666666666667_wp*d1_divV_dx_0_im1jk+&
          0.666666666666667_wp*d1_divV_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_divV_dx_0_ip2jk

d1_divV_dx_0_ijk = d1_divV_dx_0_ijk*param_float(1)

d1_divV_dy_0_ijm2k = q(i,j-2,indvars(3))

d1_divV_dy_0_ijm1k = q(i,j-1,indvars(3))

d1_divV_dy_0_ijp1k = q(i,j+1,indvars(3))

d1_divV_dy_0_ijp2k = q(i,j+2,indvars(3))

d1_divV_dy_0_ijk = 0.08333333333333333_wp*d1_divV_dy_0_ijm2k-&
          0.666666666666667_wp*d1_divV_dy_0_ijm1k+&
          0.666666666666667_wp*d1_divV_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_divV_dy_0_ijp2k

d1_divV_dy_0_ijk = d1_divV_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update stored variables    [u]_1x  + [v]_1y  *************
!                                                           
!***********************************************************


qst(i,j,indvarsst(1)) =  d1_divV_dx_0_ijk+d1_divV_dy_0_ijk

     enddo
   enddo
