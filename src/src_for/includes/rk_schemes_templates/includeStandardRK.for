 subroutine rk3_stepk(param_float,hlo,nrk,neq,neqst,ind,nx,ny,nz,sizeblck,&
                      nbc,bc,&
                      nfacei,nfacej,nfacek,&
                      nedgeij,nedgejk,nedgeik,&
                      q,q1,q2,rhs,qst,& 
                      qface_i , qface_j ,qface_k,&
                      qedge_ij, qedge_jk,qedge_ik)


#include "dtypes.h"
#include "param_fort.h"

  
  integer :: i,j,k
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk
  integer :: ibc

  real(wp), intent(in)    :: param_float(*)

  integer, intent(in) :: nx,ny,nz,hlo,nrk,neq,neqst
  integer, intent(in) :: nfacei,nfacej,nfacek,nedgeij,nedgejk,nedgeik
  integer, intent(in) :: ind(1:neq+neqst)
  integer, intent(in) :: sizeblck(3)
  

  integer, intent(in) :: bc(nbc),nbc
  real(wp) :: one_over_rho
  
#include "includeRK_coeffs.f90"
#include "includeRHS_globVar_rk3.f90"
#include "includeqbc_varrk.f90"
  
  integer :: idloop(6),idarray(6),idrhs(6),indvars(neq),indvarsst(neqst)
  integer :: nvar_f(3),nvar_e(3)



  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 


nvar_f(1) = nfacei
nvar_f(2) = nfacej
nvar_f(3) = nfacek

nvar_e(1) = nedgeij
nvar_e(2) = nedgejk
nvar_e(3) = nedgeik



!$OMP PARALLEL DEFAULT(SHARED) private(idloop)


      
  idrhs(1) = 1 
  idrhs(2) = nx

  idrhs(3) = 1 
  idrhs(4) = ny

  idrhs(5) = 1 
  idrhs(6) = nz
      
  idarray(1) = 1 -hlo
  idarray(2) = nx+hlo

  idarray(3) = 1 -hlo
  idarray(4) = ny+hlo

  idarray(5) = 1 -hlo
  idarray(6) = nz+hlo

size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2)  
  do bk=1,nz,size_bk  
    do bj=1,ny,size_bj 
      do bi=1,nx,size_bi 
    
idloop(6) = min( bk+size_bk, nz+1)-1
idloop(4) = min( bj+size_bj, ny+1)-1
idloop(2) = min( bi+size_bi, nx+1)-1

idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 

call cmprhs(param_float,ind,idloop,idarray,neq,neqst,sizeblck,q,qst,&
rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

#include "include_bcrhs.f90"
 
!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2)       
  do bk=idrhs(5),idrhs(6),size_bk  
    do bj=idrhs(3),idrhs(4),size_bj 
      do bi=idrhs(1),idrhs(2),size_bi 
 
idloop(6) = min( bk+size_bk, idrhs(6)+1)-1
idloop(4) = min( bj+size_bj, idrhs(4)+1)-1
idloop(2) = min( bi+size_bi, idrhs(2)+1)-1
    
idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 


if (nrk == 1) then

#include "LOOP_BEGIN"

#include "primitive_to_conservative.f90"

#include "LOOP_END"

endif  

#include "LOOP_BEGIN"

#include "includeRK3.f90"       
     
#include "LOOP_END"

#include "LOOP_BEGIN"

#include "includeRK3update.f90"       
     
#include "LOOP_END"

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO NOWAIT



#include "include_bcq.f90"

! print*,'rk3_stepk',q(:,ny+hlo,1)


!$OMP END PARALLEL



 end subroutine rk3_stepk
