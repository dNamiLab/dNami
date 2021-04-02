subroutine init_numa(param_float,hlo,nrk,neq,neqst,ind,nx,ny,nz,sizeblck,q,q1,rhs,qst)


#include "dtypes.h"
#include "param_fort.h"

  
  integer :: i,j,k
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk

  real(wp), intent(in)    :: param_float(*)

  integer, intent(in) :: nx,ny,nz,hlo,nrk,neq,neqst
  integer, intent(in) :: ind(1:neq+neqst)
  integer, intent(in) :: sizeblck(3)
  
#include "includeRHS_globVar_init.f90"
  
  integer :: idloop(6),idarray(6),indvars(neq),indvarsst(neqst)


  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 


!$OMP PARALLEL DEFAULT(SHARED) private(idloop)

  idarray(1) = 1 -hlo
  idarray(2) = nx+hlo

  idarray(3) = 1 -hlo
  idarray(4) = ny+hlo

  idarray(5) = 1 -hlo
  idarray(6) = nz+hlo

size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)

!$OMP DO SCHEDULE(STATIC) 
  do bk=1,nz,size_bk  
    do bj=1,ny,size_bj 
      do bi=1,nx,size_bi 
    
idloop(6) = min( bk+size_bk, nz+1)-1
idloop(4) = min( bj+size_bj, ny+1)-1
idloop(2) = min( bi+size_bi, nx+1)-1

idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 


#include "LOOP_BEGIN"

#include "init_numa.f90"

#include "LOOP_END"

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO
!$OMP END PARALLEL

end subroutine init_numa
