!====================================================================================================
!
! General Boundary Conditions: enforce physical BC and compute Eqns with Boundary Scheme
!                              Both steady and unsteady BC can be generated (i.e. on vector q or rhs)
!     
!====================================================================================================

subroutine phyBCNAME(param_float,hlo,neq,neqst,ind,idloop,idarray,nx,ny,nz,sizeblck,q,qst,q1,q2,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

  implicit none

#include "dtypes.h"
#include "param_fort.h"

  real(wp), intent(in)    :: param_float(*)

  integer, intent(in) :: neq,neqst
  integer, intent(in) :: nx,ny,nz,hlo
  integer, intent(in) :: ind(1:neq+neqst)
  integer, intent(in) :: idarray(6)
  integer, intent(inout) :: idloop(6) 
  integer, intent(in) :: sizeblck(3)
  integer, intent(in) :: nvar_f(3),nvar_e(3)
    
#include "includeRHS_globVar_rk3.f90"

#include "includeqbc_var.f90"

! LOCAL VARIABLES
 
  integer :: i,j,k
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk
  integer :: indbc(6) 

!PHYBCLOCVAR_rhs

  integer :: indvars(neq),indvarsst(neqst)
  
  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 

!f2py intent(in)    :: q1,q2,qst,nx,ny,nz,rhs,nrk
!f2py intent(inout) :: q

      
size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)

!INDBC_RANGEi
!INDBC_RANGEj
!INDBC_RANGEk

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2)  
  do bk=indbc(5),indbc(6),size_bk  
    do bj=indbc(3),indbc(4),size_bj 
      do bi=indbc(1),indbc(2),size_bi 
    
idloop(6) = min( bk+size_bk, indbc(6)+1)-1
idloop(4) = min( bj+size_bj, indbc(4)+1)-1
idloop(2) = min( bi+size_bi, indbc(2)+1)-1

idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 

!PHYBCLOOPS_rhs

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

 end subroutine phyBCNAME
