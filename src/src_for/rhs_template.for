subroutine cmprhs(param_float,ind,idloop,idarray,neq,neqst,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

use iso_c_binding, only: c_float, c_double, c_int
implicit none

#include "dtypes.h"
#include "param_fort.h"

integer, intent(in) :: neq,neqst
integer, intent(in), dimension(6)   :: idloop, idarray
integer, intent(in), dimension(neq+neqst) :: ind
integer, intent(in), dimension(3) :: sizeblck
integer, intent(in) :: nvar_f(3),nvar_e(3)

real(wp), intent(in)    :: param_float(*)

#include "includeRHS_globVar.f90"

#include "includeqbc_var.f90"

! LOCAL VARIABLES

#include "includeRHS_locVar.f90"

integer :: i,j,k,nq
integer :: bi,bj,bk
integer :: biend,bjend,bkend
integer :: size_bi,size_bj,size_bk

integer :: indvars(neq),indvarsst(neqst)


  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 

#include "include_RHSLoops.f90"

end subroutine cmprhs

subroutine cmpstored(param_float,ind,idloop,idarray,neq,neqst,sizeblck,q,qst,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

use iso_c_binding, only: c_float, c_double, c_int
implicit none

#include "dtypes.h"
#include "param_fort.h"

integer, intent(in) :: neq,neqst
integer, intent(in), dimension(6)   :: idloop, idarray
integer, intent(in), dimension(neq+neqst) :: ind
integer, intent(in), dimension(3) :: sizeblck
integer, intent(in) :: nvar_f(3),nvar_e(3)

real(wp), intent(in)    :: param_float(*)

#include "includeRHS_globVarStored.f90"

#include "includeqbc_var.f90"

! LOCAL VARIABLES

#include "includeRHS_locVarStored.f90"

integer :: i,j,k,nq
integer :: bi,bj,bk
integer :: biend,bjend,bkend
integer :: size_bi,size_bj,size_bk

integer :: indvars(neq),indvarsst(neqst)


  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 

#include "include_StoredVar.f90"

end subroutine cmpstored

subroutine cmpstoredstatic(param_float,ind,idloop,idarray,neq,neqst,sizeblck,q,qst,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

use iso_c_binding, only: c_float, c_double, c_int
implicit none

#include "dtypes.h"
#include "param_fort.h"

integer, intent(in) :: neq,neqst
integer, intent(in), dimension(6)   :: idloop, idarray
integer, intent(in), dimension(neq+neqst) :: ind
integer, intent(in), dimension(3) :: sizeblck
integer, intent(in) :: nvar_f(3),nvar_e(3)

real(wp), intent(in)    :: param_float(*)

#include "includeRHS_globVarStored.f90"

#include "includeqbc_var.f90"

! LOCAL VARIABLES

#include "includeRHS_locVarStoredStatic.f90"

integer :: i,j,k,nq
integer :: bi,bj,bk
integer :: biend,bjend,bkend
integer :: size_bi,size_bj,size_bk

integer :: indvars(neq),indvarsst(neqst)


  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 

#include "include_StoredVarStatic.f90"

end subroutine cmpstoredstatic
