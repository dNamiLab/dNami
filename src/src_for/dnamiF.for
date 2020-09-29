module dnami1
      contains
subroutine time_march(param_int,param_float,data_float) bind(c)
  use iso_c_binding, only: c_double, c_int
  implicit none  

#include "dtypes.h"
#include "param_fort.h" 

  real(c_double), intent(in)    :: param_float(*)
  real(c_double), intent(inout) ::  data_float(*)
  integer,  intent(in)    ::   param_int(*)

  ! LOCAL VARIABLES

  integer :: shiftstore
  integer :: shiftfacei,shiftfacej,shiftfacek
  integer :: shiftedgeij,shiftedgejk,shiftedgeik
  integer :: nbc
  integer :: nfacei,nfacej,nfacek
  integer :: nedgeij,nedgejk
  integer :: nvar_f(3),nvar_e(3)

  nbc       = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST)  )
  
  nvar_f(1) = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc      )
  nvar_f(2) = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 1  )
  nvar_f(3) = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 2  ) 

  nvar_e(1) = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 3      )
  nvar_e(2) = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 3 + 1  )
  nvar_e(3) = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 3 + 2  ) 
   
  nfacei    = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 6  )
  nfacej    = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 7  )
  nfacek    = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 8  )
  nedgeij   = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 9  )
  nedgejk   = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 10 )

  if ( param_int(iadrNVARST) .gt. 0 ) then
   shiftstore = param_int(iadrNVARST)*param_int(iadrNDIMTOT)
  else
   shiftstore = 1
  endif

  shiftfacei  = shiftstore  + nfacei
  shiftfacej  = shiftfacei  + nfacej
  shiftfacek  = shiftfacej  + nfacek
  shiftedgeij = shiftfacek  + nedgeij
  shiftedgejk = shiftedgeij + nedgejk

  call rk3_stepk(param_float                                                          ,&
                 param_int(iadrHLO), param_int(iadrNRK)                               ,&
                 param_int(iadrNVAR)                                                  ,&
                 param_int(iadrNVARST)                                                ,&
                 param_int(iadrVARS)                                                  ,&
                 param_int(iadrNX),param_int(iadrNY) , param_int(iadrNZ)              ,&
                 param_int(iadrCacheBLCK)                                             ,&
                 param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST))    ,&
                 param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1),&
                 nvar_f(1),nvar_f(2),nvar_f(3),&
                 nvar_e(1),nvar_e(2),nvar_e(3),&
                 data_float(1)                                                        ,& ! q
                 data_float(1+param_int(iadrNDIMTOT)  )                               ,& ! q1
                 data_float(1+param_int(iadrNDIMTOT)*2)                               ,& ! q2 
                 data_float(1+param_int(iadrNDIMTOT)*3)                               ,& ! rhs
                 data_float(1+param_int(iadrNDIMTOT)*4)                               ,& ! stored     (if any)
                 data_float(1+param_int(iadrNDIMTOT)*4 + shiftstore)                  ,& ! qbcface_i  (if any)
                 data_float(1+param_int(iadrNDIMTOT)*4 + shiftfacei)                  ,& ! qbcface_j  (if any)
                 data_float(1+param_int(iadrNDIMTOT)*4 + shiftfacej)                  ,& ! qbcface_k  (if any)
                 data_float(1+param_int(iadrNDIMTOT)*4 + shiftfacek)                  ,& ! qbcedge_ij (if any)
                 data_float(1+param_int(iadrNDIMTOT)*4 + shiftedgeij)                 ,& ! qbcedge_jk (if any)
                 data_float(1+param_int(iadrNDIMTOT)*4 + shiftedgejk)                 )  ! qbcedge_ik (if any)
      

  end subroutine time_march


  subroutine stored(param_int,param_float,data_float,type_st) bind(c)
  use iso_c_binding, only: c_double, c_int
  implicit none  

#include "dtypes.h"
#include "param_fort.h" 

  real(c_double), intent(in)    :: param_float(*)
  real(c_double), intent(inout) ::  data_float(*)
  integer,  intent(in)    ::   param_int(*)
  integer                 ::   type_st
 
 !f2py integer optional, intent(in)    :: type_st = 0

  integer :: shiftstore
  integer :: shiftfacei,shiftfacej,shiftfacek
  integer :: shiftedgeij,shiftedgejk,shiftedgeik
  integer :: nbc
  integer :: nfacei,nfacej,nfacek
  integer :: nedgeij,nedgejk
  integer :: nvar_f(3),nvar_e(3)

  nbc       = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST)  )
  
  nvar_f(1) = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc      )
  nvar_f(2) = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 1  )
  nvar_f(3) = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 2  ) 

  nvar_e(1) = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 3      )
  nvar_e(2) = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 3 + 1  )
  nvar_e(3) = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 3 + 2  ) 
   
  nfacei    = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 6  )
  nfacej    = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 7  )
  nfacek    = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 8  )
  nedgeij   = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 9  )
  nedgejk   = param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1 + nbc + 10 )

  if ( param_int(iadrNVARST) .gt. 0 ) then
   shiftstore = param_int(iadrNVARST)*param_int(iadrNDIMTOT)
  else
   shiftstore = 1
  endif

  shiftfacei  = shiftstore  + nfacei
  shiftfacej  = shiftfacei  + nfacej
  shiftfacek  = shiftfacej  + nfacek
  shiftedgeij = shiftfacek  + nedgeij
  shiftedgejk = shiftedgeij + nedgejk


  call cmp_stored(param_float                                                          ,&
                  type_st                                                              ,&
                  param_int(iadrHLO), param_int(iadrNRK)                               ,&
                  param_int(iadrNVAR)                                                  ,&
                  param_int(iadrNVARST)                                                ,&
                  param_int(iadrVARS)                                                  ,&
                  param_int(iadrNX),param_int(iadrNY) , param_int(iadrNZ)              ,&
                  param_int(iadrCacheBLCK)                                             ,&
                  param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST))    ,&
                  param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1),&                  
                  nvar_f(1),nvar_f(2),nvar_f(3),&
                  nvar_e(1),nvar_e(2),nvar_e(3),&
                  data_float(1)                                                        ,& ! q
                  data_float(1+param_int(iadrNDIMTOT)*3)                               ,& ! rhs
                  data_float(1+param_int(iadrNDIMTOT)*4)                               ,& ! stored (if any)
                  data_float(1+param_int(iadrNDIMTOT)*4 + shiftstore)                  ,& ! qbcface_i  (if any)
                  data_float(1+param_int(iadrNDIMTOT)*4 + shiftfacei)                  ,& ! qbcface_j  (if any)
                  data_float(1+param_int(iadrNDIMTOT)*4 + shiftfacej)                  ,& ! qbcface_k  (if any)
                  data_float(1+param_int(iadrNDIMTOT)*4 + shiftfacek)                  ,& ! qbcedge_ij (if any)
                  data_float(1+param_int(iadrNDIMTOT)*4 + shiftedgeij)                 ,& ! qbcedge_jk (if any)
                  data_float(1+param_int(iadrNDIMTOT)*4 + shiftedgejk)                 )  ! qbcedge_ik (if any)         


  end subroutine stored  

  subroutine cmp_stored(param_float,type_st,hlo,nrk,neq,neqst,ind,nx,ny,nz,sizeblck,&
                        nbc,bc,&
                        nfacei,nfacej,nfacek,&
                        nedgeij,nedgejk,nedgeik,&
                        q,rhs,qst,& 
                        qface_i , qface_j ,qface_k,&
                        qedge_ij, qedge_jk,qedge_ik)

  implicit none

#include "dtypes.h"
#include "param_fort.h"

  
  integer :: i,j,k
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk
  integer :: ibc

  real(wp), intent(in)    :: param_float(*)

  integer, intent(in) :: nx,ny,nz,hlo,nrk,neq,neqst,type_st
  integer, intent(in) :: nfacei,nfacej,nfacek,nedgeij,nedgejk,nedgeik
  integer, intent(in) :: ind(1:neq+neqst)
  integer, intent(in) :: sizeblck(3)
  integer, intent(in) :: bc(nbc),nbc

  real(wp),  dimension(1:3) :: rk1 = (/ 2.0_wp/ 3.0_wp, &
                                        5.0_wp/12.0_wp, &
                                        3.0_wp/ 5.0_wp /)
  real(wp),  dimension(1:3) :: rk2 = (/ 1.0_wp/ 4.0_wp, &
                                        3.0_wp/20.0_wp, &
                                        3.0_wp/ 5.0_wp /)

  real(wp) :: one_over_rho
  
#include "includeF_globVarStored.f90"

#include "includeqbc_varrk.f90"

  
  integer :: idloop(6),idarray(6),indvars(neq),indvarsst(neqst)
  integer :: nvar_f(3),nvar_e(3)

  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 


!f2py intent(in)    :: nx,ny,nz,nrk,q,nvar_f,nvar_e,neq,nbc
!f2py intent(inout) :: qst

nvar_f(1) = nfacei
nvar_f(2) = nfacej
nvar_f(3) = nfacek

nvar_e(1) = nedgeij
nvar_e(2) = nedgejk
nvar_e(3) = nedgeik

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

    SELECT CASE (type_st)

    CASE(0)

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

call cmpstored(param_float,ind,idloop,idarray,neq,neqst,sizeblck,q,qst,& 
nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO NOWAIT

#include "include_bcstored.f90"
    
    CASE(1)

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

call cmpstoredstatic(param_float,ind,idloop,idarray,neq,neqst,sizeblck,&
q,qst,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO NOWAIT

#include "include_bcstoredstatic.f90"

  END SELECT

!$OMP END PARALLEL

 end subroutine cmp_stored


 subroutine rk3_stepk(param_float,hlo,nrk,neq,neqst,ind,nx,ny,nz,sizeblck,&
                      nbc,bc,&
                      nfacei,nfacej,nfacek,&
                      nedgeij,nedgejk,nedgeik,&
                      q,q1,q2,rhs,qst,& 
                      qface_i , qface_j ,qface_k,&
                      qedge_ij, qedge_jk,qedge_ik)

  implicit none

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

  real(wp),  dimension(1:3) :: rk1 = (/ 2.0_wp/ 3.0_wp, &
                                        5.0_wp/12.0_wp, &
                                        3.0_wp/ 5.0_wp /)
  real(wp),  dimension(1:3) :: rk2 = (/ 1.0_wp/ 4.0_wp, &
                                        3.0_wp/20.0_wp, &
                                        3.0_wp/ 5.0_wp /)

  real(wp) :: one_over_rho
  
#include "includeRHS_globVar_rk3.f90"

#include "includeqbc_varrk.f90"
  
  integer :: idloop(6),idarray(6),idrhs(6),indvars(neq),indvarsst(neqst)
  integer :: nvar_f(3),nvar_e(3)



  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 

!f2py intent(in)    :: q1,q2,nx,ny,nz,rhs,nrk,qst,nvar_f,nvar_e,neq,nbc
!f2py intent(inout) :: q

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

subroutine filter(dir,param_int,param_float,data_float) bind(c)
  use iso_c_binding, only: c_double, c_int
implicit none  

#include "dtypes.h"
#include "param_fort.h" 

real(c_double), intent(in)    :: param_float(*)
real(c_double), intent(inout) ::  data_float(*)
integer,  intent(in)    :: param_int(*)

integer, intent(in)   :: dir

!f2py  intent(in)   :: param_int,param_float,dir
!f2py intent(inout) :: data_float    
 if     (dir == 1) then

    call flt_x(param_float                                                          ,&
               param_int(iadrHLO)                                                   ,&
               param_int(iadrNVAR)                                                  ,&
               param_int(iadrNX),param_int(iadrNY) , param_int(iadrNZ)              ,&
               param_int(iadrCacheBLCK)                                             ,& 
               param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) )   ,&
               param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1),&                             
               data_float(1)                                                        ,&
               data_float(1+param_int(iadrNDIMTOT)*2)                                ) ! q2 is used for filter working array
    
 elseif (dir == 2)then

    call flt_y(param_float                                                          ,&
               param_int(iadrHLO)                                                   ,&
               param_int(iadrNVAR)                                                  ,&
               param_int(iadrNX),param_int(iadrNY) , param_int(iadrNZ)              ,&
               param_int(iadrCacheBLCK)                                             ,& 
               param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) )   ,&
               param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1),&                                              
               data_float(1)                                                        ,&
               data_float(1+param_int(iadrNDIMTOT)*2)                                ) ! q2 is used for filter working array

 elseif (dir == 3)then


    call flt_z(param_float                                                          ,&
               param_int(iadrHLO)                                                   ,&
               param_int(iadrNVAR)                                                  ,&
               param_int(iadrNX),param_int(iadrNY) , param_int(iadrNZ)              ,&
               param_int(iadrCacheBLCK)                                             ,&
               param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) )   ,&
               param_int(iadrVARS + param_int(iadrNVAR) + param_int(iadrNVARST) + 1),&                                               
               data_float(1)                                                        ,&
               data_float(1+param_int(iadrNDIMTOT)*2)                                ) ! q2 is used for filter working array

 endif

 end subroutine filter


subroutine flt_x(param_float,hlo,neq,nx,ny,nz,sizeblck,nbc,bc,q,q2)

implicit none  

#include "dtypes.h"
#include "param_fort.h" 

  real(wp), intent(in)    :: param_float(*)
! LOCAL VARIABLES

  integer, intent(in) :: nx,ny,nz,hlo,neq  
  integer, intent(in) :: sizeblck(3)

  integer :: i,j,k  
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk
  integer :: ibc
  integer, intent(in) :: bc(nbc),nbc

#include "includeRHS_globVar_filter.f90"
  
  integer :: idloop(6),idarray(6),idflt(6)


!f2py intent(in)    :: nx,ny,nz,hlo,neq,indvars,param_float
!f2py intent(inout) :: q

!$OMP PARALLEL DEFAULT(SHARED) private(idloop)

idloop(1) = 1 
idloop(2) = nx
idloop(3) = 1 
idloop(4) = ny
idloop(5) = 1 
idloop(6) = nz
    
idflt(1) = 1 
idflt(2) = nx

idflt(3) = 1 
idflt(4) = ny

idflt(5) = 1 
idflt(6) = nz

idarray(1) = 1 -hlo
idarray(2) = nx+hlo

idarray(3) = 1 -hlo
idarray(4) = ny+hlo

idarray(5) = 1 -hlo
idarray(6) = nz+hlo

size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)

#include "include_bcfilterx.f90"

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=idflt(5),idflt(6),size_bk  
    do bj=idflt(3),idflt(4),size_bj 
      do bi=idflt(1),idflt(2),size_bi 
    
idloop(6) = min( bk+size_bk, idflt(6)+1)-1
idloop(4) = min( bj+size_bj, idflt(4)+1)-1
idloop(2) = min( bi+size_bi, idflt(2)+1)-1

idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 

#include "LOOP_BEGIN"

#include "Filter_x"    
     
#include "LOOP_END"

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO
 
!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=idflt(5),idflt(6),size_bk  
    do bj=idflt(3),idflt(4),size_bj 
      do bi=idflt(1),idflt(2),size_bi 
    
idloop(6) = min( bk+size_bk, idflt(6)+1)-1
idloop(4) = min( bj+size_bj, idflt(4)+1)-1
idloop(2) = min( bi+size_bi, idflt(2)+1)-1

idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 

#include "LOOP_BEGIN"  
#include "updateFilter_x"   
#include "LOOP_END"
  
     enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

idloop(1) = 1 
idloop(2) = nx
idloop(3) = 1 
idloop(4) = ny
idloop(5) = 1 
idloop(6) = nz

#include "update_bcfilterx.f90"

!$OMP END PARALLEL

 end subroutine flt_x

subroutine flt_y(param_float,hlo,neq,nx,ny,nz,sizeblck,nbc,bc,q,q2)

implicit none  

#include "dtypes.h"
#include "param_fort.h" 

  real(wp), intent(in)    :: param_float(*)
! LOCAL VARIABLES

  integer, intent(in) :: nx,ny,nz,hlo,neq  
  integer, intent(in) :: sizeblck(3)

  integer :: i,j,k  
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk
  integer :: ibc
  integer, intent(in) :: bc(nbc),nbc

#include "includeRHS_globVar_filter.f90"
  
  integer :: idloop(6),idarray(6),idflt(6)


!f2py intent(in)    :: nx,ny,nz,hlo,neq,indvars,param_float
!f2py intent(inout) :: q

!$OMP PARALLEL DEFAULT(SHARED) private(idloop)

idloop(1) = 1 
idloop(2) = nx
idloop(3) = 1 
idloop(4) = ny
idloop(5) = 1 
idloop(6) = nz
  
idflt(1) = 1 
idflt(2) = nx

idflt(3) = 1 
idflt(4) = ny

idflt(5) = 1 
idflt(6) = nz


idarray(1) = 1 -hlo
idarray(2) = nx+hlo

idarray(3) = 1 -hlo
idarray(4) = ny+hlo

idarray(5) = 1 -hlo
idarray(6) = nz+hlo


size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)

#include "include_bcfiltery.f90"

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=idflt(5),idflt(6),size_bk  
    do bj=idflt(3),idflt(4),size_bj 
      do bi=idflt(1),idflt(2),size_bi 
    
idloop(6) = min( bk+size_bk, idflt(6)+1)-1
idloop(4) = min( bj+size_bj, idflt(4)+1)-1
idloop(2) = min( bi+size_bi, idflt(2)+1)-1

  idloop(5) = bk
  idloop(3) = bj
  idloop(1) = bi

#include "LOOP_BEGIN"

#include "Filter_y"    
     
#include "LOOP_END"

     enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=idflt(5),idflt(6),size_bk  
    do bj=idflt(3),idflt(4),size_bj 
      do bi=idflt(1),idflt(2),size_bi 
    
idloop(6) = min( bk+size_bk, idflt(6)+1)-1
idloop(4) = min( bj+size_bj, idflt(4)+1)-1
idloop(2) = min( bi+size_bi, idflt(2)+1)-1

  idloop(5) = bk
  idloop(3) = bj
  idloop(1) = bi

#include "LOOP_BEGIN"  
#include "updateFilter_y"   
#include "LOOP_END"
  
     enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

idloop(1) = 1 
idloop(2) = nx
idloop(3) = 1 
idloop(4) = ny
idloop(5) = 1 
idloop(6) = nz

#include "update_bcfiltery.f90"

!$OMP END PARALLEL


 end subroutine flt_y


 subroutine flt_z(param_float,hlo,neq,nx,ny,nz,sizeblck,nbc,bc,q,q2)

implicit none  

#include "dtypes.h"
#include "param_fort.h" 

  real(wp), intent(in)    :: param_float(*)
! LOCAL VARIABLES

  integer, intent(in) :: nx,ny,nz,hlo,neq  
  integer, intent(in) :: sizeblck(3)

  integer :: i,j,k  
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk
  integer :: ibc
  integer, intent(in) :: bc(nbc),nbc

#include "includeRHS_globVar_filter.f90"
  
  integer :: idloop(6),idarray(6),idflt(6)


!f2py intent(in)    :: nx,ny,nz,hlo,neq,indvars,param_float
!f2py intent(inout) :: q

!$OMP PARALLEL DEFAULT(SHARED) private(idloop)

idloop(1) = 1 
idloop(2) = nx
idloop(3) = 1 
idloop(4) = ny
idloop(5) = 1 
idloop(6) = nz
  
idflt(1) = 1 
idflt(2) = nx

idflt(3) = 1 
idflt(4) = ny

idflt(5) = 1 
idflt(6) = nz

idarray(1) = 1 -hlo
idarray(2) = nx+hlo

idarray(3) = 1 -hlo
idarray(4) = ny+hlo

idarray(5) = 1 -hlo
idarray(6) = nz+hlo

size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)

#include "include_bcfilterz.f90"

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=idflt(5),idflt(6),size_bk  
    do bj=idflt(3),idflt(4),size_bj 
      do bi=idflt(1),idflt(2),size_bi 
    
idloop(6) = min( bk+size_bk, idflt(6)+1)-1
idloop(4) = min( bj+size_bj, idflt(4)+1)-1
idloop(2) = min( bi+size_bi, idflt(2)+1)-1

idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 


#include "LOOP_BEGIN"

#include "Filter_z"    
     
#include "LOOP_END"


     enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=idflt(5),idflt(6),size_bk  
    do bj=idflt(3),idflt(4),size_bj 
      do bi=idflt(1),idflt(2),size_bi 
    
idloop(6) = min( bk+size_bk, idflt(6)+1)-1
idloop(4) = min( bj+size_bj, idflt(4)+1)-1
idloop(2) = min( bi+size_bi, idflt(2)+1)-1

  idloop(5) = bk
  idloop(3) = bj
  idloop(1) = bi

#include "LOOP_BEGIN"  
#include "updateFilter_z"   
#include "LOOP_END"
  
     enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

idloop(1) = 1 
idloop(2) = nx
idloop(3) = 1 
idloop(4) = ny
idloop(5) = 1 
idloop(6) = nz

#include "update_bcfilterz.f90"

!$OMP END PARALLEL

 end subroutine flt_z


subroutine bc(dir,param_int,param_float,data_float)!param_float,hlo,nrk,neq,indvars,nx,ny,nz,q,q1,q2,rhs)

implicit none  

#include "dtypes.h"
#include "param_fort.h" 

real(wp), intent(in)    :: param_float(*)
real(wp), intent(inout) ::  data_float(*)
integer,  intent(in)    ::   param_int(*)

integer, intent(in)   :: dir

!f2py  intent(in)   :: param_int,param_float,dir
!f2py intent(inout) :: data_float    

 if     (dir == 1) then

    call periodic_x(param_float                                            ,&
                    param_int(iadrHLO)                                   ,&
                    param_int(iadrNVAR)                                  ,&
                    param_int(iadrNX),param_int(iadrNY) , param_int(iadrNZ),&
                    param_int(iadrCacheBLCK)                             ,&               
                    data_float(1)                                          ,&
                    data_float(1+param_int(iadrNDIMTOT)*2)                   ) ! q2 is used for filter working array
    
 elseif (dir == 2)then

    call periodic_y(param_float                                            ,&
                    param_int(iadrHLO)                                   ,&
                    param_int(iadrNVAR)                                  ,&
                    param_int(iadrNX),param_int(iadrNY) , param_int(iadrNZ),&
                    param_int(iadrCacheBLCK)                             ,&                              
                    data_float(1)                                          ,&
                    data_float(1+param_int(iadrNDIMTOT)*2)                   ) ! q2 is used for filter working array

 elseif (dir == 3)then


    call periodic_z(param_float                                            ,&
                    param_int(iadrHLO)                                   ,&
                    param_int(iadrNVAR)                                  ,&
                    param_int(iadrNX),param_int(iadrNY) , param_int(iadrNZ),&
                    param_int(iadrCacheBLCK)                             ,&                              
                    data_float(1)                                          ,&
                    data_float(1+param_int(iadrNDIMTOT)*2)                   ) ! q2 is used for filter working array

 elseif (dir == 0 )then 

  call periodic_x(param_float                                            ,&
                    param_int(iadrHLO)                                   ,&
                    param_int(iadrNVAR)                                  ,&
                    param_int(iadrNX),param_int(iadrNY) , param_int(iadrNZ),&
                    param_int(iadrCacheBLCK)                             ,&               
                    data_float(1)                                          ,&
                    data_float(1+param_int(iadrNDIMTOT)*2)                   ) ! q2 is used for filter working array
  call periodic_y(param_float                                            ,&
                    param_int(iadrHLO)                                   ,&
                    param_int(iadrNVAR)                                  ,&
                    param_int(iadrNX),param_int(iadrNY) , param_int(iadrNZ),&
                    param_int(iadrCacheBLCK)                             ,&                              
                    data_float(1)                                          ,&
                    data_float(1+param_int(iadrNDIMTOT)*2)                   ) ! q2 is used for filter working array
  call periodic_z(param_float                                            ,&
                    param_int(iadrHLO)                                   ,&
                    param_int(iadrNVAR)                                  ,&
                    param_int(iadrNX),param_int(iadrNY) , param_int(iadrNZ),&
                    param_int(iadrCacheBLCK)                             ,&                              
                    data_float(1)                                          ,&
                    data_float(1+param_int(iadrNDIMTOT)*2)                   ) ! q2 is used for filter working array

    
 endif



 end subroutine bc

subroutine periodic_x(param_float,hlo,neq,nx,ny,nz,sizeblck,q,q2)

implicit none  

#include "dtypes.h"
#include "param_fort.h" 

  real(wp), intent(in)    :: param_float(*)
! LOCAL VARIABLES

  integer, intent(in) :: nx,ny,nz,hlo,neq  
  integer, intent(in) :: sizeblck(3)

  integer :: i,j,k  
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk

#include "includeRHS_globVar_filter.f90"
  
  integer :: idloop(6),idarray(6)


!f2py intent(in)    :: nx,ny,nz,hlo,neq,indvars,param_float
!f2py intent(inout) :: q

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

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=1,nz,size_bk  
    do bj=1,ny,size_bj 
    
idloop(6) = min( bk+size_bk, nz +1)-1
idloop(4) = min( bj+size_bj, ny +1)-1
idloop(2) = 0

idloop(5) = bk
idloop(3) = bj
idloop(1) = 0 

#include "LOOP_BEGIN"

#include "periodic_x"    
     
#include "LOOP_END"

  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO
!$OMP END PARALLEL
 
 end subroutine periodic_x

subroutine periodic_y(param_float,hlo,neq,nx,ny,nz,sizeblck,q,q2)

implicit none  

#include "dtypes.h"
#include "param_fort.h" 

  real(wp), intent(in)    :: param_float(*)
! LOCAL VARIABLES

  integer, intent(in) :: nx,ny,nz,hlo,neq  
  integer, intent(in) :: sizeblck(3)

  integer :: i,j,k  
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk

#include "includeRHS_globVar_filter.f90"
  
  integer :: idloop(6),idarray(6)


!f2py intent(in)    :: nx,ny,nz,hlo,neq,indvars,param_float
!f2py intent(inout) :: q

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

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=1,nz,size_bk      
      do bi=1-hlo,nx+hlo,size_bi 
    
  idloop(6) = min( bk+size_bk, nz+1)-1
  idloop(4) = 0
  idloop(2) = min( bi+size_bi, nx+hlo+1)-1

  idloop(5) = bk
  idloop(3) = 0
  idloop(1) = bi

#include "LOOP_BEGIN"

#include "periodic_y"    
     
#include "LOOP_END"

     enddo ! END cache blocking i
enddo ! END cache blocking k
!$OMP END DO
!$OMP END PARALLEL


 end subroutine periodic_y


 subroutine periodic_z(param_float,hlo,neq,nx,ny,nz,sizeblck,q,q2)

implicit none  

#include "dtypes.h"
#include "param_fort.h" 

  real(wp), intent(in)    :: param_float(*)
! LOCAL VARIABLES

  integer, intent(in) :: nx,ny,nz,hlo,neq  
  integer, intent(in) :: sizeblck(3)

  integer :: i,j,k  
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk

#include "includeRHS_globVar_filter.f90"
  
  integer :: idloop(6),idarray(6)


!f2py intent(in)    :: nx,ny,nz,hlo,neq,indvars,param_float
!f2py intent(inout) :: q

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

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2)  
    do bj=1-hlo,ny+hlo,size_bj 
      do bi=1-hlo,nx+hlo,size_bi 
    
idloop(6) = 0
idloop(4) = min( bj+size_bj, ny+hlo+1)-1
idloop(2) = min( bi+size_bi, nx+hlo+1)-1

idloop(5) = 0
idloop(3) = bj
idloop(1) = bi 


#include "LOOP_BEGIN"

#include "periodic_z"    
     
#include "LOOP_END"


     enddo ! END cache blocking i
  enddo ! END cache blocking j
!$OMP END DO
!$OMP END PARALLEL

 end subroutine periodic_z

subroutine pack(buf,f,ibeg,iend,jbeg,jend,kbeg,kend,sizex,sizey,sizez,sizenv) bind(c)
  use iso_c_binding, only: c_double, c_int
implicit none

#include "dtypes.h"

real(c_double), intent(out):: buf(*)
real(c_double), intent(in) :: f(*)
integer(c_int),intent(in)      :: ibeg,iend,jbeg,jend,kbeg,kend,sizex,sizey,sizez,sizenv

integer :: sx,sy,sz,i,j,k,n,indbuf,indf

sx = iend - ibeg 
sy = jend - jbeg 
sz = kend - kbeg 


!$OMP PARALLEL DEFAULT(SHARED) private(indbuf,indf)

!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
do n=1,sizenv
 do k=kbeg,kend-1
    do j=jbeg,jend-1 
!$omp simd 
      do i=ibeg,iend-1 
        indbuf = i - ibeg + 1 + (j-jbeg)*sx + (k-kbeg)*sx*sy + (n-1)*sx*sy*sz
        indf   = i        + 1 + j*sizex     + k*sizex*sizey  + (n-1)*sizex*sizey*sizez
        buf(indbuf) = f(indf)
      enddo
    enddo
 enddo
enddo 
!$OMP END DO NOWAIT
!$OMP END PARALLEL

end subroutine pack  


subroutine unpack(buf,f,ibeg,iend,jbeg,jend,kbeg,kend,sizex,sizey,sizez,sizenv)bind(c)
  use iso_c_binding, only: c_double, c_int
implicit none

#include "dtypes.h"

real(c_double), intent(in) :: buf(*)
real(c_double), intent(out) :: f(*)
integer(c_int),intent(in)      :: ibeg,iend,jbeg,jend,kbeg,kend,sizex,sizey,sizez,sizenv

integer :: sx,sy,sz,i,j,k,n,indbuf,indf

sx = iend - ibeg 
sy = jend - jbeg 
sz = kend - kbeg 


!$OMP PARALLEL DEFAULT(SHARED) private(indbuf,indf)

!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
do n=1,sizenv
 do k=kbeg,kend-1
    do j=jbeg,jend-1 
!$omp simd       
      do i=ibeg,iend-1 
        indbuf = i - ibeg + 1 + (j-jbeg)*sx + (k-kbeg)*sx*sy + (n-1)*sx*sy*sz
        indf   = i        + 1 + j*sizex     + k*sizex*sizey  + (n-1)*sizex*sizey*sizez
        f(indf) = buf(indbuf)
      enddo
    enddo
 enddo
enddo 
!$OMP END DO NOWAIT
!$OMP END PARALLEL

end subroutine unpack  

 subroutine init(param_int,param_float,data_float) bind(c)
  use iso_c_binding, only: c_double, c_int
  implicit none  

#include "dtypes.h"
#include "param_fort.h" 

  real(c_double), intent(in)    :: param_float(*)
  real(c_double), intent(inout) ::  data_float(*)
  integer,  intent(in)    :: param_int(*)

  call init_numa(param_float                                            ,&
                 param_int(iadrHLO), param_int(iadrNRK)                 ,&
                 param_int(iadrNVAR)                                    ,&
                 param_int(iadrNVARST)                                  ,&
                 param_int(iadrVARS)                                    ,&
                 param_int(iadrNX),param_int(iadrNY) , param_int(iadrNZ),&
                 param_int(iadrCacheBLCK)                               ,&
                 data_float(1)                                          ,& ! q
                 data_float(1+param_int(iadrNDIMTOT)  )                 ,& ! q1
                 data_float(1+param_int(iadrNDIMTOT)*2)                 ,& ! q2 
                 data_float(1+param_int(iadrNDIMTOT)*3)                 ,& ! rhs
                 data_float(1+param_int(iadrNDIMTOT)*4)                  ) ! stored (if any)        


  end subroutine init

subroutine init_numa(param_float,hlo,nrk,neq,neqst,ind,nx,ny,nz,sizeblck,q,q1,q2,rhs,qst)

  implicit none

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
  
#include "includeRHS_globVar_rk3.f90"

  
  integer :: idloop(6),idarray(6),indvars(neq),indvarsst(neqst)


  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 

!f2py intent(in)    :: q1,q2,a,nx,ny,nz,rhs,nrk
!f2py intent(inout) :: q

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
end module
