! The following function call was inserted from the template
! file include includeWilliamsonRK_call_cmp_stored.for
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
data_float(1+param_int(iadrNDIMTOT)*2)                               ,& ! rhs
data_float(1+param_int(iadrNDIMTOT)*3)                               ,& ! stored (if any)
data_float(1+param_int(iadrNDIMTOT)*3 + shiftstore)                  ,& ! qbcface_i  (if any)
data_float(1+param_int(iadrNDIMTOT)*3 + shiftfacei)                  ,& ! qbcface_j  (if any)
data_float(1+param_int(iadrNDIMTOT)*3 + shiftfacej)                  ,& ! qbcface_k  (if any)
data_float(1+param_int(iadrNDIMTOT)*3 + shiftfacek)                  ,& ! qbcedge_ij (if any)
data_float(1+param_int(iadrNDIMTOT)*3 + shiftedgeij)                 ,& ! qbcedge_jk (if any)
data_float(1+param_int(iadrNDIMTOT)*3 + shiftedgejk)                 )  ! qbcedge_ik (if any)         