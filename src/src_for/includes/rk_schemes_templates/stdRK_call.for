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
