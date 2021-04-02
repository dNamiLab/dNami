  call init_numa(param_float                                            ,&
                 param_int(iadrHLO), param_int(iadrNRK)                 ,&
                 param_int(iadrNVAR)                                    ,&
                 param_int(iadrNVARST)                                  ,&
                 param_int(iadrVARS)                                    ,&
                 param_int(iadrNX),param_int(iadrNY) , param_int(iadrNZ),&
                 param_int(iadrCacheBLCK)                               ,&
                 data_float(1)                                          ,& ! q
                 data_float(1+param_int(iadrNDIMTOT)  )                 ,& ! q1
                 data_float(1+param_int(iadrNDIMTOT)*2)                 ,& ! rhs
                 data_float(1+param_int(iadrNDIMTOT)*3)                  ) ! stored (if any)        
