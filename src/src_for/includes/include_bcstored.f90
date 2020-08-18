

do ibc=1,nbc

    SELECT CASE (bc(ibc))

#include "selectstoredbc.f90"

    END SELECT

enddo