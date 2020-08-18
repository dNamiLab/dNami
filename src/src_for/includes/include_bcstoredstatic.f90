

do ibc=1,nbc

    SELECT CASE (bc(ibc))

#include "selectstoredstaticbc.f90"

    END SELECT

enddo