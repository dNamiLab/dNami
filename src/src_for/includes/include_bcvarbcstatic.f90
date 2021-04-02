

do ibc=1,nbc

    SELECT CASE (bc(ibc))

#include "selectvarbcstaticbc.f90"

    END SELECT

enddo