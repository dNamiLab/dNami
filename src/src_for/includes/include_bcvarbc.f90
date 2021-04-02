

do ibc=1,nbc

    SELECT CASE (bc(ibc))

#include "selectvarbcbc.f90"

    END SELECT

enddo