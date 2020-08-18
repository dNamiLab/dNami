

do ibc=1,nbc

    SELECT CASE (bc(ibc))

#include "selectupdate_filterbc_x.f90"

    END SELECT

enddo