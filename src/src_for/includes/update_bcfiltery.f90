

do ibc=1,nbc

    SELECT CASE (bc(ibc))

#include "selectupdate_filterbc_y.f90"

    END SELECT

enddo