

do ibc=1,nbc

    SELECT CASE (bc(ibc))

#include "selectfilterbc_x.f90"

    END SELECT

enddo