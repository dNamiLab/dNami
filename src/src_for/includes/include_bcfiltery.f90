

do ibc=1,nbc

    SELECT CASE (bc(ibc))

#include "selectfilterbc_y.f90"

    END SELECT

enddo