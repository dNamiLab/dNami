

do ibc=1,nbc

    SELECT CASE (bc(ibc))

#include "select_phybc_q.f90"

	END SELECT

enddo