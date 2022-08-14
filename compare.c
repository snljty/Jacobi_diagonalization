# include <stdlib.h>
/* # include "cblas.h" */
# include "lapacke.h"

# include "compare.h"

int Reference_diagonalization(int n, double const *a, double *w, double *v, char uplo)
{
    LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', n, n, a, n, v, n);
    return LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', uplo, n, v, n, w);
}
