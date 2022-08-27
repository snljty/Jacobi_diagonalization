/* Jacobi diagonalization algorithm to get the eigenvalues and eigenvectors */
/* currently there is no method for handling degenerate situation */

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <math.h>
# include <float.h>
# include <stdbool.h>

# include "compare.h"
# include "Jacobi_diagonalization.h"

/* the matrices are assumed to be column major */


int main(int argc, const char **argv)
{
    int n;
    int const n_def = 5;
    double *a = NULL, *v = NULL, *w = NULL, *vecp = NULL, *vecq = NULL;
    int i, j, k;
    int status;

    double *v_cmp = NULL, *w_cmp = NULL;
    bool const is_compare = true;

    FILE *mat_fl = NULL;
    char const mat_name[] = "mat.dat";

    /* pharse command arguments */
    if (argc - 1 == 1)
    {
        sscanf(argv[1], "%d", & n);
    }
    else
    {
        n = n_def;
    }

    /* allocate memory for arrays */
    a = (double *)malloc(n * n * sizeof(double));
    v = (double *)malloc(n * n * sizeof(double));
    w = (double *)malloc(n * sizeof(double));
    vecp = (double *)malloc(n * sizeof(double));
    vecq = (double *)malloc(n * sizeof(double));

    /* initialize a */
    /*
    memset(a, 0, n * n * sizeof(double));
    for (i = 0; i < n; ++ i)
    {
        j = (i + 1) % n;
        * Get_Mat_Ptr(a, i, j, n) = 1.0;
        * Get_Mat_Ptr(a, j, i, n) = 1.0;
    }
    * Get_Mat_Ptr(a, 0, n - 1, n) = 0.0;
    * Get_Mat_Ptr(a, n - 1, 0, n) = 0.0;
    */

    /* or */
    if (! (mat_fl = fopen(mat_name, "rb")))
    {
        fprintf(stderr, "Error! Cannot open \"%s\" for reading.\n", mat_name);
        exit(EXIT_FAILURE);
    }
    if (fread(a, sizeof(double), n * n, mat_fl) != n * n)
    {
        fprintf(stderr, "Error! Cannot read %d (%d * %d) numbers from \"%s\".\n", n * n, n, n, mat_name);
        exit(EXIT_FAILURE);
    }
    fclose(mat_fl);
    mat_fl = NULL;

    printf("%s\n", "Matrix:");
    Print_square_matrix(n, a);
    printf("\n");

    if (is_compare)
    {
        v_cmp = (double *)malloc(n * n * sizeof(double));
        w_cmp = (double *)malloc(n * sizeof(double));
        status = Reference_diagonalization(n, a, w_cmp, v_cmp, 'L');
        if (status)
        {
            fprintf(stderr, "%s%d%s\n", "Error! LAPACK dsyev returns ", status, " .");
        }
    }

    status = Jacobi_diagonalization(n, a, w, v, vecp, vecq, 'L', false);
    if (status)
    {
        Print_repeat(stderr, "#", 64);
        fprintf(stderr, "\n");
        Print_repeat(stderr, "#", 8);
        fprintf(stderr, "%s%2d%s", " Error! The algorithm returns ", status, " instead of 0 . ");
        Print_repeat(stderr, "#", 8);
        fprintf(stderr, "\n");
        Print_repeat(stderr, "#", 64);
        fprintf(stderr, "\n");
    }
    Sort_eigen(n, w, v);

    printf("%s\n", "Eigenvalues:");
    Print_vector(n, w);
    printf("\n");
    printf("%s\n", "Eigenvectors:");
    Print_square_matrix(n, v);
    printf("\n");

    free(a);
    a = NULL;
    free(w);
    w = NULL;
    free(v);
    v = NULL;
    free(vecp);
    vecp = NULL;
    free(vecq);
    vecq = NULL;

    if (is_compare)
    {
        printf("%s\n", "Reference: ");
        printf("\n");
        printf("%s\n", "Eigenvalues:");
        Print_vector(n, w_cmp);
        printf("\n");
        printf("%s\n", "Eigenvectors:");
        Print_square_matrix(n, v_cmp);
        printf("\n");
        free(w_cmp);
        w_cmp = NULL;
        free(v_cmp);
        v_cmp = NULL;
    }

    return 0;
}

inline double Get_Mat_Val(double const *a, int i, int j, int lda)
{
    return a[j * lda + i];
}

inline double *Get_Mat_Ptr(double *a, int i, int j, int lda)
{
    return a + j * lda + i;
}

void Print_interval_vector(int n, int step, double const *a)
{
    int j;

    if (step != 1)
    {
        for (j = 0; j < n; ++ j)
        {
            if (j)
            {
                printf("%2s", "");
            }
            printf("%10.7lf", a[j * step]);
        }
    }
    else
    {
        for (j = 0; j < n; ++ j)
        {
            if (j)
            {
                printf("%2s", "");
            }
            printf("%10.7lf", a[j]);
        }
    }
    printf("\n");

    return;
}

void Print_vector(int n, double const *a)
{
    Print_interval_vector(n, 1, a);

    return;
}


void Print_matrix(int m, int n, double const *a)
{
    int i;

    for (i = 0; i < m; ++ i)
    {
        Print_interval_vector(n, m, a + i);
    }

    return;
}

void Print_square_matrix(int n, double const *a)
{
    Print_matrix(n, n, a);

    return;
}

int Jacobi_diagonalization(int n, double *a, double *w, double *v, double *vecp, double *vecq, char uplo, bool destroy)
{
    int status;
    int i, j, p, q;
    double *tmp = NULL, *mat = NULL;
    double abs_pq;

    int iloop, nloops_max;
    double const tol = DBL_EPSILON;
    bool const is_greddy = true;

    /* check uplo */
    if (uplo != 'U' && uplo != 'u' && uplo != 'L' && uplo != 'l')
    {
        status = -2;
        return status;
    }

    /* check n and if we can quick return */
    status = 0;
    if (n < 1)
    {
        status = -1;
        return status;
    }
    else if (n == 1)
    {
        * v = 1.0;
        * w = * a;
        return status;
    }
    else
    {
        ;
    }

    /* allocate temporary matrix if needed, assign pointer to the actual matrix */
    if (! destroy)
    {
        tmp = (double *)malloc(n * n * sizeof(double));
        Copy_array(n * n, a, 1, tmp, 1);
    }
    mat = destroy ? a : tmp;

    /* initialize eigenvectors as unit matrix */
    memset(v, 0, n * n * sizeof(double));
    for (i = 0; i < n; ++ i)
    {
        * Get_Mat_Ptr(v, i, i, n) = 1.0;
    }

    /* maximum number of loops  */
    nloops_max = 30 * n * n; /* do not set less than 2 * n * n */
    /* 30 is from OpenCV */

    /* the main loop */
    iloop = 0;
    /* the upper triangle will be used */
    if (uplo == 'U' || uplo == 'u')
    {
        if (is_greddy)
        {
            while (iloop < nloops_max)
            {
                p = 1;
                q = 0;
                abs_pq = fabs(Get_Mat_Val(mat, p, q, n));
                for (j = 1; j < n; ++ j)
                {
                    for (i = 0; i < j; ++ i)
                    {
                        if (fabs(Get_Mat_Val(mat, i, j, n)) > abs_pq)
                        {
                            p = i;
                            q = j;
                            abs_pq = fabs(Get_Mat_Val(mat, p, q, n));
                        }
                    }
                }
                if (abs_pq <= tol)
                {
                    break;
                }
                status = Jacobi_rotate(n, p, q, mat, v, vecp, vecq);
                ++ iloop;
            }
        }
        else
        {
            while (iloop < nloops_max)
            {
                if (Is_converged(n, mat, 'U', tol))
                {
                    break;
                }
                for (q = 1; q < n; ++ q)
                {
                    for (p = 0; p < q; ++ p)
                    {
                        status = Jacobi_rotate(n, p, q, mat, v, vecp, vecq);
                        ++ iloop;
                    }
                }
            }
        }
    }
    /* the lower triangle will be used */
    else if (uplo == 'L' || uplo == 'l')
    {
        if (is_greddy)
        {
            while (iloop < nloops_max)
            {
                p = 0;
                q = 1;
                abs_pq = fabs(Get_Mat_Val(mat, p, q, n));
                for (j = 0; j < n - 1; ++ j)
                {
                    for (i = j + 1; i < n; ++ i)
                    {
                        if (fabs(Get_Mat_Val(mat, i, j, n)) > abs_pq)
                        {
                            p = i;
                            q = j;
                            abs_pq = fabs(Get_Mat_Val(mat, p, q, n));
                        }
                    }
                }
                if (abs_pq <= tol)
                {
                    break;
                }
                status = Jacobi_rotate(n, p, q, mat, v, vecp, vecq);
                ++ iloop;
            }
        }
        else
        {
            while (iloop < nloops_max)
            {
                if (Is_converged(n, mat, 'L', tol))
                {
                    break;
                }
                for (q = 0; q < n - 1; ++ q)
                {
                    for (p = q + 1; p < n; ++ p)
                    {
                        status = Jacobi_rotate(n, p, q, mat, v, vecp, vecq);
                        ++ iloop;
                    }
                }
            }
        }
    }

    /* get the eigenvalues from the diagonal of the target matrix */
    for (i = 0; i < n; ++ i)
    {
        w[i] = Get_Mat_Val(mat, i, i, n);
    }

    /* if not converged, returns 1 */
    if (iloop >= nloops_max)
    {
        status = -1;
    }

    /* nullify pointer, release memory if allocated in this subroutine */
    mat = NULL;
    if (! destroy)
    {
        free(tmp);
        tmp = NULL;
    }

    return status;
}

int Jacobi_rotate(int n, int p, int q, double *mat, double *v, double *vecp, double *vecq)
{
    int status;
    double theta;
    double c, s;
    double pp, pq, qq;
    int i;

    /* quick return if wrong in input */
    status = 0;
    if (p < 0 || p >= n || q < 0 || q >= n)
    {
        status = -3;
        return status;
    }

    /* get rotation angle and sin as well as cos values */
    pp = Get_Mat_Val(mat, p, p, n);
    pq = Get_Mat_Val(mat, p, q, n);
    qq = Get_Mat_Val(mat, q, q, n);
    if (fabs(pp - qq) > 0.0)
    {
        theta = atan(2.0 * pq / (pp - qq)) / 2.0;
    }
    else
    {
        theta = M_PI_4;
    }
    c = cos(theta);
    s = sin(theta);

    /* update eigenvectors */
    Copy_array(n, Get_Mat_Ptr(v, 0, p, n), 1, vecp, 1);
    Copy_array(n, Get_Mat_Ptr(v, 0, q, n), 1, vecq, 1);
    for (i = 0; i < n; ++ i)
    {
        * Get_Mat_Ptr(v, i, p, n) = vecq[i] * s + vecp[i] * c;
        * Get_Mat_Ptr(v, i, q, n) = vecq[i] * c - vecp[i] * s;
    }

    /* update target matrix */
    Copy_array(n, Get_Mat_Ptr(mat, p, 0, n), n, vecp, 1);
    Copy_array(n, Get_Mat_Ptr(mat, q, 0, n), n, vecq, 1);
    for (i = 0; i < n; ++ i)
    {
        * Get_Mat_Ptr(mat, p, i, n) = vecq[i] * s + vecp[i] * c;
        * Get_Mat_Ptr(mat, q, i, n) = vecq[i] * c - vecp[i] * s;
    }
    Copy_array(n, Get_Mat_Ptr(mat, 0, p, n), 1, vecp, 1);
    Copy_array(n, Get_Mat_Ptr(mat, 0, q, n), 1, vecq, 1);
    for (i = 0; i < n; ++ i)
    {
        * Get_Mat_Ptr(mat, i, p, n) = vecq[i] * s + vecp[i] * c;
        * Get_Mat_Ptr(mat, i, q, n) = vecq[i] * c - vecp[i] * s;
    }
    * Get_Mat_Ptr(mat, p, p, n) = pp * c * c + qq * s * s + 2.0 * pq * s * c;
    * Get_Mat_Ptr(mat, q, q, n) = pp * s * s + qq * c * c - 2.0 * pq * s * c;
    * Get_Mat_Ptr(mat, p, q, n) = (qq - pp) * s * c + pq * (c * c - s * s);
    /*
    * Get_Mat_Ptr(mat, p, q, n) = 0.0;
    * Get_Mat_Ptr(mat, q, p, n) = Get_Mat_Val(mat, p, q, n);
    */

    return status;
}

void Sort_eigen(int n, double *w, double *v)
{
    int i;
    int dad, son;
    int *ind = NULL;
    double *swp = NULL;
    int pos;

    /* do not use qsort because it is inconvient here */

    /* initialize the sort sequence, allocate memory for swap */
    ind = (int *)malloc(n * sizeof(int));
    for (i = 0; i < n; ++ i)
    {
        ind[i] = i;
    }
    swp = (double *)malloc(n * sizeof(double));

    /* heap sort of w, and get the sort sequence */
    for (i = n / 2 - 1; i >= 0; -- i)
    {
        dad = i;
        son = dad * 2 + 1;
        while (son <= n - 1)
        {
            if (son + 1 <= n - 1 && w[son] < w[son + 1])
            {
                ++ son;
            }
            if (w[son] <= w[dad])
            {
                break;
            }
            else
            {
                Swap_double(w + dad, w + son);
                Swap_int(ind + dad, ind + son);
                dad = son;
                son = dad * 2 + 1;
            }
        }
    }
    for (i = n - 1; i > 0; i--)
    {
        Swap_double(w, w + i);
        Swap_int(ind, ind + i);
        dad = 0;
        son = 1;
        while (son <= i - 1)
        {
            if (son + 1 <= i - 1 && w[son] < w[son + 1])
            {
                ++ son;
            }
            if (w[son] <= w[dad])
            {
                break;
            }
            else
            {
                Swap_double(w + dad, w + son);
                Swap_int(ind + dad, ind + son);
                dad = son;
                son = dad * 2 + 1;
            }
        }
    }

    /* reorder v according to ind */
    i = 0;
    while (i < n - 1)
    {
        if (ind[i] == i)
        {
            ++ i;
            continue;
        }
        pos = i;
        Copy_array(n, Get_Mat_Ptr(v, 0, i, n), 1, swp, 1);
        while (ind[i] != pos)
        {
            Copy_array(n, Get_Mat_Ptr(v, 0, ind[i], n), 1, Get_Mat_Ptr(v, 0, i, n), 1);
            Swap_int(& i, ind + i);
        }
        Copy_array(n, swp, 1, Get_Mat_Ptr(v, 0, i, n), 1);
        ind[i] = i;
        i = pos + 1;
    }

    /* release memory */
    free(ind);
    ind = NULL;
    free(swp);
    swp = NULL;

    return;
}

void Swap_double(double *a, double *b)
{
    double t;

    t = * a;
    * a = * b;
    * b = t;

    return;
}

void Swap_int(int *a, int *b)
{
    int t;

    t = * a;
    * a = * b;
    * b = t;

    return;
}

void Copy_array(int n, double const *a, int inda, double *b, int indb)
{
    int i;
    int ia, ib;

    if (inda == 1 && indb == 1)
    {
        memcpy(b, a, n * sizeof(double));
    }
    else
    {
        ia = 0;
        ib = 0;
        for (i = 0; i < n; ++ i)
        {
            b[ib] = a[ia];
            ib += indb;
            ia += inda;
        }
    }

    return;
}

void Print_repeat(FILE *fp, char const *s, int num)
{
    int i;

    for (i = 0; i < num; ++ i)
    {
        fprintf(fp, "%s", s);
    }

    return;
}

bool Is_converged(int n, double const *mat, char uplo, double tol)
{
    int i, j;

    if (uplo == 'U' || uplo == 'u')
    {
        for (j = 1; j < n; ++ j)
        {
            for (i = 0; i < j; ++ i)
            {
                if (fabs(Get_Mat_Val(mat, i, j, n)) > tol)
                {
                    return false;
                }
            }
        }
        return true;
    }
    else if (uplo == 'L' || uplo == 'l')
    {
        for (j = 0; j < n - 1; ++ j)
        {
            for (i = j + 1; i < n; ++ i)
            {
                if (fabs(Get_Mat_Val(mat, i, j, n)) > tol)
                {
                    return false;
                }
            }
        }
        return true;
    }
    else
    {
        return false;
    }

    return false; /* should never happen */
}

