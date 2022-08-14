# ifndef __JACOBI_DIAGONALIZATION_H__
# define __JACOBI_DIAGONALIZATION_H__ 1

/* the matrices are assumed to be column major */

double Get_Mat_Val(double const *a, int i, int j, int lda);

double *Get_Mat_Ptr(double *a, int i, int j, int lda);

void Print_interval_vector(int n, int step, double const *a);

void Print_vector(int n, double const *a);

void Print_matrix(int m, int n, double const *a);

void Print_square_matrix(int n, double const *a);

int Jacobi_diagonalization(int n, double *a, double *w, double *v, double *vecp, double *vecq, char uplo, bool destroy);

int Jacobi_rotate(int n, int p, int q, double *mat, double *v, double *vecp, double *vecq);

void Sort_eigen(int n, double *w, double *v);

void Swap_double(double *a, double *b);

void Swap_int(int *a, int *b);

void Copy_array(int n, double const *a, int inda, double *b, int indb);

void Print_repeat(FILE *fp, char const *s, int num);

# endif

