# ifndef __JACOBI_DIAGNOALIZATION_HPP__
# define __JACOBI_DIAGNOALIZATION_HPP__ 1

# include <cfloat>

namespace Swap_method
{
    void Swap_int(int &a, int &b);
};

class Integer_vector
{
private:
    int *vec;
public:
    int n;

    Integer_vector(int n);
    ~Integer_vector();
    int &Get(int i);
    const int *Get_vector() const;
    void Swap(int i, int j);
    void Fill_sequence_from_zero();
    void Print(int num_len = 2, int sep_len = 1);
};

class eigen;

class Vector
{
private:
    double *vec;
public:
    int n;

    Vector(int n);
    Vector(const Vector &old);
    ~Vector();
    double &Get(int i);
    const double *Get_vector() const;
    void Print(int len = 9, int prec_len = 6, int sep_len = 2);
    void Fill_zeros();
    void Swap(int i, int j);
};

class Square_matrix
{
private:
    double *mat;
    void Row_to_vector(int i, Vector &vec);
    void Column_to_vector(int j, Vector &vec);
    void Jacobi_rotate(int p, int q, Square_matrix &rotate_matrix, Vector &vecp, Vector &vecq);
public:
    int n;

    Square_matrix(int n);
    Square_matrix(const Square_matrix &old);
    ~Square_matrix();
    double &Get(int i, int j);
    const double *Get_matrix() const;
    void Print(int len = 9, int prec_len = 6, int sep_len = 2);
    void Fill_zeros();
    void Fill_unitary();
    void Fill_linear_molecule();
    void Fill_ring_molecule();
    void Fill_random();
    int Jacobi_diagonalization(eigen &result, char uplo = 'L', bool destroy = false);
    bool Is_diagonal_matrix(char uplo = 'L', double tol = DBL_EPSILON);
};

class eigen
{
public:
    int n;
    Vector eigenvalues;
    Square_matrix eigenvectors;

    eigen(int n);
    eigen(const eigen &old);

    void Sort();
    int Reference_eigen(const Square_matrix &matrix, char uplo = 'L');
};

# endif // __JACOBI_DIAGNOALIZATION_HPP__

