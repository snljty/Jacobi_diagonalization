# include <iostream>
# include <iomanip>
# include <cstring>
# include <string>
# include <cmath>
# include <cfloat>
# include <random>

# include "lapacke.h"

# include "Jacobi_diagonalization.hpp"

void Swap_method::Swap_int(int &a, int &b)
{
    int t = a;
    a = b;
    b = t;
    return;
}

Integer_vector::Integer_vector(int n) : n(n), vec(nullptr)
{
    if (n <= 0)
    {
        throw std::length_error("Error! size of a vector should be positive, but got " + std::to_string(n) + ".");
    }
    vec = new int [n];
    return;
}

Integer_vector::~Integer_vector()
{
    delete [] vec;
    vec = nullptr;
    n = 0;
    return;
}

int &Integer_vector::Get(int i)
{
    if (i < 0 || i >= n)
    {
        throw std::out_of_range("Error! Index " + std::to_string(i) + \
            "out of range [0, " + std::to_string(n) + ").");
    }
    return vec[i];
}

const int *Integer_vector::Get_vector() const
{
    return static_cast<const int *>(vec);
}

void Integer_vector::Swap(int i, int j)
{
    int t = Get(i);
    Get(i) = Get(j);
    Get(j) = t;
    return;
}

void Integer_vector::Fill_sequence_from_zero()
{
    for (int i = 0; i < n; ++ i)
    {
        Get(i) = i;
    }
    return;
}

void Integer_vector::Print(int num_len, int sep_len)
{
    for (int i = 0; i < n; ++ i)
    {
        if (! i)
        {
            std::cout << std::setw(sep_len) << "";
        }
        std::cout << std::setw(num_len) << Get(i);
    }
    std::cout << std::endl;
    return;
}

Vector::Vector(int n) : n(n), vec(nullptr)
{
    if (n <= 0)
    {
        throw std::length_error("Error! size of a vector should be positive, but got " + std::to_string(n) + ".");
    }
    vec = new double [n];
    return;
}

Vector::Vector(const Vector &old) : n(old.n), vec(new double [old.n])
{
    memcpy(vec, old.vec, n * sizeof(double));
    return;
}

Vector::~Vector()
{
    delete [] vec;
    vec = nullptr;
    n = 0;
    return;
}

double &Vector::Get(int i)
{
    if (i < 0 || i >= n)
    {
        throw std::out_of_range("Error! Index " + std::to_string(i) + \
            "out of range [0, " + std::to_string(n) + ").");
    }
    return vec[i];
}

const double *Vector::Get_vector() const
{
    return static_cast<const double *>(vec);
}

void Vector::Print(int len, int prec_len, int sep_len)
{
    std::cout << std::fixed << std::setprecision(prec_len);
    for (int i = 0; i < n; ++ i)
    {
        if (i)
        {
            std::cout << std::setw(sep_len) << "";
        }
        std::cout << std::setw(len) << Get(i);
    }
    std::cout << std::defaultfloat << std::endl;
    return;
}

void Vector::Fill_zeros()
{
    memset(vec, 0, n * sizeof(double));
    return;
}

void Vector::Swap(int i, int j)
{
    double t = Get(i);
    Get(i) = Get(j);
    Get(j) = t;
    return;
}

Square_matrix::Square_matrix(int n) : n(n), mat(nullptr)
{
    if (n <= 0)
    {
        throw std::length_error("Error! size of a square matrix should be positive, but got " + std::to_string(n) + ".");
    }
    mat = new double [n * n];
    return;
}

Square_matrix::Square_matrix(const Square_matrix &old) : n(old.n), mat(new double [old.n * old.n])
{
    memcpy(mat, old.mat, n * n * sizeof(double));
    return;
}

Square_matrix::~Square_matrix()
{
    delete [] mat;
    mat = nullptr;
    n = 0;
    return;
}

double &Square_matrix::Get(int i, int j)
{
    if (i < 0 || i >= n)
    {
        throw std::out_of_range("Error! Row index " + std::to_string(i) + \
            " out of range [0, " + std::to_string(n) + ").");
    }
    if (j < 0 || j >= n)
    {
        throw std::out_of_range("Error! Column index " + std::to_string(i) + \
            " out of range [0, " + std::to_string(n) + ").");
    }
    return mat[j * n + i];
}

const double *Square_matrix::Get_matrix() const
{
    return static_cast<const double *>(mat);
}

void Square_matrix::Print(int len, int prec_len, int sep_len)
{
    std::cout << std::fixed << std::setprecision(prec_len);
    for (int i = 0; i < n; ++ i)
    {
        for (int j = 0; j < n; ++ j)
        {
            if (j)
            {
                std::cout << std::setw(sep_len) << "";
            }
            std::cout << std::setw(len) << Get(i, j);
        }
        std::cout << std::endl;
    }
    std::cout << std::defaultfloat;
    return;
}

void Square_matrix::Fill_zeros()
{
    memset(mat, 0, n * n * sizeof(double));
    return;
}

void Square_matrix::Fill_unitary()
{
    Fill_zeros();
    for (int i = 0; i < n; ++ i)
    {
        Get(i, i) = 1.0;
    }
    return;
}

void Square_matrix::Fill_linear_molecule()
{
    Fill_zeros();
    int i, j;
    for (i = 1; i < n; ++ i)
    {
        j = i - 1;
        Get(i, j) = Get(j, i) = 1.0;
    }
    return;
}

void Square_matrix::Fill_ring_molecule()
{
    Fill_zeros();
    int i, j;
    for (i = 0; i < n; ++ i)
    {
        j = (i + 1) % n;
        Get(i, j) = Get(j, i) = 1.0;
    }
    return;
}

void Square_matrix::Fill_random()
{
    std::random_device seed;
    std::minstd_rand engine(seed());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i = 0; i < n; ++ i)
    {
        for (int j = 0; j <= i; ++ j)
        {
            Get(i, j) = distribution(engine);
            if (j != i)
            {
                Get(j, i) = Get(i, j);
            }
        }
    }
    return;
}

void Square_matrix::Row_to_vector(int i, Vector &vec)
{
    if (n != vec.n)
    {
        throw std::length_error("Error! Size of matrix (" + std::to_string(n) + \
            ") does not equal to size of vector (" + std::to_string(vec.n) + ").");
    }
    for (int j = 0; j < n; ++ j)
    {
        vec.Get(j) = Get(i, j);
    }
    return;
}

void Square_matrix::Column_to_vector(int j, Vector &vec)
{
    if (n != vec.n)
    {
        throw std::length_error("Error! Size of matrix (" + std::to_string(n) + \
            ") does not equal to size of vector (" + std::to_string(vec.n) + ").");
    }
    for (int i = 0; i < n; ++ i)
    {
        vec.Get(i) = Get(i, j);
    }
    return;
}

void Square_matrix::Jacobi_rotate(int p, int q, Square_matrix &rotate_matrix, Vector &vecp, Vector &vecq)
{
    double pp = Get(p, p), pq = Get(p, q), qq = Get(q, q);
    double theta = fabs(pp - qq) > 0.0 ? atan(2.0 * pq / (pp - qq)) / 2.0 : M_PI_4;
    double c = cos(theta), s = sin(theta);
    int i;

    // update eigenvectors
    rotate_matrix.Column_to_vector(p, vecp);
    rotate_matrix.Column_to_vector(q, vecq);
    for (i = 0; i < n; ++ i)
    {
        rotate_matrix.Get(i, p) = vecq.Get(i) * s + vecp.Get(i) * c;
        rotate_matrix.Get(i, q) = vecq.Get(i) * c - vecp.Get(i) * s;
    }

    // update target matrix
    Row_to_vector(p, vecp);
    Row_to_vector(q, vecq);
    for (i = 0; i < n; ++ i)
    {
        Get(p, i) = vecq.Get(i) * s + vecp.Get(i) * c;
        Get(q, i) = vecq.Get(i) * c - vecp.Get(i) * s;
    }
    Column_to_vector(p, vecp);
    Column_to_vector(q, vecq);
    for (i = 0; i < n; ++ i)
    {
        Get(i, p) = vecq.Get(i) * s + vecp.Get(i) * c;
        Get(i, q) = vecq.Get(i) * c - vecp.Get(i) * s;
    }
    Get(p, p) = pp * c * c + qq * s * s + 2.0 * pq * s * c;
    Get(q, q) = pp * s * s + qq * c * c - 2.0 * pq * s * c;
    Get(p, q) = (qq - pp) * s * c + pq * (c * c - s * s);
    // Get(p, q) = 0.0;
    // Get(q, p) = Get(p, q);

    return;
}

int Square_matrix::Jacobi_diagonalization(eigen &result, char uplo, bool destroy)
{
    if (toupper(uplo) != 'U' && toupper(uplo) != 'L')
    {
        throw std::invalid_argument("Error! \"uplo\" should be either \'U\' or \'L\', but got \'" + \
            std::string(1, toupper(uplo)) + "\'.");
    }

    if (n != result.n)
    {
        throw std::length_error("Error! The size of eigen result (" + std::to_string(result.n) + \
            ") should equal to the size of the matrix (" + std::to_string(n) + ").");
    }

    if (n == 1)
    {
        result.eigenvalues.Get(0) = Get(0, 0);
        result.eigenvectors.Get(0, 0) = 1.0;
        return 0;
    }

    if (! destroy)
    {
        Square_matrix mat_copy(* this);
        return mat_copy.Jacobi_diagonalization(result, uplo, true);
    }

    const bool is_greddy = false;
    const double tol = DBL_EPSILON;
    int iloop = 0;
    const int nloops_max = 30 * n * n; // 30 is from OpenCV, do not set to less than 2
    int p, q;
    int i, j;
    double abs_pq;
    Vector vecp(n), vecq(n);

    result.eigenvectors.Fill_unitary();
    if (toupper(uplo) == 'U')
    {
        if (is_greddy)
        {
            while (iloop < nloops_max)
            {
                p = 1;
                q = 0;
                abs_pq = fabs(Get(p, q));
                for (j = 1; j < n; ++ j)
                {
                    for (i = 0; i < j; ++ i)
                    {
                        if (fabs(Get(i, j)) > abs_pq)
                        {
                            p = i;
                            q = j;
                            abs_pq = Get(p, q);
                        }
                    }
                }
                if (abs_pq <= tol)
                {
                    break;
                }
                Jacobi_rotate(p, q, result.eigenvectors, vecp, vecq);
                ++ iloop;
            }
        }
        else
        {
            while (iloop < nloops_max)
            {
                if (Is_diagonal_matrix('U', tol))
                {
                    break;
                }
                for (q = 1; q < n; ++ q)
                {
                    for (p = 0; p < q; ++ p)
                    {
                        Jacobi_rotate(p, q, result.eigenvectors, vecp, vecq);
                        ++ iloop;
                    }
                }
            }
        }
    }
    /* the lower triangle will be used */
    else
    {
        if (is_greddy)
        {
            while (iloop < nloops_max)
            {
                p = 0;
                q = 1;
                abs_pq = fabs(Get(p, q));
                for (j = 0; j < n - 1; ++ j)
                {
                    for (i = j + 1; i < n; ++ i)
                    {
                        if (fabs(Get(i, j)) > abs_pq)
                        {
                            p = i;
                            q = j;
                            abs_pq = fabs(Get(p, q));;
                        }
                    }
                }
                if (abs_pq <= tol)
                {
                    break;
                }
                Jacobi_rotate(p, q, result.eigenvectors, vecp, vecq);
                ++ iloop;
            }
        }
        else
        {
            while (iloop < nloops_max)
            {
                if (Is_diagonal_matrix('L', tol))
                {
                    break;
                }
                for (q = 0; q < n - 1; ++ q)
                {
                    for (p = q + 1; p < n; ++ p)
                    {
                        Jacobi_rotate(p, q, result.eigenvectors, vecp, vecq);
                        ++ iloop;
                    }
                }
            }
        }
    }

    for (int i = 0; i < n; ++ i)
    {
        result.eigenvalues.Get(i) = Get(i, i);
    }

    if (iloop > nloops_max)
    {
        std::cerr << "Warning! Jacobi diagonalization algorithm did not convergd." << std::endl;
    }

    result.Sort();
    return iloop > nloops_max ? 1 : 0;
}

bool Square_matrix::Is_diagonal_matrix(char uplo, double tol)
{
    int i, j;
    if (toupper(uplo) == 'U')
    {
        for (j = 1; j < n; ++ j)
        {
            for (i = 0; i < j; ++ i)
            {
                if (fabs(Get(i, j)) > tol)
                {
                    return false;
                }
            }
        }
        return true;
    }
    else if (toupper(uplo) == 'L')
    {
        for (j = 0; j < n - 1; ++ j)
        {
            for (i = j + 1; i < n; ++ i)
            {
                if (fabs(Get(i, j)) > tol)
                {
                    return false;
                }
            }
        }
        return true;
    }
    else
    {
        throw std::invalid_argument("Error! \"UPLO\" should be either \'U\' or \'L\', but got \'" + \
            std::string(1, toupper(uplo)) + "\'.");
    }
}

eigen::eigen(int n) : n(n), eigenvalues(n), eigenvectors(n)
{
    return;
}

eigen::eigen(const eigen &old) : n(old.eigenvalues.n), eigenvalues(old.eigenvalues), eigenvectors(old.eigenvectors)
{
    return;
}

void eigen::Sort()
{
    using Swap_method::Swap_int;
    int i, j;
    int dad, son;
    Integer_vector ind(n);
    Vector swp(n);
    int pos;

    // do not use standard library because it is inconvient here

    // initialize the sort sequence
    ind.Fill_sequence_from_zero();

    // heap sort of eigenvalues, and get the sort sequence
    for (i = n / 2 - 1; i >= 0; -- i)
    {
        dad = i;
        son = dad * 2 + 1;
        while (son <= n - 1)
        {
            if (son + 1 <= n - 1 && eigenvalues.Get(son) < eigenvalues.Get(son + 1))
            {
                ++ son;
            }
            if (eigenvalues.Get(son) <= eigenvalues.Get(dad))
            {
                break;
            }
            else
            {
                eigenvalues.Swap(dad, son);
                ind.Swap(dad, son);
                dad = son;
                son = dad * 2 + 1;
            }
        }
    }
    for (i = n - 1; i > 0; -- i)
    {
        eigenvalues.Swap(0, i);
        ind.Swap(0, i);
        dad = 0;
        son = 1;
        while (son <= i - 1)
        {
            if (son + 1 <= i - 1 && eigenvalues.Get(son) < eigenvalues.Get(son + 1))
            {
                ++ son;
            }
            if (eigenvalues.Get(son) <= eigenvalues.Get(dad))
            {
                break;
            }
            else
            {
                eigenvalues.Swap(dad, son);
                ind.Swap(dad, son);
                dad = son;
                son = dad * 2 + 1;
            }
        }
    }

    // reorder eigenvectors according to ind
    i = 0;
    while (i < n - 1)
    {
        if (ind.Get(i) == i)
        {
            ++ i;
            continue;
        }
        pos = i;
        for (j = 0; j < n; ++ j)
        {
            swp.Get(j) = eigenvectors.Get(j, i);
        }
        while (ind.Get(i) != pos)
        {
            for (j = 0; j < n; ++ j)
            {
                eigenvectors.Get(j, i) = eigenvectors.Get(j, ind.Get(i));
            }
            Swap_int(i, ind.Get(i));
        }
        for (j = 0; j < n; ++ j)
        {
            eigenvectors.Get(j, i) = swp.Get(j);
        }
        ind.Get(i) = i;
        i = pos + 1;
    }
    return; 
}

int eigen::Reference_eigen(const Square_matrix &matrix, char uplo)
{
    if (matrix.n != n)
    {
        throw std::length_error("Error! The size of the passed matrix (" + std::to_string(matrix.n) + \
            ") should equal to the size of eigen (" + std::to_string(n) + ").");
    }
    if (toupper(uplo) != 'U' && toupper(uplo) != 'L')
    {
        throw std::invalid_argument("Error! \"uplo\" should be either \'U\' or \'L\', but got \'" + \
            std::string(1, toupper(uplo)) + "\'.");
    }
    // memcpy(& eigenvectors.Get(0, 0), & matrix.Get(0, 0), n * n * sizeof(double));
    memcpy(& eigenvectors.Get(0, 0), matrix.Get_matrix(), n * n * sizeof(double));
    double *w = & eigenvalues.Get(0);
    double *v = & eigenvectors.Get(0, 0);
    lapack_int ret;

    ret = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', uplo, n, v, n, w);

    return ret;
}

