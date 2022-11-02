# include <iostream>

# include "Jacobi_diagonalization.hpp"

int main()
{
    const int n = 5;
    Square_matrix a(n);
    a.Fill_random();
    std::cout << "Hamilton Matrix" << std::endl;
    a.Print();
    std::cout << std::endl;
    eigen res(n);
    a.Jacobi_diagonalization(res);
    std::cout << "Eigenvalues" << std::endl;
    res.eigenvalues.Print();
    std::cout << std::endl;
    std::cout << "Eigenvectors" << std::endl;
    res.eigenvectors.Print();
    std::cout << std::endl;
    std::cout << "Reference:" << std::endl << std::endl;
    eigen ref(n);
    ref.Reference_eigen(a);
    std::cout << "Eigenvalues" << std::endl;
    ref.eigenvalues.Print();
    std::cout << std::endl;
    std::cout << "Eigenvectors" << std::endl;
    ref.eigenvectors.Print();
    std::cout << std::endl;    

    return 0;
}

