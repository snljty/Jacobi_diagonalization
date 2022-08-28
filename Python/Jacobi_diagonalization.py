#! /usr/bin/env python3
# -*- Coding: UTF-8 -*-

# Jacobi diagonalization algorithm

import numpy as np

import sys
argc = len(sys.argv)

class Const(object):
    class ConstError(TypeError):
        pass
    def __setattr__(self, name, value):
        if self.__dict__.__contains__(name):
            raise self.ConstError(f"Cannot rebind const \"{name:s}\"")
        else:
            self.__dict__[name] = value

const = Const()
const.n_def = 5
const.tol_def = sys.float_info.epsilon
const.mat_name = "mat.dat"

class Square_matrix(object):
    def __init__(self, n: int):
        if (n < 1):
            raise ValueError(f"Error! n must be positive, but got {self.n:d}.")
        self.n = n
        self.mat = np.empty((n, n))
        return

    def Init_ring(self):
        self.mat[:, :] = 0.0
        for i in np.arange(self.n):
            j = (i + 1) % self.n
            self.mat[i][j] = self.mat[j][i] = 1.0
        return

    def Init_linear(self):
        self.Init_ring()
        self.mat[0][n - 1] = self.mat[n - 1][0] = 0.0
        return

    def Jacobi_diagonalization(self, uplo: str = "L", is_greddy: bool = True) -> (np.array, np.array):
        # quick return if possible
        if self.n < 1:
            raise ValueError(f"Error! n must be positive, but got {self.n:d}.")
        if uplo.upper().find("U") != 0 and uplo.upper().find("L") != 0:
            raise ValueError(f"Error! uplo should be either \"U\" or \"L\"., but got {uplo:s}.")
        mat = self.mat.copy()
        if self.n == 1:
            return np.ones(1), mat

        def Is_converged(uplo: str = "L", tol: np.double = const.tol_def) -> bool:
            if not uplo.upper().find("U"):
                for j in np.arange(1, self.n):
                    for i in np.arange(j):
                        if np.abs(mat[i, j]) > tol:
                            return False
                return True
            elif not uplo.upper().find("L"):
                for j in np.arange(n - 1):
                    for i in np.arange(j + 1, n):
                        if np.abs(mat[i, j]) > tol:
                            return False
                return True
            else:
                return False
            return False # should never happen

        v = np.identity(self.n)
        vecp, vecq = np.empty((self.n,)), np.empty((self.n,))
        def Jacobi_rotate(p: int, q: int) -> None:
            if (p < 0 or p >= n or q < 0 or q >= n):
                raise ValueError("Error! both p and q should be in [0, {:d}].".format(n - 1))

            pp = mat[p, p]
            qq = mat[q, q]
            pq = mat[p, q]
            theta = np.arctan(2.0 * pq / (pp - qq)) / 2.0 if np.abs(pp - qq) > 0.0 else np.pi / 4.0
            c = np.cos(theta)
            s = np.sin(theta)

            # update eigenvectors
            vecp[:] = v[:, p]
            vecq[:] = v[:, q]
            v[:, p] = vecq * s + vecp * c
            v[:, q] = vecq * c - vecp * s

            # update target matrix
            vecp[:] = mat[p, :]
            vecq[:] = mat[q, :]
            mat[p, :] = vecq * s + vecp * c
            mat[q, :] = vecq * c - vecp * s
            vecp[:] = mat[:, p]
            vecq[:] = mat[:, q]
            mat[:, p] = vecq * s + vecp * c
            mat[:, q] = vecq * c - vecp * s
            mat[p, p] = pp * c ** 2 + qq * s ** 2 + 2.0 * pq * s * c
            mat[q, q] = pp * s ** 2 + qq * c ** 2 - 2.0 * pq * s * c
            mat[p, q] = (qq - pp) * s * c + pq * (c ** 2 - s ** 2)

            return

        # maximum number of loops
        nloops_max = 30 * self.n ** 2

        iloop = 0
        # 30 is from OpenCV
        if not uplo.upper().find("U"):
            if is_greddy:
                while iloop < nloops_max:
                    p, q = 1, 0
                    for j in np.arange(1, self.n):
                        for i in np.arange(j):
                            if (np.abs(mat[i, j]) > np.abs(mat[p, q])):
                                p, q = i, j
                    if np.abs(mat[p, q]) <= const.tol_def:
                        break
                    Jacobi_rotate(p, q)
                    iloop += 1
            else:
                while iloop < nloops_max:
                    if Is_converged("U"):
                        break
                    for q in np.arange(1, self.n):
                        for p in np.arange(q):
                            Jacobi_rotate(p, q)
                            iloop += 1
        elif not uplo.upper().find("L"):
            if is_greddy:
                while iloop < nloops_max:
                    p, q = 1, 0
                    for j in np.arange(n - 1):
                        for i in np.arange(j + 1, n):
                            if (np.abs(mat[i, j]) > np.abs(mat[p, q])):
                                p, q = i, j
                    if np.abs(mat[p, q]) <= const.tol_def:
                        break
                    Jacobi_rotate(p, q)
                    iloop += 1
            else:
                while iloop < nloops_max:
                    if Is_converged("L"):
                        break
                    for q in np.arange(n - 1):
                        for p in np.arange(q + 1, n):
                            Jacobi_rotate(p, q)
                            iloop += 1
        else:
            raise ValueError(f"Error! uplo should be either \"U\" or \"L\"., but got {uplo:s}.") # should never happen

        print(f"iloop = {iloop:d}")
        w = mat.diagonal().copy()
        sort_arg = np.argsort(w)
        w = w[sort_arg]
        v = v[:, sort_arg]
        return w, v

if __name__ == "__main__":
    is_compare = True
    np.set_printoptions(threshold = sys.maxsize, suppress = True, formatter = {"float": "{:10.6f}".format})
    if argc - 1 == 1:
        n = int(sys.argv[1])
    else:
        n = const.n_def
    a = Square_matrix(n)
    a.mat = np.fromfile(const.mat_name, count = n * n).reshape((n, n), order = "F")
    # or
    # a.Init_linear()
    w, v = a.Jacobi_diagonalization() # is_greddy = True
    print("Matrix:")
    print(np.matrix(a.mat))
    print()
    print("Eigenvalues:")
    print(w)
    print()
    print("Eigenvectors:")
    print(v)
    print()
    if is_compare:
        w_ref, v_ref = np.linalg.eigh(a.mat)
        print("Reference:")
        print()
        print("Eigenvalues:")
        print(w_ref)
        print()
        print("Eigenvectors:")
        print(v_ref)
        print()
    
