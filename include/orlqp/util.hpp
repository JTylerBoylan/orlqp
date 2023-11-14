#ifndef ORLQP_UTIL_HPP_
#define ORLQP_UTIL_HPP_

#include <iostream>

#include "orlqp/types.hpp"
#include "osqp/osqp.h"

namespace orlqp
{

    void deleteOSQPCscMatrix(OSQPCscMatrix *M)
    {
        delete[] M->p;
        delete[] M->i;
        delete[] M->x;
        delete M;
    }

    void convertEigenSparseToCSC(const EigenSparseMatrix &matrix, OSQPCscMatrix *&M)
    {
        if (M)
            deleteOSQPCscMatrix(M);

        M = new OSQPCscMatrix;
        M->nzmax = matrix.nonZeros();
        M->m = matrix.rows();
        M->n = matrix.cols();
        M->p = new OSQPInt[M->n + 1];
        M->i = new OSQPInt[M->nzmax];
        M->x = new OSQPFloat[M->nzmax];
        M->nz = -1;

        int k = 0;
        M->p[0] = 0;
        for (int j = 0; j < matrix.outerSize(); ++j)
        {
            for (EigenSparseMatrix::InnerIterator it(matrix, j); it; ++it)
            {
                M->x[k] = it.value();
                M->i[k] = it.row();
                ++k;
            }
            M->p[j + 1] = k;
        }
    }

    void printOSQPCscMatrix(const OSQPCscMatrix *mat)
    {
        if (mat != nullptr)
        {
            std::cout << "Matrix dimensions: " << mat->m << " x " << mat->n << std::endl;
            std::cout << "Number of non-zero values: " << mat->nzmax << std::endl;
            std::cout << "Column pointers (p): ";
            for (OSQPInt col = 0; col <= mat->n; ++col)
            {
                std::cout << mat->p[col] << " ";
            }
            std::cout << "\nRow indices (i): ";
            for (OSQPInt nz = 0; nz < mat->nzmax; ++nz)
            {
                std::cout << mat->i[nz] << " ";
            }
            std::cout << "\nNon-zero values (x): ";
            for (OSQPInt nz = 0; nz < mat->nzmax; ++nz)
            {
                std::cout << mat->x[nz] << " ";
            }
            std::cout << std::endl;
        }
        else
        {
            std::cout << "Matrix is nullptr." << std::endl;
        }
    }

}

#endif