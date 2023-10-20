#ifndef ORLQP_OSQP_SOLVER_HPP_
#define ORLQP_OSQP_SOLVER_HPP_

#include <execution>

#include "orlqp/types.hpp"
#include "osqp/osqp.h"
#include "orlqp/QPProblem.hpp"

namespace orlqp
{

    class OSQP
    {
    public:
        using Ptr = std::shared_ptr<OSQP>;

        OSQP();

        ~OSQP();

        bool isOK() { return ok; }

        OSQPInt solve();

        bool isSetup() { return is_setup; }

        OSQPInt setup(QPProblem::Ptr qp);

        OSQPInt update();

        QPSolution::Ptr getQPSolution();

        OSQPSolver *getSolver()
        {
            return solver;
        }

        OSQPSettings *getSettings()
        {
            this->update_settings = true;
            return settings;
        }

    private:
        QPProblem::Ptr QP;

        OSQPInt n, m;

        OSQPSolver *solver = nullptr;
        OSQPSettings *settings = nullptr;

        bool is_setup = false;
        bool ok = true;
        bool update_settings = false;

        OSQPFloat *q = nullptr;
        OSQPFloat *l = nullptr;
        OSQPFloat *u = nullptr;

        OSQPCscMatrix *P = nullptr;
        OSQPInt Pnnz;
        OSQPFloat *Px = nullptr;
        OSQPInt *Pi = nullptr;
        OSQPInt *Pp = nullptr;

        OSQPCscMatrix *A = nullptr;
        OSQPInt Annz;
        OSQPFloat *Ax = nullptr;
        OSQPInt *Ai = nullptr;
        OSQPInt *Ap = nullptr;

        void convertEigenSparseToCSC(const EigenSparseMatrix &matrix,
                                     OSQPCscMatrix *&M, OSQPInt &Mnnz, OSQPFloat *&Mx, OSQPInt *&Mi, OSQPInt *&Mp);
    };

}

#endif