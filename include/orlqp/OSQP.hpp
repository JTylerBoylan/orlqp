#ifndef ORLQP_OSQP_SOLVER_HPP_
#define ORLQP_OSQP_SOLVER_HPP_

#include <execution>

#include "orlqp/types.hpp"
#include "osqp/osqp.h"
#include "orlqp/QPProblem.hpp"
#include "orlqp/util.hpp"

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

        void printPMatrix();
        void printAMatrix();

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
        OSQPCscMatrix *A = nullptr;
    };

}

#endif