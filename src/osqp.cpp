#include "orlqp/OSQP.hpp"

namespace orlqp
{

    OSQP::OSQP()
    {
        this->settings = new OSQPSettings;
        osqp_set_default_settings(this->settings);
    }

    OSQP::~OSQP()
    {
        osqp_cleanup(solver);
        delete P;
        delete A;
        delete settings;
    }

    OSQPInt OSQP::solve()
    {
        if (!this->is_setup)
            return -1;

        return osqp_solve(this->solver);
    }

    OSQPInt OSQP::setup(QPProblem::Ptr qp)
    {

        this->QP = qp;

        if (this->is_setup)
            osqp_cleanup(this->solver);
        else
            this->is_setup = true;

        this->n = qp->num_variables;
        this->m = qp->num_constraints;

        this->q = qp->gradient.data();
        this->l = qp->lower_bound.data();
        this->u = qp->upper_bound.data();

        convertEigenSparseToCSC(qp->hessian, this->P, this->Pnnz, this->Px, this->Pi, this->Pp);
        convertEigenSparseToCSC(qp->linear_constraint, this->A, this->Annz, this->Ax, this->Ai, this->Ap);

        return osqp_setup(&this->solver, this->P, this->q, this->A, this->l, this->u, this->m, this->n, this->settings);
    }

    OSQPInt OSQP::update()
    {
        OSQPInt flag;
        if (QP->update.hessian)
        {
            QP->update.hessian = false;
            convertEigenSparseToCSC(QP->hessian, this->P, this->Pnnz, this->Px, this->Pi, this->Pp);
            flag = osqp_update_data_mat(this->solver, this->Px, this->Pi, this->Pnnz, OSQP_NULL, OSQP_NULL, OSQP_NULL);
        }
        if (QP->update.gradient)
        {
            QP->update.gradient = false;
            this->q = QP->gradient.data();
            flag = osqp_update_data_vec(this->solver, this->q, OSQP_NULL, OSQP_NULL);
        }
        if (QP->update.linear_constraint)
        {
            QP->update.linear_constraint = false;
            convertEigenSparseToCSC(QP->linear_constraint, this->A, this->Annz, this->Ax, this->Ai, this->Ap);
            flag = osqp_update_data_mat(this->solver, this->Px, this->Pi, this->Pnnz, OSQP_NULL, OSQP_NULL, OSQP_NULL);
        }
        if (QP->update.lower_bound && QP->update.upper_bound)
        {
            QP->update.lower_bound = false;
            QP->update.upper_bound = false;
            this->u = QP->upper_bound.data();
            this->l = QP->lower_bound.data();
            flag = osqp_update_data_vec(this->solver, OSQP_NULL, this->l, this->u);
        }
        if (QP->update.lower_bound)
        {
            QP->update.lower_bound = false;
            this->l = QP->lower_bound.data();
            flag = osqp_update_data_vec(this->solver, OSQP_NULL, this->l, OSQP_NULL);
        }
        if (QP->update.upper_bound)
        {
            QP->update.upper_bound = false;
            this->u = QP->upper_bound.data();
            flag = osqp_update_data_vec(this->solver, OSQP_NULL, OSQP_NULL, this->u);
        }
        if (this->update_settings)
        {
            this->update_settings = false;
            flag = osqp_update_settings(this->solver, this->settings);
        }
        return flag;
    }

    QPSolution::Ptr OSQP::getQPSolution()
    {
        QPSolution::Ptr qp_solution = std::make_shared<QPSolution>();
        qp_solution->xstar = Eigen::Map<EigenVector>(this->solver->solution->x, this->n);
        qp_solution->run_time_s = this->solver->info->run_time;
        qp_solution->setup_time_s = this->solver->info->setup_time;
        qp_solution->solve_time_s = this->solver->info->solve_time;
        return qp_solution;
    }

    void OSQP::convertEigenSparseToCSC(const EigenSparseMatrix &matrix,
                                       OSQPCscMatrix *&M, OSQPInt &Mnnz, OSQPFloat *&Mx, OSQPInt *&Mi, OSQPInt *&Mp)
    {
        M = new OSQPCscMatrix;
        Mnnz = matrix.nonZeros();
        Mx = new OSQPFloat[Mnnz];
        Mi = new OSQPInt[Mnnz];
        Mp = new OSQPInt[matrix.cols() + 1];

        int k = 0;
        Mp[0] = 0;
        for (int j = 0; j < matrix.outerSize(); ++j)
        {
            for (EigenSparseMatrix::InnerIterator it(matrix, j); it; ++it)
            {
                Mx[k] = it.value();
                Mi[k] = it.row();
                ++k;
            }
            Mp[j + 1] = k;
        }
        csc_set_data(M, matrix.rows(), matrix.cols(), Mnnz, Mx, Mi, Mp);
    }

}