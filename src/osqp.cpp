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
        deleteOSQPCscMatrix(P);
        deleteOSQPCscMatrix(A);
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

        auto ut_qp = this->QP->hessian.triangularView<Eigen::Upper>();
        convertEigenSparseToCSC(ut_qp, this->P);
        convertEigenSparseToCSC(qp->linear_constraint, this->A);

        return osqp_setup(&this->solver, this->P, this->q, this->A, this->l, this->u, this->m, this->n, this->settings);
    }

    OSQPInt OSQP::update()
    {
        OSQPInt flag;
        if (QP->update.hessian)
        {
            auto ut_qp = this->QP->hessian.triangularView<Eigen::Upper>();
            convertEigenSparseToCSC(ut_qp, this->P);
            flag = osqp_update_data_mat(this->solver, this->P->x, this->P->i, this->P->nzmax, OSQP_NULL, OSQP_NULL, OSQP_NULL);
            if (flag)
                return flag;
            QP->update.hessian = false;
        }
        if (QP->update.gradient)
        {
            this->q = QP->gradient.data();
            flag = osqp_update_data_vec(this->solver, this->q, OSQP_NULL, OSQP_NULL);
            if (flag)
                return flag;
            QP->update.gradient = false;
        }
        if (QP->update.linear_constraint)
        {
            convertEigenSparseToCSC(QP->linear_constraint, this->A);
            flag = osqp_update_data_mat(this->solver, OSQP_NULL, OSQP_NULL, OSQP_NULL, this->A->x, this->A->i, this->A->nzmax);
            if (flag)
                return flag;
            QP->update.linear_constraint = false;
        }
        if (QP->update.lower_bound && QP->update.upper_bound)
        {
            this->u = QP->upper_bound.data();
            this->l = QP->lower_bound.data();
            flag = osqp_update_data_vec(this->solver, OSQP_NULL, this->l, this->u);
            if (flag)
                return flag;
            QP->update.lower_bound = false;
            QP->update.upper_bound = false;
        }
        if (QP->update.lower_bound)
        {
            QP->update.lower_bound = false;
            this->l = QP->lower_bound.data();
            flag = osqp_update_data_vec(this->solver, OSQP_NULL, this->l, OSQP_NULL);
            if (flag)
                return flag;
        }
        if (QP->update.upper_bound)
        {
            this->u = QP->upper_bound.data();
            flag = osqp_update_data_vec(this->solver, OSQP_NULL, OSQP_NULL, this->u);
            if (flag)
                return flag;
            QP->update.upper_bound = false;
        }
        if (this->update_settings)
        {
            flag = osqp_update_settings(this->solver, this->settings);
            if (flag)
                return flag;
            this->update_settings = false;
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

    void OSQP::printPMatrix()
    {
        printOSQPCscMatrix(this->P);
    }
    void OSQP::printAMatrix()
    {
        printOSQPCscMatrix(this->A);
    }

}