#include "orlqp/QPArrayProblem.hpp"

namespace orlqp
{

    QPArrayProblem::QPArrayProblem(const std::vector<QPProblem::Ptr> &qp_array)
        : num_problems(qp_array.size()), qp_problems(std::move(qp_array))
    {
    }

    void QPArrayProblem::update()
    {
        for (int p = 0; p < this->num_problems; p++)
        {
            const QPProblem::Ptr QP_i = this->qp_problems[p];
            if (QP_i->update.hessian)
            {
                updateQPArrayHessian(p);
                QP_i->update.hessian = false;
            }
            if (QP_i->update.gradient)
            {
                updateQPArrayGradient(p);
                QP_i->update.gradient = false;
            }
            if (QP_i->update.linear_constraint)
            {
                updateQPArrayLinearConstraint(p);
                QP_i->update.linear_constraint = false;
            }
            if (QP_i->update.lower_bound)
            {
                updateQPArrayLowerBound(p);
                QP_i->update.lower_bound = false;
            }
            if (QP_i->update.upper_bound)
            {
                updateQPArrayUpperBound(p);
                QP_i->update.upper_bound = false;
            }
        }
    }

    QPProblem::Ptr QPArrayProblem::getQP()
    {
        if (!this->QP)
        {
            const int Np = this->num_problems;

            int n = 0, m = 0;
            this->variable_idx_map.resize(Np);
            this->constraint_idx_map.resize(Np);
            for (int p = 0; p < Np; p++)
            {
                this->variable_idx_map[p] = n;
                this->constraint_idx_map[p] = m;
                n += this->qp_problems[p]->num_variables;
                m += this->qp_problems[p]->num_constraints;
            }

            this->QP = std::make_shared<QPProblem>(n, m);

            calculateQPArrayHessian();
            calculateQPArrayGradient();
            calculateQPArrayLinearConstraint();
            calculateQPArrayLowerBound();
            calculateQPArrayUpperBound();
        }
        return this->QP;
    }

    std::vector<QPSolution::Ptr> QPArrayProblem::splitQPSolution(const QPSolution::Ptr qp_array_solution)
    {
        const int Np = this->num_problems;
        std::vector<QPSolution::Ptr> qp_solutions(this->num_problems);
        for (int p = 0; p < this->num_problems; p++)
        {
            QPSolution::Ptr qp_solution_i = std::make_shared<QPSolution>();

            qp_solution_i->run_time_s = qp_array_solution->run_time_s;
            qp_solution_i->setup_time_s = qp_array_solution->setup_time_s;
            qp_solution_i->solve_time_s = qp_array_solution->solve_time_s;

            const int n_i = this->qp_problems[p]->num_variables;
            const int var_idx = this->variable_idx_map[p];
            qp_solution_i->xstar.resize(n_i);
            qp_solution_i->xstar = qp_array_solution->xstar.block(var_idx, 0, n_i, 1);

            qp_solutions[p] = qp_solution_i;
        }
        return qp_solutions;
    }

    void QPArrayProblem::calculateQPArrayHessian()
    {
        std::vector<EigenTriplet> triplets;
        for (int p = 0; p < this->num_problems; p++)
        {
            const int var_idx = this->variable_idx_map[p];
            const EigenSparseMatrix &hessian_p = qp_problems[p]->hessian;
            triplets.reserve(hessian_p.nonZeros());
            for (int k = 0; k < hessian_p.outerSize(); k++)
                for (EigenSparseMatrix::InnerIterator it(hessian_p, k); it; ++it)
                    triplets.push_back(EigenTriplet(var_idx + it.row(), var_idx + it.col(), it.value()));
        }
        this->QP->hessian.setFromTriplets(triplets.begin(), triplets.end());
    }

    void QPArrayProblem::calculateQPArrayGradient()
    {
        for (int p = 0; p < this->num_problems; p++)
        {
            const int var_idx = this->variable_idx_map[p];
            this->QP->gradient.block(var_idx, 0, qp_problems[p]->num_variables, 1) = qp_problems[p]->gradient;
        }
    }

    void QPArrayProblem::calculateQPArrayLinearConstraint()
    {
        std::vector<EigenTriplet> triplets;
        for (int p = 0; p < this->num_problems; p++)
        {
            const int var_idx = this->variable_idx_map[p];
            const int con_idx = this->constraint_idx_map[p];
            const EigenSparseMatrix &linear_constraint_p = qp_problems[p]->linear_constraint;
            triplets.reserve(linear_constraint_p.nonZeros());
            for (int k = 0; k < linear_constraint_p.outerSize(); k++)
                for (EigenSparseMatrix::InnerIterator it(linear_constraint_p, k); it; ++it)
                    triplets.push_back(EigenTriplet(con_idx + it.row(), var_idx + it.col(), it.value()));
        }
        this->QP->linear_constraint.setFromTriplets(triplets.begin(), triplets.end());
    }

    void QPArrayProblem::calculateQPArrayLowerBound()
    {
        for (int p = 0; p < this->num_problems; p++)
        {
            const int con_idx = this->constraint_idx_map[p];
            this->QP->lower_bound.block(con_idx, 0, qp_problems[p]->num_constraints, 1) = qp_problems[p]->lower_bound;
        }
    }

    void QPArrayProblem::calculateQPArrayUpperBound()
    {
        for (int p = 0; p < this->num_problems; p++)
        {
            const int con_idx = this->constraint_idx_map[p];
            this->QP->upper_bound.block(con_idx, 0, qp_problems[p]->num_constraints, 1) = qp_problems[p]->upper_bound;
        }
    }

    void QPArrayProblem::updateQPArrayHessian(const int index)
    {
        const int var_idx = this->variable_idx_map[index];
        const EigenSparseMatrix &hessian_i = this->qp_problems[index]->hessian;
        for (int k = 0; k < this->QP->hessian.outerSize(); ++k)
            for (EigenSparseMatrix::InnerIterator it(this->QP->hessian, k); it; ++it)
                if (it.row() >= var_idx && it.col() >= var_idx &&
                    it.row() < var_idx + hessian_i.rows() && it.col() < var_idx + hessian_i.cols())
                    it.valueRef() = 0;

        for (int k = 0; k < hessian_i.outerSize(); ++k)
            for (EigenSparseMatrix::InnerIterator it(hessian_i, k); it; ++it)
                this->QP->hessian.coeffRef(it.row() + var_idx, it.col() + var_idx) = it.value();

        this->QP->hessian.prune([](const int &, const int &, const double &value)
                                { return value != 0.0; });
        this->QP->update.hessian = true;
    }

    void QPArrayProblem::updateQPArrayGradient(const int index)
    {
        const int var_idx = this->variable_idx_map[index];
        const EigenVector &gradient_i = this->qp_problems[index]->gradient;
        this->QP->gradient.block(var_idx, 0, gradient_i.rows(), 1) = gradient_i;
        this->QP->update.gradient = true;
    }

    void QPArrayProblem::updateQPArrayLinearConstraint(const int index)
    {
        const int var_idx = this->variable_idx_map[index];
        const int con_idx = this->constraint_idx_map[index];
        const EigenSparseMatrix &linear_constraint_i = this->qp_problems[index]->linear_constraint;
        for (int k = 0; k < this->QP->linear_constraint.outerSize(); ++k)
            for (EigenSparseMatrix::InnerIterator it(this->QP->linear_constraint, k); it; ++it)
                if (it.row() >= var_idx && it.col() >= var_idx &&
                    it.row() < var_idx + linear_constraint_i.rows() && it.col() < var_idx + linear_constraint_i.cols())
                    it.valueRef() = 0;

        for (int k = 0; k < linear_constraint_i.outerSize(); ++k)
            for (EigenSparseMatrix::InnerIterator it(linear_constraint_i, k); it; ++it)
                this->QP->linear_constraint.coeffRef(it.row() + var_idx, it.col() + var_idx) = it.value();

        this->QP->linear_constraint.prune([](const int &, const int &, const double &value)
                                          { return value != 0.0; });
        this->QP->update.linear_constraint = true;
    }

    void QPArrayProblem::updateQPArrayLowerBound(const int index)
    {
        const int con_idx = this->constraint_idx_map[index];
        const EigenVector &lower_bound_i = this->qp_problems[index]->lower_bound;
        this->QP->lower_bound.block(con_idx, 0, lower_bound_i.rows(), 1) = lower_bound_i;
        this->QP->update.lower_bound = true;
    }

    void QPArrayProblem::updateQPArrayUpperBound(const int index)
    {
        const int con_idx = this->constraint_idx_map[index];
        const EigenVector &upper_bound_i = this->qp_problems[index]->upper_bound;
        this->QP->upper_bound.block(con_idx, 0, upper_bound_i.rows(), 1) = upper_bound_i;
        this->QP->update.upper_bound = true;
    }

}
