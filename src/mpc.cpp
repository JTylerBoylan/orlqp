#include "orlqp/MPCProblem.hpp"

namespace orlqp
{

    MPCProblem::MPCProblem(const int Nx, const int Nu, const int N)
        : num_states(Nx), num_controls(Nu), num_nodes(N)
    {
        this->x0.resize(Nx);
        this->xf.resize(Nx);
        this->uf.resize(Nu);
        this->state_objective.resize(Nx, Nx);
        this->control_objective.resize(Nu, Nu);
        this->state_dynamics.resize(Nx, Nx);
        this->control_dynamics.resize(Nx, Nu);
        this->x_min.resize(Nx);
        this->x_max.resize(Nx);
        this->u_min.resize(Nu);
        this->u_max.resize(Nu);
    }

    QPProblem::Ptr MPCProblem::getQP()
    {
        if (!this->QP)
        {
            const int Nx = this->num_states;
            const int Nu = this->num_controls;
            const int Nn = this->num_nodes;
            const int n = Nx * (Nn + 1) + Nu * Nn;
            const int m = 2 * Nx * (Nn + 1) + Nu * Nn;
            this->QP = std::make_shared<QPProblem>(n, m);
            this->calculateQPHessian();
            this->calculateQPGradient();
            this->calculateQPLinearConstraint();
            this->calculateQPLowerBound();
            this->calculateQPUpperBound();
        }
        return this->QP;
    }

    void MPCProblem::updateQP()
    {
        if (this->update.state_objective || this->update.control_objective)
        {
            this->update.control_objective = false;
            this->calculateQPHessian();
            this->QP->update.hessian = true;
        }
        if (this->update.state_objective || this->update.xf)
        {
            this->update.state_objective = false;
            this->update.xf = false;
            this->calculateQPGradient();
            this->QP->update.gradient = true;
        }
        if (this->update.state_dynamics || this->update.control_dynamics)
        {
            this->update.state_dynamics = false;
            this->update.control_dynamics = false;
            this->calculateQPLinearConstraint();
            this->QP->update.linear_constraint = true;
        }
        if (this->update.x_bounds || this->update.u_bounds || this->update.x0)
        {
            this->update.x_bounds = false;
            this->update.u_bounds = false;
            this->update.x0 = false;
            this->calculateQPLowerBound();
            this->calculateQPUpperBound();
            this->QP->update.lower_bound = true;
            this->QP->update.upper_bound = true;
        }
    }

    MPCSolution::Ptr MPCProblem::getMPCSolution(QPSolution::Ptr qp_solution)
    {
        const int Nx = this->num_states;
        const int Nu = this->num_controls;
        const int Nn = this->num_nodes;
        auto mpc_solution = std::make_shared<MPCSolution>();
        mpc_solution->run_time_s = qp_solution->run_time_s;
        mpc_solution->setup_time_s = qp_solution->setup_time_s;
        mpc_solution->solve_time_s = qp_solution->solve_time_s;
        mpc_solution->xstar = qp_solution->xstar.block(0, 0, Nx * (Nn + 1), 1);
        mpc_solution->ustar = qp_solution->xstar.block(Nx * (Nn + 1), 0, Nu * Nn, 1);
        return mpc_solution;
    }

    void MPCProblem::setInitialState(const EigenVector &x0)
    {
        this->x0 = x0;
        const int Nx = this->num_states;
        this->QP->lower_bound.block(0, 0, Nx, 1) = -x0;
        this->QP->upper_bound.block(0, 0, Nx, 1) = -x0;
        this->QP->update.lower_bound = true;
        this->QP->update.upper_bound = true;
    }

    void MPCProblem::setDesiredState(const EigenVector &xf)
    {
        this->xf = xf;
        this->calculateQPGradient();
        this->QP->update.gradient = true;
    }

    void MPCProblem::setDesiredControl(const EigenVector &uf)
    {
        this->uf = uf;
        this->calculateQPGradient();
        this->QP->update.gradient = true;
    }

    void MPCProblem::calculateQPHessian()
    {
        const int Nx = this->num_states;
        const int Nu = this->num_controls;
        const int Nn = this->num_nodes;
        std::vector<EigenTriplet> triplets;
        for (int i = 0; i < Nn + 1; i++)
        {
            const int idx = Nx * i;
            for (int j = 0; j < Nx; j++)
                for (int k = 0; k < Nx; k++)
                {
                    const Float value = this->state_objective(j, k);
                    if (value != 0)
                        triplets.push_back(EigenTriplet(idx + k, idx + k, value));
                }
        }
        for (int i = 0; i < Nn; i++)
        {
            const int idx = Nx * (Nn + 1) + Nu * i;
            for (int j = 0; j < Nu; j++)
                for (int k = 0; k < Nu; k++)
                {
                    const Float value = this->control_objective(j, k);
                    if (value != 0)
                        triplets.push_back(EigenTriplet(idx + k, idx + k, value));
                }
        }
        this->QP->hessian.setFromTriplets(triplets.begin(), triplets.end());
    }

    void MPCProblem::calculateQPGradient()
    {
        const int Nx = this->num_states;
        const int Nu = this->num_controls;
        const int Nn = this->num_nodes;
        const int Nnx = Nx * (Nn + 1);
        const EigenVector Gx = this->state_objective * -(this->xf);
        const EigenVector Gu = this->control_objective * -(this->uf);
        for (int xi = 0; xi < Nnx; xi++)
        {
            int gIdx = xi % Nx;
            this->QP->gradient(xi, 0) = Gx(gIdx);
        }
        for (int ui = 0; ui < Nu * Nn; ui++)
        {
            int gIdx = ui % Nu;
            this->QP->gradient(Nnx + ui, 0) = Gu(gIdx);
        }
    }

    void MPCProblem::calculateQPLinearConstraint()
    {
        const int n = this->QP->num_variables;
        const int Nx = this->num_states;
        const int Nu = this->num_controls;
        const int Nn = this->num_nodes;
        std::vector<EigenTriplet> triplets;
        for (int i = 0; i < Nx * (Nn + 1); i++)
            triplets.push_back(EigenTriplet(i, i, -1));

        for (int i = 0; i < Nn; i++)
            for (int j = 0; j < Nx; j++)
                for (int k = 0; k < Nx; k++)
                {
                    const Float value = this->state_dynamics(j, k);
                    if (value != 0)
                        triplets.push_back(EigenTriplet(
                            Nx * (i + 1) + j,
                            Nx * i + k,
                            value));
                }

        for (int i = 0; i < Nn; i++)
            for (int j = 0; j < Nx; j++)
                for (int k = 0; k < Nu; k++)
                {
                    const Float value = this->control_dynamics(j, k);
                    if (value != 0)
                        triplets.push_back(EigenTriplet(
                            Nx * (i + 1) + j,
                            Nu * i + k + Nx * (Nn + 1),
                            value));
                }

        for (int i = 0; i < n; i++)
            triplets.push_back(EigenTriplet(i + (Nn + 1) * Nx, i, 1));

        this->QP->linear_constraint.setFromTriplets(triplets.begin(), triplets.end());
    }

    void MPCProblem::calculateQPLowerBound()
    {
        const int Nx = this->num_states;
        const int Nu = this->num_controls;
        const int Nn = this->num_nodes;
        const int Neq = Nx * (Nn + 1);
        const int Nineq = Nx * (Nn + 1) + Nu * Nn;

        EigenVector lower_inequality = EigenVector::Zero(Nineq);
        for (int i = 0; i < Nn + 1; i++)
            lower_inequality.block(Nx * i, 0, Nx, 1) = x_min;
        for (int i = 0; i < Nn; i++)
            lower_inequality.block(Nu * i + Nx * (Nn + 1), 0, Nu, 1) = u_min;

        EigenVector lower_equality = EigenVector::Zero(Neq);
        lower_equality.block(0, 0, Nx, 1) = -x0;

        this->QP->lower_bound.block(0, 0, Neq, 1) = lower_equality;
        this->QP->lower_bound.block(Neq, 0, Nineq, 1) = lower_inequality;
    }

    void MPCProblem::calculateQPUpperBound()
    {
        const int Nx = this->num_states;
        const int Nu = this->num_controls;
        const int Nn = this->num_nodes;
        const int Neq = Nx * (Nn + 1);
        const int Nineq = Nx * (Nn + 1) + Nu * Nn;

        EigenVector upper_inequality = EigenVector::Zero(Nineq);
        for (int i = 0; i < Nn + 1; i++)
            upper_inequality.block(Nx * i, 0, Nx, 1) = x_max;
        for (int i = 0; i < Nn; i++)
            upper_inequality.block(Nu * i + Nx * (Nn + 1), 0, Nu, 1) = u_max;

        EigenVector upper_equality = EigenVector::Zero(Neq);
        upper_equality.block(0, 0, Nx, 1) = -x0;

        this->QP->upper_bound.block(0, 0, Neq, 1) = upper_equality;
        this->QP->upper_bound.block(Neq, 0, Nineq, 1) = upper_inequality;
    }

}
