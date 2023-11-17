#ifndef ORLQP_MPC_PROBLEM_HPP_
#define ORLQP_MPC_PROBLEM_HPP_

#include "orlqp/types.hpp"
#include "orlqp/QPProblem.hpp"

namespace orlqp
{

    struct MPCSolution
    {
        using Ptr = std::shared_ptr<MPCSolution>;

        float run_time_s;
        float setup_time_s;
        float solve_time_s;
        EigenVector xstar;
        EigenVector ustar;
    };

    class ExpandedMPCProblem
    {
        /*************************************

           Expanded MPC Problem:
           minimize x(k)'*Q(k)*x(k) + u(k)'*R(k)*u(k)
           subject to x(k+1) = A(k)*x(k) + B(k)*u(k)
                      x_min(k) <= x(k) <= x_max(k)
                      u_min(k) <= u(k) <= u_max(k)

           Q : State objective
           R : Control objective
           A : State dynamics
           B : Control dynamics

       *************************************/

    public:
        using Ptr = std::shared_ptr<ExpandedMPCProblem>;

        const int num_states;
        const int num_controls;
        const int num_nodes;
        EigenVector x0, xf, uf;
        std::vector<EigenMatrix> state_objective;
        std::vector<EigenMatrix> control_objective;
        std::vector<EigenMatrix> state_dynamics;
        std::vector<EigenMatrix> control_dynamics;
        std::vector<EigenVector> x_min, x_max;
        std::vector<EigenVector> u_min, u_max;

        struct
        {
            bool state_objective = false,
                 control_objective = false,
                 state_dynamics = false,
                 control_dynamics = false,
                 x_bounds = false,
                 u_bounds = false,
                 x0 = false,
                 xf = false;
        } update;

        ExpandedMPCProblem(const int Nx, const int Nu, const int N);

        QPProblem::Ptr getQP();

        void updateQP();

        MPCSolution::Ptr getMPCSolution(QPSolution::Ptr qp_solution);

        void setInitialState(const EigenVector &x0);

        void setDesiredState(const EigenVector &xf);

        void setDesiredControl(const EigenVector &uf);

    private:
        QPProblem::Ptr QP;

        void calculateQPHessian();
        void calculateQPGradient();
        void calculateQPLinearConstraint();
        void calculateQPLowerBound();
        void calculateQPUpperBound();
    };

}

#endif