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

    class MPCProblem
    {
        /*************************************

           MPC Problem:
           minimize x'*Q*x + u'*R*u
           subject to x(k+1) = A*x(k) + B*u(k)
                      x_min <= x <= x_max
                      u_min <= u <= u_max

           Q : State objective
           R : Control objective
           A : State dynamics
           B : Control dynamics

       *************************************/

    public:
        using Ptr = std::shared_ptr<MPCProblem>;

        const int num_states;
        const int num_controls;
        const int num_nodes;
        EigenVector x0, xf;
        EigenMatrix state_objective;
        EigenMatrix control_objective;
        EigenMatrix state_dynamics;
        EigenMatrix control_dynamics;
        EigenVector x_min, x_max;
        EigenVector u_min, u_max;

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

        MPCProblem(const int Nx, const int Nu, const int N);

        QPProblem::Ptr getQP();

        void updateQP();

        MPCSolution::Ptr getMPCSolution(QPSolution::Ptr qp_solution);

        void setInitialState(const EigenVector &x0);

        void setDesiredState(const EigenVector &xf);

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