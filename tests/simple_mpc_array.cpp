#include <iostream>
#include <chrono>

#include "orlqp/MPCProblem.hpp"
#include "orlqp/OSQP.hpp"
#include "orlqp/QPArrayProblem.hpp"

#define NUMBER_OF_PROBLEMS 10

#define NUMBER_OF_STATES 2
#define NUMBER_OF_CONTROLS 1

#define NUMBER_OF_NODES 11

#define POSITION_ERROR_COST_WEIGHT 10.0
#define VELOCITY_ERROR_COST_WEIGHT 1.0
#define FORCE_COST_WEIGHT 0.0

#define HORIZON_TIME 1.0
#define MASS 1.0

#define MIN_POSITION -10.0
#define MAX_POSITION +10.0
#define MIN_VELOCITY -10.0
#define MAX_VELOCITY +10.0
#define MIN_FORCE -5.0
#define MAX_FORCE +5.0

#define RANDOM_NOISE_GAIN 0.5
#define NUMBER_OF_MPC_ITERATIONS 10000

using namespace orlqp;

int main()
{

    std::vector<MPCProblem::Ptr> mpcs(NUMBER_OF_PROBLEMS);
    std::vector<QPProblem::Ptr> qps(NUMBER_OF_PROBLEMS);
    for (int p = 0; p < NUMBER_OF_PROBLEMS; p++)
    {
        EigenVector x0(NUMBER_OF_STATES), xf(NUMBER_OF_STATES);
        x0 << 5.0, 0.0;
        xf << 0.0, 0.0;

        MPCProblem::Ptr mpc = std::make_shared<MPCProblem>(NUMBER_OF_STATES, NUMBER_OF_CONTROLS, NUMBER_OF_NODES);
        mpc->x0 = x0;
        mpc->xf = xf;
        mpc->state_objective << POSITION_ERROR_COST_WEIGHT, 0.0, 0.0, VELOCITY_ERROR_COST_WEIGHT;
        mpc->control_objective << FORCE_COST_WEIGHT;
        mpc->state_dynamics << 1.0, (HORIZON_TIME / NUMBER_OF_NODES), 0.0, 1.0;
        mpc->control_dynamics << 0.0, (HORIZON_TIME / NUMBER_OF_NODES) / MASS;
        mpc->x_min << MIN_POSITION, MIN_VELOCITY;
        mpc->x_max << MAX_POSITION, MAX_VELOCITY;
        mpc->u_min << MIN_FORCE;
        mpc->u_max << MAX_FORCE;

        mpcs[p] = mpc;
        qps[p] = mpc->getQP();
    }

    QPArrayProblem::Ptr qp_array = std::make_shared<QPArrayProblem>(qps);

    OSQP::Ptr osqp = std::make_shared<OSQP>();
    osqp->getSettings()->verbose = false;
    osqp->getSettings()->warm_starting = true;
    osqp->getSettings()->polishing = true;

    osqp->setup(qp_array->getQP());

    const auto cstart = std::chrono::high_resolution_clock::now();
    for (int k = 1; k <= NUMBER_OF_MPC_ITERATIONS; k++)
    {

        osqp->solve();

        const QPSolution::Ptr qp_solution = osqp->getQPSolution();
        const std::vector<QPSolution::Ptr> qp_solution_array = qp_array->splitQPSolution(qp_solution);
        for (int p = 0; p < NUMBER_OF_PROBLEMS; p++)
        {
            const QPSolution::Ptr qp_solution_i = qp_solution_array[p];
            const MPCSolution::Ptr mpc_solution_i = mpcs[p]->getMPCSolution(qp_solution_i);
            const MPCProblem::Ptr mpc_i = mpcs[p];

            const float rand_float = (float)(rand()) / (float)(RAND_MAX);
            const float rand_force = RANDOM_NOISE_GAIN * (rand_float - 0.5F);
            const Float u0 = mpc_solution_i->ustar(0, 0);
            const EigenVector x0 = mpc_i->state_dynamics * mpc_i->x0 +
                                   mpc_i->control_dynamics * (u0 + rand_force);

            if (k % 10000 == 0)
            {
                std::cout << "x0 = [" << x0.transpose() << "]\n";
            }

            mpc_i->setInitialState(x0);
        }

        qp_array->update();
        osqp->update();

        const auto cend = std::chrono::high_resolution_clock::now();

        if (k % 10000 == 0)
        {
            const time_t duration = std::chrono::duration_cast<std::chrono::microseconds>(cend - cstart).count();
            const double kHz = (double)(k * 1E3) / (double)(duration);
            std::cout << "f = " << kHz << " kHz\n";
            std::cout << "----------------------\n";
        }
    }

    return 0;
}