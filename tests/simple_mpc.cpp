#include <iostream>
#include <chrono>

#include "orlqp/types.hpp"
#include "orlqp/OSQP.hpp"
#include "orlqp/MPCProblem.hpp"

using namespace orlqp;

#define NUM_PROBLEMS 3

#define NUM_STATES 2
#define NUM_CONTROLS 1

#define NUM_NODES 11

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
#define NUMBER_OF_MPC_ITERATIONS 1000000

int main()
{

    srand(time(NULL));

    EigenVector x0(NUM_STATES), xf(NUM_STATES);
    x0 << 5.0, 0.0;
    xf << 0.0, 0.0;

    MPCProblem::Ptr double_integrator_mpc = std::make_shared<MPCProblem>(NUM_STATES, NUM_CONTROLS, NUM_NODES);
    double_integrator_mpc->x0 = x0;
    double_integrator_mpc->xf = xf;
    double_integrator_mpc->state_objective << POSITION_ERROR_COST_WEIGHT, 0.0, 0.0, VELOCITY_ERROR_COST_WEIGHT;
    double_integrator_mpc->control_objective << FORCE_COST_WEIGHT;
    double_integrator_mpc->state_dynamics << 1.0, (HORIZON_TIME / NUM_NODES), 0.0, 1.0;
    double_integrator_mpc->control_dynamics << 0.0, (HORIZON_TIME / NUM_NODES) / MASS;
    double_integrator_mpc->x_min << MIN_POSITION, MIN_VELOCITY;
    double_integrator_mpc->x_max << MAX_POSITION, MAX_VELOCITY;
    double_integrator_mpc->u_min << MIN_FORCE;
    double_integrator_mpc->u_max << MAX_FORCE;

    OSQP::Ptr osqp = std::make_shared<OSQP>();
    osqp->getSettings()->verbose = false;
    osqp->getSettings()->warm_starting = true;
    osqp->getSettings()->polishing = true;

    osqp->setup(double_integrator_mpc->getQP());

    const auto cstart = std::chrono::high_resolution_clock::now();
    for (int k = 1; k <= NUMBER_OF_MPC_ITERATIONS; k++)
    {

        osqp->solve();

        QPSolution::Ptr qp_solution = osqp->getQPSolution();
        MPCSolution::Ptr mpc_solution = double_integrator_mpc->getMPCSolution(qp_solution);

        const float rand_float = (float)(rand()) / (float)(RAND_MAX);
        const float rand_force = RANDOM_NOISE_GAIN * (rand_float - 0.5F);

        const Float u0 = mpc_solution->ustar(0, 0);
        x0 = double_integrator_mpc->state_dynamics * x0 + double_integrator_mpc->control_dynamics * (u0 + rand_force);

        double_integrator_mpc->setInitialState(x0);

        osqp->update();

        const auto cend = std::chrono::high_resolution_clock::now();

        if (k % 100000 == 0)
        {
            const time_t duration = std::chrono::duration_cast<std::chrono::microseconds>(cend - cstart).count();
            const double kHz = (double)(k * 1E3) / (double)(duration);
            std::cout << "----------------------\n";
            std::cout << "f = " << kHz << " kHz\n";
            std::cout << "x0: [" << x0.transpose() << "]\n"
                      << "u0: [" << u0 << "]\n";
            std::cout << "----------------------\n";
        }
    }

    return 0;
}
