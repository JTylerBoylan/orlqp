#include <iostream>

#include "orlqp/orlqp.hpp"
#include "orlqp/symbolic/util.hpp"
#include "orlqp/SymbolicQPProblem.hpp"

#include <chrono>
#include <time.h>

#define INF 1E8

using namespace orlqp;

int main(void)
{

    // Params
    const int nodes_per_step = 2;
    const int num_steps = 2;
    const int num_nodes = nodes_per_step * num_steps;

    /* Objective Function */

    const int num_states = num_nodes * 4;
    auto q = createSymbolVector(num_states, "q");
    auto q_ref = createSymbolVector(num_states, "Qr");      // const
    auto w_target = createSymbolVector(num_states, "Wtar"); // const
    GinacEx J_target;
    for (int i = 0; i < num_states; i++)
    {
        J_target += w_target[i] * (q[i] - q_ref[i]) * (q[i] - q_ref[i]);
    }

    auto sx = createSymbolVector(num_steps, "sx");
    auto sx_ref = createSymbolVector(num_steps, "Sxr");    // const
    auto w_stepx = createSymbolVector(num_steps, "Wstpx"); // const
    GinacEx J_stepx;
    for (int i = 0; i < num_steps; i++)
    {
        J_stepx += w_stepx[i] * (sx[i] - sx_ref[i]) * (sx[i] - sx_ref[i]);
    }

    auto sy = createSymbolVector(num_steps, "sy");
    auto rho_step = createSymbolVector(num_nodes, "rho_s");
    auto sy_ref = createSymbolVector(num_steps, "Syr");    // const
    auto w_stepy = createSymbolVector(num_nodes, "Wstpy"); // const
    GinacEx J_stepy;
    for (int i = 0; i < num_nodes; i++)
    {
        J_stepy += w_stepy[i] * rho_step[i] * rho_step[i];
    }

    const int num_controls = num_nodes * 2;
    auto u = createSymbolVector(num_controls, "u");
    auto p_foot_x = GinacSymbol("Pftx"); // const
    auto p_foot_y = GinacSymbol("Pfty"); // const
    GinacMatrix u_ref(num_controls, 1);
    for (int i = 0; i < num_nodes; i++)
    {
        const int step_i = std::floor(i / nodes_per_step);
        GinacEx sum_sx, sum_sy;
        for (int j = 0; j < step_i; j++)
        {
            sum_sx += sx[j];
            sum_sy += sy[j];
        }
        u_ref[2 * i] = p_foot_x + sum_sx;
        u_ref[2 * i + 1] = p_foot_y + sum_sy;
    }
    auto w_effort = createSymbolVector(num_controls, "Weff"); // const
    GinacEx J_effort;
    for (int i = 0; i < num_controls; i++)
    {
        J_effort += w_effort[i] * (u[i] - u_ref[i]) * (u[i] - u_ref[i]);
    }

    auto rho_input = createSymbolVector(num_nodes, "rho_i");
    auto w_input = createSymbolVector(num_nodes, "Winp"); // const
    GinacEx J_input;
    for (int i = 0; i < num_nodes; i++)
    {
        J_input += w_input[i] * rho_input[i] * rho_input[i];
    }

    GinacEx J = J_target + J_stepx + J_stepy + J_effort + J_input;

    /* Constraints + Bounds */

    const int num_constraints = 4 * num_nodes;
    GinacMatrix constraints(num_constraints, 1);
    GinacMatrix lower_bound(num_constraints, 1);
    GinacMatrix upper_bound(num_constraints, 1);

    for (int i = 0; i < num_nodes; i++)
    {
        const int step_i = std::floor(i / nodes_per_step);
        if (step_i % 2)
        {
            // right step
            lower_bound(i, 0) = -INF;
            upper_bound(i, 0) = -sy_ref[step_i];
        }
        else
        {
            // left step
            lower_bound(i, 0) = sy_ref[step_i];
            upper_bound(i, 0) = INF;
        }
        constraints(i, 0) = sy[step_i] - rho_step[i];
    }

    GinacSymbol r_max("r_max");
    for (int i = 0; i < num_nodes; i++)
    {
        const int ci = i + num_nodes;
        lower_bound(ci, 0) = -r_max;
        upper_bound(ci, 0) = +r_max;
        constraints(ci, 0) = u[i] - q[i];
    }

    GinacSymbol u_off("u_off");
    for (int i = 0; i < num_nodes; i++)
    {
        const int ci = i + 2 * num_nodes;
        lower_bound(ci, 0) = -INF;
        upper_bound(ci, 0) = u_off;
        constraints(ci, 0) = u[i] - u_ref[i] - rho_input[i];
    }
    for (int i = 0; i < num_nodes; i++)
    {
        const int ci = i + 3 * num_nodes;
        lower_bound(ci, 0) = -INF;
        upper_bound(ci, 0) = u_off;
        constraints(ci, 0) = -u[i] + u_ref[i] - rho_input[i];
    }

    SymbolVector con_params = {p_foot_x, p_foot_y, r_max, u_off};

    /* Symbolic QP Problem */

    SymbolVector x = combineSymbolVectors({q, u, sx, sy, rho_step, rho_input});
    SymbolVector c = combineSymbolVectors({q_ref, sx_ref, sy_ref, w_target, w_stepx, w_stepy, w_effort, w_input, con_params});

    SymbolicQPProblem::Ptr sym_qp = std::make_shared<SymbolicQPProblem>(x, c);
    sym_qp->objective = J;
    sym_qp->constraints = constraints;
    sym_qp->lower_bound = lower_bound;
    sym_qp->upper_bound = upper_bound;

    /* Constants Eval */

    std::vector<Float> const_eval(c.size(), 0.0);
    int ci = 0;
    // q_ref
    for (int i = 0; i < num_nodes; i++)
    {
        const_eval[ci++] = 0.0; // x
        const_eval[ci++] = 0.0; // y
        const_eval[ci++] = 0.0; // vx
        const_eval[ci++] = 0.0; // vy
    }
    // sx_ref
    const Float stride_length = 0.0;
    for (int i = 0; i < num_steps; i++)
        const_eval[ci++] = stride_length * i;
    // sy_ref
    const Float stride_width = 0.1;
    for (int i = 0; i < num_steps; i++)
        const_eval[ci++] = i % 2 ? -stride_width : stride_width;
    // W_target
    for (int i = 0; i < num_states; i++)
        const_eval[ci++] = 1.0;
    // W_stepx
    for (int i = 0; i < num_steps; i++)
        const_eval[ci++] = 1.0;
    // W_stepy
    for (int i = 0; i < num_steps; i++)
        const_eval[ci++] = 1.0;
    // W_effort
    for (int i = 0; i < num_controls; i++)
        const_eval[ci++] = 1.0;
    // W_input
    for (int i = 0; i < num_nodes; i++)
        const_eval[ci++] = 1.0;
    // P_footx, P_footy, r_max, u_off
    const_eval[ci++] = 0.0;
    const_eval[ci++] = 0.0;
    const_eval[ci++] = 1.0;
    const_eval[ci++] = 0.0;

    sym_qp->evaluateConstants(const_eval);

    QPProblem::Ptr qp = sym_qp->getQP();

    std::cout << "Hessian:\n"
              << qp->hessian << std::endl;
    std::cout << "Gradient:\n"
              << qp->gradient << std::endl;
    std::cout << "Lin. Constraint:\n"
              << qp->linear_constraint << std::endl;
    std::cout << "Lower Bound:\n"
              << qp->lower_bound << std::endl;
    std::cout << "Upper Bound:\n"
              << qp->upper_bound << std::endl;

    OSQP::Ptr osqp = std::make_shared<OSQP>();
    osqp->getSettings()->verbose = true;
    osqp->getSettings()->warm_starting = true;
    osqp->getSettings()->polishing = true;

    osqp->setup(sym_qp->getQP());

    osqp->solve();

    const int runs = 1;
    const auto cstart = std::chrono::high_resolution_clock::now();
    for (int r = 0; r < runs; r++)
    {
        sym_qp->evaluateConstants(const_eval);
        osqp->update();
        osqp->solve();
    }
    const auto cend = std::chrono::high_resolution_clock::now();
    const long dur = std::chrono::duration_cast<std::chrono::microseconds>(cend - cstart).count();

    std::cout << "Run time: " << dur << "us\n";
    std::cout << "Frequency: " << 1.0 / (dur * 1E-6) << "Hz\n";

    return 0;
}