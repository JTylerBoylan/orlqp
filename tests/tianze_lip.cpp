#include <iostream>

#include "orlqp/symbolic/util.hpp"
#include "orlqp/SymbolicQPProblem.hpp"

using namespace orlqp;

int main(void)
{

    // Params
    const int nodes_per_step = 4;
    const int num_steps = 3;
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
    auto sx_ref = createSymbolVector(num_steps, "Sxr");     // const
    auto w_stepx = createSymbolVector(num_states, "Wstpx"); // const
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
    auto p_foot = createSymbolVector(2, "Pft"); // const
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
        u_ref[2 * i - 1] = p_foot[0] + sum_sx;
        u_ref[2 * i] = p_foot[1] + sum_sy;
    }
    auto w_effort = createSymbolVector(num_controls, "Weff"); // const
    GinacEx J_effort;
    for (int i = 0; i < num_controls; i++)
    {
        J_effort += w_effort[i] * (u[i] - u_ref[i]) * (u[i] - u_ref[i]);
    }

    auto rho_avoid = createSymbolVector(num_nodes, "rho_a");
    auto w_avoid = createSymbolVector(num_nodes, "Wavd"); // const
    GinacEx J_avoid;
    for (int i = 0; i < num_nodes; i++)
    {
        J_avoid += w_avoid[i] * rho_avoid[i] * rho_avoid[i];
    }

    auto rho_input = createSymbolVector(num_nodes, "rho_i");
    auto w_input = createSymbolVector(num_nodes, "Winp"); // const
    GinacEx J_input;
    for (int i = 0; i < num_nodes; i++)
    {
        J_input += w_input[i] * rho_input[i] * rho_input[i];
    }

    GinacEx J = J_target + J_stepx + J_stepy + J_effort + J_avoid + J_input;

    /* Constraints + Bounds */
    /* TODO */

    SymbolVector x = combineSymbolVectors({q, u, sx, sy, rho_step, rho_avoid, rho_input});
    SymbolVector c = combineSymbolVectors({q_ref, sx_ref, sy_ref, p_foot, w_target, w_stepx, w_stepy, w_effort, w_avoid, w_input});

    /* Symbolic QP Problem */
    /* TODO */

    GinacMatrix hessian = calculateExpressionHessian(J, x);
    GinacMatrix gradient = calculateExpressionGradient(J, x);

    return 0;
}