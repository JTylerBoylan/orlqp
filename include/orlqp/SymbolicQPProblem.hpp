#ifndef ORLQP_SYMBOLIC_QP_PROBLEM_HPP_
#define ORLQP_SYMBOLIC_QP_PROBLEM_HPP_

#include "orlqp/symbolic/types.hpp"
#include "orlqp/QPProblem.hpp"

namespace orlqp
{

    class SymbolicQPProblem
    {
        /*************************************

           Symbolic QP Problem

           minimize f(x,c)
           subject to lb <= h(x,c) <= ub

           x : Decision variables
           c : Constants
           f : Objective function
           h : Inequality constraints
           lb : Constraint lower bound
           ub : Constraint upper bound

       *************************************/

    public:
        using Ptr = std::shared_ptr<SymbolicQPProblem>;

        const SymbolVector x;
        const SymbolVector c;

        GinacEx objective;
        GinacMatrix constraints;
        GinacMatrix lower_bound;
        GinacMatrix upper_bound;

        struct
        {
            bool objective = false,
                 constraints = false,
                 lower_bound = false,
                 upper_bound = false;
        } update;

        SymbolicQPProblem(const SymbolVector &x, const SymbolVector &c);

        QPProblem::Ptr getQP();

        void updateQP();

        void evaluateConstants(const std::vector<Float> &);

    private:
        GinacMatrix sym_hessian;
        GinacMatrix sym_gradient;
        GinacMatrix sym_lin_constraints;
        GinacMatrix sym_lower_bound;
        GinacMatrix sym_upper_bound;
        void calculateSymbolicHessian();
        void calculateSymbolicGradient();
        void calculateSymbolicLinearConstraint();
        void calculateSymbolicLowerBound();
        void calculateSymbolicUpperBound();

        QPProblem::Ptr QP;
        void calculateQPHessian();
        void calculateQPGradient();
        void calculateQPLinearConstraint();
        void calculateQPLowerBound();
        void calculateQPUpperBound();

        std::vector<Float> ceval;
    };

}

#endif