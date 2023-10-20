#ifndef ORLQP_QP_PROBLEM_ARRAY_HPP_
#define ORLQP_QP_PROBLEM_ARRAY_HPP_

#include "orlqp/types.hpp"
#include "orlqp/QPProblem.hpp"

namespace orlqp
{

    class QPArrayProblem
    {

    public:
        using Ptr = std::shared_ptr<QPArrayProblem>;

        const int num_problems;
        const std::vector<QPProblem::Ptr> qp_problems;

        QPArrayProblem(const std::vector<QPProblem::Ptr> &qp_array);

        void update();

        std::vector<QPSolution::Ptr> splitQPSolution(const QPSolution::Ptr qp_array_solution);

        QPProblem::Ptr getQP();

    private:
        QPProblem::Ptr QP;

        std::vector<int> variable_idx_map;
        std::vector<int> constraint_idx_map;

        void calculateQPArrayHessian();
        void calculateQPArrayGradient();
        void calculateQPArrayLinearConstraint();
        void calculateQPArrayLowerBound();
        void calculateQPArrayUpperBound();

        void updateQPArrayHessian(const int index);
        void updateQPArrayGradient(const int index);
        void updateQPArrayLinearConstraint(const int index);
        void updateQPArrayLowerBound(const int index);
        void updateQPArrayUpperBound(const int index);
    };

}

#endif