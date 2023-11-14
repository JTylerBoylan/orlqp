#ifndef ORLQP_SYMBOLIC_EXPORT_HPP_
#define ORLQP_SYMBOLIC_EXPORT_HPP_

#include "orlqp/symbolic/util.hpp"
#include "orlqp/SymbolicQPProblem.hpp"

namespace orlqp
{

    class ExportSymbolicQPProblem
    {
    public:
        ExportSymbolicQPProblem(const SymbolicQPProblem::Ptr sym_qp,
                                const std::string &name,
                                const std::string &path);

    private:
        const std::string ClassName;
        const SymbolicQPProblem::Ptr QP;

        std::string generateHeader();
        std::string generateHessian();
        std::string generateGradient();
        std::string generateLinearConstraint();
        std::string generateLowerBound();
        std::string generateUpperBound();
        std::string generateEnder();
    };

}

#endif