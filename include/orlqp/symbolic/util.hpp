#ifndef ORLQP_SYMBOLIC_UTIL_HPP_
#define ORLQP_SYMBOLIC_UTIL_HPP_

#include <string>

#include "orlqp/symbolic/types.hpp"

namespace orlqp
{
    SymbolVector createSymbolVector(const int, const std::string &);

    GinacMatrix calculateExpressionGradient(const GinacEx &, const SymbolVector &);

    GinacMatrix calculateVectorJacobian(const GinacMatrix &, const SymbolVector &);

    GinacMatrix calculateExpressionHessian(const GinacEx &, const SymbolVector &);

    SymbolVector combineSymbolVectors(const std::vector<SymbolVector> &);

    GinacEx evaluateExpression(const GinacEx &, const SymbolVector &, const std::vector<Float> &);

    GinacMatrix evaluateMatrix(const GinacMatrix &, const SymbolVector &, const std::vector<Float> &);

    void printSymbolicMatrix(const GinacMatrix &);
}

#endif