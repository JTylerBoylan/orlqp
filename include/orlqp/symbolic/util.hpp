#ifndef ORLQP_SYMBOLIC_UTIL_HPP_
#define ORLQP_SYMBOLIC_UTIL_HPP_

#include <string>

#include "orlqp/symbolic/types.hpp"

namespace orlqp
{
    SymbolVector createSymbolVector(const int, const std::string &);

    GinacMatrix calculateExpressionGradient(GinacEx &, SymbolVector &);

    GinacMatrix calculateVectorJacobian(GinacMatrix &, SymbolVector &);

    GinacMatrix calculateExpressionHessian(GinacEx &, SymbolVector &);

    SymbolVector combineSymbolVectors(const std::vector<SymbolVector> &);
}

#endif