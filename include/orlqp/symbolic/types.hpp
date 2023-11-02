#ifndef ORLQP_SYMBOLIC_TYPES_HPP_
#define ORLQP_SYMBOLIC_TYPES_HPP_

#include "orlqp/types.hpp"

#include "ginac/ginac.h"

namespace orlqp
{

    using GinacEx = GiNaC::ex;
    using GinacMatrix = GiNaC::matrix;
    using GinacSymbol = GiNaC::symbol;
    using SymbolVector = std::vector<GinacSymbol>;

}

#endif