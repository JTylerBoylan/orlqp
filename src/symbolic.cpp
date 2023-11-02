/* SYMBOLIC UTILITY */
#include "orlqp/symbolic/util.hpp"
namespace orlqp
{

    SymbolVector createSymbolVector(const int size, const std::string &var)
    {
        SymbolVector vec(size);
        for (int i = 0; i < size; i++)
        {
            const std::string var_i = var + "[" + std::to_string(i) + "]";
            vec[i] = GinacSymbol(var_i);
        }
        return vec;
    }

    GinacMatrix calculateExpressionGradient(GinacEx &f, SymbolVector &x)
    {
        const int N = x.size();
        GinacMatrix G(N, 1);
        for (int i = 0; i < N; i++)
        {
            G(i, 0) = f.diff(x[i]);
        }
        return G;
    }

    GinacMatrix calculateVectorJacobian(GinacMatrix &G, SymbolVector &x)
    {
        const int N1 = G.rows();
        const int N2 = x.size();
        GinacMatrix J(N1, N2);
        for (int i = 0; i < N1; i++)
        {
            const GinacEx &expr = G[i];
            for (int j = 0; j < N2; j++)
            {
                const GinacSymbol sym = x[j];
                J(i, j) = expr.diff(sym);
            }
        }
        return J;
    }

    GinacMatrix calculateExpressionHessian(GinacEx &f, SymbolVector &x)
    {
        auto G = calculateExpressionGradient(f, x);
        auto H = calculateVectorJacobian(G, x);
        return H;
    }

    SymbolVector combineSymbolVectors(const std::vector<SymbolVector> &vecs)
    {
        size_t totalSize = 0;
        for (const auto &vec : vecs)
            totalSize += vec.size();
        SymbolVector result;
        result.reserve(totalSize);
        for (const auto &vec : vecs)
            result.insert(result.cend(), vec.cbegin(), vec.cend());
        return result;
    }

}

/* SYMBOLIC QP PROBLEM */
#include "orlqp/SymbolicQPProblem.hpp"
namespace orlqp
{

    SymbolicQPProblem::SymbolicQPProblem()
    {
        /* TODO */
    }

    QPProblem::Ptr SymbolicQPProblem::getQP()
    {
        /* TODO */
    }

    void SymbolicQPProblem::calculateSymbolicHessian()
    {
        /* TODO */
    }
    void SymbolicQPProblem::calculateSymbolicGradient()
    {
        /* TODO */
    }
    void SymbolicQPProblem::calculateSymbolicLinearConstraint()
    {
        /* TODO */
    }
    void SymbolicQPProblem::calculateSymbolicLowerBound()
    {
        /* TODO */
    }
    void SymbolicQPProblem::calculateSymbolicUpperBound()
    {
        /* TODO */
    }

    void SymbolicQPProblem::calculateQPHessian()
    {
        /* TODO */
    }
    void SymbolicQPProblem::calculateQPGradient()
    {
        /* TODO */
    }
    void SymbolicQPProblem::calculateQPLinearConstraint()
    {
        /* TODO */
    }
    void SymbolicQPProblem::calculateQPLowerBound()
    {
        /* TODO */
    }
    void SymbolicQPProblem::calculateQPUpperBound()
    {
        /* TODO */
    }

}