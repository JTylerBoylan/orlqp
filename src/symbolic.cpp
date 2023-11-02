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

    GinacMatrix calculateExpressionGradient(const GinacEx &f, const SymbolVector &x)
    {
        const int N = x.size();
        GinacMatrix G(N, 1);
        for (int i = 0; i < N; i++)
        {
            G(i, 0) = f.diff(x[i]);
        }
        return G;
    }

    GinacMatrix calculateVectorJacobian(const GinacMatrix &G, const SymbolVector &x)
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

    GinacMatrix calculateExpressionHessian(const GinacEx &f, const SymbolVector &x)
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

    GinacEx evaluateExpression(const GinacEx &ex, const SymbolVector &vars, const std::vector<Float> &cnsts)
    {
        GinacEx ex_out;
        for (int i = 0; i < vars.size(); i++)
            ex_out = ex.subs({{vars[i], cnsts[i]}});
        return ex_out;
    }
    GinacMatrix evaluateMatrix(const GinacMatrix &matrix, const SymbolVector &vars, const std::vector<Float> &cnsts)
    {
        GinacMatrix matrix_out(matrix.rows(), matrix.cols());
        for (int r = 0; r < matrix.rows(); r++)
            for (int c = 0; c < matrix.cols(); c++)
            {
                matrix_out(r,c) = matrix(r,c);
                for (int i = 0; i < vars.size(); i++)
                {
                    matrix_out(r, c) = matrix_out(r,c).subs({{vars[i], cnsts[i]}});
                }
            }
        return matrix_out;
    }

}

/* SYMBOLIC QP PROBLEM */
#include "orlqp/SymbolicQPProblem.hpp"
namespace orlqp
{

    SymbolicQPProblem::SymbolicQPProblem(const SymbolVector &vars, const SymbolVector &cons)
        : x(vars), c(cons)
    {
    }

    QPProblem::Ptr SymbolicQPProblem::getQP()
    {
        if (!this->QP)
        {
            calculateSymbolicHessian();
            calculateSymbolicGradient();
            calculateSymbolicLinearConstraint();
            calculateSymbolicLowerBound();
            calculateSymbolicUpperBound();

            this->QP = std::make_shared<QPProblem>(this->x.size(), this->constraints.rows());
            calculateQPHessian();
            calculateQPGradient();
            calculateQPLinearConstraint();
            calculateQPLowerBound();
            calculateQPUpperBound();
        }
        return this->QP;
    }

    void SymbolicQPProblem::updateQP()
    {
        if (this->update.objective)
        {
            calculateSymbolicHessian();
            calculateSymbolicGradient();
            calculateQPHessian();
            calculateQPGradient();
            this->update.objective = false;
        }
        if (this->update.constraints)
        {
            calculateSymbolicLinearConstraint();
            calculateQPLinearConstraint();
            this->update.constraints = false;
        }
        if (this->update.lower_bound)
        {
            calculateSymbolicLowerBound();
            calculateQPLowerBound();
            this->update.lower_bound = false;
        }
        if (this->update.upper_bound)
        {
            calculateSymbolicUpperBound();
            calculateQPUpperBound();
            this->update.upper_bound = false;
        }
    }

    void SymbolicQPProblem::evaluateConstants(const std::vector<Float> &constants)
    {
        this->ceval = constants;
        calculateQPHessian();
        calculateQPGradient();
        calculateQPLinearConstraint();
        calculateQPLowerBound();
        calculateQPUpperBound();
    }

    void SymbolicQPProblem::calculateSymbolicHessian()
    {
        this->sym_hessian = calculateExpressionHessian(this->objective, this->x);
    }
    void SymbolicQPProblem::calculateSymbolicGradient()
    {
        auto gradient = calculateExpressionGradient(this->objective, this->x);
        gradient = evaluateMatrix(gradient, this->x, std::vector<Float>(this->x.size(), 0.0));
        this->sym_gradient = gradient;
    }
    void SymbolicQPProblem::calculateSymbolicLinearConstraint()
    {
        this->sym_lin_constraints = calculateVectorJacobian(this->constraints, this->x);
    }
    void SymbolicQPProblem::calculateSymbolicLowerBound()
    {
        this->sym_lower_bound = this->lower_bound;
    }
    void SymbolicQPProblem::calculateSymbolicUpperBound()
    {
        this->sym_upper_bound = this->upper_bound;
    }

    void SymbolicQPProblem::calculateQPHessian()
    {
        const GinacMatrix hessian_eval = evaluateMatrix(this->sym_hessian, this->c, this->ceval);
        std::vector<EigenTriplet> triplets;
        for (int r = 0; r < hessian_eval.rows(); r++)
            for (int c = 0; c < hessian_eval.cols(); c++)
            {
                const Float val = GiNaC::ex_to<GiNaC::numeric>(hessian_eval(r, c)).to_double();
                if (val != 0)
                {
                    triplets.push_back(EigenTriplet(r, c, val));
                }
            }
        this->QP->hessian.setFromTriplets(triplets.begin(), triplets.end());
    }
    void SymbolicQPProblem::calculateQPGradient()
    {
        const GinacMatrix gradient_eval = evaluateMatrix(this->sym_gradient, this->c, this->ceval);
        for (int r = 0; r < gradient_eval.rows(); r++)
        {
            const Float val = GiNaC::ex_to<GiNaC::numeric>(gradient_eval(r, 0)).to_double();
            this->QP->gradient(r, 0) = val;
        }
    }
    void SymbolicQPProblem::calculateQPLinearConstraint()
    {
        const GinacMatrix lin_constraints_eval = evaluateMatrix(this->sym_lin_constraints, this->c, this->ceval);
        std::vector<EigenTriplet> triplets;
        for (int r = 0; r < lin_constraints_eval.rows(); r++)
            for (int c = 0; c < lin_constraints_eval.cols(); c++)
            {
                const Float val = GiNaC::ex_to<GiNaC::numeric>(lin_constraints_eval(r, c)).to_double();
                if (val != 0)
                {
                    triplets.push_back(EigenTriplet(r, c, val));
                }
            }
        this->QP->linear_constraint.setFromTriplets(triplets.begin(), triplets.end());
    }
    void SymbolicQPProblem::calculateQPLowerBound()
    {
        const GinacMatrix lower_bound_eval = evaluateMatrix(this->sym_lower_bound, this->c, this->ceval);
        for (int r = 0; r < lower_bound_eval.rows(); r++)
        {
            const Float val = GiNaC::ex_to<GiNaC::numeric>(lower_bound_eval(r, 0)).to_double();
            this->QP->lower_bound(r, 0) = val;
        }
    }
    void SymbolicQPProblem::calculateQPUpperBound()
    {
        const GinacMatrix upper_bound_eval = evaluateMatrix(this->sym_upper_bound, this->c, this->ceval);
        for (int r = 0; r < upper_bound_eval.rows(); r++)
        {
            const Float val = GiNaC::ex_to<GiNaC::numeric>(upper_bound_eval(r, 0)).to_double();
            this->QP->upper_bound(r, 0) = val;
        }
    }

}