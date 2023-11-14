#ifndef ORLQP_EXPORT_QP_TIANZELIP_HPP_
#define ORLQP_EXPORT_QP_TIANZELIP_HPP_

#include "orlqp/orlqp.hpp"

namespace orlqp
{
    class TianzeLIP
    {
    public:
        const int num_variables = 100;
        const int num_constraints = 100;
        using Ptr = std::shared_ptr<TianzeLIP>;

        TianzeLIP() {}

        QPProblem::Ptr getQP(const EigenVector &Qr, const EigenVector &Sxr, const EigenVector &Syr)
        {
            QP = std::make_shared<QPProblem>(num_variables, num_constraints);
            QP->hessian.setFromTriplets(hessian_triplets.begin(), hessian_triplets.end());
            QP->gradient << 0.0, 0.0;
            QP->linear_constraint.setFromTriplets(lin_constraint_triplets.begin(), lin_constraint_triplets.end());
            QP->upper_bound << 0.0, 0.0;
            QP->lower_bound << 0.0, 0.0;
            return QP;
        }

    private:
        QPProblem::Ptr QP;
        std::vector<EigenTriplet> hessian_triplets;
        std::vector<EigenTriplet> lin_constraint_triplets;
    };
}

#endif