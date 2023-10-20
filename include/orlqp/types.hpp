#ifndef ORLQP_TYPES_HPP_
#define ORLQP_TYPES_HPP_

#include <memory>

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include <osqp/osqp.h>

namespace orlqp
{

#ifdef OSQP_USE_FLOAT
    using Float = float;
    using EigenVector = Eigen::VectorXf;
    using EigenMatrix = Eigen::MatrixXf;
#else
    using Float = double;
    using EigenVector = Eigen::VectorXd;
    using EigenMatrix = Eigen::MatrixXd;
#endif
    using EigenSparseMatrix = Eigen::SparseMatrix<Float>;
    using EigenTriplet = Eigen::Triplet<Float>;

}

#endif