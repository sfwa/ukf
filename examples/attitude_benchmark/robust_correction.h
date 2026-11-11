#pragma once

#include <Eigen/Cholesky>
#include <stdexcept>
#include "UKF/Types.h"

namespace robust_detail {
struct GravityUpdate {
    UKF::Vector<3> delta;
    UKF::Matrix<3, 3> posteriorP;
};

// Gravity cannot observe rotation about world +Z. Constrain that row of the
// gain in LEFT/world attitude-error coordinates. The usual optimal-gain
// subtraction P-KSK^T no longer applies after imposing this constraint.
inline GravityUpdate constrained_gravity_update(
        const UKF::Matrix<3, 3>& priorP,
        const UKF::Matrix<3, 2>& cross_covariance,
        const UKF::Matrix<2, 2>& innovation_covariance,
        const UKF::Vector<2>& innovation) {
    if (!priorP.allFinite() || !cross_covariance.allFinite() ||
        !innovation_covariance.allFinite() || !innovation.allFinite()) {
        throw std::domain_error("Gravity correction requires finite moments");
    }
    const Eigen::LLT<UKF::Matrix<2, 2>> factor(innovation_covariance);
    if (factor.info() != Eigen::Success) {
        throw std::domain_error("Gravity innovation covariance must be positive definite");
    }
    UKF::Matrix<3, 2> gain = factor.solve(cross_covariance.transpose()).transpose();
    gain.row(2).setZero();
    const UKF::Matrix<3, 3> posterior = priorP
        - gain * cross_covariance.transpose()
        - cross_covariance * gain.transpose()
        + gain * innovation_covariance * gain.transpose();
    GravityUpdate result;
    result.delta = gain * innovation;
    result.posteriorP = ((posterior + posterior.transpose()) * real_t(0.5)).eval();
    if (!result.delta.allFinite() || !result.posteriorP.allFinite()) {
        throw std::domain_error("Gravity correction produced non-finite values");
    }
    return result;
}
}
