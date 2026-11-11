#pragma once

#include <cmath>
#include "UKF/Types.h"
#include "UKF/Covariance.h"

namespace robust_detail {
// Covariance reset for left/world Gibbs errors (MRP_A=0, MRP_F=2).
// For correction d, J = d.w^2 I + d.w [d.vec]_cross. This is the
// derivative of Gibbs(delta(r) * delta(r_hat)^-1) at r = r_hat.
// Return false, leaving root unchanged, unless both quaternions are finite
// unit quaternions, the correction stays inside the Gibbs chart, and the
// transported covariance is positive semidefinite.
inline bool recenter_attitude_root(
        const UKF::Quaternion& predicted, const UKF::Quaternion& corrected,
        const UKF::Matrix<3, 3>& posterior_root, UKF::Matrix<3, 3>& root) {
    if (!predicted.coeffs().allFinite() || !corrected.coeffs().allFinite() ||
        std::abs(predicted.norm() - 1) > 1e-6 || std::abs(corrected.norm() - 1) > 1e-6)
        return false;
    const UKF::Quaternion d = (corrected * predicted.conjugate()).normalized();
    if (std::abs(d.w()) < 1e-6)
        return false;
    UKF::Matrix<3, 3> cross;
    cross << 0, -d.z(), d.y(), d.z(), 0, -d.x(), -d.y(), d.x(), 0;
    const UKF::Matrix<3, 3> jacobian = d.w() * d.w() * UKF::Matrix<3, 3>::Identity() + d.w() * cross;
    const UKF::Matrix<3, 3> transformed = jacobian * posterior_root;
    const UKF::Matrix<3, 3> covariance = transformed * transformed.transpose();
    return UKF::Detail::covariance_root(covariance, root);
}
}
