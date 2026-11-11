// MIT license; see the repository LICENSE.
#pragma once

#include <Eigen/Cholesky>
#include <cmath>
#include <stdexcept>
#include "UKF/Types.h"
#include "UKF/Covariance.h"

namespace mav_detail {
template <int N> struct Correction {
    UKF::Vector<N> delta;
    UKF::Matrix<N, N> posterior;
};

// State error order is left/world attitude, body velocity, then optional body
// gyro bias. World yaw is a gauge: these IMU-only observations cannot fix it.
// A constrained gain requires the general covariance update, not P-K*S*K^T.
template <int N>
Correction<N> constrained_update(
        const UKF::Matrix<N, N>& prior,
        const UKF::Matrix<N, Eigen::Dynamic>& cross,
        const UKF::Matrix<Eigen::Dynamic, Eigen::Dynamic>& innovation_covariance,
        const UKF::Vector<Eigen::Dynamic>& innovation,
        bool freeze_bias = false) {
    static_assert(N >= 3, "Attitude must occupy the first three error coordinates");
    const int m = innovation.size();
    if (m == 0 || cross.cols() != m || innovation_covariance.rows() != m ||
            innovation_covariance.cols() != m || !prior.allFinite() ||
            !cross.allFinite() || !innovation_covariance.allFinite() || !innovation.allFinite()) {
        throw std::domain_error("MAV correction requires compatible finite moments");
    }
    const Eigen::LLT<UKF::Matrix<Eigen::Dynamic, Eigen::Dynamic>> factor(innovation_covariance);
    if (factor.info() != Eigen::Success)
        throw std::domain_error("MAV innovation covariance must be positive definite");
    UKF::Matrix<N, Eigen::Dynamic> gain = factor.solve(cross.transpose()).transpose();
    gain.row(2).setZero();
    if (freeze_bias) {
        if constexpr (N >= 9) gain.template bottomRows<3>().setZero();
        else throw std::domain_error("MAV bias constraint requires a bias state");
    }
    const UKF::Matrix<N, N> posterior = prior - gain * cross.transpose()
        - cross * gain.transpose() + gain * innovation_covariance * gain.transpose();
    Correction<N> result;
    result.delta = gain * innovation;
    result.posterior = ((posterior + posterior.transpose()) * real_t(0.5)).eval();
    if (!result.delta.allFinite() || !result.posterior.allFinite())
        throw std::domain_error("MAV correction produced non-finite values");
    return result;
}

// Recenter the complete joint covariance. Velocity and bias are additive in
// the fixed sensor frame, so only the attitude block of the chart Jacobian
// changes. Multiplying the complete covariance also transports cross blocks.
template <int N>
UKF::Matrix<N, N> recenter_root(const UKF::Quaternion& predicted,
        const UKF::Quaternion& corrected, const UKF::Matrix<N, N>& posterior,
        real_t reference_scale = 0) {
    static_assert(N >= 3, "Attitude must occupy the first three error coordinates");
    if (!predicted.coeffs().allFinite() || !corrected.coeffs().allFinite() ||
            std::abs(predicted.norm() - 1) > 1e-6 || std::abs(corrected.norm() - 1) > 1e-6 ||
            !std::isfinite(reference_scale) || reference_scale < 0)
        throw std::domain_error("MAV covariance reset requires unit quaternions and a finite scale");
    const UKF::Quaternion correction = (corrected * predicted.conjugate()).normalized();
    if (std::abs(correction.w()) < 1e-6)
        throw std::domain_error("MAV attitude correction reached a Gibbs-chart singularity");
    UKF::Matrix<3, 3> cross;
    cross << 0, -correction.z(), correction.y(),
             correction.z(), 0, -correction.x(),
             -correction.y(), correction.x(), 0;
    UKF::Matrix<N, N> jacobian = UKF::Matrix<N, N>::Identity();
    jacobian.template topLeftCorner<3, 3>() = correction.w() * correction.w()
        * UKF::Matrix<3, 3>::Identity() + correction.w() * cross;
    const UKF::Matrix<N, N> transported = jacobian
        * ((posterior + posterior.transpose()) * real_t(0.5)) * jacobian.transpose();
    return UKF::Detail::covariance_root(transported, reference_scale);
}
}
