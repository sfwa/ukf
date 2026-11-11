// MIT license; see the repository LICENSE.
#pragma once

#include <Eigen/Cholesky>
#include <cmath>
#include "UKF/Types.h"
#include "UKF/Covariance.h"

namespace mav_detail {
template <int N> struct Correction {
    UKF::Vector<N> delta;
    UKF::Matrix<N, N> posterior;
};

// Measurement-sized storage is bounded at compile time by MaxM rows, so the
// correction path keeps every Eigen temporary on the stack.
template <int N, int MaxM>
using CrossMatrix = Eigen::Matrix<real_t, N, Eigen::Dynamic, Eigen::ColMajor, N, MaxM>;
template <int MaxM>
using SquareMatrix = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, MaxM, MaxM>;
template <int MaxM>
using MeasurementVector = Eigen::Matrix<real_t, Eigen::Dynamic, 1, Eigen::ColMajor, MaxM, 1>;

// State error order is left/world attitude, body velocity, then optional body
// gyro bias. World yaw is a gauge: these IMU-only observations cannot fix it.
// A constrained gain requires the general covariance update, not P-K*S*K^T.
// Return false, leaving result unspecified, for incompatible or non-finite
// moments, an innovation covariance that is not positive definite, a bias
// constraint without a bias state, or a non-finite result.
template <int N, int MaxM>
bool constrained_update(
        const UKF::Matrix<N, N>& prior,
        const CrossMatrix<N, MaxM>& cross,
        const SquareMatrix<MaxM>& innovation_covariance,
        const MeasurementVector<MaxM>& innovation,
        bool freeze_bias,
        Correction<N>& result) {
    static_assert(N >= 3, "Attitude must occupy the first three error coordinates");
    const int m = innovation.size();
    if (m == 0 || cross.cols() != m || innovation_covariance.rows() != m ||
            innovation_covariance.cols() != m || !prior.allFinite() ||
            !cross.allFinite() || !innovation_covariance.allFinite() || !innovation.allFinite()) {
        return false;
    }
    const Eigen::LLT<SquareMatrix<MaxM>> factor(innovation_covariance);
    if (factor.info() != Eigen::Success)
        return false;
    CrossMatrix<N, MaxM> gain = factor.solve(cross.transpose()).transpose();
    gain.row(2).setZero();
    if (freeze_bias) {
        if constexpr (N >= 9) gain.template bottomRows<3>().setZero();
        else return false;
    }
    const UKF::Matrix<N, N> posterior = prior - gain * cross.transpose()
        - cross * gain.transpose() + gain * innovation_covariance * gain.transpose();
    result.delta = gain * innovation;
    result.posterior = ((posterior + posterior.transpose()) * real_t(0.5)).eval();
    return result.delta.allFinite() && result.posterior.allFinite();
}

// Recenter the complete joint covariance. Velocity and bias are additive in
// the fixed sensor frame, so only the attitude block of the chart Jacobian
// changes. Multiplying the complete covariance also transports cross blocks.
// Return false, leaving root unchanged, unless both quaternions are finite
// unit quaternions, the correction stays inside the Gibbs chart, the scale is
// finite and nonnegative, and the transported covariance is semidefinite.
template <int N>
bool recenter_root(const UKF::Quaternion& predicted,
        const UKF::Quaternion& corrected, const UKF::Matrix<N, N>& posterior,
        UKF::Matrix<N, N>& root, real_t reference_scale = 0) {
    static_assert(N >= 3, "Attitude must occupy the first three error coordinates");
    if (!predicted.coeffs().allFinite() || !corrected.coeffs().allFinite() ||
            std::abs(predicted.norm() - 1) > 1e-6 || std::abs(corrected.norm() - 1) > 1e-6 ||
            !std::isfinite(reference_scale) || reference_scale < 0)
        return false;
    const UKF::Quaternion correction = (corrected * predicted.conjugate()).normalized();
    if (std::abs(correction.w()) < 1e-6)
        return false;
    UKF::Matrix<3, 3> cross;
    cross << 0, -correction.z(), correction.y(),
             correction.z(), 0, -correction.x(),
             -correction.y(), correction.x(), 0;
    UKF::Matrix<N, N> jacobian = UKF::Matrix<N, N>::Identity();
    jacobian.template topLeftCorner<3, 3>() = correction.w() * correction.w()
        * UKF::Matrix<3, 3>::Identity() + correction.w() * cross;
    const UKF::Matrix<N, N> transported = jacobian
        * ((posterior + posterior.transpose()) * real_t(0.5)) * jacobian.transpose();
    return UKF::Detail::covariance_root(transported, root, reference_scale);
}
}
