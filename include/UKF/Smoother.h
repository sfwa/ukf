#ifndef UKF_SMOOTHER_H
#define UKF_SMOOTHER_H

#include <Eigen/LU>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

#include "UKF/Covariance.h"
#include "UKF/StateVector.h"

/*
Unscented Rauch-Tung-Striebel smoothing: unscented_rts_step for offline use
and FixedLagSmoother for a bounded delay. Experimental in 1.0: these
interfaces and a_priori_step_with_cross_covariance may change in a minor
release (see RELEASING.md).
*/

namespace UKF {

/* A state and its full covariance in the tangent coordinates at that state. */
template <typename State>
struct SmoothingEstimate {
  State state = State::Zero();
  typename State::CovarianceMatrix covariance = State::CovarianceMatrix::Zero();
};

namespace Detail {

template <typename State>
bool valid_smoothing_estimate(const SmoothingEstimate<State>& estimate) {
  if (!estimate.state.allFinite() ||
      !estimate.covariance.isApprox(estimate.covariance.transpose())) {
    return false;
  }
  typename State::CovarianceMatrix root;
  return covariance_root(estimate.covariance, root);
}

/*
Solve G P_pred = C in unit-variance coordinates. A deterministic direction
can have zero variance, including after an exact measurement, so use the
pseudoinverse of the correlation matrix when a Cholesky pivot is within
roundoff of zero. An exactly dependent row can acquire a small positive
pivot through rounding.
*/
template <typename MatrixType>
bool smoothing_gain(const MatrixType& predicted, const MatrixType& cross,
                    MatrixType& gain) {
  if (!cross.allFinite() || !predicted.allFinite() ||
      !predicted.isApprox(predicted.transpose())) {
    return false;
  }
  const auto scale = predicted.diagonal()
                         .unaryExpr([](real_t v) {
                           return v > 0 ? real_t(1) / std::sqrt(v) : real_t(1);
                         })
                         .eval();
  const MatrixType correlation =
      scale.asDiagonal() * predicted * scale.asDiagonal();
  const MatrixType right = scale.asDiagonal() * cross.transpose();
  if (!correlation.allFinite() || !right.allFinite()) return false;
  const real_t tolerance = real_t(16) * std::numeric_limits<real_t>::epsilon() *
                           predicted.rows() * correlation.cwiseAbs().maxCoeff();
  const Eigen::LLT<MatrixType> llt(correlation);
  if (llt.info() == Eigen::Success &&
      llt.matrixLLT().diagonal().minCoeff() > std::sqrt(tolerance)) {
    gain = llt.solve(right).transpose() * scale.asDiagonal();
  } else {
    const Eigen::SelfAdjointEigenSolver<MatrixType> eig(correlation);
    if (eig.info() != Eigen::Success) return false;
    if (eig.eigenvalues().minCoeff() < -tolerance) return false;
    const auto reciprocal =
        eig.eigenvalues()
            .unaryExpr([&](real_t value) {
              return value > tolerance ? real_t(1) / value : real_t(0);
            })
            .eval();
    const auto& basis = eig.eigenvectors();
    gain = right.transpose() * basis * reciprocal.asDiagonal() *
           basis.transpose() * scale.asDiagonal();
  }
  return gain.allFinite();
}

/* Transport full covariance, including cross terms with Euclidean fields. */
template <typename State, bool Inverse>
void smoothing_transport(typename State::CovarianceMatrix& covariance,
                         const typename State::StateVectorDelta& delta) {
  if constexpr (AttitudeReset<State>::has_quaternion) {
    AttitudeReset<State>::for_each(delta, [&](std::size_t offset,
                                              const Matrix<3, 3>& reset) {
      Matrix<3, 3> transform = reset;
      if constexpr (Inverse) transform = reset.inverse();
      covariance.template middleRows<3>(offset) =
          (transform * covariance.template middleRows<3>(offset)).eval();
      covariance.template middleCols<3>(offset) =
          (covariance.template middleCols<3>(offset) * transform.transpose())
              .eval();
    });
  }
}

/* One RTS recursion with a gain already computed from the forward pass. */
template <typename State>
bool smooth_with_gain(const SmoothingEstimate<State>& filtered,
                      const SmoothingEstimate<State>& predicted,
                      const typename State::CovarianceMatrix& gain,
                      const SmoothingEstimate<State>& next,
                      SmoothingEstimate<State>& output) {
  const auto next_delta = next.state.delta_from(predicted.state);
  const typename State::StateVectorDelta delta = gain * next_delta;
  SmoothingEstimate<State> result = filtered;
  result.state.apply_delta(delta);
  auto future_cov = next.covariance;
  smoothing_transport<State, true>(future_cov, next_delta);
  result.covariance.noalias() +=
      gain * (future_cov - predicted.covariance) * gain.transpose();
  smoothing_transport<State, false>(result.covariance, delta);
  result.covariance =
      ((result.covariance + result.covariance.transpose()) * real_t(0.5))
          .eval();
  typename State::CovarianceMatrix root;
  if (!result.state.allFinite() ||
      !covariance_root(result.covariance, root,
                       filtered.covariance.cwiseAbs().maxCoeff(),
                       &result.covariance)) {
    return false;
  }
  output = result;
  return true;
}
}  // namespace Detail

/*
One offline unscented Rauch-Tung-Striebel backward step. filtered is k|k,
predicted is k+1|k, cross is Cov(x_k,x_{k+1}|z_{1:k}), and next is k+1|T.
On success output holds k|T. On failure it is unchanged.
*/
template <typename State>
bool unscented_rts_step(const SmoothingEstimate<State>& filtered,
                        const SmoothingEstimate<State>& predicted,
                        const typename State::CovarianceMatrix& cross,
                        const SmoothingEstimate<State>& next,
                        SmoothingEstimate<State>& output) {
  typename State::CovarianceMatrix gain;
  return Detail::valid_smoothing_estimate(filtered) &&
         Detail::valid_smoothing_estimate(next) &&
         predicted.state.allFinite() &&
         Detail::smoothing_gain(predicted.covariance, cross, gain) &&
         Detail::smooth_with_gain(filtered, predicted, gain, next, output);
}

/*
A bounded ring of Lag transitions and the latest filtered estimate. Call
reset once, push one completed forward update per transition, and query any
retained time by its age. A zero-lag window stores just the current estimate.

Storage is fixed: Lag records of two states and three n-by-n matrices (the
filtered and predicted covariances and the gain), plus the latest estimate,
for a state uncertainty of dimension n. Each push solves one n-by-n gain;
estimate(k) applies k backward steps, each a few n-by-n matrix products.
*/
template <typename State, std::size_t Lag>
class FixedLagSmoother {
 public:
  using Estimate = SmoothingEstimate<State>;
  using Covariance = typename State::CovarianceMatrix;

  bool reset(const Estimate& filtered) {
    if (!Detail::valid_smoothing_estimate(filtered)) return false;
    latest_ = filtered;
    newest_ = 0;
    count_ = 1;
    return true;
  }

  bool push(const Estimate& predicted, const Covariance& cross,
            const Estimate& filtered) {
    Covariance gain;
    if (count_ == 0 || !Detail::valid_smoothing_estimate(filtered) ||
        !predicted.state.allFinite() ||
        !Detail::smoothing_gain(predicted.covariance, cross, gain)) {
      return false;
    }
    if constexpr (Lag > 0) {
      history_[newest_] = {latest_, predicted, gain};
      newest_ = (newest_ + 1) % Lag;
      if (count_ <= Lag) ++count_;
    }
    latest_ = filtered;
    return true;
  }

  std::size_t size() const { return count_; }

  bool estimate(std::size_t steps_back, Estimate& output) const {
    if (steps_back >= count_) return false;
    Estimate result = latest_;
    if constexpr (Lag > 0) {
      auto index = newest_;
      for (std::size_t i = 0; i < steps_back; ++i) {
        index = index == 0 ? Lag - 1 : index - 1;
        const auto& record = history_[index];
        if (!Detail::smooth_with_gain(record.filtered, record.predicted,
                                      record.gain, result, result))
          return false;
      }
    }
    output = result;
    return true;
  }

 private:
  struct Record {
    Estimate filtered;
    Estimate predicted;
    Covariance gain = Covariance::Zero();
  };
  std::array<Record, Lag> history_;
  Estimate latest_;
  std::size_t newest_ = 0;
  std::size_t count_ = 0;
};
}  // namespace UKF
#endif
