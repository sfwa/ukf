/*
Copyright (C) 2016 Thiemar Pty Ltd

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef CORE_H
#define CORE_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

#include "UKF/Covariance.h"
#include "UKF/Integrator.h"
#include "UKF/RankUpdate.h"
#include "UKF/StateVector.h"
#include "UKF/Types.h"

namespace UKF {

/* A full gated step retains its prediction on Rejected or Empty and restores
   its pre-step values on Failed. A split posterior only changes the state on
   Applied. Failed also covers a split posterior without a valid innovation.
   A result is not to be ignored: the compiler warns when one is discarded. */
enum class [[nodiscard]] GateResult { Applied, Rejected, Empty, Failed };

namespace Detail {

/*
Propagate every sigma point through the process model: either the state
vector's direct process_model(delta, input...) specialisation, when the core
uses NotUsedIntegrator, or the integrated derivative.
*/
template <typename StateVectorType, typename IntegratorType, typename... U>
void propagate_sigma_points(
    typename StateVectorType::SigmaPointDistribution& sigma_points,
    real_t delta, const U&... input) {
  for (std::size_t i = 0; i < StateVectorType::num_sigma(); i++) {
    if constexpr (std::is_same_v<IntegratorType, NotUsedIntegrator>) {
      sigma_points.col(i) << StateVectorType(sigma_points.col(i))
                                 .process_model(delta, input...);
    } else {
      sigma_points.col(i) << StateVectorType(sigma_points.col(i))
                                 .template process_model<IntegratorType>(
                                     delta, input...);
    }
  }
}

/*
Calculate the cross-correlation matrix described in equations 70 and 71 from
the Kraft paper, as one weighted product over the outer sigma points plus the
central point's rank-one term.
*/
template <typename StateVectorType, typename ResultType, typename StateDeltas,
          typename MeasurementDeltas>
void calculate_cross_correlation(ResultType& result, const StateDeltas& w_prime,
                                 const MeasurementDeltas& z_prime) {
  constexpr int outer = StateVectorType::num_sigma() - 1;
  result.noalias() = (Parameters::Sigma_WCI<StateVectorType> *
                      w_prime.template rightCols<outer>()) *
                     z_prime.template rightCols<outer>().transpose();
  result.noalias() +=
      (Parameters::Sigma_WC0<StateVectorType> * w_prime.col(0)) *
      z_prime.col(0).transpose();
}

/*
From the lower triangular factor L of A, compute a factor of
A + weight * sum(v_i v_i^T) over the columns v_i of vectors. Rank-one
updates or downdates are tried first; if a downdate would lose positive
definiteness, the factor is rebuilt from the updated covariance, which
accepts a semidefinite result within roundoff. Return false, leaving the
input factor unchanged, if the result is indefinite or not finite.
*/
template <typename MatrixType, typename VectorsType>
bool calculate_updated_root(const MatrixType& factor, MatrixType& updated,
                            const VectorsType& vectors,
                            typename MatrixType::Scalar weight) {
  updated = factor;
  if (rank_update<Eigen::Lower>(updated, vectors, weight) >= 0) {
    const MatrixType prior = factor * factor.transpose();
    MatrixType covariance = prior;
    covariance.noalias() += (weight * vectors) * vectors.transpose();
    return covariance_root(covariance, updated, prior.cwiseAbs().maxCoeff());
  }
  return updated.allFinite();
}

template <typename MatrixType, typename VectorsType>
bool update_root(MatrixType& factor, const VectorsType& vectors,
                 typename MatrixType::Scalar weight) {
  MatrixType updated;
  if (!calculate_updated_root(factor, updated, vectors, weight)) {
    return false;
  }
  factor = updated;
  return true;
}

/*
Calculate the lower triangular root L of the covariance of a set of sigma
point deltas d_i about their mean, plus additive noise with root covariance N:

    L L^T = WC0 d_0 d_0^T + WCI sum(d_i d_i^T) + N N^T

This is the QR decomposition and rank-one Cholesky update of the square-root
UKF paper. The central delta d_0 joins the QR decomposition unless its weight
is negative, as it is for small alpha; it is then a downdate, which fails if
the covariance is indefinite. Return false, leaving the root unchanged, in
that case or if the result is not finite.
*/
template <typename StateVectorType, typename RootType, typename DeltasType,
          typename NoiseType>
bool sigma_point_root(RootType& root, const DeltasType& deltas,
                      const NoiseType& noise_root) {
  constexpr real_t central = Parameters::Sigma_WC0<StateVectorType>;
  constexpr real_t outer = Parameters::Sigma_WCI<StateVectorType>;

  /* One row for each weighted delta, followed by one for each noise vector. */
  Eigen::Matrix<real_t, Eigen::Dynamic, DeltasType::RowsAtCompileTime,
                Eigen::ColMajor,
                int(DeltasType::MaxColsAtCompileTime) +
                    int(NoiseType::MaxColsAtCompileTime),
                DeltasType::MaxRowsAtCompileTime>
      vectors(deltas.cols() + noise_root.cols(), deltas.rows());
  vectors.topRows(deltas.cols()) = std::sqrt(outer) * deltas.transpose();
  vectors.row(0) *= std::sqrt(std::max(central, real_t(0)) / outer);
  vectors.bottomRows(noise_root.cols()) = noise_root.transpose();

  RootType factor = lower_factor(vectors);
  if (!factor.allFinite() ||
      (central < 0 && !update_root(factor, deltas.col(0), central))) {
    return false;
  }
  root = factor;
  return true;
}

/*
Whiten the cross-correlation P and the innovation v by the lower triangular
innovation root S: replace P with U = P S^-T and v with w = S^-1 v, by
forward substitution. The Kalman gain K = P (S S^T)^-1 itself is not needed,
because the state changes by K v = U w and its covariance by
-K S S^T K^T = -U U^T.

An innovation which cannot vary in some direction has a pivot that is
negligible beside the rest of its row: a field vector's rotation, for one,
has no component along the predicted field. That column of S is rotated into
the later columns, which leaves S S^T unchanged, and its direction, in which
there is nothing to learn, is left out of the correction.
*/
template <typename CrossType, typename InnovationType, typename RootType>
void whiten(CrossType& cross, InnovationType& innovation,
            const RootType& innovation_root) {
  using Scalar = typename RootType::Scalar;
  using Index = typename RootType::Index;
  const Scalar tolerance = std::sqrt(std::numeric_limits<Scalar>::epsilon());
  RootType root = innovation_root;
  for (Index k = 0; k < root.cols(); ++k) {
    const Index rest = root.cols() - k - 1;
    if (std::abs(root(k, k)) > tolerance * root.row(k).head(k).norm()) {
      cross.col(k) /= root(k, k);
      innovation(k) /= root(k, k);
      cross.rightCols(rest).noalias() -=
          cross.col(k) * root.col(k).tail(rest).transpose();
      innovation.tail(rest) -= innovation(k) * root.col(k).tail(rest);
    } else {
      /* The last column has nothing after it. A root of one column never
         has a later corner, and leaving it out at compile time keeps GCC's
         array bounds analysis quiet at -O2 and above. */
      if constexpr (RootType::MaxColsAtCompileTime != 1) {
        if (rest > 0) {
          auto later = root.bottomRightCorner(rest, rest);
          rank_update_lower(later, root.col(k).tail(rest), Scalar(1));
          root.col(k).tail(rest).setZero();
        }
      }
      cross.col(k).setZero();
      innovation(k) = Scalar(0);
    }
  }
}

/*
An innovation covariance S, checked and factored for the Kalman gain. S is
scaled to unit diagonal first, a zero or negative diagonal being left as it
is, so that the units of each sensor do not decide its numerical rank. The
scaled matrix must be symmetric, and is factored as P^T S P = L D L^T by a
symmetrically pivoted LDL^T decomposition, which eliminates the largest
remaining diagonal at each stage and stops at the first that is negligible.
What remains must then be negligible too, because no coefficient of a
semidefinite matrix exceeds its diagonal: a negative variance, or a
correlation that no covariance can have, leaves a remainder that is not.
Both tests take roundoff against the largest scaled coefficient.

The directions left at that point are those in which the innovation cannot
vary, such as a field vector's rotation about the predicted field, or which
are negative only by roundoff. The gain leaves out exactly these, so the
check and the solve agree on what is negligible. Plain loops keep GCC's
array bounds analysis quiet for a single row.
*/
template <typename MatrixType>
struct SemidefiniteFactor {
  using Scalar = typename MatrixType::Scalar;
  using Index = typename MatrixType::Index;
  using Vector = typename MatrixType::ColXpr::PlainObject;
  using Order =
      Eigen::Matrix<Index, MatrixType::RowsAtCompileTime, 1, Eigen::ColMajor,
                    MatrixType::MaxRowsAtCompileTime, 1>;

  /* L D below and on the diagonal of the scaled, permuted matrix. */
  MatrixType factor;
  Vector scale;
  /* The row of S eliminated at each stage. */
  Order order;
  Index rank = 0;

  /* Returns false if S is not finite, symmetric and semidefinite. */
  bool compute(const MatrixType& covariance) {
    const Index n = covariance.rows();
    if (!covariance.allFinite()) {
      return false;
    }
    scale.resize(n);
    order.resize(n);
    for (Index i = 0; i < n; ++i) {
      scale(i) = covariance(i, i) > Scalar(0)
                     ? Scalar(1) / std::sqrt(covariance(i, i))
                     : Scalar(1);
      order(i) = i;
    }
    factor = scale.asDiagonal() * covariance * scale.asDiagonal();
    Scalar largest = Scalar(0);
    for (Index j = 0; j < n; ++j) {
      for (Index i = 0; i < n; ++i) {
        largest = std::max(largest, std::abs(factor(i, j)));
      }
    }
    const Scalar tolerance = Scalar(16) *
                             std::numeric_limits<Scalar>::epsilon() *
                             Scalar(n) * largest;
    for (Index j = 0; j < n; ++j) {
      for (Index i = j + 1; i < n; ++i) {
        if (std::abs(factor(i, j) - factor(j, i)) > tolerance) {
          return false;
        }
      }
    }

    rank = 0;
    for (Index k = 0; k < n; ++k) {
      Index pivot = k;
      for (Index i = k + 1; i < n; ++i) {
        if (factor(i, i) > factor(pivot, pivot)) {
          pivot = i;
        }
      }
      if (factor(pivot, pivot) <= tolerance) {
        for (Index j = k; j < n; ++j) {
          for (Index i = j; i < n; ++i) {
            if (std::abs(factor(i, j)) > tolerance) {
              return false;
            }
          }
        }
        return true;
      }
      /* Exchange rows and columns k and pivot within the lower triangle. */
      if (pivot != k) {
        for (Index j = 0; j < k; ++j) {
          std::swap(factor(k, j), factor(pivot, j));
        }
        std::swap(factor(k, k), factor(pivot, pivot));
        for (Index i = k + 1; i < pivot; ++i) {
          std::swap(factor(i, k), factor(pivot, i));
        }
        for (Index i = pivot + 1; i < n; ++i) {
          std::swap(factor(i, k), factor(i, pivot));
        }
        std::swap(order(k), order(pivot));
      }
      const Scalar inverse = Scalar(1) / factor(k, k);
      for (Index j = k + 1; j < n; ++j) {
        const Scalar multiple = factor(j, k) * inverse;
        for (Index i = j; i < n; ++i) {
          factor(i, j) -= factor(i, k) * multiple;
        }
      }
      rank = k + 1;
    }
    return true;
  }

  /*
  Solve K S = P for the gain K over the first rank pivots, through the
  generalised inverse P [L11^-T D^-1 L11^-1, 0; 0, 0] P^T of the factored
  S, and give the remaining directions no gain.
  */
  template <typename CrossType>
  void solve(const CrossType& cross, CrossType& gain) const {
    const Index rows = cross.rows();
    const Index n = factor.rows();
    CrossType work(rows, n);
    for (Index k = 0; k < n; ++k) {
      for (Index r = 0; r < rows; ++r) {
        work(r, k) = cross(r, order(k)) * scale(order(k));
      }
    }
    for (Index k = 0; k < rank; ++k) {
      const Scalar inverse = Scalar(1) / factor(k, k);
      for (Index j = k + 1; j < rank; ++j) {
        const Scalar multiple = factor(j, k) * inverse;
        for (Index r = 0; r < rows; ++r) {
          work(r, j) -= work(r, k) * multiple;
        }
      }
      for (Index r = 0; r < rows; ++r) {
        work(r, k) *= inverse;
      }
    }
    /* A single row has no back substitution; leaving it out at compile time
       keeps GCC's array bounds analysis quiet. */
    if constexpr (MatrixType::MaxColsAtCompileTime != 1) {
      for (Index i = rank - 1; i > 0; --i) {
        for (Index j = 0; j < i; ++j) {
          const Scalar multiple = factor(i, j) / factor(j, j);
          for (Index r = 0; r < rows; ++r) {
            work(r, j) -= work(r, i) * multiple;
          }
        }
      }
    }
    gain.setZero(rows, n);
    for (Index k = 0; k < rank; ++k) {
      for (Index r = 0; r < rows; ++r) {
        gain(r, order(k)) = work(r, k) * scale(order(k));
      }
    }
  }
};

/*
Solve K S = P for the Kalman gain K, given the innovation covariance S and
the cross-correlation P, leaving out the directions in which S cannot vary.
An exactly zero S leaves nothing to learn and gives a zero gain. An S that is
not a covariance, which a split-step caller can write into the public
innovation covariance, is refused rather than turned into a gain.
*/
template <typename MatrixType, typename CrossType>
bool kalman_gain(const MatrixType& covariance, const CrossType& cross,
                 CrossType& gain) {
  if (!cross.allFinite()) {
    return false;
  }
  SemidefiniteFactor<MatrixType> factor;
  if (!factor.compute(covariance)) {
    return false;
  }
  factor.solve(cross, gain);
  return gain.allFinite();
}

/* Make a computed covariance exactly symmetric. */
template <typename MatrixType>
void symmetrise(MatrixType& covariance) {
  using Index = typename MatrixType::Index;
  for (Index j = 0; j < covariance.cols(); ++j) {
    for (Index i = j + 1; i < covariance.rows(); ++i) {
      const typename MatrixType::Scalar mean =
          (covariance(i, j) + covariance(j, i)) / 2;
      covariance(i, j) = mean;
      covariance(j, i) = mean;
    }
  }
}

/*
A successful Cholesky factor already preserves the units of each observation.
Equilibrate a singular innovation before the spectral fallback, so the size
of an independent sensor's variance in its chosen units does not determine
its numerical rank. Restoring the row scales gives a root in the original
units; whiten tests each pivot against its own row.
*/
template <typename MatrixType>
bool innovation_root(const MatrixType& covariance, MatrixType& root) {
  using Scalar = typename MatrixType::Scalar;
  using Index = typename MatrixType::Index;
  if (!covariance.allFinite()) {
    return false;
  }
  const Eigen::LLT<MatrixType> llt(covariance);
  if (llt.info() == Eigen::Success) {
    root = llt.matrixL();
    return root.allFinite();
  }

  typename MatrixType::ColXpr::PlainObject scale(covariance.rows());
  MatrixType correlation = covariance;
  for (Index i = 0; i < covariance.rows(); ++i) {
    scale(i) =
        covariance(i, i) > Scalar(0) ? std::sqrt(covariance(i, i)) : Scalar(1);
    correlation.row(i) /= scale(i);
  }
  for (Index i = 0; i < covariance.cols(); ++i) {
    correlation.col(i) /= scale(i);
  }
  if (!covariance_root(correlation, root)) {
    return false;
  }
  for (Index i = 0; i < root.rows(); ++i) {
    if (covariance(i, i) == Scalar(0)) {
      /* A zero variance in a PSD matrix means an exactly zero row.
         The eigensolver may leave roundoff here, which the subsequent
         row-relative whitening would mistake for an observable pivot. */
      root.row(i).setZero();
    } else {
      root.row(i) *= scale(i);
    }
  }
  return root.allFinite();
}

/* Draw a fresh measurement prior in the state's local error coordinates. */
template <typename S>
bool draw_measurement_prior(const S& state,
                            const typename S::CovarianceMatrix& root,
                            typename S::SigmaPointDistribution& points,
                            typename S::SigmaPointDeltas& deltas) {
  if (!state.allFinite() || !root.allFinite()) {
    return false;
  }
  points = state.calculate_sigma_point_distribution(
      root * std::sqrt(Parameters::ScaledDimension<S>));
  deltas = state.calculate_sigma_point_deltas(points);
  return points.allFinite() && deltas.allFinite();
}

/* The posterior commits state and uncertainty together. Only the innovation
   diagnostics need to be restored if a valid innovation cannot be applied. */
template <typename Innovation, typename Covariance, typename Correction>
bool correct_measurement(Innovation& innovation,
                         Covariance& innovation_covariance,
                         Correction&& correct) {
  const Innovation previous_innovation = innovation;
  const Covariance previous_covariance = innovation_covariance;
  if (correct()) {
    return true;
  }
  innovation = previous_innovation;
  innovation_covariance = previous_covariance;
  return false;
}

}  // namespace Detail

/*
UKF core class. This class contains the intermediate values used in the filter
iteration, and the methods used to carry out a filter iteration itself.

The algorithm used here is based largely on the work of Edgar Kraft, described
in the paper "A Quaternion-based Unscented Kalman Filter for Orientation
Tracking", retrieved from:
http://kodlab.seas.upenn.edu/uploads/Arun/UKFpaper.pdf
Comments in the code will occasionally refer to equations or sections in the
paper.

The attitude-related code makes use of the MRP method described in the paper
"Unscented Filtering for Spacecraft Attitude Estimation" by John L. Crassidis
and F. Landis Markley.

Every step returns true on success. On failure it returns false and leaves
the public state, covariance, innovation and innovation covariance exactly as
they were; later steps that depend on it also return false until the next
successful a_priori_step. step() restores them even when only its correction
fails. A step fails when a covariance is not finite or is indefinite beyond
roundoff. The cores never throw.
*/
template <typename StateVectorType, typename MeasurementVectorType,
          typename IntegratorType = NotUsedIntegrator,
          typename MeasurementNoiseType =
              typename MeasurementVectorType::CovarianceVector>
class Core {
 public:
  /*
  The innovation has one entry for each component of the measurement's
  uncertainty: a quaternion field's is its three-dimensional rotation.
  Without quaternion fields it is the measurement vector type.
  */
  using InnovationVectorType = typename MeasurementVectorType::InnovationVector;

 protected:
  typename StateVectorType::SigmaPointDistribution sigma_points;
  typename StateVectorType::SigmaPointDeltas w_prime;
  typename MeasurementVectorType::template SigmaPointDeltas<StateVectorType>
      z_prime;

  /* Whether the sigma points and innovation are ready for the next step. */
  bool prior_valid = false;
  bool innovation_valid = false;

  /* Aliases needed during filter iteration. */
  using CrossCorrelation =
      Eigen::Matrix<real_t, StateVectorType::covariance_size(),
                    InnovationVectorType::RowsAtCompileTime,
                    StateVectorType::covariance_size() == 1 ? Eigen::RowMajor
                                                            : Eigen::ColMajor,
                    StateVectorType::covariance_size(),
                    InnovationVectorType::MaxRowsAtCompileTime>;

 public:
  /* State and covariance are public for simplicity of initialisation. */
  StateVectorType state;
  typename StateVectorType::CovarianceMatrix covariance;
  typename StateVectorType::CovarianceMatrix process_noise_covariance;

  /*
  Measurement noise covariance is an attribute of the filter instance,
  rather than the MeasurementVector specialisation.
  Defaults to a vector of diagonal variances. Use DenseCovarianceMatrix as
  MeasurementNoiseType to retain correlations between sensor components.
  */
  MeasurementNoiseType measurement_covariance;

  /*
  Innovation and innovation covariance are public to ease implementation of
  filter health monitoring.
  */
  InnovationVectorType innovation;
  typename MeasurementVectorType::CovarianceMatrix innovation_covariance;

  /*
  The sigma points drawn for the most recent prediction or measurement
  preparation, one per column in state layout (quaternions as x, y, z, w).
  Prediction includes the process noise. Measurement preparation draws from
  the current uncertainty. The innovation step passes these points through
  the measurement model and the posterior leaves them available for inspection.
  They are initially zero and are unspecified after a failed step.
  */
  const typename StateVectorType::SigmaPointDistribution& get_sigma_points()
      const {
    return sigma_points;
  }

  Core() {
    /*
    Give the fixed-size values which the steps fill in a defined state
    before any step, so that a filter can be copied or inspected at any
    time.
    */
    sigma_points.setZero();
    w_prime.setZero();
    if constexpr (MeasurementVectorType::RowsAtCompileTime != Eigen::Dynamic) {
      innovation.setZero();
      innovation_covariance.setZero();
    }
  }

  /* Prepare the current state for a split measurement update, without
     running the process model or adding process noise. The sigma points are
     drawn afresh from the stored state and uncertainty, so this replaces
     the points of any earlier prediction and factors the uncertainty once
     more. */
  bool prepare_measurement_step() {
    innovation_valid = false;
    typename StateVectorType::CovarianceMatrix root;
    prior_valid =
        Detail::covariance_root(covariance, root) &&
        Detail::draw_measurement_prior(state, root, sigma_points, w_prime);
    return prior_valid;
  }

  /* Fuse an independent observation at the current time. On numerical failure
     the public state, uncertainty and innovation diagnostics are preserved.
     Straight after a_priori_step, innovation_step and a_posteriori_step use
     the predicted points; this call draws new ones. */
  template <typename... U>
  bool measurement_step(const MeasurementVectorType& z, const U&... input) {
    const bool ok =
        prepare_measurement_step() &&
        Detail::correct_measurement(innovation, innovation_covariance, [&] {
          return innovation_step(z, input...) && a_posteriori_step();
        });
    prior_valid = false;
    innovation_valid = false;
    return ok;
  }

  /*
  Top-level function used to carry out a filter step. It is all-or-nothing:
  if the correction fails after a successful prediction, every public value
  is restored. Use the split steps to keep a prediction whose measurement
  cannot be applied.
  */
  template <typename... U>
  bool step(real_t delta, const MeasurementVectorType& z, const U&... input) {
    const StateVectorType previous_state = state;
    const typename StateVectorType::CovarianceMatrix previous_covariance =
        covariance;
    const InnovationVectorType previous_innovation = innovation;
    const typename MeasurementVectorType::CovarianceMatrix
        previous_innovation_covariance = innovation_covariance;
    if (a_priori_step(delta, input...) && innovation_step(z, input...) &&
        a_posteriori_step()) {
      return true;
    }
    state = previous_state;
    covariance = previous_covariance;
    innovation = previous_innovation;
    innovation_covariance = previous_innovation_covariance;
    prior_valid = false;
    innovation_valid = false;
    return false;
  }

  /*
  Apply a correction only when a gate accepts its valid innovation. A
  rejected or empty observation keeps the prediction and diagnostics;
  numerical failure restores the public values from before this call,
  matching step(). False can therefore mean rejection, empty input or
  failure. Application callbacks are invoked synchronously and not copied.
  */
  template <typename Gate, typename... U>
  bool step_with_gate(real_t delta, const MeasurementVectorType& z, Gate&& gate,
                      const U&... input) {
    return step_with_gate_result(delta, z, std::forward<Gate>(gate),
                                 input...) == GateResult::Applied;
  }

  /* Report the outcome of the prediction, gate and correction. */
  template <typename Gate, typename... U>
  GateResult step_with_gate_result(real_t delta, const MeasurementVectorType& z,
                                   Gate&& gate, const U&... input) {
    const StateVectorType previous_state = state;
    const typename StateVectorType::CovarianceMatrix previous_covariance =
        covariance;
    const InnovationVectorType previous_innovation = innovation;
    const typename MeasurementVectorType::CovarianceMatrix
        previous_innovation_covariance = innovation_covariance;
    if (a_priori_step(delta, input...) && innovation_step(z, input...)) {
      if (innovation.size() == 0) {
        return GateResult::Empty;
      }
      if (!std::forward<Gate>(gate)(static_cast<const Core&>(*this), z,
                                    input...)) {
        return GateResult::Rejected;
      }
      if (a_posteriori_step()) {
        return GateResult::Applied;
      }
    }
    state = previous_state;
    covariance = previous_covariance;
    innovation = previous_innovation;
    innovation_covariance = previous_innovation_covariance;
    prior_valid = false;
    innovation_valid = false;
    return GateResult::Failed;
  }

  /*
  The a priori step calculates the sigma point distribution using the
  previous state covariance and previous state estimate, and then
  propagates all sigma points through the process model.
  From the transformed sigma point distribution, the a priori mean and
  covariance are calculated.

  Returns false, without changing the state or covariance, if either the
  prior or the predicted covariance is not finite or is indefinite.
  */
  template <typename... U>
  bool a_priori_step(real_t delta, const U&... input) {
    return predict(nullptr, delta, input...);
  }

  /*
  Also return Cov(x_k, x_{k+1}|z_{1:k}) for an RTS smoother. The two
  deviations in each outer product come from the same process trajectory.
  Independent additive process noise contributes to the prediction but
  not to the cross covariance. On failure the output is unchanged.
  Experimental in 1.0, with the smoother: see RELEASING.md.
  */
  template <typename... U>
  bool a_priori_step_with_cross_covariance(
      typename StateVectorType::CovarianceMatrix& cross, real_t delta,
      const U&... input) {
    return predict(&cross, delta, input...);
  }

 private:
  template <typename Cross, typename... U>
  bool predict(Cross cross, real_t delta, const U&... input) {
    prior_valid = false;
    innovation_valid = false;

    /*
    Generate sigma points from the prior covariance. Process noise is
    additive in the predicted state: x_next = f(x, input) + noise.
    */
    const real_t sigma_scale =
        std::sqrt(Parameters::ScaledDimension<StateVectorType>);
    typename StateVectorType::CovarianceMatrix root;
    if (!Detail::covariance_root(covariance, root)) {
      return false;
    }
    sigma_points = state.calculate_sigma_point_distribution(root * sigma_scale);

    constexpr bool trace = !std::is_null_pointer_v<Cross>;
    [[maybe_unused]] std::conditional_t<
        trace, typename StateVectorType::SigmaPointDeltas, bool> prior_deltas;
    if constexpr (trace) {
      prior_deltas = state.calculate_sigma_point_deltas(sigma_points);
    }

    /* Propagate the sigma points through the process model. */
    Detail::propagate_sigma_points<StateVectorType, IntegratorType>(
        sigma_points, delta, input...);

    /* Calculate the a priori estimate mean, deltas and covariance. */
    const StateVectorType predicted_state =
        StateVectorType::calculate_sigma_point_mean(sigma_points);
    w_prime = predicted_state.calculate_sigma_point_deltas(sigma_points);
    typename StateVectorType::CovarianceMatrix predicted_covariance =
        StateVectorType::calculate_sigma_point_covariance(w_prime);
    predicted_covariance += process_noise_covariance;
    if (!predicted_state.allFinite() ||
        !Detail::covariance_root(predicted_covariance, root, real_t(0),
                                 &predicted_covariance)) {
      return false;
    }

    if constexpr (trace) {
      typename StateVectorType::CovarianceMatrix candidate;
      Detail::calculate_cross_correlation<StateVectorType>(
          candidate, prior_deltas, w_prime);
      if (!candidate.allFinite()) return false;
      *cross = candidate;
    }

    /* Commit the prediction only once it is known to be valid. */
    state = predicted_state;
    covariance = predicted_covariance;

    /* Include the additive noise in measurement and cross covariance. */
    sigma_points = state.calculate_sigma_point_distribution(root * sigma_scale);
    w_prime = state.calculate_sigma_point_deltas(sigma_points);
    prior_valid = true;
    return true;
  }

 public:
  /*
  In the innovation step, the a priori sigma point distribution is further
  propagated using the measurement model. For the measurement sigma point
  distribution, the estimated measurement mean and covariance are
  calculated. A measurement vector is then supplied, and the innovation
  and innovation covariance are calculated.

  Returns false if there is no valid a priori step to build on, or if the
  innovation is not finite.
  */
  template <typename... U>
  bool innovation_step(const MeasurementVectorType& z, const U&... input) {
    innovation_valid = false;
    if (!prior_valid) {
      return false;
    }

    /* Propagate the sigma points through the measurement model. */
    typename MeasurementVectorType::template SigmaPointDistribution<
        StateVectorType>
        measurement_sigma_points =
            z.template calculate_sigma_point_distribution<StateVectorType>(
                sigma_points, input...);

    /* Calculate the measurement prediction mean and deltas. */
    MeasurementVectorType z_pred =
        z.template calculate_sigma_point_mean<StateVectorType>(
            measurement_sigma_points);
    z_prime = z_pred.template calculate_sigma_point_deltas<StateVectorType>(
        measurement_sigma_points);

    /*
    Calculate innovation and innovation covariance.
    See equations 44 and 45 from the Kraft paper for details.
    */
    const InnovationVectorType predicted_innovation =
        z_pred.calculate_innovation(z);
    const typename MeasurementVectorType::CovarianceMatrix spread =
        z_pred.template calculate_sigma_point_covariance<StateVectorType>(
            z_prime);
    const typename MeasurementVectorType::CovarianceMatrix noise =
        z.calculate_measurement_covariance(measurement_covariance, z_pred);
    typename MeasurementVectorType::CovarianceMatrix
        predicted_innovation_covariance = spread + noise;
    Detail::symmetrise(predicted_innovation_covariance);
    /* The innovation covariance must be positive semidefinite. */
    typename MeasurementVectorType::CovarianceMatrix innovation_root;
    if (!predicted_innovation.allFinite() ||
        (predicted_innovation_covariance.size() != 0 &&
         !Detail::innovation_root(predicted_innovation_covariance,
                                  innovation_root))) {
      return false;
    }

    innovation = predicted_innovation;
    innovation_covariance = predicted_innovation_covariance;
    innovation_valid = true;
    return true;
  }

  /*
  The a posteriori step solves for the Kalman gain and applies the
  correction. It is split from the innovation step to allow consistency
  checks before applying the observation. A successful, nonempty correction
  consumes its prediction and innovation.

  Returns false if there is no valid innovation step to build on, or if the
  corrected state or covariance is not finite.
  */
  bool a_posteriori_step() {
    if (!innovation_valid) {
      return false;
    }

    /* An empty dynamic measurement leaves the prediction available. */
    if (innovation.size() == 0) {
      return true;
    }

    /*
    Calculate the cross-correlation matrix described in equations 70 and
    71 from from the Kraft paper.
    */
    CrossCorrelation cross_correlation(StateVectorType::covariance_size(),
                                       innovation.size());
    Detail::calculate_cross_correlation<StateVectorType>(cross_correlation,
                                                         w_prime, z_prime);

    /*
    Calculate the Kalman gain of equation 72 from the Kraft paper, by
    solving K * innovation_covariance = cross_correlation rather than
    forming an inverse; a singular innovation is left out in the directions
    in which it cannot vary.
    */
    CrossCorrelation kalman_gain(StateVectorType::covariance_size(),
                                 innovation.size());
    if (!Detail::kalman_gain(innovation_covariance, cross_correlation,
                             kalman_gain)) {
      return false;
    }
    StateVectorType updated_state = state;
    const typename StateVectorType::StateVectorDelta update_delta =
        kalman_gain * innovation;
    updated_state.apply_delta(update_delta);

    /* Update the covariance using equation 75 from the Kraft paper. */
    typename StateVectorType::CovarianceMatrix updated_covariance = covariance;
    updated_covariance.noalias() -=
        kalman_gain * innovation_covariance * kalman_gain.transpose();
    /*
    Reject an indefinite posterior, allowing for roundoff in the update
    against the prior's scale. A posterior that is semidefinite within
    that roundoff is stored as repaired, so that the next prediction can
    take its root.
    */
    typename StateVectorType::CovarianceMatrix posterior_root;
    if (!updated_state.allFinite() ||
        !Detail::covariance_root(updated_covariance, posterior_root,
                                 covariance.cwiseAbs().maxCoeff(),
                                 &updated_covariance)) {
      return false;
    }

    state = updated_state;
    covariance = updated_covariance;
    prior_valid = false;
    innovation_valid = false;
    return true;
  }

  /*
  Split-step gating hook. A missing or failed innovation never calls the
  gate. Rejection retains prediction and diagnostics; an accepted gate
  returns the actual correction result. Context is passed unchanged.
  */
  template <typename Gate, typename... U>
  bool a_posteriori_step(Gate&& gate, const U&... context) {
    return a_posteriori_step_with_gate_result(
               std::forward<Gate>(gate), context...) == GateResult::Applied;
  }

  template <typename Gate, typename... U>
  GateResult a_posteriori_step_with_gate_result(Gate&& gate,
                                                const U&... context) {
    if (!innovation_valid) {
      return GateResult::Failed;
    }
    if (innovation.size() == 0) {
      return GateResult::Empty;
    }
    if (!std::forward<Gate>(gate)(static_cast<const Core&>(*this),
                                  context...)) {
      return GateResult::Rejected;
    }
    return a_posteriori_step() ? GateResult::Applied : GateResult::Failed;
  }
};

/*
The class implements the square-root form of the UKF, described in the paper
"The Square-Root Unscented Kalman Filter for State and Parameter-Estimation" by
Rudolph van der Merwe and Eric A. Wan, retrieved from:
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.80.1421&rep=rep1&type=pdf

The SR-UKF has a number of numerical and stability advantages over the
standard UKF, and a lower computational cost. Steps report failure in the
same way as Core.
*/
template <typename StateVectorType, typename MeasurementVectorType,
          typename IntegratorType = NotUsedIntegrator,
          typename MeasurementNoiseType =
              typename MeasurementVectorType::CovarianceVector>
class SquareRootCore {
 public:
  /*
  The innovation has one entry for each component of the measurement's
  uncertainty: a quaternion field's is its three-dimensional rotation.
  Without quaternion fields it is the measurement vector type.
  */
  using InnovationVectorType = typename MeasurementVectorType::InnovationVector;

 protected:
  typename StateVectorType::SigmaPointDistribution sigma_points;
  typename StateVectorType::SigmaPointDeltas w_prime;
  typename MeasurementVectorType::template SigmaPointDeltas<StateVectorType>
      z_prime;

  /* Whether the sigma points and innovation are ready for the next step. */
  bool prior_valid = false;
  bool innovation_valid = false;

  /* Aliases needed during filter iteration. */
  using CrossCorrelation =
      Eigen::Matrix<real_t, StateVectorType::covariance_size(),
                    InnovationVectorType::RowsAtCompileTime,
                    StateVectorType::covariance_size() == 1 ? Eigen::RowMajor
                                                            : Eigen::ColMajor,
                    StateVectorType::covariance_size(),
                    InnovationVectorType::MaxRowsAtCompileTime>;

 public:
  /*
  For the SR-UKF, covariance is stored in square-root form. The variable
  is named accordingly to avoid mix-ups.
  */
  StateVectorType state;
  typename StateVectorType::CovarianceMatrix root_covariance;
  typename StateVectorType::CovarianceMatrix process_noise_root_covariance;

  /*
  Measurement noise root covariance is an attribute of the filter instance,
  rather than the MeasurementVector specialisation.
  Defaults to a vector of standard deviations. Use DenseCovarianceMatrix
  as MeasurementNoiseType for a full factor L with noise covariance L*L^T.
  */
  MeasurementNoiseType measurement_root_covariance;

  /*
  Innovation covariance is also stored in square-root form as described
  above.
  */
  InnovationVectorType innovation;
  typename MeasurementVectorType::CovarianceMatrix innovation_root_covariance;

  /*
  The sigma points drawn for the most recent prediction or measurement
  preparation, one per column in state layout (quaternions as x, y, z, w).
  Prediction includes the process noise. Measurement preparation draws from
  the current uncertainty. The innovation step passes these points through
  the measurement model and the posterior leaves them available for inspection.
  They are initially zero and are unspecified after a failed step.
  */
  const typename StateVectorType::SigmaPointDistribution& get_sigma_points()
      const {
    return sigma_points;
  }

  SquareRootCore() {
    /*
    Give the fixed-size values which the steps fill in a defined state
    before any step, so that a filter can be copied at any time.
    */
    sigma_points.setZero();
    w_prime.setZero();
    if constexpr (MeasurementVectorType::RowsAtCompileTime != Eigen::Dynamic) {
      z_prime.setZero();
      innovation.setZero();
      innovation_root_covariance.setZero();
    }
  }

  /* Prepare the current state for a split measurement update, without
     running the process model or adding process noise. The sigma points are
     drawn afresh from the stored state and root covariance, so this replaces
     the points of any earlier prediction; the stored root is used as it is,
     with no factorisation. */
  bool prepare_measurement_step() {
    innovation_valid = false;
    prior_valid = Detail::draw_measurement_prior(state, root_covariance,
                                                 sigma_points, w_prime);
    return prior_valid;
  }

  /* Fuse an independent observation at the current time. On numerical failure
     the public state, uncertainty and innovation diagnostics are preserved.
     Straight after a_priori_step, innovation_step and a_posteriori_step use
     the predicted points; this call draws new ones. */
  template <typename... U>
  bool measurement_step(const MeasurementVectorType& z, const U&... input) {
    const bool ok =
        prepare_measurement_step() &&
        Detail::correct_measurement(
            innovation, innovation_root_covariance, [&] {
              return innovation_step(z, input...) && a_posteriori_step();
            });
    prior_valid = false;
    innovation_valid = false;
    return ok;
  }

  /*
  Top-level function used to carry out a filter step. Like Core::step, it
  restores every public value if any stage fails.
  */
  template <typename... U>
  bool step(real_t delta, const MeasurementVectorType& z, const U&... input) {
    const StateVectorType previous_state = state;
    const typename StateVectorType::CovarianceMatrix previous_root =
        root_covariance;
    const InnovationVectorType previous_innovation = innovation;
    const typename MeasurementVectorType::CovarianceMatrix
        previous_innovation_root = innovation_root_covariance;
    if (a_priori_step(delta, input...) && innovation_step(z, input...) &&
        a_posteriori_step()) {
      return true;
    }
    state = previous_state;
    root_covariance = previous_root;
    innovation = previous_innovation;
    innovation_root_covariance = previous_innovation_root;
    prior_valid = false;
    innovation_valid = false;
    return false;
  }

  /*
  Apply a correction only when a gate accepts its valid innovation. A
  rejected or empty observation keeps the prediction and diagnostics;
  numerical failure restores the public values from before this call,
  matching step(). False can therefore mean rejection, empty input or
  failure. Application callbacks are invoked synchronously and not copied.
  */
  template <typename Gate, typename... U>
  bool step_with_gate(real_t delta, const MeasurementVectorType& z, Gate&& gate,
                      const U&... input) {
    return step_with_gate_result(delta, z, std::forward<Gate>(gate),
                                 input...) == GateResult::Applied;
  }

  /* Report the outcome of the prediction, gate and correction. */
  template <typename Gate, typename... U>
  GateResult step_with_gate_result(real_t delta, const MeasurementVectorType& z,
                                   Gate&& gate, const U&... input) {
    const StateVectorType previous_state = state;
    const typename StateVectorType::CovarianceMatrix previous_covariance =
        root_covariance;
    const InnovationVectorType previous_innovation = innovation;
    const typename MeasurementVectorType::CovarianceMatrix
        previous_innovation_covariance = innovation_root_covariance;
    if (a_priori_step(delta, input...) && innovation_step(z, input...)) {
      if (innovation.size() == 0) {
        return GateResult::Empty;
      }
      if (!std::forward<Gate>(gate)(static_cast<const SquareRootCore&>(*this),
                                    z, input...)) {
        return GateResult::Rejected;
      }
      if (a_posteriori_step()) {
        return GateResult::Applied;
      }
    }
    state = previous_state;
    root_covariance = previous_covariance;
    innovation = previous_innovation;
    innovation_root_covariance = previous_innovation_covariance;
    prior_valid = false;
    innovation_valid = false;
    return GateResult::Failed;
  }

  /*
  The a priori step differs from the standard UKF in that it doesn't
  compute the LLT decomposition of the state covariance matrix.

  In addition, the a priori state root covariance is recovered from the
  sigma point distribution using a QR decomposition. The central sigma
  point is part of it unless its weight is negative (small alpha); it is
  then applied as a rank-one Cholesky downdate. If that cannot be applied,
  the factor is rebuilt from the covariance, and the step fails without
  changing the filter if the predicted covariance is indefinite.
  */
  template <typename... U>
  bool a_priori_step(real_t delta, const U&... input) {
    return predict(nullptr, delta, input...);
  }

  /*
  Also return Cov(x_k, x_{k+1}|z_{1:k}) for an RTS smoother. The two
  deviations in each outer product come from the same process trajectory.
  Independent additive process noise contributes to the prediction but
  not to the cross covariance. On failure the output is unchanged.
  Experimental in 1.0, with the smoother: see RELEASING.md.
  */
  template <typename... U>
  bool a_priori_step_with_cross_covariance(
      typename StateVectorType::CovarianceMatrix& cross, real_t delta,
      const U&... input) {
    return predict(&cross, delta, input...);
  }

 private:
  template <typename Cross, typename... U>
  bool predict(Cross cross, real_t delta, const U&... input) {
    prior_valid = false;
    innovation_valid = false;

    const real_t sigma_scale =
        std::sqrt(Parameters::ScaledDimension<StateVectorType>);
    sigma_points =
        state.calculate_sigma_point_distribution(root_covariance * sigma_scale);

    constexpr bool trace = !std::is_null_pointer_v<Cross>;
    [[maybe_unused]] std::conditional_t<
        trace, typename StateVectorType::SigmaPointDeltas, bool> prior_deltas;
    if constexpr (trace) {
      prior_deltas = state.calculate_sigma_point_deltas(sigma_points);
    }

    /* Propagate the sigma points through the process model. */
    Detail::propagate_sigma_points<StateVectorType, IntegratorType>(
        sigma_points, delta, input...);

    /*
    Calculate the a priori estimate mean and deltas, and from the deltas
    and the process noise the a priori root covariance.
    */
    const StateVectorType predicted_state =
        StateVectorType::calculate_sigma_point_mean(sigma_points);
    w_prime = predicted_state.calculate_sigma_point_deltas(sigma_points);
    typename StateVectorType::CovarianceMatrix predicted_root;
    if (!predicted_state.allFinite() ||
        !Detail::sigma_point_root<StateVectorType>(
            predicted_root, w_prime, process_noise_root_covariance)) {
      return false;
    }

    if constexpr (trace) {
      typename StateVectorType::CovarianceMatrix candidate;
      Detail::calculate_cross_correlation<StateVectorType>(
          candidate, prior_deltas, w_prime);
      if (!candidate.allFinite()) return false;
      *cross = candidate;
    }

    /* Commit the prediction only once it is known to be valid. */
    state = predicted_state;
    root_covariance = predicted_root;

    /* Recalculate the sigma points using the a priori covariance. */
    sigma_points =
        state.calculate_sigma_point_distribution(root_covariance * sigma_scale);
    w_prime = state.calculate_sigma_point_deltas(sigma_points);
    prior_valid = true;
    return true;
  }

 public:
  /*
  The innovation step differs from the standard UKF in the same way as the
  a priori step; it recovers the innovation root covariance using a QR
  decomposition and, when the central weight is negative, a rank-one
  Cholesky downdate.
  */
  template <typename... U>
  bool innovation_step(const MeasurementVectorType& z, const U&... input) {
    innovation_valid = false;
    if (!prior_valid) {
      return false;
    }

    /* An empty dynamic measurement means prediction only. */
    if constexpr (MeasurementVectorType::RowsAtCompileTime == Eigen::Dynamic) {
      if (z.size() == 0) {
        innovation = InnovationVectorType();
        innovation_root_covariance.resize(0, 0);
        z_prime.resize(0, StateVectorType::num_sigma());
        innovation_valid = true;
        return true;
      }
    }

    /* Propagate the sigma points through the measurement model. */
    typename MeasurementVectorType::template SigmaPointDistribution<
        StateVectorType>
        measurement_sigma_points =
            z.template calculate_sigma_point_distribution<StateVectorType>(
                sigma_points, input...);

    /* Calculate the measurement prediction mean and deltas. */
    MeasurementVectorType z_pred =
        z.template calculate_sigma_point_mean<StateVectorType>(
            measurement_sigma_points);
    z_prime = z_pred.template calculate_sigma_point_deltas<StateVectorType>(
        measurement_sigma_points);

    /*
    Calculate the innovation, and from the innovation deltas and the
    measurement noise the innovation root covariance.
    */
    const InnovationVectorType predicted_innovation =
        z_pred.calculate_innovation(z);
    typename MeasurementVectorType::CovarianceMatrix predicted_root;
    const auto noise_root = z.calculate_measurement_root_covariance(
        measurement_root_covariance, z_pred);
    if (!predicted_innovation.allFinite() ||
        !Detail::sigma_point_root<StateVectorType>(predicted_root, z_prime,
                                                   noise_root)) {
      return false;
    }

    innovation = predicted_innovation;
    innovation_root_covariance = predicted_root;
    innovation_valid = true;
    return true;
  }

  /*
  Compared to the standard UKF, the a posteriori step whitens the
  cross-correlation with a triangular solve rather than calculating the
  Kalman gain, and then updates the root covariance using a series of
  rank-one Cholesky downdates. Nothing is changed unless every downdate
  succeeds.
  */
  bool a_posteriori_step() {
    if (!innovation_valid) {
      return false;
    }
    if (innovation.size() == 0) {
      return true;
    }

    /*
    Calculate the cross-correlation matrix described in equations 70 and
    71 from from the Kraft paper.
    */
    CrossCorrelation cross_correlation(StateVectorType::covariance_size(),
                                       innovation.size());
    Detail::calculate_cross_correlation<StateVectorType>(cross_correlation,
                                                         w_prime, z_prime);

    /*
    Whiten the cross-correlation and the innovation by the innovation
    root covariance. Their product is the update delta vector, which is
    applied to a copy of the a priori estimate.
    */
    InnovationVectorType whitened_innovation = innovation;
    Detail::whiten(cross_correlation, whitened_innovation,
                   innovation_root_covariance);
    StateVectorType updated_state = state;
    const typename StateVectorType::StateVectorDelta update_delta =
        cross_correlation * whitened_innovation;
    updated_state.apply_delta(update_delta);

    /*
    Update the root covariance using a series of rank-one Cholesky
    downdates, one for each column of the whitened cross-correlation.
    LLT downdates require strict positivity; an exact observation can
    instead leave a semidefinite covariance, which is rebuilt from the
    original factor rather than the partial downdate. The root covariance
    is left as it was unless every downdate succeeds.
    */
    if (!updated_state.allFinite() ||
        !Detail::update_root(root_covariance, cross_correlation, real_t(-1))) {
      return false;
    }
    state = updated_state;
    prior_valid = false;
    innovation_valid = false;
    return true;
  }

  /*
  Split-step gating hook. A missing or failed innovation never calls the
  gate. Rejection retains prediction and diagnostics; an accepted gate
  returns the actual correction result. Context is passed unchanged.
  */
  template <typename Gate, typename... U>
  bool a_posteriori_step(Gate&& gate, const U&... context) {
    return a_posteriori_step_with_gate_result(
               std::forward<Gate>(gate), context...) == GateResult::Applied;
  }

  template <typename Gate, typename... U>
  GateResult a_posteriori_step_with_gate_result(Gate&& gate,
                                                const U&... context) {
    if (!innovation_valid) {
      return GateResult::Failed;
    }
    if (innovation.size() == 0) {
      return GateResult::Empty;
    }
    if (!std::forward<Gate>(gate)(static_cast<const SquareRootCore&>(*this),
                                  context...)) {
      return GateResult::Rejected;
    }
    return a_posteriori_step() ? GateResult::Applied : GateResult::Failed;
  }
};

/*
The class implements the optimised form of the square-root UKF for parameter
estimation, described in the paper "The Square-Root Unscented Kalman Filter for
State and Parameter-Estimation" by Rudolph van der Merwe and Eric A. Wan,
retrieved from:
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.80.1421&rep=rep1&type=pdf

Steps report failure in the same way as Core.
*/
template <typename StateVectorType, typename MeasurementVectorType,
          typename MeasurementNoiseType =
              typename MeasurementVectorType::CovarianceVector>
class SquareRootParameterEstimationCore {
 public:
  /*
  The innovation has one entry for each component of the measurement's
  uncertainty: a quaternion field's is its three-dimensional rotation.
  Without quaternion fields it is the measurement vector type.
  */
  using InnovationVectorType = typename MeasurementVectorType::InnovationVector;

 protected:
  typename StateVectorType::SigmaPointDistribution sigma_points;
  typename StateVectorType::SigmaPointDeltas w_prime;
  typename MeasurementVectorType::template SigmaPointDeltas<StateVectorType>
      z_prime;

  /* Whether the sigma points and innovation are ready for the next step. */
  bool prior_valid = false;
  bool innovation_valid = false;

  /* Aliases needed during filter iteration. */
  using CrossCorrelation =
      Eigen::Matrix<real_t, StateVectorType::covariance_size(),
                    InnovationVectorType::RowsAtCompileTime,
                    StateVectorType::covariance_size() == 1 ? Eigen::RowMajor
                                                            : Eigen::ColMajor,
                    StateVectorType::covariance_size(),
                    InnovationVectorType::MaxRowsAtCompileTime>;

 public:
  /*
  For the SR-UKF, covariance is stored in square-root form. The variable
  is named accordingly to avoid mix-ups.
  */
  StateVectorType state;
  typename StateVectorType::CovarianceMatrix root_covariance;
  typename StateVectorType::CovarianceMatrix process_noise_root_covariance;

  /*
  Measurement noise root covariance is an attribute of the filter instance,
  rather than the MeasurementVector specialisation.
  Defaults to a vector of standard deviations. Use DenseCovarianceMatrix
  as MeasurementNoiseType for a full factor L with noise covariance L*L^T.
  */
  MeasurementNoiseType measurement_root_covariance;

  /*
  Innovation covariance is also stored in square-root form as described
  above.
  */
  InnovationVectorType innovation;
  typename MeasurementVectorType::CovarianceMatrix innovation_root_covariance;

  /*
  The sigma points drawn for the most recent prediction or measurement
  preparation, one per column in state layout (quaternions as x, y, z, w).
  Prediction includes the process noise. Measurement preparation draws from
  the current uncertainty. The innovation step passes these points through
  the measurement model and the posterior leaves them available for inspection.
  They are initially zero and are unspecified after a failed step.
  */
  const typename StateVectorType::SigmaPointDistribution& get_sigma_points()
      const {
    return sigma_points;
  }

  SquareRootParameterEstimationCore() {
    /*
    Give the fixed-size values which the steps fill in a defined state
    before any step, so that a filter can be copied at any time.
    */
    sigma_points.setZero();
    w_prime.setZero();
    if constexpr (MeasurementVectorType::RowsAtCompileTime != Eigen::Dynamic) {
      z_prime.setZero();
      innovation.setZero();
      innovation_root_covariance.setZero();
    }
  }

  /* Prepare the current state for a split measurement update, without
     running the process model or adding process noise. The sigma points are
     drawn afresh from the stored state and root covariance, so this replaces
     the points of any earlier prediction; the stored root is used as it is,
     with no factorisation. */
  bool prepare_measurement_step() {
    innovation_valid = false;
    prior_valid = Detail::draw_measurement_prior(state, root_covariance,
                                                 sigma_points, w_prime);
    return prior_valid;
  }

  /* Fuse an independent observation at the current time. On numerical failure
     the public state, uncertainty and innovation diagnostics are preserved.
     Straight after a_priori_step, innovation_step and a_posteriori_step use
     the predicted points; this call draws new ones. */
  template <typename... U>
  bool measurement_step(const MeasurementVectorType& z, const U&... input) {
    const bool ok =
        prepare_measurement_step() &&
        Detail::correct_measurement(
            innovation, innovation_root_covariance, [&] {
              return innovation_step(z, input...) && a_posteriori_step();
            });
    prior_valid = false;
    innovation_valid = false;
    return ok;
  }

  /*
  Top-level function used to carry out a filter step. Like Core::step, it
  restores every public value if any stage fails.
  */
  template <typename... U>
  bool step(const MeasurementVectorType& z, const U&... input) {
    const StateVectorType previous_state = state;
    const typename StateVectorType::CovarianceMatrix previous_root =
        root_covariance;
    const InnovationVectorType previous_innovation = innovation;
    const typename MeasurementVectorType::CovarianceMatrix
        previous_innovation_root = innovation_root_covariance;
    if (a_priori_step() && innovation_step(z, input...) &&
        a_posteriori_step()) {
      return true;
    }
    state = previous_state;
    root_covariance = previous_root;
    innovation = previous_innovation;
    innovation_root_covariance = previous_innovation_root;
    prior_valid = false;
    innovation_valid = false;
    return false;
  }

  /*
  Apply a correction only when a gate accepts its valid innovation. A
  rejected or empty observation keeps the prediction and diagnostics;
  numerical failure restores the public values from before this call,
  matching step(). False can therefore mean rejection, empty input or
  failure. Application callbacks are invoked synchronously and not copied.
  */
  template <typename Gate, typename... U>
  bool step_with_gate(const MeasurementVectorType& z, Gate&& gate,
                      const U&... input) {
    return step_with_gate_result(z, std::forward<Gate>(gate), input...) ==
           GateResult::Applied;
  }

  /* Report the outcome of the prediction, gate and correction. */
  template <typename Gate, typename... U>
  GateResult step_with_gate_result(const MeasurementVectorType& z, Gate&& gate,
                                   const U&... input) {
    const StateVectorType previous_state = state;
    const typename StateVectorType::CovarianceMatrix previous_covariance =
        root_covariance;
    const InnovationVectorType previous_innovation = innovation;
    const typename MeasurementVectorType::CovarianceMatrix
        previous_innovation_covariance = innovation_root_covariance;
    if (a_priori_step() && innovation_step(z, input...)) {
      if (innovation.size() == 0) {
        return GateResult::Empty;
      }
      if (!std::forward<Gate>(gate)(
              static_cast<const SquareRootParameterEstimationCore&>(*this), z,
              input...)) {
        return GateResult::Rejected;
      }
      if (a_posteriori_step()) {
        return GateResult::Applied;
      }
    }
    state = previous_state;
    root_covariance = previous_covariance;
    innovation = previous_innovation;
    innovation_root_covariance = previous_innovation_covariance;
    prior_valid = false;
    innovation_valid = false;
    return GateResult::Failed;
  }

  /*
  For the parameter estimation filter, there is no need for a process
  model, which allows the a priori step to be optimised significantly.
  There's no need to propagate the sigma points through a process model,
  and we only have to add the process noise to the root covariance.
  */
  bool a_priori_step() {
    prior_valid = false;
    innovation_valid = false;

    /*
    Recover a factor of S*S' + Q*Q' from the columns of S and Q.
    Updating only S's diagonal would also change the off-diagonal
    covariance whenever the parameters are correlated.
    */
    Matrix<2 * StateVectorType::covariance_size(),
           StateVectorType::covariance_size()>
        vectors;
    vectors << root_covariance.transpose(),
        process_noise_root_covariance.transpose();
    const typename StateVectorType::CovarianceMatrix predicted_root =
        Detail::lower_factor(vectors);
    if (!predicted_root.allFinite() || !state.allFinite()) {
      return false;
    }
    root_covariance = predicted_root;

    sigma_points = state.calculate_sigma_point_distribution(
        root_covariance *
        std::sqrt(Parameters::ScaledDimension<StateVectorType>));

    /* Calculate the sigma point deltas. */
    w_prime = state.calculate_sigma_point_deltas(sigma_points);
    prior_valid = true;
    return true;
  }

  /*
  The innovation step differs from the standard UKF in the same way as the
  a priori step; it recovers the innovation root covariance using a QR
  decomposition and, when the central weight is negative, a rank-one
  Cholesky downdate.
  */
  template <typename... U>
  bool innovation_step(const MeasurementVectorType& z, const U&... input) {
    innovation_valid = false;
    if (!prior_valid) {
      return false;
    }

    /* An empty dynamic measurement means no correction. */
    if constexpr (MeasurementVectorType::RowsAtCompileTime == Eigen::Dynamic) {
      if (z.size() == 0) {
        innovation = InnovationVectorType();
        innovation_root_covariance.resize(0, 0);
        z_prime.resize(0, StateVectorType::num_sigma());
        innovation_valid = true;
        return true;
      }
    }

    /* Propagate the sigma points through the measurement model. */
    typename MeasurementVectorType::template SigmaPointDistribution<
        StateVectorType>
        measurement_sigma_points =
            z.template calculate_sigma_point_distribution<StateVectorType>(
                sigma_points, input...);

    /* Calculate the measurement prediction mean and deltas. */
    MeasurementVectorType z_pred =
        z.template calculate_sigma_point_mean<StateVectorType>(
            measurement_sigma_points);
    z_prime = z_pred.template calculate_sigma_point_deltas<StateVectorType>(
        measurement_sigma_points);

    /*
    Calculate the innovation, and from the innovation deltas and the
    measurement noise the innovation root covariance.
    */
    const InnovationVectorType predicted_innovation =
        z_pred.calculate_innovation(z);
    typename MeasurementVectorType::CovarianceMatrix predicted_root;
    const auto noise_root = z.calculate_measurement_root_covariance(
        measurement_root_covariance, z_pred);
    if (!predicted_innovation.allFinite() ||
        !Detail::sigma_point_root<StateVectorType>(predicted_root, z_prime,
                                                   noise_root)) {
      return false;
    }

    innovation = predicted_innovation;
    innovation_root_covariance = predicted_root;
    innovation_valid = true;
    return true;
  }

  /*
  Compared to the standard UKF, the a posteriori step whitens the
  cross-correlation with a triangular solve rather than calculating the
  Kalman gain, and then updates the root covariance using a series of
  rank-one Cholesky downdates. Nothing is changed unless every downdate
  succeeds.
  */
  bool a_posteriori_step() {
    if (!innovation_valid) {
      return false;
    }
    if (innovation.size() == 0) {
      return true;
    }

    /*
    Calculate the cross-correlation matrix described in equations 70 and
    71 from from the Kraft paper.
    */
    CrossCorrelation cross_correlation(StateVectorType::covariance_size(),
                                       innovation.size());
    Detail::calculate_cross_correlation<StateVectorType>(cross_correlation,
                                                         w_prime, z_prime);

    /*
    Whiten the cross-correlation and the innovation by the innovation
    root covariance. Their product is the update delta vector, which is
    applied to a copy of the a priori estimate.
    */
    InnovationVectorType whitened_innovation = innovation;
    Detail::whiten(cross_correlation, whitened_innovation,
                   innovation_root_covariance);
    StateVectorType updated_state = state;
    const typename StateVectorType::StateVectorDelta update_delta =
        cross_correlation * whitened_innovation;
    updated_state.apply_delta(update_delta);

    /*
    Update the root covariance using a series of rank-one Cholesky
    downdates, rebuilding it from the original factor if the posterior
    is only semidefinite. It is left as it was unless every downdate
    succeeds.
    */
    if (!updated_state.allFinite() ||
        !Detail::update_root(root_covariance, cross_correlation, real_t(-1))) {
      return false;
    }
    state = updated_state;
    prior_valid = false;
    innovation_valid = false;
    return true;
  }

  /*
  Split-step gating hook. A missing or failed innovation never calls the
  gate. Rejection retains prediction and diagnostics; an accepted gate
  returns the actual correction result. Context is passed unchanged.
  */
  template <typename Gate, typename... U>
  bool a_posteriori_step(Gate&& gate, const U&... context) {
    return a_posteriori_step_with_gate_result(
               std::forward<Gate>(gate), context...) == GateResult::Applied;
  }

  template <typename Gate, typename... U>
  GateResult a_posteriori_step_with_gate_result(Gate&& gate,
                                                const U&... context) {
    if (!innovation_valid) {
      return GateResult::Failed;
    }
    if (innovation.size() == 0) {
      return GateResult::Empty;
    }
    if (!std::forward<Gate>(gate)(
            static_cast<const SquareRootParameterEstimationCore&>(*this),
            context...)) {
      return GateResult::Rejected;
    }
    return a_posteriori_step() ? GateResult::Applied : GateResult::Failed;
  }
};

}  // namespace UKF

#endif
