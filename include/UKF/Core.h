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
#include <cmath>
#include <limits>
#include "UKF/Types.h"
#include "UKF/Integrator.h"
#include "UKF/StateVector.h"
#include "UKF/Covariance.h"
#include "UKF/RankUpdate.h"
#include <type_traits>

namespace UKF {

namespace Detail {

/*
Propagate every sigma point through the process model: either the state
vector's direct process_model(delta, input...) specialisation, when the core
uses NotUsedIntegrator, or the integrated derivative.
*/
template <typename StateVectorType, typename IntegratorType, typename... U>
void propagate_sigma_points(typename StateVectorType::SigmaPointDistribution& sigma_points,
        real_t delta, const U&... input) {
    for(std::size_t i = 0; i < StateVectorType::num_sigma(); i++) {
        if constexpr (std::is_same_v<IntegratorType, NotUsedIntegrator>) {
            sigma_points.col(i) << StateVectorType(sigma_points.col(i)).process_model(delta, input...);
        } else {
            sigma_points.col(i) << StateVectorType(sigma_points.col(i))
                .template process_model<IntegratorType>(delta, input...);
        }
    }
}

/*
Calculate the cross-correlation matrix described in equations 70 and 71 from
the Kraft paper, as one weighted product over the outer sigma points plus the
central point's rank-one term.
*/
template <typename StateVectorType, typename ResultType, typename StateDeltas, typename MeasurementDeltas>
void calculate_cross_correlation(ResultType& result, const StateDeltas& w_prime,
        const MeasurementDeltas& z_prime) {
    constexpr int outer = StateVectorType::num_sigma() - 1;
    result.noalias() = (Parameters::Sigma_WCI<StateVectorType> * w_prime.template rightCols<outer>()) *
        z_prime.template rightCols<outer>().transpose();
    result.noalias() += (Parameters::Sigma_WC0<StateVectorType> * w_prime.col(0)) * z_prime.col(0).transpose();
}

/*
Replace the lower triangular factor L of A with a factor of
A + weight * sum(v_i v_i^T) over the columns v_i of vectors. Rank-one
updates or downdates are tried first; if a downdate would lose positive
definiteness, the factor is rebuilt from the updated covariance, which
accepts a semidefinite result within roundoff. Return false, leaving the
factor unchanged, if the result is indefinite or not finite.
*/
template <typename MatrixType, typename VectorsType>
bool update_root(MatrixType& factor, const VectorsType& vectors, typename MatrixType::Scalar weight) {
    MatrixType updated = factor;
    for(typename MatrixType::Index i = 0; i < vectors.cols(); i++) {
        if (rank_update<Eigen::Lower>(updated, vectors.col(i), weight) >= 0) {
            const MatrixType prior = factor * factor.transpose();
            MatrixType covariance = prior;
            covariance.noalias() += (weight * vectors) * vectors.transpose();
            return covariance_root(covariance, factor, prior.cwiseAbs().maxCoeff());
        }
    }
    if (!updated.allFinite()) {
        return false;
    }
    factor = updated;
    return true;
}

/*
Solve K (S S^T) = P for the Kalman gain K, given the cross-correlation P and
the lower triangular innovation root S. This is (S'\(S\P'))', equivalent to
(P/S')/S in the literature. When S is well conditioned, two triangular solves
are used. Otherwise, as for field-vector measurements, whose innovation
covariance has no component along the predicted field, rank-revealing QR
solves are used instead.
*/
template <typename GainType, typename CrossType, typename RootType>
void square_root_gain(GainType& gain, const CrossType& cross, const RootType& root) {
    using Scalar = typename RootType::Scalar;
    const Scalar tolerance = std::sqrt(std::numeric_limits<Scalar>::epsilon());
    const Scalar largest = root.diagonal().cwiseAbs().maxCoeff();
    if (root.diagonal().cwiseAbs().minCoeff() > tolerance * largest) {
        gain = root.transpose().template triangularView<Eigen::Upper>().solve(
            root.template triangularView<Eigen::Lower>().solve(cross.transpose())).transpose();
    } else {
        gain = root.transpose().fullPivHouseholderQr().solve(
            root.fullPivHouseholderQr().solve(cross.transpose())).transpose();
    }
}

}

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
successful a_priori_step. The cores never throw.
*/
template <typename StateVectorType, typename MeasurementVectorType, typename IntegratorType=NotUsedIntegrator>
class Core {
protected:
    typename StateVectorType::SigmaPointDistribution sigma_points;
    typename StateVectorType::SigmaPointDeltas w_prime;
    typename MeasurementVectorType::template SigmaPointDeltas<StateVectorType> z_prime;

    /* Whether the sigma points and innovation are ready for the next step. */
    bool prior_valid = false;
    bool innovation_valid = false;

    /* Aliases needed during filter iteration. */
    using CrossCorrelation = Eigen::Matrix<
        real_t,
        StateVectorType::covariance_size(),
        MeasurementVectorType::RowsAtCompileTime,
        StateVectorType::covariance_size() == 1 ? Eigen::RowMajor : Eigen::ColMajor,
        StateVectorType::covariance_size(),
        MeasurementVectorType::MaxRowsAtCompileTime>;

public:
    /* State and covariance are public for simplicity of initialisation. */
    StateVectorType state;
    typename StateVectorType::CovarianceMatrix covariance;
    typename StateVectorType::CovarianceMatrix process_noise_covariance;

    /*
    Measurement noise covariance is an attribute of the filter instance,
    rather than the MeasurementVector specialisation.
    Only diagonal measurement noise covariance is supported at the moment.
    */
    typename MeasurementVectorType::CovarianceVector measurement_covariance;

    /*
    Innovation and innovation covariance are public to ease implementation of
    filter health monitoring.
    */
    MeasurementVectorType innovation;
    typename MeasurementVectorType::CovarianceMatrix innovation_covariance;

    /* Top-level function used to carry out a filter step. */
    template <typename... U>
    bool step(real_t delta, const MeasurementVectorType& z, const U&... input) {
        return a_priori_step(delta, input...) && innovation_step(z, input...) && a_posteriori_step();
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
        prior_valid = false;
        innovation_valid = false;

        /*
        Generate sigma points from the prior covariance. Process noise is
        additive in the predicted state: x_next = f(x, input) + noise.
        */
        const real_t sigma_scale = std::sqrt(Parameters::ScaledDimension<StateVectorType>);
        typename StateVectorType::CovarianceMatrix root;
        if (!Detail::covariance_root(covariance, root)) {
            return false;
        }
        sigma_points = state.calculate_sigma_point_distribution(root * sigma_scale);

        /* Propagate the sigma points through the process model. */
        Detail::propagate_sigma_points<StateVectorType, IntegratorType>(sigma_points, delta, input...);

        /* Calculate the a priori estimate mean, deltas and covariance. */
        const StateVectorType predicted_state = StateVectorType::calculate_sigma_point_mean(sigma_points);
        w_prime = predicted_state.calculate_sigma_point_deltas(sigma_points);
        typename StateVectorType::CovarianceMatrix predicted_covariance =
            StateVectorType::calculate_sigma_point_covariance(w_prime);
        predicted_covariance += process_noise_covariance;
        if (!predicted_state.allFinite() || !Detail::covariance_root(predicted_covariance, root)) {
            return false;
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
        typename MeasurementVectorType::template SigmaPointDistribution<StateVectorType>
            measurement_sigma_points = z.template calculate_sigma_point_distribution<StateVectorType>(
                sigma_points, input...);

        /* Calculate the measurement prediction mean and deltas. */
        MeasurementVectorType z_pred =
            z.template calculate_sigma_point_mean<StateVectorType>(measurement_sigma_points);
        z_prime = z_pred.template calculate_sigma_point_deltas<StateVectorType>(measurement_sigma_points);

        /*
        Calculate innovation and innovation covariance.
        See equations 44 and 45 from the Kraft paper for details.
        */
        const MeasurementVectorType predicted_innovation = z_pred.calculate_innovation(z);
        typename MeasurementVectorType::CovarianceMatrix predicted_innovation_covariance =
            z_pred.template calculate_sigma_point_covariance<StateVectorType>(z_prime);
        predicted_innovation_covariance += z.calculate_measurement_covariance(measurement_covariance, z_pred);
        if (!predicted_innovation.allFinite() || !predicted_innovation_covariance.allFinite()) {
            return false;
        }

        innovation = predicted_innovation;
        innovation_covariance = predicted_innovation_covariance;
        innovation_valid = true;
        return true;
    }

    /*
    In the a posteriori step, the cross-correlation matrix and innovation
    covariance inverse are calculated and used to calculate the Kalman gain.
    This step is split from the innovation step to allow the user to carry
    out innovation consistency checks and identify failed sensors, if
    desired.

    Returns false if there is no valid innovation step to build on, or if the
    corrected state or covariance is not finite.
    */
    bool a_posteriori_step() {
        if (!innovation_valid) {
            return false;
        }

        /* No observations means prediction only; Eigen's QR requires a
        nonempty matrix. This also supports the split-step interface. */
        if (innovation.size() == 0) {
            return true;
        }

        /*
        Calculate the cross-correlation matrix described in equations 70 and
        71 from from the Kraft paper.
        */
        CrossCorrelation cross_correlation(StateVectorType::covariance_size(), innovation.size());
        Detail::calculate_cross_correlation<StateVectorType>(cross_correlation, w_prime, z_prime);

        /*
        Calculate the Kalman gain as described in equation 72 from the Kraft
        paper, by solving K * innovation_covariance = cross_correlation with a
        column-pivoting QR decomposition rather than forming an inverse. The
        rank-revealing decomposition allows the innovation covariance to be
        positive semidefinite.
        */
        CrossCorrelation kalman_gain = innovation_covariance.colPivHouseholderQr().solve(
            cross_correlation.transpose()).transpose();

        /*
        Calculate the update delta vector, to be applied to the a priori
        estimate.
        */
        typename StateVectorType::StateVectorDelta update_delta = kalman_gain * innovation;

        /* Apply the update delta to a copy of the state vector. */
        StateVectorType updated_state = state;
        updated_state.apply_delta(update_delta);

        /* Update the covariance using equation 75 from the Kraft paper. */
        typename StateVectorType::CovarianceMatrix updated_covariance = covariance;
        updated_covariance.noalias() -= (kalman_gain * innovation_covariance * kalman_gain.transpose());
        if (!updated_state.allFinite() || !updated_covariance.allFinite()) {
            return false;
        }

        state = updated_state;
        covariance = updated_covariance;
        return true;
    }
};

/*
The class implements the square-root form of the UKF, described in the paper
"The Square-Root Unscented Kalman Filter for State and Parameter-Estimation"
by Rudolph van der Merwe and Eric A. Wan, retrieved from:
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.80.1421&rep=rep1&type=pdf

The SR-UKF has a number of numerical and stability advantages over the
standard UKF, and a lower computational cost. Steps report failure in the
same way as Core.
*/
template <typename StateVectorType, typename MeasurementVectorType, typename IntegratorType=NotUsedIntegrator>
class SquareRootCore {
protected:
    typename StateVectorType::SigmaPointDistribution sigma_points;
    typename StateVectorType::SigmaPointDeltas w_prime;
    typename MeasurementVectorType::template SigmaPointDeltas<StateVectorType> z_prime;

    /* Whether the sigma points and innovation are ready for the next step. */
    bool prior_valid = false;
    bool innovation_valid = false;

    /* Aliases needed during filter iteration. */
    using CrossCorrelation = Eigen::Matrix<
        real_t,
        StateVectorType::covariance_size(),
        MeasurementVectorType::RowsAtCompileTime,
        StateVectorType::covariance_size() == 1 ? Eigen::RowMajor : Eigen::ColMajor,
        StateVectorType::covariance_size(),
        MeasurementVectorType::MaxRowsAtCompileTime>;

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
    Only diagonal measurement noise root covariance is supported at the
    moment.
    */
    typename MeasurementVectorType::CovarianceVector measurement_root_covariance;

    /*
    Innovation covariance is also stored in square-root form as described
    above.
    */
    MeasurementVectorType innovation;
    typename MeasurementVectorType::CovarianceMatrix innovation_root_covariance;

    /* Top-level function used to carry out a filter step. */
    template <typename... U>
    bool step(real_t delta, const MeasurementVectorType& z, const U&... input) {
        return a_priori_step(delta, input...) && innovation_step(z, input...) && a_posteriori_step();
    }

    /*
    The a priori step differs from the standard UKF in that it doesn't
    compute the LLT decomposition of the state covariance matrix.

    In addition, the a priori state root covariance is recovered from the
    sigma point distribution using a QR decomposition and a rank-one Cholesky
    update. That update is a downdate when the central weight is negative
    (small alpha); if it cannot be applied, the factor is rebuilt from the
    covariance, and the step fails without changing the filter if the
    predicted covariance is indefinite.
    */
    template <typename... U>
    bool a_priori_step(real_t delta, const U&... input) {
        prior_valid = false;
        innovation_valid = false;

        const real_t sigma_scale = std::sqrt(Parameters::ScaledDimension<StateVectorType>);
        sigma_points = state.calculate_sigma_point_distribution(root_covariance * sigma_scale);

        /* Propagate the sigma points through the process model. */
        Detail::propagate_sigma_points<StateVectorType, IntegratorType>(sigma_points, delta, input...);

        /* Calculate the a priori estimate mean, deltas and root covariance. */
        const StateVectorType predicted_state = StateVectorType::calculate_sigma_point_mean(sigma_points);
        w_prime = predicted_state.calculate_sigma_point_deltas(sigma_points);

        /*
        Create an augmented matrix containing all but the centre sigma point
        delta, with the process noise root covariance to the right.
        */
        using AugmentedSigmaPointDeltas = Matrix<
            StateVectorType::covariance_size(),
            StateVectorType::num_sigma() - 1 + StateVectorType::covariance_size()>;

        AugmentedSigmaPointDeltas augmented_w_prime;
        augmented_w_prime <<
            std::sqrt(Parameters::Sigma_WCI<StateVectorType>) * w_prime.rightCols(StateVectorType::num_sigma() - 1),
            process_noise_root_covariance;

        /*
        Calculate the QR decomposition of the augmented sigma point deltas,
        and transpose its triangular factor to a lower root.
        */
        typename StateVectorType::CovarianceMatrix predicted_root =
            augmented_w_prime.transpose().householderQr().matrixQR().topLeftCorner(
                StateVectorType::covariance_size(),
                StateVectorType::covariance_size()).template triangularView<Eigen::Upper>();
        predicted_root.transposeInPlace();

        /*
        Do a rank-one Cholesky update of the root covariance matrix using the
        central sigma point delta.
        */
        if (!predicted_state.allFinite() || !Detail::update_root(
                predicted_root, w_prime.col(0), Parameters::Sigma_WC0<StateVectorType>)) {
            return false;
        }

        /* Commit the prediction only once it is known to be valid. */
        state = predicted_state;
        root_covariance = predicted_root;

        /* Recalculate the sigma points using the a priori covariance. */
        sigma_points = state.calculate_sigma_point_distribution(root_covariance * sigma_scale);
        w_prime = state.calculate_sigma_point_deltas(sigma_points);
        prior_valid = true;
        return true;
    }

    /*
    The innovation step differs from the standard UKF in the same way as the
    a priori step; it recovers the innovation root covariance using a QR
    decomposition and a rank-one Cholesky update.
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
                innovation = z;
                innovation_root_covariance.resize(0, 0);
                z_prime.resize(0, StateVectorType::num_sigma());
                innovation_valid = true;
                return true;
            }
        }

        /* Propagate the sigma points through the measurement model. */
        typename MeasurementVectorType::template SigmaPointDistribution<StateVectorType>
            measurement_sigma_points = z.template calculate_sigma_point_distribution<StateVectorType>(
                sigma_points, input...);

        /* Calculate the measurement prediction mean and deltas. */
        MeasurementVectorType z_pred =
            z.template calculate_sigma_point_mean<StateVectorType>(measurement_sigma_points);
        z_prime = z_pred.template calculate_sigma_point_deltas<StateVectorType>(measurement_sigma_points);

        /* Calculate innovation and innovation root covariance. */
        const MeasurementVectorType predicted_innovation = z_pred.calculate_innovation(z);

        /*
        Create an augmented matrix containing all but the centre innovation
        delta, with the measurement root covariance to the right.
        */
        constexpr std::ptrdiff_t augmented_innovation_cols =
            (MeasurementVectorType::RowsAtCompileTime != Eigen::Dynamic) ?
            StateVectorType::num_sigma() - 1 + MeasurementVectorType::RowsAtCompileTime : Eigen::Dynamic;

        using AugmentedInnovationDeltas = Eigen::Matrix<
            real_t,
            MeasurementVectorType::RowsAtCompileTime,
            augmented_innovation_cols,
            MeasurementVectorType::MaxRowsAtCompileTime == 1 ? Eigen::RowMajor : Eigen::ColMajor,
            MeasurementVectorType::MaxRowsAtCompileTime,
            StateVectorType::num_sigma() - 1 + MeasurementVectorType::MaxRowsAtCompileTime>;

        AugmentedInnovationDeltas augmented_z_prime(z.size(), StateVectorType::num_sigma() - 1 + z.size());
        augmented_z_prime.block(0, 0, z.size(), StateVectorType::num_sigma() - 1) =
            std::sqrt(Parameters::Sigma_WCI<StateVectorType>) * z_prime.rightCols(StateVectorType::num_sigma() - 1);
        augmented_z_prime.block(0, StateVectorType::num_sigma() - 1, z.size(), z.size()) =
            z.calculate_measurement_root_covariance(measurement_root_covariance, z_pred);

        /*
        Calculate the QR decomposition of the augmented innovation deltas,
        and transpose its triangular factor to a lower root.
        */
        typename MeasurementVectorType::CovarianceMatrix predicted_root =
            augmented_z_prime.transpose().householderQr().matrixQR().topLeftCorner(
                z.size(), z.size()).template triangularView<Eigen::Upper>();
        predicted_root.transposeInPlace();

        /*
        Do a rank-one Cholesky update of the innovation root covariance
        matrix using the central innovation delta.
        */
        if (!predicted_innovation.allFinite() || !Detail::update_root(
                predicted_root, z_prime.col(0), Parameters::Sigma_WC0<StateVectorType>)) {
            return false;
        }

        innovation = predicted_innovation;
        innovation_root_covariance = predicted_root;
        innovation_valid = true;
        return true;
    }

    /*
    Compared to the standard UKF, the a posteriori step calculates the Kalman
    gain using triangular solves, and then updates the root covariance using
    a series of rank-one Cholesky downdates. Nothing is changed unless every
    update succeeds.
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
        CrossCorrelation cross_correlation(StateVectorType::covariance_size(), innovation.size());
        Detail::calculate_cross_correlation<StateVectorType>(cross_correlation, w_prime, z_prime);

        /* Calculate the Kalman gain from the innovation root covariance. */
        CrossCorrelation kalman_gain(StateVectorType::covariance_size(), innovation.size());
        Detail::square_root_gain(kalman_gain, cross_correlation, innovation_root_covariance);

        /*
        Calculate the update delta vector, to be applied to the a priori
        estimate.
        */
        typename StateVectorType::StateVectorDelta update_delta = kalman_gain * innovation;

        /*
        Calculate the Cholesky update matrix. Reuse the cross-correlation
        variable, since we don't need it again.
        */
        cross_correlation.noalias() = kalman_gain * innovation_root_covariance;

        /*
        Update the root covariance using a series of rank-one Cholesky
        downdates. LLT downdates require strict positivity; an exact
        observation can instead leave a semidefinite covariance, which is
        rebuilt from the original factor rather than the partial downdate.
        */
        typename StateVectorType::CovarianceMatrix updated_root = root_covariance;
        StateVectorType updated_state = state;
        updated_state.apply_delta(update_delta);
        if (!updated_state.allFinite() || !Detail::update_root(updated_root, cross_correlation, real_t(-1.0))) {
            return false;
        }

        /* Commit the correction only after every covariance update succeeds. */
        state = updated_state;
        root_covariance = updated_root;
        return true;
    }
};

/*
The class implements the optimised form of the square-root UKF for parameter
estimation, described in the paper
"The Square-Root Unscented Kalman Filter for State and Parameter-Estimation"
by Rudolph van der Merwe and Eric A. Wan, retrieved from:
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.80.1421&rep=rep1&type=pdf

Steps report failure in the same way as Core.
*/
template <typename StateVectorType, typename MeasurementVectorType>
class SquareRootParameterEstimationCore {
protected:
    typename StateVectorType::SigmaPointDistribution sigma_points;
    typename StateVectorType::SigmaPointDeltas w_prime;
    typename MeasurementVectorType::template SigmaPointDeltas<StateVectorType> z_prime;

    /* Whether the sigma points and innovation are ready for the next step. */
    bool prior_valid = false;
    bool innovation_valid = false;

    /* Aliases needed during filter iteration. */
    using CrossCorrelation = Eigen::Matrix<
        real_t,
        StateVectorType::covariance_size(),
        MeasurementVectorType::RowsAtCompileTime,
        StateVectorType::covariance_size() == 1 ? Eigen::RowMajor : Eigen::ColMajor,
        StateVectorType::covariance_size(),
        MeasurementVectorType::MaxRowsAtCompileTime>;

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
    Only diagonal measurement noise root covariance is supported at the
    moment.
    */
    typename MeasurementVectorType::CovarianceVector measurement_root_covariance;

    /*
    Innovation covariance is also stored in square-root form as described
    above.
    */
    MeasurementVectorType innovation;
    typename MeasurementVectorType::CovarianceMatrix innovation_root_covariance;

    /* Top-level function used to carry out a filter step. */
    template <typename... U>
    bool step(const MeasurementVectorType& z, const U&... input) {
        return a_priori_step() && innovation_step(z, input...) && a_posteriori_step();
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
        Recover a factor of S*S' + Q*Q' from the augmented matrix [S Q].
        Updating only S's diagonal would also change the off-diagonal
        covariance whenever the parameters are correlated.
        */
        Matrix<StateVectorType::covariance_size(), 2 * StateVectorType::covariance_size()> augmented;
        augmented << root_covariance, process_noise_root_covariance;
        typename StateVectorType::CovarianceMatrix predicted_root =
            augmented.transpose().householderQr().matrixQR().topLeftCorner(
                StateVectorType::covariance_size(), StateVectorType::covariance_size())
            .template triangularView<Eigen::Upper>();
        predicted_root.transposeInPlace();
        if (!predicted_root.allFinite() || !state.allFinite()) {
            return false;
        }
        root_covariance = predicted_root;

        sigma_points = state.calculate_sigma_point_distribution(root_covariance *
            std::sqrt(Parameters::ScaledDimension<StateVectorType>));

        /* Calculate the sigma point deltas. */
        w_prime = state.calculate_sigma_point_deltas(sigma_points);
        prior_valid = true;
        return true;
    }

    /*
    The innovation step differs from the standard UKF in the same way as the
    a priori step; it recovers the innovation root covariance using a QR
    decomposition and a rank-one Cholesky update.
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
                innovation = z;
                innovation_root_covariance.resize(0, 0);
                z_prime.resize(0, StateVectorType::num_sigma());
                innovation_valid = true;
                return true;
            }
        }

        /* Propagate the sigma points through the measurement model. */
        typename MeasurementVectorType::template SigmaPointDistribution<StateVectorType>
            measurement_sigma_points = z.template calculate_sigma_point_distribution<StateVectorType>(
                sigma_points, input...);

        /* Calculate the measurement prediction mean and deltas. */
        MeasurementVectorType z_pred =
            z.template calculate_sigma_point_mean<StateVectorType>(measurement_sigma_points);
        z_prime = z_pred.template calculate_sigma_point_deltas<StateVectorType>(measurement_sigma_points);

        /* Calculate innovation and innovation root covariance. */
        const MeasurementVectorType predicted_innovation = z_pred.calculate_innovation(z);

        /*
        Create an augmented matrix containing all but the centre innovation
        delta, with the measurement root covariance to the right.
        */
        constexpr std::ptrdiff_t augmented_innovation_cols =
            (MeasurementVectorType::RowsAtCompileTime != Eigen::Dynamic) ?
            StateVectorType::num_sigma() - 1 + MeasurementVectorType::RowsAtCompileTime : Eigen::Dynamic;

        using AugmentedInnovationDeltas = Eigen::Matrix<
            real_t,
            MeasurementVectorType::RowsAtCompileTime,
            augmented_innovation_cols,
            MeasurementVectorType::MaxRowsAtCompileTime == 1 ? Eigen::RowMajor : Eigen::ColMajor,
            MeasurementVectorType::MaxRowsAtCompileTime,
            StateVectorType::num_sigma() - 1 + MeasurementVectorType::MaxRowsAtCompileTime>;

        AugmentedInnovationDeltas augmented_z_prime(z.size(), StateVectorType::num_sigma() - 1 + z.size());
        augmented_z_prime.block(0, 0, z.size(), StateVectorType::num_sigma() - 1) =
            std::sqrt(Parameters::Sigma_WCI<StateVectorType>) * z_prime.rightCols(StateVectorType::num_sigma() - 1);
        augmented_z_prime.block(0, StateVectorType::num_sigma() - 1, z.size(), z.size()) =
            z.calculate_measurement_root_covariance(measurement_root_covariance, z_pred);

        /*
        Calculate the QR decomposition of the augmented innovation deltas,
        and transpose its triangular factor to a lower root.
        */
        typename MeasurementVectorType::CovarianceMatrix predicted_root =
            augmented_z_prime.transpose().householderQr().matrixQR().topLeftCorner(
                z.size(), z.size()).template triangularView<Eigen::Upper>();
        predicted_root.transposeInPlace();

        /*
        Do a rank-one Cholesky update of the innovation root covariance
        matrix using the central innovation delta.
        */
        if (!predicted_innovation.allFinite() || !Detail::update_root(
                predicted_root, z_prime.col(0), Parameters::Sigma_WC0<StateVectorType>)) {
            return false;
        }

        innovation = predicted_innovation;
        innovation_root_covariance = predicted_root;
        innovation_valid = true;
        return true;
    }

    /*
    Compared to the standard UKF, the a posteriori step calculates the Kalman
    gain using triangular solves, and then updates the root covariance using
    a series of rank-one Cholesky downdates. Nothing is changed unless every
    update succeeds.
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
        CrossCorrelation cross_correlation(StateVectorType::covariance_size(), innovation.size());
        Detail::calculate_cross_correlation<StateVectorType>(cross_correlation, w_prime, z_prime);

        /* Calculate the Kalman gain from the innovation root covariance. */
        CrossCorrelation kalman_gain(StateVectorType::covariance_size(), innovation.size());
        Detail::square_root_gain(kalman_gain, cross_correlation, innovation_root_covariance);

        /*
        Calculate the update delta vector, to be applied to the a priori
        estimate.
        */
        typename StateVectorType::StateVectorDelta update_delta = kalman_gain * innovation;

        /*
        Calculate the Cholesky update matrix. Reuse the cross-correlation
        variable, since we don't need it again.
        */
        cross_correlation.noalias() = kalman_gain * innovation_root_covariance;

        /*
        Update the root covariance using a series of rank-one Cholesky
        downdates, rebuilding it from the original factor if the posterior
        is only semidefinite.
        */
        typename StateVectorType::CovarianceMatrix updated_root = root_covariance;
        StateVectorType updated_state = state;
        updated_state.apply_delta(update_delta);
        if (!updated_state.allFinite() || !Detail::update_root(updated_root, cross_correlation, real_t(-1.0))) {
            return false;
        }

        /* Commit the correction only after every covariance update succeeds. */
        state = updated_state;
        root_covariance = updated_root;
        return true;
    }
};

}

#endif
