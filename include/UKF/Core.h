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
#include "UKF/Types.h"
#include "UKF/StateVector.h"

namespace UKF {

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
*/
template <typename StateVectorType, typename MeasurementVectorType, typename IntegratorType>
class Core {
protected:
    typename StateVectorType::SigmaPointDistribution sigma_points;
    typename StateVectorType::SigmaPointDeltas w_prime;
    typename MeasurementVectorType::template SigmaPointDeltas<StateVectorType> z_prime;

    /* Aliases needed during filter iteration. */
    using CrossCorrelation = Eigen::Matrix<
        real_t,
        StateVectorType::covariance_size(),
        MeasurementVectorType::RowsAtCompileTime,
        0,
        StateVectorType::covariance_size(),
        MeasurementVectorType::MaxRowsAtCompileTime>;

public:
    /* State and covariance are public for simplicity of initialisation. */
    StateVectorType state;
    typename StateVectorType::CovarianceMatrix covariance;
    typename StateVectorType::CovarianceMatrix process_noise_covariance;

    /*
    Innovation and innovation covariance are public to ease implementation of
    filter health monitoring.
    */
    MeasurementVectorType innovation;
    typename MeasurementVectorType::CovarianceMatrix innovation_covariance;

    /* Top-level function used to carry out a filter step. */
    template <typename... U>
    void step(real_t delta, const MeasurementVectorType& z, const U&... input) {
        a_priori_step(delta, input...);
        innovation_step(z, input...);
        a_posteriori_step();
    }

    /*
    The a priori step calculates the sigma point distribution using the
    previous state covariance and previous state estimate, and then
    propagates all sigma points through the process model.
    From the transformed sigma point distribution, the a priori mean and
    covariance are calculated.
    */
    template <typename... U>
    void a_priori_step(real_t delta, const U&... input) {
        /*
        Add process noise covariance to the state covariance and calculate
        the LLT decomposition of the scaled covariance matrix.
        */
        assert(((covariance * (StateVectorType::covariance_size() +
            Parameters::Lambda<StateVectorType>)).llt().info() == Eigen::Success)
            && "Covariance matrix is not positive definite");

        sigma_points = state.calculate_sigma_point_distribution(((covariance + process_noise_covariance) *
            (StateVectorType::covariance_size() + Parameters::Lambda<StateVectorType>)).llt().matrixL());

        /* Propagate the sigma points through the process model. */
        for(std::size_t i = 0; i < StateVectorType::num_sigma(); i++) {
            sigma_points.col(i) <<
                StateVectorType(sigma_points.col(i)).template process_model<IntegratorType>(delta, input...);
        }

        /* Calculate the a priori estimate mean, deltas and covariance. */
        state = StateVectorType::calculate_sigma_point_mean(sigma_points);
        w_prime = state.calculate_sigma_point_deltas(sigma_points);
        covariance = StateVectorType::calculate_sigma_point_covariance(w_prime);
    }

    /*
    In the innovation step, the a priori sigma point distribution is further
    propagated using the measurement model. For the measurement sigma point
    distribution, the estimated measurement mean and covariance are
    calculated. A measurement vector is then supplied, and the innovation
    and innovation covariance are calculated.
    */
    template <typename... U>
    void innovation_step(const MeasurementVectorType& z, const U&... input) {
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
        innovation = z_pred.template calculate_innovation(z);
        innovation_covariance = z_pred.template calculate_sigma_point_covariance<StateVectorType>(z_prime);
        innovation_covariance += z.template calculate_measurement_covariance(z_pred);
    }

    /*
    In the a posteriori step, the cross-correlation matrix and innovation
    covariance inverse are calculated and used to calculate the Kalman gain.
    This step is split from the innovation step to allow the user to carry
    out innovation consistency checks and identify failed sensors, if
    desired.
    */
    void a_posteriori_step() {
        /*
        Calculate the cross-correlation matrix described in equations 70 and
        71 from from the Kraft paper.
        */
        CrossCorrelation cross_correlation = CrossCorrelation::Zero(
            StateVectorType::covariance_size(), innovation.size());
        for(std::size_t i = 1; i < StateVectorType::num_sigma(); i++) {
            cross_correlation.noalias() += Parameters::Sigma_WCI<StateVectorType> *
                (w_prime.col(i) * z_prime.col(i).transpose());
        }
        cross_correlation.noalias() += Parameters::Sigma_WC0<StateVectorType> *
            (w_prime.col(0) * z_prime.col(0).transpose());

        /*
        Calculate the Kalman gain as described in equation 72 from the Kraft
        paper, with a slight change in that we're using a column-pivoting QR
        decomposition to calculate the Moore-Penrose pseudoinverse rather than
        the inverse. This allows the innovation covariance to be positive semi-
        definite.
        */
        CrossCorrelation kalman_gain = cross_correlation * innovation_covariance.colPivHouseholderQr().solve(
            MeasurementVectorType::CovarianceMatrix::Identity(innovation.size(), innovation.size()));

        /*
        Calculate the update delta vector, to be applied to the a priori
        estimate.
        */
        typename StateVectorType::StateVectorDelta update_delta = kalman_gain * innovation;

        /* Apply the update delta to the state vector. */
        state.apply_delta(update_delta);

        /* Update the covariance using equation 75 from the Kraft paper. */
        covariance.noalias() -= (kalman_gain * innovation_covariance * kalman_gain.transpose());
    }
};

/*
The class implements the square-root form of the UKF, described in the paper
"The Square-Root Unscented Kalman Filter for State and Parameter-Estimation"
by Rudolph van der Merwe and Eric A. Wan, retrieved from:
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.80.1421&rep=rep1&type=pdf

The SR-UKF has a number of numerical and stability advantages over the
standard UKF, and a lower computational cost.
*/
template <typename StateVectorType, typename MeasurementVectorType, typename IntegratorType>
class SquareRootCore {
protected:
    typename StateVectorType::SigmaPointDistribution sigma_points;
    typename StateVectorType::SigmaPointDeltas w_prime;
    typename MeasurementVectorType::template SigmaPointDeltas<StateVectorType> z_prime;

    /* Aliases needed during filter iteration. */
    using CrossCorrelation = Eigen::Matrix<
        real_t,
        StateVectorType::covariance_size(),
        MeasurementVectorType::RowsAtCompileTime,
        0,
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
    Innovation covariance is also stored in square-root form as described
    above.
    */
    MeasurementVectorType innovation;
    typename MeasurementVectorType::CovarianceMatrix innovation_root_covariance;

    /* Top-level function used to carry out a filter step. */
    template <typename... U>
    void step(real_t delta, const MeasurementVectorType& z, const U&... input) {
        a_priori_step(delta, input...);
        innovation_step(z, input...);
        a_posteriori_step();
    }

    /*
    The a priori step differs from the standard UKF in that it doesn't
    compute the LLT decomposition of the state covariance matrix.

    In addition, the a priori state root covariance is recovered from the
    sigma point distribution using a QR decomposition and a rank-one Cholesky
    update.
    */
    template <typename... U>
    void a_priori_step(real_t delta, const U&... input) {
        sigma_points = state.calculate_sigma_point_distribution(root_covariance *
            std::sqrt(StateVectorType::covariance_size() + Parameters::Lambda<StateVectorType>));

        /* Propagate the sigma points through the process model. */
        for(std::size_t i = 0; i < StateVectorType::num_sigma(); i++) {
            sigma_points.col(i) <<
                StateVectorType(sigma_points.col(i)).template process_model<IntegratorType>(delta, input...);
        }

        /* Calculate the a priori estimate mean, deltas and root covariance. */
        state = StateVectorType::calculate_sigma_point_mean(sigma_points);
        w_prime = state.calculate_sigma_point_deltas(sigma_points);

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
        Calculate the QR decomposition of the augmented sigma point deltas.
        */
        root_covariance = augmented_w_prime.transpose().householderQr().matrixQR().topLeftCorner(
            StateVectorType::covariance_size(),
            StateVectorType::covariance_size()).template triangularView<Eigen::Upper>();

        /*
        Do a rank-one Cholesky update of the root covariance matrix using the
        central sigma point delta.
        */
        Eigen::internal::llt_inplace<real_t, Eigen::Upper>::rankUpdate(
            root_covariance, w_prime.col(0), Parameters::Sigma_WC0<StateVectorType>);
        root_covariance.transposeInPlace();

        /* Recalculate the sigma points using the a priori covariance. */
        sigma_points = state.calculate_sigma_point_distribution(root_covariance *
            std::sqrt(StateVectorType::covariance_size() + Parameters::Lambda<StateVectorType>));
        w_prime = state.calculate_sigma_point_deltas(sigma_points);
    }

    /*
    The innovation step differs from the standard UKF in the same way as the
    a priori step; it recovers the innovation root covariance using a QR
    decomposition and a rank-one Cholesky update.
    */
    template <typename... U>
    void innovation_step(const MeasurementVectorType& z, const U&... input) {
        /* Propagate the sigma points through the measurement model. */
        typename MeasurementVectorType::template SigmaPointDistribution<StateVectorType>
            measurement_sigma_points = z.template calculate_sigma_point_distribution<StateVectorType>(
                sigma_points, input...);

        /* Calculate the measurement prediction mean and deltas. */
        MeasurementVectorType z_pred =
            z.template calculate_sigma_point_mean<StateVectorType>(measurement_sigma_points);
        z_prime = z_pred.template calculate_sigma_point_deltas<StateVectorType>(measurement_sigma_points);

        /* Calculate innovation and innovation root covariance. */
        innovation = z_pred.template calculate_innovation(z);

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
            0,
            MeasurementVectorType::MaxRowsAtCompileTime,
            StateVectorType::num_sigma() - 1 + MeasurementVectorType::MaxRowsAtCompileTime>;

        AugmentedInnovationDeltas augmented_z_prime(z.size(), StateVectorType::num_sigma() - 1 + z.size());
        augmented_z_prime.block(0, 0, z.size(), StateVectorType::num_sigma() - 1) =
            std::sqrt(Parameters::Sigma_WCI<StateVectorType>) * z_prime.rightCols(StateVectorType::num_sigma() - 1);
        augmented_z_prime.block(0, StateVectorType::num_sigma() - 1, z.size(), z.size()) =
            z.template calculate_measurement_root_covariance(z_pred);

        /*
        Calculate the QR decomposition of the augmented innovation deltas.
        */
        innovation_root_covariance = augmented_z_prime.transpose().householderQr().matrixQR().topLeftCorner(
            z.size(), z.size()).template triangularView<Eigen::Upper>();

        /*
        Do a rank-one Cholesky update of the innovation root covariance
        matrix using the central innovation delta.
        */
        Eigen::internal::llt_inplace<real_t, Eigen::Upper>::rankUpdate(
            innovation_root_covariance, z_prime.col(0), Parameters::Sigma_WC0<StateVectorType>);
        innovation_root_covariance.transposeInPlace();
    }

    /*
    Compared to the standard UKF, the a posteriori step calculates the Kalman
    gain using a QR decomposition, and then updates the root covariance using
    a series of rank-one Cholesky updates.
    */
    void a_posteriori_step() {
        /*
        Calculate the cross-correlation matrix described in equations 70 and
        71 from from the Kraft paper.
        */
        CrossCorrelation cross_correlation = CrossCorrelation::Zero(
            StateVectorType::covariance_size(), innovation.size());
        for(std::size_t i = 1; i < StateVectorType::num_sigma(); i++) {
            cross_correlation.noalias() += Parameters::Sigma_WCI<StateVectorType> *
                (w_prime.col(i) * z_prime.col(i).transpose());
        }
        cross_correlation.noalias() += Parameters::Sigma_WC0<StateVectorType> *
            (w_prime.col(0) * z_prime.col(0).transpose());

        /*
        Calculate the Kalman gain using QR decomposition. This expression
        implements (S’\(S\P’))’, which is equivalent to (P/S’)/S given in
        literature. Eigen's QR decomposition implements a left-division,
        rather than the right-division assumed in the literature.
        */
        CrossCorrelation kalman_gain = innovation_root_covariance.transpose().fullPivHouseholderQr().solve(
            innovation_root_covariance.fullPivHouseholderQr().solve(
                cross_correlation.transpose())).transpose();

        /*
        Calculate the update delta vector, to be applied to the a priori
        estimate.
        */
        typename StateVectorType::StateVectorDelta update_delta = kalman_gain * innovation;

        /* Apply the update delta to the state vector. */
        state.apply_delta(update_delta);

        /*
        Calculate the Cholesky update matrix. Reuse the cross-correlation
        variable, since we don't need it again.
        */
        cross_correlation.noalias() = kalman_gain * innovation_root_covariance;

        /*
        Update the root covariance using a series of rank-one Cholesky
        downdates.
        */
        for(std::ptrdiff_t i = 0; i < cross_correlation.cols(); i++) {
            Eigen::internal::llt_inplace<real_t, Eigen::Lower>::rankUpdate(
                root_covariance, cross_correlation.col(i), real_t(-1.0));
        }
    }
};

/*
The class implements the optimised form of the square-root UKF for parameter
estimation, described in the paper
"The Square-Root Unscented Kalman Filter for State and Parameter-Estimation"
by Rudolph van der Merwe and Eric A. Wan, retrieved from:
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.80.1421&rep=rep1&type=pdf
*/
template <typename StateVectorType, typename MeasurementVectorType>
class SquareRootParameterEstimationCore {
protected:
    typename StateVectorType::SigmaPointDistribution sigma_points;
    typename StateVectorType::SigmaPointDeltas w_prime;
    typename MeasurementVectorType::template SigmaPointDeltas<StateVectorType> z_prime;

    /* Aliases needed during filter iteration. */
    using CrossCorrelation = Eigen::Matrix<
        real_t,
        StateVectorType::covariance_size(),
        MeasurementVectorType::RowsAtCompileTime,
        0,
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
    Innovation covariance is also stored in square-root form as described
    above.
    */
    MeasurementVectorType innovation;
    typename MeasurementVectorType::CovarianceMatrix innovation_root_covariance;

    /* Top-level function used to carry out a filter step. */
    template <typename... U>
    void step(const MeasurementVectorType& z, const U&... input) {
        a_priori_step();
        innovation_step(z, input...);
        a_posteriori_step();
    }

    /*
    For the parameter estimation filter, there is no need for a process
    model, which allows the a priori step to be optimised significantly.
    There's no need to propagate the sigma points through a process model,
    and we only have to add the process noise to the root covariance.
    */
    void a_priori_step() {
        /*
        Add the process noise to the root covariance. Note that there are two
        methods for doing this; this one has been chosen as it still allows
        individual fields to have different process noise values.
        */
        typename StateVectorType::StateVectorDelta R_diag = process_noise_root_covariance.diagonal();
        typename StateVectorType::StateVectorDelta S_diag = root_covariance.diagonal();
        root_covariance.diagonal() = (S_diag.cwiseProduct(S_diag) + R_diag.cwiseProduct(R_diag)).cwiseSqrt();

        sigma_points = state.calculate_sigma_point_distribution(root_covariance *
            std::sqrt(StateVectorType::covariance_size() + Parameters::Lambda<StateVectorType>));

        /* Calculate the sigma point deltas. */
        w_prime = state.calculate_sigma_point_deltas(sigma_points);
    }

    /*
    The innovation step differs from the standard UKF in the same way as the
    a priori step; it recovers the innovation root covariance using a QR
    decomposition and a rank-one Cholesky update.
    */
    template <typename... U>
    void innovation_step(const MeasurementVectorType& z, const U&... input) {
        /* Propagate the sigma points through the measurement model. */
        typename MeasurementVectorType::template SigmaPointDistribution<StateVectorType>
            measurement_sigma_points = z.template calculate_sigma_point_distribution<StateVectorType>(
                sigma_points, input...);

        /* Calculate the measurement prediction mean and deltas. */
        MeasurementVectorType z_pred =
            z.template calculate_sigma_point_mean<StateVectorType>(measurement_sigma_points);
        z_prime = z_pred.template calculate_sigma_point_deltas<StateVectorType>(measurement_sigma_points);

        /* Calculate innovation and innovation root covariance. */
        innovation = z_pred.template calculate_innovation(z);

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
            0,
            MeasurementVectorType::MaxRowsAtCompileTime,
            StateVectorType::num_sigma() - 1 + MeasurementVectorType::MaxRowsAtCompileTime>;

        AugmentedInnovationDeltas augmented_z_prime(z.size(), StateVectorType::num_sigma() - 1 + z.size());
        augmented_z_prime.block(0, 0, z.size(), StateVectorType::num_sigma() - 1) =
            std::sqrt(Parameters::Sigma_WCI<StateVectorType>) * z_prime.rightCols(StateVectorType::num_sigma() - 1);
        augmented_z_prime.block(0, StateVectorType::num_sigma() - 1, z.size(), z.size()) =
            z.template calculate_measurement_root_covariance(z_pred);

        /*
        Calculate the QR decomposition of the augmented innovation deltas.
        */
        innovation_root_covariance = augmented_z_prime.transpose().householderQr().matrixQR().topLeftCorner(
            z.size(), z.size()).template triangularView<Eigen::Upper>();

        /*
        Do a rank-one Cholesky update of the innovation root covariance
        matrix using the central innovation delta.
        */
        Eigen::internal::llt_inplace<real_t, Eigen::Upper>::rankUpdate(
            innovation_root_covariance, z_prime.col(0), Parameters::Sigma_WC0<StateVectorType>);
        innovation_root_covariance.transposeInPlace();
    }

    /*
    Compared to the standard UKF, the a posteriori step calculates the Kalman
    gain using a QR decomposition, and then updates the root covariance using
    a series of rank-one Cholesky updates.
    */
    void a_posteriori_step() {
        /*
        Calculate the cross-correlation matrix described in equations 70 and
        71 from from the Kraft paper.
        */
        CrossCorrelation cross_correlation = CrossCorrelation::Zero(
            StateVectorType::covariance_size(), innovation.size());
        for(std::size_t i = 1; i < StateVectorType::num_sigma(); i++) {
            cross_correlation.noalias() += Parameters::Sigma_WCI<StateVectorType> *
                (w_prime.col(i) * z_prime.col(i).transpose());
        }
        cross_correlation.noalias() += Parameters::Sigma_WC0<StateVectorType> *
            (w_prime.col(0) * z_prime.col(0).transpose());

        /*
        Calculate the Kalman gain using QR decomposition. This expression
        implements (S’\(S\P’))’, which is equivalent to (P/S’)/S given in
        literature. Eigen's QR decomposition implements a left-division,
        rather than the right-division assumed in the literature.
        */
        CrossCorrelation kalman_gain = innovation_root_covariance.transpose().fullPivHouseholderQr().solve(
            innovation_root_covariance.fullPivHouseholderQr().solve(
                cross_correlation.transpose())).transpose();

        /*
        Calculate the update delta vector, to be applied to the a priori
        estimate.
        */
        typename StateVectorType::StateVectorDelta update_delta = kalman_gain * innovation;

        /* Apply the update delta to the state vector. */
        state.apply_delta(update_delta);

        /*
        Calculate the Cholesky update matrix. Reuse the cross-correlation
        variable, since we don't need it again.
        */
        cross_correlation.noalias() = kalman_gain * innovation_root_covariance;

        /*
        Update the root covariance using a series of rank-one Cholesky
        downdates.
        */
        for(std::ptrdiff_t i = 0; i < cross_correlation.cols(); i++) {
            Eigen::internal::llt_inplace<real_t, Eigen::Lower>::rankUpdate(
                root_covariance, cross_correlation.col(i), real_t(-1.0));
        }
    }
};

}

#endif
