// Independent finite-difference checks of the robust estimator's covariance reset.
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>
#include "robust_geometry.h"
#include "robust_correction.h"
#include "UKF/StateVector.h"

namespace {
// Checks report the first failure through this flag; no exceptions are used.
const char* failure = nullptr;

bool check(bool condition, const char* message) {
    if (!condition && !failure) failure = message;
    return condition;
}

using Matrix = Eigen::Matrix<double, 3, 3>;
using Vector = Eigen::Matrix<double, 3, 1>;
using Quaternion = Eigen::Quaternion<double>;

Quaternion from_gibbs(const Vector& r) {
    return Quaternion(1.0, r.x() / 2, r.y() / 2, r.z() / 2).normalized();
}

Vector recentered_error(const Vector& r, const Quaternion& correction) {
    const Quaternion relative = from_gibbs(r) * correction.conjugate();
    return 2 * relative.vec() / relative.w();
}

Matrix finite_difference_reset(const Quaternion& correction) {
    const Vector center = 2 * correction.vec() / correction.w();
    constexpr double epsilon = 1e-6;
    Matrix result;
    for (int axis = 0; axis < 3; ++axis) {
        Vector plus = center, minus = center;
        plus[axis] += epsilon;
        minus[axis] -= epsilon;
        result.col(axis) = (recentered_error(plus, correction)
                          - recentered_error(minus, correction)) / (2 * epsilon);
    }
    return result;
}

void close(const Matrix& a, const Matrix& b, const char* message) {
    check(a.allFinite() && (a - b).norm() <= 2e-8 * b.norm() + 1e-14, message);
}

void check_reset(double degrees) {
    const Quaternion predicted = Quaternion(Eigen::AngleAxis<double>(0.37, Vector::UnitX()))
        * Quaternion(Eigen::AngleAxis<double>(-0.81, Vector::UnitY()));
    const Quaternion correction(Eigen::AngleAxis<double>(degrees * std::acos(-1.0) / 180,
                                                         Vector(.3, -.4, .5).normalized()));
    const Quaternion corrected = correction * predicted;
    Matrix old_root;
    old_root << .1, 0, 0, .02, .05, 0, -.03, .01, .08;
    const Matrix covariance = old_root * old_root.transpose();
    Matrix root;
    if (!check(robust_detail::recenter_attitude_root(predicted, corrected, old_root, root),
               "reset rejected a valid correction")) return;
    const Matrix J = finite_difference_reset(correction);
    close(root * root.transpose(), J * covariance * J.transpose(),
          "reset covariance disagrees with finite difference of nonlinear error coordinates");
    check(root.template triangularView<Eigen::StrictlyUpper>().toDenseMatrix().norm() <= 1e-14,
          "reset did not return a lower triangular covariance factor");

    Quaternion opposite_predicted = predicted;
    opposite_predicted.coeffs() *= -1;
    Matrix opposite;
    if (!check(robust_detail::recenter_attitude_root(opposite_predicted, corrected, old_root, opposite),
               "reset rejected a sign-flipped prediction")) return;
    close(opposite * opposite.transpose(), root * root.transpose(),
          "quaternion sign changes reset covariance");

    const Quaternion yaw(Eigen::AngleAxis<double>(0.73, Vector::UnitZ()));
    const Matrix rotation = yaw.toRotationMatrix();
    const Matrix yaw_covariance = rotation * covariance * rotation.transpose();
    const Matrix yaw_root = yaw_covariance.llt().matrixL();
    Matrix transformed;
    if (!check(robust_detail::recenter_attitude_root(yaw * predicted, yaw * corrected, yaw_root, transformed),
               "reset rejected a yaw-rotated correction")) return;
    close(transformed * transformed.transpose(), rotation * root * root.transpose() * rotation.transpose(),
          "world yaw frame changes the physical reset covariance");
}

void check_constrained_gravity() {
    using Measurement = Eigen::Matrix<double, 2, 1>;
    using MeasurementMatrix = Eigen::Matrix<double, 2, 2>;
    using Cross = Eigen::Matrix<double, 3, 2>;
    using TestState = UKF::StateVector<UKF::Field<701, UKF::Quaternion>>;
    Matrix prior_root;
    prior_root << .1, 0, 0, .02, .05, 0, -.03, .01, .08;
    const Matrix prior = prior_root * prior_root.transpose();
    Eigen::Matrix<double, 2, 3> H;
    H << .7, -.2, 0, .1, .8, 0;
    MeasurementMatrix noise_root;
    noise_root << .03, 0, .01, .04;
    const Cross cross = prior * H.transpose();
    const MeasurementMatrix innovation_covariance = H * prior * H.transpose()
        + noise_root * noise_root.transpose();
    const Measurement innovation(.3, -.2);
    robust_detail::GravityUpdate result;
    if (!check(robust_detail::constrained_gravity_update(prior, cross, innovation_covariance, innovation, result),
               "constrained gravity update rejected valid moments")) return;
    check(result.delta.z() == 0 && result.delta.head<2>().norm() >= .01,
          "gravity correction must update tilt but not world yaw");
    check(std::abs(result.posteriorP(2, 2) - prior(2, 2)) <= 1e-15,
          "constrained gravity update changed prior yaw variance");
    check(result.posteriorP.topLeftCorner<2, 2>().trace() < prior.topLeftCorner<2, 2>().trace(),
          "gravity correction did not reduce tilt uncertainty");
    check(Eigen::SelfAdjointEigenSolver<Matrix>(result.posteriorP).eigenvalues().minCoeff() >= -1e-13,
          "constrained gravity posterior is indefinite");
    MeasurementMatrix indefinite = innovation_covariance;
    indefinite(1, 1) = -1;
    robust_detail::GravityUpdate rejected;
    check(!robust_detail::constrained_gravity_update(prior, cross, indefinite, innovation, rejected),
          "indefinite gravity innovation covariance was accepted");

    // Recover the gain actually used from unit innovations. Independently form
    // residual errors for joint state/measurement sigma samples, then compute
    // their sample covariance. This catches the incorrect P-KSK^T shortcut.
    Cross gain;
    for (int column = 0; column < 2; ++column) {
        Measurement unit = Measurement::Zero();
        unit[column] = 1;
        robust_detail::GravityUpdate unit_result;
        if (!check(robust_detail::constrained_gravity_update(prior, cross, innovation_covariance, unit, unit_result),
                   "constrained gravity update rejected a unit innovation")) return;
        gain.col(column) = unit_result.delta;
    }
    Matrix sample_covariance = Matrix::Zero();
    for (int component = 0; component < 5; ++component) {
        for (int sign : {-1, 1}) {
            Vector state_error = Vector::Zero();
            Measurement noise = Measurement::Zero();
            if (component < 3) state_error = sign * std::sqrt(5.0) * prior_root.col(component);
            else noise = sign * std::sqrt(5.0) * noise_root.col(component - 3);
            const Vector residual = state_error - gain * (H * state_error + noise);
            sample_covariance += residual * residual.transpose() / 10;
        }
    }
    close(result.posteriorP, sample_covariance,
          "constrained covariance disagrees with actual error-sample covariance");

    const Quaternion predicted = Quaternion(Eigen::AngleAxis<double>(.6, Vector::UnitX()))
        * Quaternion(Eigen::AngleAxis<double>(-.4, Vector::UnitY()));
    TestState state;
    state.set_field<701>(predicted);
    state.apply_delta(result.delta);
    const Quaternion corrected = state.get_field<701>();
    const Quaternion correction = corrected * predicted.conjugate();
    check(std::abs(correction.z()) <= 1e-14, "actual quaternion correction contains a yaw twist");
    TestState opposite;
    Quaternion negated = predicted;
    negated.coeffs() *= -1;
    opposite.set_field<701>(negated);
    opposite.apply_delta(result.delta);
    check((opposite.get_field<701>().toRotationMatrix() - corrected.toRotationMatrix()).norm() <= 1e-13,
          "quaternion sign changes constrained correction");

    const Quaternion yaw(Eigen::AngleAxis<double>(.73, Vector::UnitZ()));
    const Matrix rotation = yaw.toRotationMatrix();
    robust_detail::GravityUpdate transformed;
    if (!check(robust_detail::constrained_gravity_update(
            (rotation * prior * rotation.transpose()).eval(), (rotation * cross).eval(),
            innovation_covariance, innovation, transformed),
            "constrained gravity update rejected yaw-rotated moments")) return;
    close(transformed.posteriorP, rotation * result.posteriorP * rotation.transpose(),
          "world yaw frame changes constrained covariance");
    check((transformed.delta - rotation * result.delta).norm() <= 1e-13,
          "world yaw frame changes constrained delta");
    TestState rotated;
    rotated.set_field<701>(yaw * predicted);
    rotated.apply_delta(transformed.delta);
    check((rotated.get_field<701>().toRotationMatrix() - (yaw * corrected).toRotationMatrix()).norm() <= 1e-13,
          "world yaw frame changes physical corrected attitude");
}
}

int main() {
    check_constrained_gravity();
    for (double correction : {0.0, 5.0, 30.0, 90.0, 150.0}) check_reset(correction);
    const Quaternion half_turn(Eigen::AngleAxis<double>(std::acos(-1.0), Vector::UnitZ()));
    Matrix unchanged = Matrix::Constant(7);
    check(!robust_detail::recenter_attitude_root(Quaternion::Identity(), half_turn, Matrix::Identity(), unchanged),
          "Gibbs chart singularity was silently accepted");
    check(unchanged == Matrix::Constant(7), "rejected reset modified its output");
    if (failure) {
        std::cerr << failure << '\n';
        return 1;
    }
    std::cout << "Constrained gravity and covariance reset: sample covariances, finite differences, signs and frame invariance passed\n";
    return 0;
}
