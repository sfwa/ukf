// MIT license; see the repository LICENSE.
// Independent error-sample and nonlinear-coordinate checks of joint covariance.
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "mav_correction.h"
#include "UKF/StateVector.h"

namespace {
constexpr int N = 9;
using Matrix = UKF::Matrix<N, N>;
using Vector = UKF::Vector<N>;
using V3 = UKF::Vector<3>;
using Q = UKF::Quaternion;
using Cross = UKF::Matrix<N, Eigen::Dynamic>;
using DynamicMatrix = UKF::Matrix<Eigen::Dynamic, Eigen::Dynamic>;
using DynamicVector = UKF::Vector<Eigen::Dynamic>;
using State = UKF::StateVector<UKF::Field<901, Q>, UKF::Field<902, V3>, UKF::Field<903, V3>>;

void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void close(const Matrix& actual, const Matrix& expected, const char* message) {
    check(actual.allFinite() && (actual - expected).norm() <= 3e-8 * expected.norm() + 1e-13, message);
}

Matrix prior_root() {
    Matrix root = Matrix::Zero();
    for (int row = 0; row < N; ++row) {
        root(row, row) = .04 + .015 * row;
        for (int column = 0; column < row; ++column)
            root(row, column) = .012 * std::sin(1.3 * row + .7 * column);
    }
    return root;
}

Q from_gibbs(const V3& error) {
    return Q(1, error.x()/2, error.y()/2, error.z()/2).normalized();
}

Vector new_coordinates(const Vector& error, const Vector& center, const Q& correction) {
    Vector output = error - center;
    const Q relative = from_gibbs(error.head<3>()) * correction.conjugate();
    output.head<3>() = 2 * relative.vec() / relative.w();
    return output;
}

void check_transport(double degrees) {
    const Q predicted = Q(Eigen::AngleAxis<real_t>(.47, V3::UnitX()))
        * Q(Eigen::AngleAxis<real_t>(-.29, V3::UnitY()));
    const Q correction(Eigen::AngleAxis<real_t>(degrees * std::acos(-1.)/180,
                                               V3(.3, -.4, .5).normalized()));
    const Q corrected = correction * predicted;
    Vector center = Vector::LinSpaced(N, -.12, .09);
    center.head<3>() = 2 * correction.vec() / correction.w();
    Matrix numerical;
    constexpr double epsilon = 1e-6;
    for (int column = 0; column < N; ++column) {
        Vector plus = center, minus = center;
        plus[column] += epsilon;
        minus[column] -= epsilon;
        numerical.col(column) = (new_coordinates(plus, center, correction)
            - new_coordinates(minus, center, correction)) / (2 * epsilon);
    }
    const Matrix old_root = prior_root();
    const Matrix posterior = old_root * old_root.transpose();
    check(posterior.topRightCorner<3, 6>().norm() > 1e-4,
          "fixture must contain attitude/velocity/bias correlations");
    const Matrix root = mav_detail::recenter_root<N>(predicted, corrected, posterior);
    close(root * root.transpose(), numerical * posterior * numerical.transpose(),
          "joint covariance transport disagrees with nonlinear finite differences");
    check(root.triangularView<Eigen::StrictlyUpper>().toDenseMatrix().norm() < 1e-13,
          "joint covariance root is not lower triangular");
    Q opposite = predicted;
    opposite.coeffs() *= -1;
    const Matrix opposite_root = mav_detail::recenter_root<N>(opposite, corrected, posterior);
    close(opposite_root * opposite_root.transpose(), root * root.transpose(),
          "quaternion sign changed joint covariance transport");
}

void check_constrained_update(bool freeze_bias) {
    const Matrix state_root = prior_root();
    const Matrix prior = state_root * state_root.transpose();
    UKF::Matrix<3, N> H = UKF::Matrix<3, N>::Zero();
    H(0, 3) = .7;
    H(1, 4) = -.4;
    H(2, 5) = .5;
    UKF::Matrix<3, 3> noise_root;
    noise_root << .03, 0, 0, .01, .04, 0, -.005, .01, .02;
    const Cross cross = prior * H.transpose();
    const DynamicMatrix S = H * prior * H.transpose() + noise_root * noise_root.transpose();
    const DynamicVector innovation = V3(.3, -.2, .1);
    const auto result = mav_detail::constrained_update<N>(prior, cross, S, innovation, freeze_bias);
    check(result.delta[2] == 0, "IMU model directly corrected world yaw");
    check(std::abs(result.posterior(2, 2) - prior(2, 2)) < 1e-14,
          "constrained update changed world-yaw variance before recentering");
    check(Eigen::SelfAdjointEigenSolver<Matrix>(result.posterior).eigenvalues().minCoeff() >= -1e-13,
          "constrained joint posterior is indefinite");
    if (freeze_bias) {
        check(result.delta.tail<3>().norm() == 0, "frozen gyro bias changed");
        check((result.posterior.bottomRightCorner<3, 3>() - prior.bottomRightCorner<3, 3>()).norm() < 1e-14,
              "frozen gyro-bias covariance changed before recentering");
    } else {
        check(result.delta.tail<3>().norm() > 1e-4, "fixture does not exercise motion-bias correction");
    }

    // Recover the implemented gain using unit residuals, then form actual
    // state errors after correction for independent state/noise sigma samples.
    Cross gain(N, 3);
    for (int column = 0; column < 3; ++column) {
        DynamicVector unit = DynamicVector::Zero(3);
        unit[column] = 1;
        gain.col(column) = mav_detail::constrained_update<N>(prior, cross, S, unit, freeze_bias).delta;
    }
    Matrix sample_covariance = Matrix::Zero();
    for (int component = 0; component < N + 3; ++component) {
        for (int sign : {-1, 1}) {
            Vector error = Vector::Zero();
            V3 noise = V3::Zero();
            if (component < N) error = sign * std::sqrt(double(N + 3)) * state_root.col(component);
            else noise = sign * std::sqrt(double(N + 3)) * noise_root.col(component - N);
            const Vector residual = error - gain * (H * error + noise);
            sample_covariance += residual * residual.transpose() / (2 * (N + 3));
        }
    }
    close(result.posterior, sample_covariance,
          "joint posterior disagrees with independent corrected error samples");

    // World-yaw gauge rotates only the left attitude errors; body-frame
    // velocity and gyro bias retain their original coordinates.
    const Q yaw(Eigen::AngleAxis<real_t>(.73, V3::UnitZ()));
    Matrix transform = Matrix::Identity();
    transform.topLeftCorner<3, 3>() = yaw.toRotationMatrix();
    const auto transformed = mav_detail::constrained_update<N>(
        (transform * prior * transform.transpose()).eval(), (transform * cross).eval(),
        S, innovation, freeze_bias);
    close(transformed.posterior, transform * result.posterior * transform.transpose(),
          "world-yaw frame changed the physical joint posterior");
    check((transformed.delta - transform * result.delta).norm() < 1e-13,
          "world-yaw frame changed the physical joint correction");

    const Q predicted = Q(Eigen::AngleAxis<real_t>(.4, V3::UnitX()));
    State state;
    state.set_field<901>(predicted);
    state.set_field<902>(V3(.2, -.3, .4));
    state.set_field<903>(V3(.01, -.02, .03));
    State rotated = state;
    rotated.set_field<901>(yaw * predicted);
    state.apply_delta(result.delta);
    rotated.apply_delta(transformed.delta);
    const Q corrected = state.get_field<901>();
    check(std::abs((corrected * predicted.conjugate()).z()) < 1e-13,
          "actual UKF quaternion correction contains a world-yaw twist");
    check((rotated.get_field<901>().toRotationMatrix() - (yaw * corrected).toRotationMatrix()).norm() < 1e-13,
          "actual UKF state correction is not world-yaw equivariant");
    check((rotated.get_field<902>() - state.get_field<902>()).norm() < 1e-13 &&
          (rotated.get_field<903>() - state.get_field<903>()).norm() < 1e-13,
          "world-yaw gauge changed body-frame velocity or bias correction");
}
}

int main() {
    try {
        check_constrained_update(false);
        check_constrained_update(true);
        for (double degrees : {0., 5., 30., 90., 150.}) check_transport(degrees);
        bool rejected = false;
        try {
            const Q singular(Eigen::AngleAxis<real_t>(std::acos(-1.), V3::UnitZ()));
            mav_detail::recenter_root<N>(Q::Identity(), singular, Matrix::Identity());
        } catch (const std::domain_error&) { rejected = true; }
        check(rejected, "Gibbs-chart singularity was accepted");
        std::cout << "Joint MAV covariance: independent error samples, bias constraints, cross-block transport and yaw gauge passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
