// MIT license; see the repository LICENSE.
// Independent error-sample and nonlinear-coordinate checks of joint covariance.
// Unless the build forbids Eigen heap allocation outright, enable Eigen's
// runtime check so the correction path can be shown not to allocate.
#if !defined(EIGEN_NO_MALLOC) && !defined(EIGEN_RUNTIME_NO_MALLOC)
#define EIGEN_RUNTIME_NO_MALLOC
#endif
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>
#include "mav_correction.h"
#include "UKF/StateVector.h"

namespace {
constexpr int N = 9;
using Matrix = UKF::Matrix<N, N>;
using Vector = UKF::Vector<N>;
using V3 = UKF::Vector<3>;
using Q = UKF::Quaternion;
constexpr int MaxM = 7;
using Cross = mav_detail::CrossMatrix<N, MaxM>;
using DynamicMatrix = mav_detail::SquareMatrix<MaxM>;
using DynamicVector = mav_detail::MeasurementVector<MaxM>;
using State = UKF::StateVector<UKF::Field<901, Q>, UKF::Field<902, V3>, UKF::Field<903, V3>>;

// Checks report the first failure through this flag; no exceptions are used.
const char* failure = nullptr;

bool check(bool condition, const char* message) {
    if (!condition && !failure) failure = message;
    return condition;
}

// While alive, any Eigen heap allocation fails an assertion (in builds with
// assertions enabled and EIGEN_RUNTIME_NO_MALLOC available).
struct NoHeapAllocation {
#if defined(EIGEN_RUNTIME_NO_MALLOC) && !defined(EIGEN_NO_MALLOC)
    NoHeapAllocation() { Eigen::internal::set_is_malloc_allowed(false); }
    ~NoHeapAllocation() { Eigen::internal::set_is_malloc_allowed(true); }
#endif
};

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
    Matrix root;
    bool transported = false;
    {
        [[maybe_unused]] NoHeapAllocation guard;
        transported = mav_detail::recenter_root<N>(predicted, corrected, posterior, root);
    }
    if (!check(transported, "joint covariance transport rejected a valid correction")) return;
    close(root * root.transpose(), numerical * posterior * numerical.transpose(),
          "joint covariance transport disagrees with nonlinear finite differences");
    check(root.triangularView<Eigen::StrictlyUpper>().toDenseMatrix().norm() < 1e-13,
          "joint covariance root is not lower triangular");
    Q opposite = predicted;
    opposite.coeffs() *= -1;
    Matrix opposite_root;
    if (!check(mav_detail::recenter_root<N>(opposite, corrected, posterior, opposite_root),
               "joint covariance transport rejected a sign-flipped prediction")) return;
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
    mav_detail::Correction<N> result;
    bool updated = false;
    {
        [[maybe_unused]] NoHeapAllocation guard;
        updated = mav_detail::constrained_update<N, MaxM>(prior, cross, S, innovation, freeze_bias, result);
    }
    if (!check(updated, "constrained update rejected valid moments")) return;
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
        mav_detail::Correction<N> unit_result;
        if (!check(mav_detail::constrained_update<N, MaxM>(prior, cross, S, unit, freeze_bias, unit_result),
                   "constrained update rejected a unit innovation")) return;
        gain.col(column) = unit_result.delta;
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
    mav_detail::Correction<N> transformed;
    if (!check(mav_detail::constrained_update<N, MaxM>(
            (transform * prior * transform.transpose()).eval(), Cross(transform * cross),
            S, innovation, freeze_bias, transformed),
            "constrained update rejected yaw-rotated moments")) return;
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

// The bias constraint needs a bias block; a smaller state must be rejected.
void check_bias_constraint_requires_bias() {
    constexpr int Small = 6;
    const UKF::Matrix<Small, Small> prior = UKF::Matrix<Small, Small>::Identity();
    mav_detail::CrossMatrix<Small, MaxM> cross = mav_detail::CrossMatrix<Small, MaxM>::Zero(Small, 1);
    cross(3, 0) = .5;
    DynamicMatrix S = DynamicMatrix::Identity(1, 1);
    DynamicVector innovation = DynamicVector::Constant(1, .1);
    mav_detail::Correction<Small> result;
    check(mav_detail::constrained_update<Small, MaxM>(prior, cross, S, innovation, false, result),
          "constrained update rejected a state without a bias block");
    check(!mav_detail::constrained_update<Small, MaxM>(prior, cross, S, innovation, true, result),
          "bias constraint was accepted without a bias state");
}
}

int main() {
    check_constrained_update(false);
    check_constrained_update(true);
    check_bias_constraint_requires_bias();
    for (double degrees : {0., 5., 30., 90., 150.}) check_transport(degrees);
    const Q singular(Eigen::AngleAxis<real_t>(std::acos(-1.), V3::UnitZ()));
    Matrix unchanged = Matrix::Constant(7);
    check(!mav_detail::recenter_root<N>(Q::Identity(), singular, Matrix::Identity(), unchanged),
          "Gibbs-chart singularity was accepted");
    check(unchanged == Matrix::Constant(7), "rejected transport modified its output");
    if (failure) {
        std::cerr << failure << '\n';
        return 1;
    }
    std::cout << "Joint MAV covariance: independent error samples, bias constraints, cross-block transport and yaw gauge passed\n";
    return 0;
}
