#include <Eigen/Eigenvalues>
#include <cmath>
#include <limits>

/*
A bound for each precision. The tests build in double precision (unittest)
and in single precision (unittest_single), which rounds about 5e8 times as
coarsely, so a check that is tight in double needs its own bound in single.
*/
constexpr real_t precision_tolerance(double double_tolerance,
                                     double single_tolerance) {
  return real_t(sizeof(real_t) == sizeof(float) ? single_tolerance
                                                : double_tolerance);
}

/*
A symmetric matrix is positive semidefinite to within roundoff: no
eigenvalue is below -16 n epsilon times its largest entry, the bound the
cores use. A determinant cannot show this for a singular covariance, whose
roundoff-sized eigenvalue takes either sign.
*/
template <typename MatrixType>
inline ::testing::AssertionResult CmpHelperPositiveSemidefinite(
    const char* expression, const MatrixType& m) {
  using Scalar = typename MatrixType::Scalar;
  const Eigen::SelfAdjointEigenSolver<MatrixType> eig(m);
  const Scalar tolerance = Scalar(16) * std::numeric_limits<Scalar>::epsilon() *
                           Scalar(m.rows()) * m.cwiseAbs().maxCoeff();
  if (eig.info() == Eigen::Success &&
      eig.eigenvalues().minCoeff() >= -tolerance) {
    return ::testing::AssertionSuccess();
  }
  return ::testing::AssertionFailure()
         << expression << " has eigenvalue " << eig.eigenvalues().minCoeff()
         << ", below -" << tolerance;
}

#define EXPECT_POSITIVE_SEMIDEFINITE(matrix) \
  EXPECT_PRED_FORMAT1(CmpHelperPositiveSemidefinite, matrix)

template <typename Q1, typename Q2>
inline ::testing::AssertionResult CmpHelperQuaternionEq(
    const char* expected_expression, const char* actual_expression, Q1 a, Q2 b,
    real_t tol = 0.02) {
  a.normalize();
  b.normalize();

  real_t angle = (real_t)2.0 * acos(std::abs(a.dot(b)));

  if (angle < tol) {
    return ::testing::AssertionSuccess();
  } else {
    return ::testing::AssertionFailure()
           << "Value of: " << actual_expression
           << "\n"
              "Expected: "
           << expected_expression << "\n"
           << "  Actual: Quaternion(" << b.w() << ", " << b.x() << ", " << b.y()
           << ", " << b.z() << ")\n"
           << "   Error: " << angle << " rad\n";
  }
}

template <typename T1, typename T2>
inline ::testing::AssertionResult CmpHelperVectorEq(
    const char* expected_expression, const char* actual_expression, T1 a, T2 b,
    real_t tol = 0.0001) {
  real_t d = (a - b).norm(), an = a.norm();

  Eigen::IOFormat CmpFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ",
                            ", ", "", "", "Vector(", ")");

  if (d / std::max(an, (real_t)1e-4) < tol) {
    return ::testing::AssertionSuccess();
  } else {
    return ::testing::AssertionFailure()
           << "Value of: " << actual_expression
           << "\n"
              "Expected: "
           << expected_expression
           << "\n"
              "Which is: "
           << a.format(CmpFormat) << "\n"
           << "  Actual: " << b.format(CmpFormat) << "\n";
  }
}

#define EXPECT_QUATERNION_EQ(expected, actual) \
  EXPECT_PRED_FORMAT2(CmpHelperQuaternionEq, expected, actual)

#define EXPECT_VECTOR_EQ(expected, actual) \
  EXPECT_PRED_FORMAT2(CmpHelperVectorEq, expected, actual)
