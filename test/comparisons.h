#include <cmath>

template <typename Q1, typename Q2>
inline ::testing::AssertionResult CmpHelperQuaternionEq(
const char* expected_expression, const char* actual_expression,
Q1 a, Q2 b, real_t tol = 0.02) {
    a.normalize();
    b.normalize();

    real_t angle = (real_t)2.0 * acos(std::abs(a.dot(b)));

    if (angle < tol) {
        return ::testing::AssertionSuccess();
    } else {
        return ::testing::AssertionFailure() <<
            "Value of: " << actual_expression << "\n"
            "Expected: " << expected_expression << "\n" <<
            "  Actual: Quaternion(" << b.w() << ", " << b.x() << ", " <<
                b.y() << ", " << b.z() << ")\n" <<
            "   Error: " << angle << " rad\n";
    }
}

template <typename T1, typename T2>
inline ::testing::AssertionResult CmpHelperVectorEq(
const char* expected_expression, const char* actual_expression,
T1 a, T2 b, real_t tol = 0.0001) {
    real_t d = (a - b).norm(), an = a.norm();

    Eigen::IOFormat CmpFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "Vector(", ")");

    if (d / std::max(an, (real_t)1e-4) < tol) {
        return ::testing::AssertionSuccess();
    } else {
        return ::testing::AssertionFailure() <<
            "Value of: " << actual_expression << "\n"
            "Expected: " << a.format(CmpFormat) << "\n" <<
            "  Actual: " << b.format(CmpFormat) << "\n";
    }
}

#define EXPECT_QUATERNION_EQ(expected, actual)\
    EXPECT_PRED_FORMAT2(CmpHelperQuaternionEq, expected, actual)

#define EXPECT_VECTOR_EQ(expected, actual)\
    EXPECT_PRED_FORMAT2(CmpHelperVectorEq, expected, actual)
