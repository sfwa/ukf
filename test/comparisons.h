#include <cmath>

template <typename Q1, typename Q2>
inline ::testing::AssertionResult CmpHelperQuaternionEq(
const char* expected_expression, const char* actual_expression,
Q1 a, Q2 b) {
    a.normalize();
    b.normalize();

    real_t angle = (real_t)2.0 * acos(std::abs(a.dot(b)));

    if (angle < 0.02) {
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
T1 a, T2 b) {
    real_t d = (a - b).norm(), an = a.norm();

    if (d / std::max(an, (real_t)1e-3) < 0.01) {
        return ::testing::AssertionSuccess();
    } else {
        return ::testing::AssertionFailure() <<
            "Value of: " << actual_expression << "\n"
            "Expected: " << expected_expression << "\n" <<
            "  Actual: Vector(" << b.x() << ", " << b.y() << ", " << b.z()
                << ")\n";
    }
}

#define EXPECT_QUATERNION_EQ(expected, actual)\
    EXPECT_PRED_FORMAT2(CmpHelperQuaternionEq, expected, actual)

#define EXPECT_VECTOR_EQ(expected, actual)\
    EXPECT_PRED_FORMAT2(CmpHelperVectorEq, expected, actual)
