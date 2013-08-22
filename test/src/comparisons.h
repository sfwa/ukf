#include <cmath>

inline ::testing::AssertionResult CmpHelperQuaternionEq(
const char* expected_expression, const char* actual_expression,
Quaternionr a, Quaternionr b) {
    a.normalize();
    b.normalize();

    real_t angle = 2.0 * acos(std::abs(a.dot(b)));

    if (angle < 0.02) {
        return ::testing::AssertionSuccess();
    } else {
        return ::testing::AssertionFailure() <<
            "Value of: " << actual_expression << "\n"
            "Expected: " << expected_expression << "\n" <<
            "  Actual: Quaternionr(" << b.w() << ", " << b.x() << ", " <<
                b.y() << ", " << b.z() << ")\n" <<
            "   Error: " << angle << " rad\n";
    }
}

inline ::testing::AssertionResult CmpHelperVectorEq(
const char* expected_expression, const char* actual_expression,
Vector3r a, Vector3r b) {
    real_t d = (a - b).norm(), an = a.norm();

    if (d / std::max(an, 1e-3) < 0.01) {
        return ::testing::AssertionSuccess();
    } else {
        return ::testing::AssertionFailure() <<
            "Value of: " << actual_expression << "\n"
            "Expected: " << expected_expression << "\n" <<
            "  Actual: Vector3r(" << b.x() << ", " << b.y() << ", " << b.z()
                << ")\n";
    }
}

inline ::testing::AssertionResult CmpHelperMeasurementEq(
const char* expected_expression, const char* actual_expression,
MeasurementVector a, MeasurementVector b) {
    bool valid = true;

    if (a.rows() != b.rows()) {
        valid = false;
    } else {
        for (size_t i = 0, l = (size_t)a.rows(); i < l && valid; i++) {
            if (std::abs(a[i] - b[i]) > 0.001) {
                valid = false;
            }
        }
    }


    if (valid) {
        return ::testing::AssertionSuccess();
    } else {
        return ::testing::AssertionFailure() <<
            "Value of: " << actual_expression << "\n"
            "Expected: " << expected_expression << "\n" <<
            "  Actual: MeasurementVector(" << b.transpose() << ")\n";

    }
}

#define EXPECT_QUATERNION_EQ(expected, actual)\
    EXPECT_PRED_FORMAT2(CmpHelperQuaternionEq, expected, actual)

#define EXPECT_VECTOR_EQ(expected, actual)\
    EXPECT_PRED_FORMAT2(CmpHelperVectorEq, expected, actual)

#define EXPECT_MEASUREMENT_EQ(expected, actual)\
    EXPECT_PRED_FORMAT2(CmpHelperMeasurementEq, expected, actual)
