#include <gtest/gtest.h>
#include "types.h"
#include "state.h"
#include "sensors.h"
#include "dynamics.h"
#include "ukf.h"
#include "comparisons.h"

#include <Eigen/Core>
#include <xmmintrin.h>

TEST(UKFTest, Instantiation) {
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    UnscentedKalmanFilter ukf = UnscentedKalmanFilter(test);
    CentripetalModel model = CentripetalModel();
    ukf.set_dynamics_model(&model);
}

TEST(UKFTest, NoSensorConstantVelocity) {
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    UnscentedKalmanFilter ukf = UnscentedKalmanFilter(test);
    State test_state;
    test_state << 0, 0, 0,
                  10, 0, 0,
                  0, 0, 0,
                  0, 0, 0, 1,
                  0, 0, 1,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0;
    IntegratorRK4 test_integrator = IntegratorRK4();
    ukf.set_state(test_state);

    for(real_t i = 0; i < 10; i += 0.001) {
        test_state = test_integrator.integrate(test_state, 0.001);
        ukf.iterate(0.001, ControlVector());
    }

    EXPECT_VECTOR_EQ(Vector3r(1.5785803e-05, 0, 0),
        ukf.get_state().position());
    EXPECT_VECTOR_EQ(Vector3r(10, 0, 0), ukf.get_state().velocity());
    EXPECT_VECTOR_EQ(Vector3r(0, 0, 0), ukf.get_state().acceleration());
    EXPECT_QUATERNION_EQ(Quaternionr(0.283662, 0, 0, 0.958924),
        Quaternionr(ukf.get_state().attitude()));
    EXPECT_VECTOR_EQ(Vector3r(0, 0, 1), ukf.get_state().angular_velocity());
    EXPECT_GE(0.001, ukf.get_state().angular_acceleration().norm());
}

TEST(UKFTest, NoSensorCircularMotion) {
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0));
    UnscentedKalmanFilter ukf = UnscentedKalmanFilter(test);
    State test_state;
    test_state << 0, 0, 0,
                  10, 0, 0,
                  0, 6.283, 0,
                  0, 0, 0, 1,
                  0, 0, 0.6283,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0;
    IntegratorRK4 test_integrator = IntegratorRK4();
    ukf.set_state(test_state);

    for(real_t i = 0; i < 10.001; i += 0.001) {
        test_state = test_integrator.integrate(test_state, 0.001);
        ukf.iterate(0.001, ControlVector());
    }

    EXPECT_VECTOR_EQ(Vector3r(0, 0, 0), ukf.get_state().position());
    EXPECT_VECTOR_EQ(Vector3r(10, 0, 0), ukf.get_state().velocity());
    EXPECT_VECTOR_EQ(Vector3r(0, 6.283, 0), ukf.get_state().acceleration());
    EXPECT_QUATERNION_EQ(Quaternionr(1, 0, 0, 0),
        Quaternionr(ukf.get_state().attitude()));
    EXPECT_VECTOR_EQ(Vector3r(0, 0, 0.6283),
        ukf.get_state().angular_velocity());
    EXPECT_GE(0.001, ukf.get_state().angular_acceleration().norm());
}

TEST(UKFTest, AccelerometerConstantAngularVelocity) {
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    IOBoardModel test_model = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(1, 0, 0));
    UnscentedKalmanFilter ukf = UnscentedKalmanFilter(test_model);
    State test_state;
    test_state << 0, 0, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0, 1,
                  0, 1, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0;
    ukf.set_state(test_state);
    MeasurementVector sensor_covariance(18);
    sensor_covariance <<
        10, 10, 10,
        0.01, 0.01, 0.01,
        0.01, 0.01, 0.01,
        5, 5, 5,
        5, 5, 5,
        1,
        1,
        1;
    test_model.set_covariance(sensor_covariance);
    CentripetalModel model = CentripetalModel();
    ukf.set_dynamics_model(&model);
    Eigen::Matrix<real_t, 4, 1> ref = Eigen::Matrix<real_t, 4, 1>(0, 0, 0, 1);

    for(real_t i = 0; i < 10; i += 0.001) {
        Quaternionr temp_q =
            (Quaternionr(0, 0, 1, 0).conjugate()*Quaternionr(ref));
        Eigen::Matrix<real_t, 4, 1> temp_v;
        temp_v << temp_q.vec(), temp_q.w();
        ref += 0.001 * 0.5 * temp_v;
        temp_q = Quaternionr(ref).normalized();
        ref << temp_q.vec(), temp_q.w();

        test_model.clear();
        test_model.set_accelerometer(
            Quaternionr(ref) *
            Vector3r(0, 0, -G_ACCEL));
        test_model.set_gps_velocity(Vector3r(0, 0, 0));
        ukf.iterate(0.001, ControlVector());
    }

    EXPECT_QUATERNION_EQ(Quaternionr(0.283662, 0, 0.958924, 0),
        Quaternionr(ukf.get_state().attitude()));
}

TEST(UKFTest, GyroscopeConstantAngularVelocity) {
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    IOBoardModel test_model = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(1, 0, 0));
    UnscentedKalmanFilter ukf = UnscentedKalmanFilter(test_model);
    State test_state;
    test_state << 0, 0, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0, 1,
                  0, 1, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0;
    ukf.set_state(test_state);
    MeasurementVector sensor_covariance(18);
    sensor_covariance <<
        0.001, 0.001, 0.001,
        0.01, 0.01, 0.01,
        0.01, 0.01, 0.01,
        5, 5, 5,
        5, 5, 5,
        1,
        1,
        1;
    test_model.set_covariance(sensor_covariance);

    for(real_t i = 0; i < 10; i += 0.01) {
        test_model.clear();
        test_model.set_gyroscope(Vector3r(0, 1, 0));
        ukf.iterate(0.01, ControlVector());
    }

    EXPECT_QUATERNION_EQ(Quaternionr(0.283662, 0, 0.958924, 0),
        Quaternionr(ukf.get_state().attitude()));
}

TEST(UKFTest, AngularSensorsConstantAngularVelocity) {
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    IOBoardModel test_model = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(1, 0, 0));
    UnscentedKalmanFilter ukf = UnscentedKalmanFilter(test_model);
    State test_state;
    test_state << 0, 0, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0, 1,
                  0, 1, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0;
    ukf.set_state(test_state);
    MeasurementVector sensor_covariance(18);
    sensor_covariance <<
        0.001, 0.001, 0.001,
        0.01, 0.01, 0.01,
        0.01, 0.01, 0.01,
        5, 5, 5,
        5, 5, 5,
        1,
        1,
        1;
    test_model.set_covariance(sensor_covariance);
    Eigen::Matrix<real_t, 4, 1> ref = Eigen::Matrix<real_t, 4, 1>(0, 0, 0, 1);

    for(real_t i = 0; i < 10; i += 0.01) {
        Quaternionr temp_q =
            (Quaternionr(0, 0, 1, 0).conjugate()*Quaternionr(ref));
        Eigen::Matrix<real_t, 4, 1> temp_v;
        temp_v << temp_q.vec(), temp_q.w();
        ref += 0.01 * 0.5 * temp_v;
        temp_q = Quaternionr(ref).normalized();
        ref << temp_q.vec(), temp_q.w();

        test_model.clear();
        test_model.set_accelerometer(
            Quaternionr(ref) *
            Vector3r(0, 0, -G_ACCEL));
        test_model.set_magnetometer(
            Quaternionr(ref) *
            Vector3r(1, 0, 0));
        test_model.set_gyroscope(Vector3r(0, 1, 0));
        ukf.iterate(0.01, ControlVector());
    }

    EXPECT_QUATERNION_EQ(Quaternionr(0.283662, 0, 0.958924, 0),
        Quaternionr(ukf.get_state().attitude()));
}

TEST(UKFTest, MagnetometerAccelerometerAtRest) {
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    IOBoardModel test = IOBoardModel(
        Quaternionr(1, 0, 0, 0),
        Vector3r(0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Quaternionr(1, 0, 0, 0),
        Vector3r(21.2578, 4.4132, -55.9578));
    UnscentedKalmanFilter ukf = UnscentedKalmanFilter(test);
    State test_state;
    test_state << 0, 0, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0, 1,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0;
    IntegratorRK4 test_integrator = IntegratorRK4();
    ukf.set_state(test_state);
    MeasurementVector sensor_covariance(18);
    sensor_covariance <<
        10, 10, 10,
        0.01, 0.01, 0.01,
        25, 25, 25,
        5, 5, 5,
        5, 5, 5,
        1,
        1,
        1;
    test.set_covariance(sensor_covariance);
    CentripetalModel model = CentripetalModel();
    ukf.set_dynamics_model(&model);

    for(real_t i = 0; i < 10; i += 0.001) {
        test.clear();
        test.set_magnetometer(
            Vector3r(-4.4132, 21.2578, -55.9578));
        test.set_accelerometer(
            Vector3r(0, 0, -G_ACCEL));
        test.set_gps_velocity(Vector3r(0, 0, 0));
        ukf.iterate(0.001, ControlVector());
    }

    EXPECT_QUATERNION_EQ(Quaternionr(0.707107, 0, 0, 0.707107),
        Quaternionr(ukf.get_state().attitude()));
}
