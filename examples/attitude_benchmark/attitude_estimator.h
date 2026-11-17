/*
Attitude estimation from a gyroscope, an accelerometer and (optionally) a
magnetometer, as used for the BROAD and Sassari benchmark recordings.

The filter is a small strapdown inertial navigator. The gyroscope turns the
attitude, and the accelerometer, rotated into the world frame with gravity
removed, changes the velocity. Nothing measures that velocity, but a sensor
moved by hand never gets far, so the filter is simply told that its velocity
is zero, give or take. That is enough to find which way is down: if the
attitude estimate is tilted, some of gravity is mistaken for acceleration,
the velocity estimate runs away, and the "measurement" of zero pulls both
the velocity and the tilt back. Real accelerations do little harm, because
they average out within a few seconds while the effect of a tilt error
keeps growing.

The gyroscope bias is part of the same state vector, so the filter learns it
from those same velocity errors while the sensor moves. Whenever the sensor
is at rest, the gyroscope reading itself is a direct measurement of the bias.

None of this can observe heading. If there is a magnetometer, a second,
one-state filter estimates the rotation about the vertical which points the
attitude at magnetic north. Keeping it separate means that a disturbed
magnetic field can spoil the heading, but never the tilt.
*/

#ifndef ATTITUDE_ESTIMATOR_H
#define ATTITUDE_ESTIMATOR_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include "UKF/Types.h"
#include "UKF/Integrator.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"
#include "imu_sample.h"
#include "sensor_checks.h"

/* Value of g in m/s^2. */
#define G_ACCEL (9.80665)

/* Convert degrees to radians. */
#define RADIANS(degrees) ((degrees) * 3.14159265358979323846 / 180.0)

enum Attitude_Keys {
    /* Attitude filter state vector fields. */
    Attitude,
    Velocity,
    GyroBias,

    /* Attitude filter measurement vector fields. */
    ZeroVelocity,
    RestGyroscope,

    /* Heading filter state and measurement vector fields. */
    HeadingOffset,
    MagneticHeading
};

/*
The attitude filter state vector contains the following:
- Attitude as a quaternion (body frame to world frame; world +Z is up)
- Velocity (world frame, m/s)
- Gyroscope bias (body frame, rad/s)
*/
using Attitude_StateVector = UKF::StateVector<
    UKF::Field<Attitude, UKF::Quaternion>,
    UKF::Field<Velocity, UKF::Vector<3>>,
    UKF::Field<GyroBias, UKF::Vector<3>>
>;

/*
The measurement vector is dynamic, because the gyroscope bias can only be
measured directly while the sensor is at rest:
- Velocity pseudo-measurement, always zero (world frame, m/s)
- Gyroscope reading at rest (body frame, rad/s)
*/
using Attitude_MeasurementVector = UKF::DynamicMeasurementVector<
    UKF::Field<ZeroVelocity, UKF::Vector<3>>,
    UKF::Field<RestGyroscope, UKF::Vector<3>>
>;

/* The process model input: the IMU readings to integrate over one step. */
struct Attitude_Input {
    UKF::Vector<3> gyroscope;       /* Angular velocity (body frame, rad/s) */
    UKF::Vector<3> accelerometer;   /* Specific force (body frame, m/s^2) */
};

/*
The heading filter estimates a single angle: the rotation about world +Z
which turns the attitude filter's arbitrary heading to magnetic north. Its
measurement is that same angle, as seen by the magnetometer.
*/
using Heading_StateVector = UKF::StateVector<
    UKF::Field<HeadingOffset, real_t>
>;

using Heading_MeasurementVector = UKF::FixedMeasurementVector<
    UKF::Field<MagneticHeading, real_t>
>;

namespace UKF {

/* Attitude filter process model. */
template <> template <>
inline Attitude_StateVector Attitude_StateVector::derivative<Attitude_Input>(
        const Attitude_Input& input) const {
    Attitude_StateVector output;

    /* Calculate change in attitude from the bias-corrected angular velocity. */
    UKF::Quaternion omega_q;
    omega_q.vec() = (input.gyroscope - get_field<GyroBias>()) * real_t(0.5);
    omega_q.w() = 0;
    output.set_field<Attitude>(get_field<Attitude>() * omega_q);

    /*
    Calculate change in velocity. An accelerometer measures acceleration
    less gravity, so gravity is added back once the reading has been rotated
    into the world frame.
    */
    output.set_field<Velocity>(
        get_field<Attitude>() * input.accelerometer + UKF::Vector<3>(0, 0, -G_ACCEL));

    /* Assume constant gyroscope bias; the process noise lets it wander. */
    output.set_field<GyroBias>(UKF::Vector<3>(0, 0, 0));

    return output;
}

/*
Attitude filter measurement model. The velocity is "measured" as zero, so
the expected measurement is just the velocity state.
*/
template <> template <>
inline UKF::Vector<3> Attitude_MeasurementVector::expected_measurement
<Attitude_StateVector, ZeroVelocity>(const Attitude_StateVector& state) {
    return state.get_field<Velocity>();
}

/* At rest the true angular velocity is zero, leaving only the bias. */
template <> template <>
inline UKF::Vector<3> Attitude_MeasurementVector::expected_measurement
<Attitude_StateVector, RestGyroscope>(const Attitude_StateVector& state) {
    return state.get_field<GyroBias>();
}

/*
Heading filter process model. The offset does not change by itself; its
process noise accounts for the attitude filter's heading drifting.
*/
template <> template <>
inline Heading_StateVector Heading_StateVector::process_model<>(real_t) const {
    return *this;
}

/* Heading filter measurement model. */
template <> template <>
inline real_t Heading_MeasurementVector::expected_measurement
<Heading_StateVector, MagneticHeading>(const Heading_StateVector& state) {
    return state.get_field<HeadingOffset>();
}

}

/*
The same models run on either core. The standard UKF propagates covariance
matrices; the square-root UKF propagates their square roots instead, which
is numerically more robust.
*/
using Attitude_UKF = UKF::Core<
    Attitude_StateVector,
    Attitude_MeasurementVector,
    UKF::IntegratorRK4
>;

using Attitude_SRUKF = UKF::SquareRootCore<
    Attitude_StateVector,
    Attitude_MeasurementVector,
    UKF::IntegratorRK4
>;

using Heading_UKF = UKF::Core<
    Heading_StateVector,
    Heading_MeasurementVector
>;

using Heading_SRUKF = UKF::SquareRootCore<
    Heading_StateVector,
    Heading_MeasurementVector
>;

/*
The two cores are used in exactly the same way, except that UKF::Core takes
covariances (variances on the diagonal) where UKF::SquareRootCore takes root
covariances (standard deviations on the diagonal). These helpers accept
standard deviations and fill in whichever the filter uses, so that the
estimator below can be written once for both.
*/
template <typename S, typename M, typename I>
void set_state_sigma(UKF::Core<S, M, I>& filter, const typename S::StateVectorDelta& sigma) {
    filter.covariance = sigma.cwiseAbs2().asDiagonal();
}

template <typename S, typename M, typename I>
void set_state_sigma(UKF::SquareRootCore<S, M, I>& filter, const typename S::StateVectorDelta& sigma) {
    filter.root_covariance = sigma.asDiagonal();
}

template <typename S, typename M, typename I>
void set_process_noise_sigma(UKF::Core<S, M, I>& filter, const typename S::StateVectorDelta& sigma) {
    filter.process_noise_covariance = sigma.cwiseAbs2().asDiagonal();
}

template <typename S, typename M, typename I>
void set_process_noise_sigma(UKF::SquareRootCore<S, M, I>& filter, const typename S::StateVectorDelta& sigma) {
    filter.process_noise_root_covariance = sigma.asDiagonal();
}

template <typename S, typename M, typename I>
void set_measurement_sigma(UKF::Core<S, M, I>& filter, const typename M::CovarianceVector& sigma) {
    filter.measurement_covariance = sigma.cwiseAbs2();
}

template <typename S, typename M, typename I>
void set_measurement_sigma(UKF::SquareRootCore<S, M, I>& filter, const typename M::CovarianceVector& sigma) {
    filter.measurement_root_covariance = sigma;
}

/*
Nothing the attitude filter measures depends on heading, so its uncertainty
about the heading can only grow. It must not grow without limit: the
unscented transform spreads its sigma points in proportion to the
uncertainty, and attitude errors approaching half a turn cannot be
represented. So if the heading uncertainty exceeds a limit, scale it back.

The attitude error is expressed in the world frame, so its third component
(index 2 of the covariance, the attitude being the first field) is the
heading error. Scaling that row and column keeps the heading error's
correlation with every other state, and leaves the rest of the covariance
untouched.
*/
template <typename S, typename M, typename I>
void limit_heading_sigma(UKF::Core<S, M, I>& filter, real_t limit) {
    real_t sigma = std::sqrt(filter.covariance(2, 2));
    if (sigma > limit) {
        filter.covariance.row(2) *= limit / sigma;
        filter.covariance.col(2) *= limit / sigma;
    }
}

template <typename S, typename M, typename I>
void limit_heading_sigma(UKF::SquareRootCore<S, M, I>& filter, real_t limit) {
    real_t sigma = filter.root_covariance.row(2).norm();
    if (sigma > limit) {
        filter.root_covariance.row(2) *= limit / sigma;
    }
}

/*
Filter tuning, in one place. Sensor noise is given as a density: the
standard deviation which builds up over one second. The estimator scales it
to the actual time step, so the same values work at any sample rate.
*/
struct Attitude_Tuning {
    /* Initial uncertainty (standard deviations). */
    real_t initial_tilt = 0.05;                         /* rad */
    real_t initial_velocity = 0.1;                      /* m/s */
    real_t initial_gyro_bias = RADIANS(0.5);            /* rad/s */

    /* Process noise. */
    real_t gyroscope_noise = 0.003;                     /* rad/sqrt(s) */
    real_t gyroscope_scale_error = 0.002;               /* rad/sqrt(s) per rad/s of rotation */
    real_t accelerometer_noise = 0.05;                  /* m/s/sqrt(s) */
    real_t gyro_bias_walk = RADIANS(0.001);             /* rad/s/sqrt(s) */

    /* Measurement noise. */
    real_t zero_velocity_noise = 0.09;                  /* m/s * sqrt(s) */
    real_t rest_gyroscope_noise = RADIANS(0.012);       /* rad/s * sqrt(s) */

    /* Largest heading uncertainty kept by the attitude filter. */
    real_t heading_sigma_limit = 0.3;                   /* rad */

    /* When is the sensor at rest? */
    RestDetector::Limits rest = {
        0.5,                                            /* averaging_time (s) */
        RADIANS(2.0),                                   /* gyroscope_change (rad/s) */
        0.5,                                            /* accelerometer_change (m/s^2) */
        0.7 * G_ACCEL,                                  /* accelerometer_minimum (m/s^2) */
        1.3 * G_ACCEL,                                  /* accelerometer_maximum (m/s^2) */
        RADIANS(2.0),                                   /* gyroscope_bias (rad/s) */
        1.5                                             /* minimum_time (s) */
    };

    /* Heading filter. */
    real_t initial_heading = RADIANS(10.0);             /* rad */
    real_t heading_walk = RADIANS(0.17);                /* rad/sqrt(s) */
    real_t magnetic_heading_noise = RADIANS(3.0);       /* rad * sqrt(s) */

    /* When is the magnetic field disturbed? */
    MagneticFieldCheck::Limits magnetic_field = {
        0.1,                                            /* strength_change (fraction) */
        RADIANS(10.0),                                  /* dip_change (rad) */
        0.05,                                           /* smoothing_time (s) */
        0.5,                                            /* settle_time (s) */
        20.0,                                           /* reference_time (s) */
        20.0,                                           /* new_field_time (s) */
        RADIANS(20.0)                                   /* new_field_rotation (rad/s) */
    };
};

/*
The estimator: an attitude filter, a heading filter, and the two sensor
checks which decide what each filter gets to measure. It works with either
pair of filter types defined above.

Because the filter types are template parameters, the methods below first
name the filter's state vector with its actual type (which is the same for
both cores); its fields can then be accessed in the usual way.
*/
template <typename AttitudeFilter, typename HeadingFilter>
class Attitude_Estimator {
public:
    explicit Attitude_Estimator(const Attitude_Tuning& tuning = Attitude_Tuning()) :
        tuning(tuning), rest_detector(tuning.rest), field_check(tuning.magnetic_field) {}

    /*
    Start from the first sample. At rest the accelerometer points straight
    up, which fixes the initial tilt; the initial heading is arbitrary.
    Returns false if the accelerometer reading is zero.
    */
    bool initialise(const ImuSample& sample, bool with_magnetometer) {
        if (sample.accelerometer.norm() < real_t(1e-6)) {
            return false;
        }
        use_magnetometer = with_magnetometer;

        /* Initialise state vector and covariance. */
        Attitude_StateVector& state = attitude_filter.state;
        state.set_field<Attitude>(UKF::Quaternion::FromTwoVectors(
            sample.accelerometer, UKF::Vector<3>(0, 0, 1)));
        state.set_field<Velocity>(UKF::Vector<3>(0, 0, 0));
        state.set_field<GyroBias>(UKF::Vector<3>(0, 0, 0));

        Attitude_StateVector::StateVectorDelta sigma;
        sigma <<
            tuning.initial_tilt * UKF::Vector<3>::Ones(),
            tuning.initial_velocity * UKF::Vector<3>::Ones(),
            tuning.initial_gyro_bias * UKF::Vector<3>::Ones();
        set_state_sigma(attitude_filter, sigma);

        at_rest = rest_detector.update(sample.gyroscope, sample.accelerometer, 0);

        return !use_magnetometer || update_heading(sample, 0);
    }

    /*
    Advance the estimate from the previous sample to this one. Returns false
    if a filter step fails.
    */
    bool update(const ImuSample& previous, const ImuSample& sample) {
        const Attitude_StateVector& state = attitude_filter.state;
        real_t dt = sample.time - previous.time;

        /* Integrate the average of the two IMU readings over the step. */
        Attitude_Input input;
        input.gyroscope = (previous.gyroscope + sample.gyroscope) * real_t(0.5);
        input.accelerometer = (previous.accelerometer + sample.accelerometer) * real_t(0.5);
        angular_velocity = (input.gyroscope - state.get_field<GyroBias>()).norm();

        /*
        Set the process noise which builds up over this step. A gyroscope's
        scale error grows with the rate of rotation, so fast turns are
        trusted a little less.
        */
        Attitude_StateVector::StateVectorDelta process_noise;
        process_noise <<
            std::hypot(tuning.gyroscope_noise, tuning.gyroscope_scale_error * angular_velocity) *
                UKF::Vector<3>::Ones(),
            tuning.accelerometer_noise * UKF::Vector<3>::Ones(),
            tuning.gyro_bias_walk * UKF::Vector<3>::Ones();
        set_process_noise_sigma(attitude_filter, process_noise * std::sqrt(dt));

        /*
        Set the measurement noise. Each measurement must count for less at a
        higher sample rate, where there are more of them every second.
        */
        Attitude_MeasurementVector::CovarianceVector measurement_noise;
        measurement_noise.set_field<ZeroVelocity>(tuning.zero_velocity_noise * UKF::Vector<3>::Ones());
        measurement_noise.set_field<RestGyroscope>(tuning.rest_gyroscope_noise * UKF::Vector<3>::Ones());
        set_measurement_sigma(attitude_filter, measurement_noise / std::sqrt(dt));

        /*
        The velocity pseudo-measurement applies at every step. The gyroscope
        measures its own bias only while the sensor is at rest.
        */
        at_rest = rest_detector.update(sample.gyroscope, sample.accelerometer, dt);
        Attitude_MeasurementVector measurement;
        measurement.set_field<ZeroVelocity>(UKF::Vector<3>(0, 0, 0));
        if (at_rest) {
            measurement.set_field<RestGyroscope>(sample.gyroscope);
        }

        if (!attitude_filter.a_priori_step(dt, input) ||
                !attitude_filter.innovation_step(measurement) ||
                !attitude_filter.a_posteriori_step()) {
            return false;
        }

        /* The attitude filter cannot observe heading; see above. */
        limit_heading_sigma(attitude_filter, tuning.heading_sigma_limit);

        return !use_magnetometer || update_heading(sample, dt);
    }

    /* The estimated attitude (body frame to world frame), heading included. */
    UKF::Quaternion attitude() const {
        const Attitude_StateVector& state = attitude_filter.state;
        UKF::Quaternion heading(Eigen::AngleAxis<real_t>(heading_offset(), UKF::Vector<3>(0, 0, 1)));
        return (heading * state.get_field<Attitude>()).normalized();
    }

    /* The rotation about world +Z applied by the heading filter (rad). */
    real_t heading_offset() const {
        const Heading_StateVector& heading_state = heading_filter.state;
        return heading_started ? heading_state.get_field<HeadingOffset>() : real_t(0.0);
    }

    /* Other estimates and decisions, for diagnostics. */
    UKF::Vector<3> velocity() const {
        const Attitude_StateVector& state = attitude_filter.state;
        return state.get_field<Velocity>();
    }

    UKF::Vector<3> gyro_bias() const {
        const Attitude_StateVector& state = attitude_filter.state;
        return state.get_field<GyroBias>();
    }

    bool resting() const { return at_rest; }
    bool magnetometer_rejected() const { return field_rejected; }

private:
    /* Wrap an angle into the range -pi to pi. */
    static real_t wrap(real_t angle) {
        return std::remainder(angle, real_t(2.0 * 3.14159265358979323846));
    }

    bool update_heading(const ImuSample& sample, real_t dt) {
        const Attitude_StateVector& state = attitude_filter.state;
        Heading_StateVector& heading_state = heading_filter.state;

        /*
        Rotate the magnetic field into the attitude filter's world frame.
        Its horizontal part points to magnetic north; the heading offset is
        the rotation about +Z which would bring that onto the +X axis. A
        vertical (or missing) field says nothing about heading, and a
        disturbed one must not be believed.
        */
        UKF::Vector<3> field = state.get_field<Attitude>() * sample.magnetometer;
        bool usable = field.head<2>().norm() > real_t(1e-6) * field.norm();
        field_rejected = !usable || field_check.update(field, dt, angular_velocity);
        real_t measured_offset = -std::atan2(field.y(), field.x());

        /* The first usable reading sets the initial heading. */
        if (!heading_started) {
            if (!usable) {
                return true;
            }
            heading_state.set_field<HeadingOffset>(measured_offset);

            Heading_StateVector::StateVectorDelta sigma;
            sigma << tuning.initial_heading;
            set_state_sigma(heading_filter, sigma);

            heading_started = true;
            return true;
        }

        /*
        Predict at every step, whether or not a measurement follows, so
        that the offset becomes less certain for as long as the
        magnetometer cannot be used.
        */
        Heading_StateVector::StateVectorDelta process_noise;
        process_noise << tuning.heading_walk * std::sqrt(dt);
        set_process_noise_sigma(heading_filter, process_noise);

        if (!heading_filter.a_priori_step(dt)) {
            return false;
        }

        /* Without a trustworthy measurement, the prediction is all there is. */
        if (field_rejected) {
            return true;
        }

        /*
        Angles wrap around, but the filter subtracts them like any other
        number. So express the measurement as the equivalent angle nearest
        to the current estimate.
        */
        real_t offset = heading_state.get_field<HeadingOffset>();
        Heading_MeasurementVector measurement;
        measurement.set_field<MagneticHeading>(offset + wrap(measured_offset - offset));

        Heading_MeasurementVector::CovarianceVector measurement_noise;
        measurement_noise.set_field<MagneticHeading>(tuning.magnetic_heading_noise / std::sqrt(dt));
        set_measurement_sigma(heading_filter, measurement_noise);

        return heading_filter.innovation_step(measurement) && heading_filter.a_posteriori_step();
    }

    Attitude_Tuning tuning;
    AttitudeFilter attitude_filter;
    HeadingFilter heading_filter;
    RestDetector rest_detector;
    MagneticFieldCheck field_check;
    real_t angular_velocity = 0.0;
    bool use_magnetometer = false;
    bool heading_started = false;
    bool at_rest = false;
    bool field_rejected = false;
};

#endif
