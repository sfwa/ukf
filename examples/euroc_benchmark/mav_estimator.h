/*
Attitude estimation for a multirotor, as used for the EuRoC MAV benchmark
recordings (an AscTec Firefly hexrotor carrying a six-axis IMU).

This is the strapdown filter of the attitude benchmark (see
../attitude_benchmark/attitude_estimator.h) with one addition: a model of
how a multirotor flies. Its rotors can only push along one axis, the thrust
axis. The only other force on the airframe is rotor drag, which opposes the
velocity across the rotor plane and is proportional to it. So in powered
flight

    lateral accelerometer reading = -drag rate * lateral velocity

which turns the accelerometer into a velocity sensor. That is far more
informative than the "velocity is roughly zero" assumption made for a
hand-held sensor, and through the velocity it pins down the tilt. See
Leishman et al., "Quadrotors and Accelerometers", IEEE Control Systems
Magazine 34(1), 2014.

The rotor model applies only while the rotors are turning, which is detected
from the vibration they cause. Otherwise (while the vehicle is carried, or
standing on the ground) the filter behaves like the hand-held one.

Unlike the hand-held filter, the velocity is expressed in the IMU frame.
Rotor drag acts in that frame, so the rotor model then never involves the
heading, which a six-axis IMU cannot observe.
*/

#ifndef MAV_ESTIMATOR_H
#define MAV_ESTIMATOR_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include "UKF/Types.h"
#include "UKF/Integrator.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"
#include "imu_csv.h"
#include "sensor_checks.h"
#include "rotor_vibration.h"

/* Value of g in m/s^2. */
#define G_ACCEL (9.80665)

/* Convert degrees to radians. */
#define RADIANS(degrees) ((degrees) * 3.14159265358979323846 / 180.0)

/*
Rotor drag: the lateral deceleration per unit of lateral velocity, in 1/s.
*/
#define MAV_DRAG_RATE (0.25)

/*
The IMU is not mounted square to the airframe: the rotor thrust is tilted
by this angle from the IMU's +X axis towards its -Z axis. The value comes
from calibrate_mount.py.
*/
#define MAV_THRUST_TILT RADIANS(20.0)

/*
The direction of rotor thrust, and two perpendicular directions which span
the rotor plane, as unit vectors in the IMU frame.
*/
static const UKF::Vector<3> MAV_ThrustAxis(std::cos(MAV_THRUST_TILT), 0, -std::sin(MAV_THRUST_TILT));
static const UKF::Vector<3> MAV_LateralAxis1(0, 1, 0);
static const UKF::Vector<3> MAV_LateralAxis2(std::sin(MAV_THRUST_TILT), 0, std::cos(MAV_THRUST_TILT));

/* The two components of an IMU-frame vector which lie in the rotor plane. */
inline UKF::Vector<2> mav_lateral(const UKF::Vector<3>& in) {
    return UKF::Vector<2>(MAV_LateralAxis1.dot(in), MAV_LateralAxis2.dot(in));
}

enum MAV_Keys {
    /* State vector fields. */
    Attitude,
    Velocity,
    GyroBias,

    /* Measurement vector fields. */
    RotorDrag,
    ZeroVelocity,
    RestGyroscope
};

/*
The state vector contains the following:
- Attitude as a quaternion (IMU frame to world frame; world +Z is up)
- Velocity (IMU frame, m/s)
- Gyroscope bias (IMU frame, rad/s)
*/
using MAV_StateVector = UKF::StateVector<
    UKF::Field<Attitude, UKF::Quaternion>,
    UKF::Field<Velocity, UKF::Vector<3>>,
    UKF::Field<GyroBias, UKF::Vector<3>>
>;

/*
The measurement vector is dynamic, because what can be measured depends on
what the vehicle is doing:
- Lateral specific force due to rotor drag (rotor plane, m/s^2), in flight
- Velocity pseudo-measurement, always zero (IMU frame, m/s), otherwise
- Gyroscope reading (IMU frame, rad/s), at rest
*/
using MAV_MeasurementVector = UKF::DynamicMeasurementVector<
    UKF::Field<RotorDrag, UKF::Vector<2>>,
    UKF::Field<ZeroVelocity, UKF::Vector<3>>,
    UKF::Field<RestGyroscope, UKF::Vector<3>>
>;

/* The process model input: the IMU readings to integrate over one step. */
struct MAV_Input {
    UKF::Vector<3> gyroscope;       /* Angular velocity (IMU frame, rad/s) */
    UKF::Vector<3> accelerometer;   /* Specific force (IMU frame, m/s^2) */
    bool powered;                   /* True if the rotor model applies */
};

namespace UKF {

/* MAV process model. */
template <> template <>
inline MAV_StateVector MAV_StateVector::derivative<MAV_Input>(const MAV_Input& input) const {
    MAV_StateVector output;

    /* Calculate change in attitude from the bias-corrected angular velocity. */
    UKF::Vector<3> angular_velocity = input.gyroscope - get_field<GyroBias>();
    UKF::Quaternion omega_q;
    omega_q.vec() = angular_velocity * 0.5;
    omega_q.w() = 0;
    output.set_field<Attitude>(get_field<Attitude>() * omega_q);

    /*
    Work out the specific force on the vehicle. Along the thrust axis it is
    whatever the accelerometer measures. Across it, in powered flight, it is
    rotor drag, which follows from the velocity; the lateral accelerometer
    reading is then used as a measurement instead (see below).
    */
    UKF::Vector<3> velocity = get_field<Velocity>();
    UKF::Vector<3> specific_force = input.accelerometer;
    if (input.powered) {
        specific_force =
            MAV_ThrustAxis * MAV_ThrustAxis.dot(input.accelerometer) -
            MAV_DRAG_RATE * (velocity - MAV_ThrustAxis * MAV_ThrustAxis.dot(velocity));
    }

    /*
    Calculate change in velocity: specific force plus gravity. The velocity
    is expressed in the IMU frame, which is rotating; hence the last term.
    */
    output.set_field<Velocity>(
        specific_force + get_field<Attitude>().conjugate() * UKF::Vector<3>(0, 0, -G_ACCEL) -
        angular_velocity.cross(velocity));

    /* Assume constant gyroscope bias; the process noise lets it wander. */
    output.set_field<GyroBias>(UKF::Vector<3>(0, 0, 0));

    return output;
}

/* MAV measurement model. In flight, the lateral accelerometer axes measure rotor drag. */
template <> template <>
inline UKF::Vector<2> MAV_MeasurementVector::expected_measurement
<MAV_StateVector, RotorDrag>(const MAV_StateVector& state) {
    return -MAV_DRAG_RATE * mav_lateral(state.get_field<Velocity>());
}

/* The velocity is "measured" as zero, so the expected measurement is the state. */
template <> template <>
inline UKF::Vector<3> MAV_MeasurementVector::expected_measurement
<MAV_StateVector, ZeroVelocity>(const MAV_StateVector& state) {
    return state.get_field<Velocity>();
}

/* At rest the true angular velocity is zero, leaving only the bias. */
template <> template <>
inline UKF::Vector<3> MAV_MeasurementVector::expected_measurement
<MAV_StateVector, RestGyroscope>(const MAV_StateVector& state) {
    return state.get_field<GyroBias>();
}

}

using MAV_Filter = UKF::SquareRootCore<
    MAV_StateVector,
    MAV_MeasurementVector,
    UKF::IntegratorRK4
>;

/*
Filter tuning, in one place. Sensor noise is given as a density: the
standard deviation which builds up over one second. The estimator scales it
to the actual time step, so the same values work at any sample rate.
*/
struct MAV_Tuning {
    /* Initial uncertainty (standard deviations). */
    real_t initial_tilt = 0.05;                         /* rad */
    real_t initial_velocity = 0.1;                      /* m/s */
    real_t initial_gyro_bias = RADIANS(5.0);            /* rad/s */

    /* Process noise. */
    real_t gyroscope_noise = 0.003;                     /* rad/sqrt(s) */
    real_t gyroscope_scale_error = 0.002;               /* rad/sqrt(s) per rad/s of rotation */
    real_t accelerometer_noise = 0.05;                  /* m/s/sqrt(s) */
    real_t rotor_model_noise = 0.2;                     /* m/s/sqrt(s), in powered flight */
    real_t gyro_bias_walk = RADIANS(0.01);              /* rad/s/sqrt(s) */

    /* Measurement noise. */
    real_t rotor_drag_noise = 0.067;                    /* m/s^2 * sqrt(s) */
    real_t zero_velocity_noise = 0.09;                  /* m/s * sqrt(s) */
    real_t rest_gyroscope_noise = RADIANS(0.09);        /* rad/s * sqrt(s) */

    /* Largest heading uncertainty kept by the filter. */
    real_t heading_sigma_limit = 0.3;                   /* rad */

    /*
    When is the vehicle at rest? The rest detector is given the IMU signals
    with everything above the cutoff frequency filtered out.
    */
    real_t rest_cutoff_frequency = 2.0;                 /* Hz */
    RestDetector::Limits rest = {
        0.5,                                            /* averaging_time (s) */
        RADIANS(2.0),                                   /* gyroscope_change (rad/s) */
        0.5,                                            /* accelerometer_change (m/s^2) */
        0.7 * G_ACCEL,                                  /* accelerometer_minimum (m/s^2) */
        1.3 * G_ACCEL,                                  /* accelerometer_maximum (m/s^2) */
        RADIANS(10.0),                                  /* gyroscope_bias (rad/s) */
        1.5                                             /* minimum_time (s) */
    };

    /* When are the rotors turning? */
    RotorDetector::Limits rotors = {
        20.0,                                           /* cutoff_frequency (Hz) */
        0.5,                                            /* averaging_time (s) */
        0.8,                                            /* on_vibration (m/s^2 RMS) */
        0.5                                             /* off_vibration (m/s^2 RMS) */
    };

    /* The rotor model is trusted while the thrust is in this range. */
    real_t minimum_thrust = 0.5 * G_ACCEL;              /* m/s^2 */
    real_t maximum_thrust = 2.0 * G_ACCEL;              /* m/s^2 */

    /*
    Set to false to ignore the rotor model altogether. This is the control
    which shows how much of the accuracy is due to the model.
    */
    bool use_rotor_model = true;
};

class MAV_Estimator {
public:
    explicit MAV_Estimator(const MAV_Tuning& tuning = MAV_Tuning()) :
        tuning(tuning),
        rest_detector(tuning.rest),
        rotor_detector(tuning.rotors),
        gyroscope_filter(tuning.rest_cutoff_frequency),
        accelerometer_filter(tuning.rest_cutoff_frequency) {}

    /*
    Start from the first sample. At rest the accelerometer points straight
    up, which fixes the initial tilt; the initial heading is arbitrary.
    Returns false if the accelerometer reading is zero.
    */
    bool initialise(const ImuSample& sample) {
        if (sample.accelerometer.norm() < 1e-6) {
            return false;
        }

        /* Initialise state vector and covariance. */
        filter.state.set_field<Attitude>(UKF::Quaternion::FromTwoVectors(
            sample.accelerometer, UKF::Vector<3>(0, 0, 1)));
        filter.state.set_field<Velocity>(UKF::Vector<3>(0, 0, 0));
        filter.state.set_field<GyroBias>(UKF::Vector<3>(0, 0, 0));

        filter.root_covariance = MAV_StateVector::CovarianceMatrix::Zero();
        filter.root_covariance.diagonal() <<
            tuning.initial_tilt * UKF::Vector<3>::Ones(),
            tuning.initial_velocity * UKF::Vector<3>::Ones(),
            tuning.initial_gyro_bias * UKF::Vector<3>::Ones();

        update_detectors(sample, 0);

        return true;
    }

    /*
    Advance the estimate from the previous sample to this one. Returns false
    if a filter step fails.
    */
    bool update(const ImuSample& previous, const ImuSample& sample) {
        real_t dt = sample.time - previous.time;

        /*
        Decide which model applies. The rotor model needs the rotors to be
        running and producing a plausible amount of thrust.
        */
        update_detectors(sample, dt);
        real_t thrust = MAV_ThrustAxis.dot(sample.accelerometer);
        powered = tuning.use_rotor_model && rotors_running &&
            thrust > tuning.minimum_thrust && thrust < tuning.maximum_thrust;

        /* Integrate the average of the two IMU readings over the step. */
        MAV_Input input;
        input.gyroscope = (previous.gyroscope + sample.gyroscope) * 0.5;
        input.accelerometer = (previous.accelerometer + sample.accelerometer) * 0.5;
        input.powered = powered;

        /*
        Set the process noise which builds up over this step. A gyroscope's
        scale error grows with the rate of rotation, so fast turns are
        trusted a little less. The rotor model is a rougher guide to the
        velocity than the accelerometer reading which it replaces.
        */
        real_t angular_velocity = (input.gyroscope - filter.state.get_field<GyroBias>()).norm();
        filter.process_noise_root_covariance = MAV_StateVector::CovarianceMatrix::Zero();
        filter.process_noise_root_covariance.diagonal() <<
            std::hypot(tuning.gyroscope_noise, tuning.gyroscope_scale_error * angular_velocity) *
                UKF::Vector<3>::Ones(),
            (powered ? tuning.rotor_model_noise : tuning.accelerometer_noise) * UKF::Vector<3>::Ones(),
            tuning.gyro_bias_walk * UKF::Vector<3>::Ones();
        filter.process_noise_root_covariance *= std::sqrt(dt);

        /*
        Set the measurement noise. Each measurement must count for less at a
        higher sample rate, where there are more of them every second.
        */
        filter.measurement_root_covariance.set_field<RotorDrag>(
            tuning.rotor_drag_noise / std::sqrt(dt) * UKF::Vector<2>::Ones());
        filter.measurement_root_covariance.set_field<ZeroVelocity>(
            tuning.zero_velocity_noise / std::sqrt(dt) * UKF::Vector<3>::Ones());
        filter.measurement_root_covariance.set_field<RestGyroscope>(
            tuning.rest_gyroscope_noise / std::sqrt(dt) * UKF::Vector<3>::Ones());

        /* Set only the measurements which apply at this step. */
        MAV_MeasurementVector measurement;
        if (powered) {
            measurement.set_field<RotorDrag>(mav_lateral(sample.accelerometer));
        } else {
            measurement.set_field<ZeroVelocity>(UKF::Vector<3>(0, 0, 0));
        }
        if (at_rest) {
            measurement.set_field<RestGyroscope>(sample.gyroscope);
        }

        if (!filter.a_priori_step(dt, input) ||
                !filter.innovation_step(measurement) ||
                !filter.a_posteriori_step()) {
            return false;
        }

        /*
        Nothing measured here depends on heading, so the filter's
        uncertainty about it can only grow. It must not grow without limit:
        the unscented transform spreads its sigma points in proportion to
        the uncertainty, and attitude errors approaching half a turn cannot
        be represented. The attitude error is expressed in the world frame,
        so the third row of the root covariance belongs to the heading
        error; scaling it back keeps the heading error's correlation with
        every other state, and leaves the rest of the covariance untouched.
        */
        real_t heading_sigma = filter.root_covariance.row(2).norm();
        if (heading_sigma > tuning.heading_sigma_limit) {
            filter.root_covariance.row(2) *= tuning.heading_sigma_limit / heading_sigma;
        }

        return true;
    }

    /* The estimated attitude (IMU frame to world frame). */
    UKF::Quaternion attitude() const { return filter.state.get_field<Attitude>().normalized(); }

    /* Other estimates and decisions, for diagnostics. */
    UKF::Vector<3> velocity() const { return filter.state.get_field<Velocity>(); }
    UKF::Vector<3> gyro_bias() const { return filter.state.get_field<GyroBias>(); }
    bool resting() const { return at_rest; }
    bool rotor_model_active() const { return powered; }
    real_t rotor_vibration() const { return rotor_detector.vibration_rms(); }

    /*
    The filter itself, read-only: its state, root covariance, noise settings
    and latest innovation. The explainer (see explainer/) uses this to show
    how sure the filter is; the benchmark does not need it.
    */
    const MAV_Filter& ukf() const { return filter; }

private:
    /*
    Work out whether the rotors are running, and whether the vehicle is at
    rest. The gyroscope can only measure its own bias with the rotors off.
    */
    void update_detectors(const ImuSample& sample, real_t dt) {
        rotors_running = rotor_detector.update(sample.accelerometer, dt);
        at_rest = rest_detector.update(
            gyroscope_filter.update(sample.gyroscope, dt),
            accelerometer_filter.update(sample.accelerometer, dt), dt) && !rotors_running;
    }

    MAV_Tuning tuning;
    MAV_Filter filter;
    RestDetector rest_detector;
    RotorDetector rotor_detector;
    LowPassFilter gyroscope_filter;
    LowPassFilter accelerometer_filter;
    bool rotors_running = false;
    bool powered = false;
    bool at_rest = false;
};

#endif
