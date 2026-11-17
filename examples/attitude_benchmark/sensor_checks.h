/*
Two small checks on the raw sensor readings that tell the attitude filter
which measurements it may use at each step:

  - RestDetector: is the sensor lying still? If so, the gyroscope is reading
    nothing but its own bias.
  - MagneticFieldCheck: does the magnetic field look like the Earth's? If
    not, the magnetometer must not be used for heading.

Neither is part of the UKF. They only decide which fields are set in the
measurement vector.
*/

#ifndef SENSOR_CHECKS_H
#define SENSOR_CHECKS_H

#include <algorithm>
#include <cmath>
#include "UKF/Types.h"

/*
Move "average" towards "sample" as an exponential moving average with the
given time constant. Using the elapsed time rather than a fixed gain keeps
the behaviour the same at any sample rate.
*/
template <typename T>
void smooth(T& average, const T& sample, real_t dt, real_t time_constant) {
    average += (real_t(1.0) - std::exp(-dt / time_constant)) * (sample - average);
}

/*
The sensor is at rest once the gyroscope and accelerometer have both stayed
close to their own recent averages for a while. Two more conditions exclude
steady readings which are not rest. The accelerometer must be measuring
gravity: in free fall, or if its readings drop out, it measures nothing at
all. And the average angular velocity must be small: a slow, steady turn
would otherwise look just like a gyroscope bias.
*/
class RestDetector {
public:
    struct Limits {
        real_t averaging_time;          /* Time constant of the running averages (s). */
        real_t gyroscope_change;        /* Largest deviation from the average (rad/s). */
        real_t accelerometer_change;    /* Largest deviation from the average (m/s^2). */
        real_t accelerometer_minimum;   /* Smallest accelerometer magnitude at rest (m/s^2). */
        real_t accelerometer_maximum;   /* Largest accelerometer magnitude at rest (m/s^2). */
        real_t gyroscope_bias;          /* Largest average angular velocity (rad/s). */
        real_t minimum_time;            /* How long the readings must stay quiet (s). */
    };

    explicit RestDetector(const Limits& limits) : limits(limits) {}

    /* Feed one sample; returns true while the sensor is at rest. */
    bool update(const UKF::Vector<3>& gyroscope, const UKF::Vector<3>& accelerometer, real_t dt) {
        if (!started) {
            average_gyroscope = gyroscope;
            average_accelerometer = accelerometer;
            started = true;
        }
        smooth(average_gyroscope, gyroscope, dt, limits.averaging_time);
        smooth(average_accelerometer, accelerometer, dt, limits.averaging_time);

        const real_t magnitude = accelerometer.norm();
        const bool quiet =
            (gyroscope - average_gyroscope).norm() <= limits.gyroscope_change &&
            (accelerometer - average_accelerometer).norm() <= limits.accelerometer_change &&
            magnitude >= limits.accelerometer_minimum &&
            magnitude <= limits.accelerometer_maximum &&
            average_gyroscope.cwiseAbs().maxCoeff() <= limits.gyroscope_bias;
        quiet_time = quiet ? quiet_time + dt : real_t(0.0);

        return quiet_time >= limits.minimum_time;
    }

private:
    Limits limits;
    UKF::Vector<3> average_gyroscope = UKF::Vector<3>::Zero();
    UKF::Vector<3> average_accelerometer = UKF::Vector<3>::Zero();
    real_t quiet_time = 0.0;
    bool started = false;
};

/*
The Earth's magnetic field has the same strength and the same dip angle
(its slope below the horizontal) wherever the sensor is turned. Steel
furniture, magnets and electronics change both, so a reading whose strength
or dip has strayed from the remembered reference is treated as disturbed.

The reference follows the measured field slowly while it is trusted. If the
sensor is moved into a different but steady field (another room, say), that
field is adopted as the new reference once it has stayed consistent through
enough rotation; a field that only looks steady because the sensor is not
turning proves nothing.
*/
class MagneticFieldCheck {
public:
    struct Limits {
        real_t strength_change;     /* Allowed change in field strength (fraction). */
        real_t dip_change;          /* Allowed change in dip angle (rad). */
        real_t smoothing_time;      /* Smoothing of the measured strength and dip (s). */
        real_t settle_time;         /* Undisturbed time before readings are trusted again (s). */
        real_t reference_time;      /* Time constant with which the reference adapts (s). */
        real_t new_field_time;      /* Rotation time needed to adopt a new field (s). */
        real_t new_field_rotation;  /* Angular velocity that counts as rotation (rad/s). */
    };

    explicit MagneticFieldCheck(const Limits& limits) : limits(limits) {}

    /*
    Feed one magnetometer reading, already rotated into the world frame, and
    the current angular velocity magnitude. Returns true if the reading is
    disturbed and should not be used.
    */
    bool update(const UKF::Vector<3>& field, real_t dt, real_t angular_velocity) {
        const real_t measured_strength = field.norm();
        const real_t measured_dip = std::asin(std::clamp(field.z() / measured_strength, real_t(-1.0), real_t(1.0)));
        if (!started) {
            strength = reference_strength = candidate_strength = measured_strength;
            dip = reference_dip = candidate_dip = measured_dip;
            started = true;
            return false;
        }
        smooth(strength, measured_strength, dt, limits.smoothing_time);
        smooth(dip, measured_dip, dt, limits.smoothing_time);

        /* Compare with the reference; trust resumes after a quiet spell. */
        if (matches(reference_strength, reference_dip)) {
            undisturbed_time += dt;
            if (undisturbed_time >= limits.settle_time) {
                disturbed = false;
                smooth(reference_strength, strength, dt, limits.reference_time);
                smooth(reference_dip, dip, dt, limits.reference_time);
            }
        } else {
            undisturbed_time = 0.0;
            disturbed = true;
        }

        /* Track a possible new field, and adopt it after enough rotation. */
        if (matches(candidate_strength, candidate_dip)) {
            if (angular_velocity >= limits.new_field_rotation) {
                candidate_time += dt;
            }
            smooth(candidate_strength, strength, dt, limits.reference_time);
            smooth(candidate_dip, dip, dt, limits.reference_time);
        } else {
            candidate_strength = strength;
            candidate_dip = dip;
            candidate_time = 0.0;
        }
        if (candidate_time >= limits.new_field_time) {
            reference_strength = candidate_strength;
            reference_dip = candidate_dip;
            undisturbed_time = limits.settle_time;
            candidate_time = 0.0;
            disturbed = false;
        }

        return disturbed;
    }

private:
    bool matches(real_t other_strength, real_t other_dip) const {
        return std::abs(strength - other_strength) <= limits.strength_change * other_strength &&
               std::abs(dip - other_dip) <= limits.dip_change;
    }

    Limits limits;
    real_t strength = 0.0, dip = 0.0;
    real_t reference_strength = 0.0, reference_dip = 0.0;
    real_t candidate_strength = 0.0, candidate_dip = 0.0;
    real_t undisturbed_time = 0.0, candidate_time = 0.0;
    bool disturbed = false;
    bool started = false;
};

#endif
