/*
Detecting whether a multirotor's rotors are turning, from the vibration they
cause in the accelerometer signal. Like the checks in sensor_checks.h, this
is not part of the UKF: it only decides which model and measurements apply.
*/

#ifndef ROTOR_VIBRATION_H
#define ROTOR_VIBRATION_H

#include <cmath>
#include "UKF/Types.h"
#include "sensor_checks.h"

/*
Second-order (Butterworth) low-pass filter for a three-axis signal, which
removes vibration above the cutoff frequency. The update is the exact
solution of the filter's differential equation over one time step, so it
behaves identically at any sample rate.
*/
class LowPassFilter {
public:
    explicit LowPassFilter(real_t cutoff_frequency) :
        time_constant(std::sqrt(real_t(2.0)) / (real_t(2.0 * 3.14159265358979323846) * cutoff_frequency)) {}

    UKF::Vector<3> update(const UKF::Vector<3>& sample, real_t dt) {
        if (!started) {
            value = sample;
            started = true;
            return value;
        }

        real_t step = dt / time_constant;
        real_t decay = std::exp(-step), cosine = std::cos(step), sine = std::sin(step);
        UKF::Vector<3> error = value - sample;
        value = sample + decay * ((cosine + sine) * error + time_constant * sine * rate);
        rate = decay * ((real_t(-2.0) * sine / time_constant) * error + (cosine - sine) * rate);

        return value;
    }

private:
    real_t time_constant;
    UKF::Vector<3> value = UKF::Vector<3>::Zero();
    UKF::Vector<3> rate = UKF::Vector<3>::Zero();
    bool started = false;
};

/*
The rotors are running when the accelerometer vibrates strongly. Vibration
is what is left of the signal after its slowly varying part is removed; its
RMS level switches the detector on above one threshold and off below a lower
one, so that it does not flicker in between.
*/
class RotorDetector {
public:
    struct Limits {
        real_t cutoff_frequency;    /* Vibration is the signal above this frequency (Hz). */
        real_t averaging_time;      /* Time constant of the RMS average (s). */
        real_t on_vibration;        /* RMS vibration above which the rotors are on (m/s^2). */
        real_t off_vibration;       /* RMS vibration below which the rotors are off (m/s^2). */
    };

    explicit RotorDetector(const Limits& limits) : limits(limits), low_pass(limits.cutoff_frequency) {}

    /* Feed one accelerometer sample; returns true while the rotors are running. */
    bool update(const UKF::Vector<3>& accelerometer, real_t dt) {
        UKF::Vector<3> vibration = accelerometer - low_pass.update(accelerometer, dt);
        smooth(mean_square, vibration.squaredNorm(), dt, limits.averaging_time);

        if (vibration_rms() > limits.on_vibration) {
            running = true;
        } else if (vibration_rms() < limits.off_vibration) {
            running = false;
        }

        return running;
    }

    real_t vibration_rms() const { return std::sqrt(mean_square); }

private:
    Limits limits;
    LowPassFilter low_pass;
    real_t mean_square = 0.0;
    bool running = false;
};

#endif
