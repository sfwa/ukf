/* Sensor input shared by the estimator and the host-side CSV tools. */
#ifndef IMU_SAMPLE_H
#define IMU_SAMPLE_H

#include "UKF/Types.h"

/* One row of sensor readings, all in the sensor's own (body) frame. */
struct ImuSample {
    /* Keep timestamps in double precision even when filter arithmetic is float. */
    double time;                    /* Seconds; must increase from row to row. */
    UKF::Vector<3> gyroscope;       /* Angular velocity (rad/s). */
    UKF::Vector<3> accelerometer;   /* Specific force (m/s^2): +g along "up" at rest. */
    UKF::Vector<3> magnetometer;    /* Magnetic field, in any consistent unit. */
};

#endif
