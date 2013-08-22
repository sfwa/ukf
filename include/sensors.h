/*
Copyright (C) 2013 Daniel Dyer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef SENSORS_H
#define SENSORS_H

#include "types.h"
#include "state.h"

/*
Sensor model base class. The public interface consists of the following:
    - The collate() method, which returns the actual measurement vector,
      containing only the sensor data that has been supplied so far (since
      the last initialise() call).
    - The predict() method, which returns the predicted measurement vector
      for the supplied state vector, including only those sensor values which
      have been supplied already.
    - A bunch of methods specific to the sensor model implementation which
      set the actual sensor readings (eg. set_accelerometer()).
*/
class SensorModel {
public:
    virtual ~SensorModel();
    virtual size_t size() const = 0;
    virtual MeasurementVector collate() const = 0;
    virtual MeasurementVector predict(const State &in) const = 0;
    virtual MeasurementVector get_covariance() const = 0;
    virtual void set_covariance(const MeasurementVector &in) = 0;
    virtual MeasurementVector calculate_mean(
        const MatrixXr &in,
        const VectorXr &weights);
    virtual MatrixXr calculate_deltas(
        const MatrixXr &in,
        const MeasurementVector &mean);
};

/*
Sensor model including accelerometer, magnetometer, gyro, gps, pitot,
barometric pressure and derivative of barometric pressure.
Contents of returned vector are as follows (when all sensor data is
populated):
    - Accelerometer (3-vector, m/s^2)
    - Gyroscope (3-vector, rad/s)
    - Magnetometer (3-vector, µT, body frame)
    - GPS position (3-vector, latitude (rad), longitude (rad), altitude (m))
    - GPS NED velocity (3-vector, m/s)
    - Pitot true airspeed (scalar, m/s)
    - Barometer altitude AMSL (scalar, m)

NOTE: The 'field' vector in the magnetometer struct should be set to the
magnetic field vector at the current position, in the ECEF reference frame.
*/
class IOBoardModel: public SensorModel {
    /* Sensor data and parameters. */
    struct {
        Quaternionr orientation;
        Vector3r data;
        Vector3r offset;
    } accelerometer;

    struct {
        Quaternionr orientation;
        Vector3r data;
    } gyroscope;

    struct {
        Quaternionr orientation;
        Vector3r data;
        Vector3r field;
    } magnetometer;

    Vector3r gps_position;
    Vector3r gps_velocity;
    real_t pitot_tas;
    real_t barometer_amsl;

    /*
    Stores sensor noise covariance, in the following order:
        - Accelerometer (3-vector)
        - Gyroscope (3-vector)
        - Magnetometer (3-vector)
        - GPS Position (3-vector)
        - GPS Velocity (3-vector)
        - Pitot TAS (scalar)
        - Barometer AMSL (scalar)
    */
    MeasurementVector covariance;

    /* Sensor flags. */
    struct {
        bool accelerometer : 1;
        bool gyroscope : 1;
        bool magnetometer : 1;
        bool gps_position : 1;
        bool gps_velocity : 1;
        bool pitot_tas : 1;
        bool barometer_amsl : 1;
    } flags;

public:
    IOBoardModel(Quaternionr accelerometer_orientation,
                 Vector3r accelerometer_offset,
                 Quaternionr gyroscope_orientation,
                 Quaternionr magnetometer_orientation,
                 Vector3r magnetic_field) {
        accelerometer.orientation = accelerometer_orientation;
        accelerometer.offset = accelerometer_offset;
        gyroscope.orientation = gyroscope_orientation;
        magnetometer.orientation = magnetometer_orientation;
        magnetometer.field = magnetic_field;
        covariance = MeasurementVector::Constant(UKF_MEASUREMENT_DIM, 1, 10e100);
        clear();
    }
    void clear() { memset(&flags, 0, sizeof(flags)); }
    size_t size() const;
    MeasurementVector collate() const;
    MeasurementVector predict(const State &in) const;
    MeasurementVector get_covariance() const;
    void set_covariance(const MeasurementVector &in) {
        covariance = in;
    }
    MeasurementVector calculate_mean(
        const MatrixXr &in,
        const VectorXr &weights);
    MeasurementVector calculate_innovation(const MeasurementVector &in);
    MatrixXr calculate_deltas(
        const MatrixXr &in,
        const MeasurementVector &mean);
    void set_accelerometer(Vector3r data) {
        accelerometer.data = data;
        flags.accelerometer = true;
    }
    void set_gyroscope(Vector3r data) {
        gyroscope.data = data;
        flags.gyroscope = true;
    }
    void set_magnetometer(Vector3r data) {
        magnetometer.data = data;
        flags.magnetometer = true;
    }

    void set_pitot_tas(real_t tas) {
        pitot_tas = tas;
        flags.pitot_tas = true;
    }
    void set_barometer_amsl(real_t amsl) {
        barometer_amsl = amsl;
        flags.barometer_amsl = true;
    }
    void set_magnetic_field(Vector3r field) {
        magnetometer.field = field;
    }
    void set_gps_position(Vector3r data) {
        gps_position = data;
        flags.gps_position = true;
    }
    void set_gps_velocity(Vector3r data) {
        gps_velocity = data;
        flags.gps_velocity = true;
    }
};

#endif
