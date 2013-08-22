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

#include <cmath>

#include "types.h"
#include "sensors.h"
#include "state.h"
#include "debug.h"

SensorModel::~SensorModel() {}


/*
Takes the mean of a sigma point distribution given a column vector of weights.
*/
MeasurementVector SensorModel::calculate_mean(
const MatrixXr &in, const VectorXr &weights) {
    return in * weights;
}


/*
Calculates the difference between a sigma point distribution and the supplied
mean vector.
*/
MatrixXr SensorModel::calculate_deltas(
const MatrixXr &in,
const MeasurementVector &mean) {
    return in.colwise() - mean;
}


/* Calculates the current size of the measurement vector. */
size_t IOBoardModel::size() const {
    size_t i = 0;

    i += flags.accelerometer ? 6 : 0;
    i += flags.gyroscope ? 3 : 0;
    i += flags.magnetometer ? 3 : 0;
    i += flags.gps_position ? 3 : 0;
    i += flags.gps_velocity ? 3 : 0;
    i += flags.pitot_tas ? 1 : 0;
    i += flags.barometer_amsl ? 1 : 0;

    return i;
}


/*
Takes all populated sensor data and puts it into a vector for use by the UKF.
*/
MeasurementVector IOBoardModel::collate() const {
    MeasurementVector::Index max_size = (MeasurementVector::Index)size();
    MeasurementVector::Index i = 0;
    if(flags.accelerometer) {
        max_size -= 3;
    }
    MeasurementVector measurement(max_size);

    if(flags.accelerometer) {
        measurement.segment<3>(i) << accelerometer.data;
        i += 3;
    }

    if(flags.gyroscope) {
        measurement.segment<3>(i) << gyroscope.data;
        i += 3;
    }

    if(flags.magnetometer) {
        measurement.segment<3>(i) << magnetometer.data;
        i += 3;
    }

    if(flags.gps_position) {
        measurement.segment<3>(i) << gps_position;
        i += 3;
    }

    if(flags.gps_velocity) {
        measurement.segment<3>(i) << gps_velocity;
        i += 3;
    }

    if(flags.pitot_tas) {
        measurement[i] = pitot_tas;
        i++;
    }

    if(flags.barometer_amsl) {
        measurement[i] = barometer_amsl;
        i++;
    }

    return measurement;
}


/*
Takes a state vector and creates a predicted measurement vector containing
only the sensor values which have been supplied.
*/
MeasurementVector IOBoardModel::predict(const State &in) const {
    MeasurementVector::Index max_size = (MeasurementVector::Index)size();
    MeasurementVector::Index i = 0;
    MeasurementVector predicted(max_size);
    Quaternionr attitude = Quaternionr(in.attitude());

    AssertNormalized(attitude);
    AssertNormalized(accelerometer.orientation);
    AssertNormalized(gyroscope.orientation);
    AssertNormalized(magnetometer.orientation);

    /*
    The accelerometer reading is predicted by adding the expected gravity
    vector (transformed by orientation) to the state vector acceleration.
    Also, if the offset vector is non-zero, need to include the angular
    acceleration term and the centripetal acceleration term.
    */
    if(flags.accelerometer) {
        predicted.segment<3>(i) << accelerometer.orientation *
            (in.acceleration() +
             in.angular_acceleration().cross(accelerometer.offset) +
             in.angular_velocity().cross(in.angular_velocity().cross(
                accelerometer.offset)));
        predicted.segment<3>(i+3) << accelerometer.orientation *
            (attitude * Vector3r(0, 0, -G_ACCEL));
        i += 6;
    }

    /*
    The gyroscope reading is simply the angular velocity transformed by the
    sensor orientation.
    */
    if(flags.gyroscope) {
        predicted.segment<3>(i) << gyroscope.orientation *
            (in.angular_velocity() + in.gyro_bias());
        i += 3;
    }

    /*
    Predicted magnetometer measurement is easy, as long as the supplied
    magnetic field vector is in the NED frame of reference. Just need to
    convert the magnetic field vector into the body frame, and also adjust by
    the orientation vector.
    */
    if(flags.magnetometer) {
        predicted.segment<3>(i) << magnetometer.orientation *
            (attitude * magnetometer.field);
        i += 3;
    }

    /*
    GPS position prediction is the previous estimate's lat/lon/alt.
    */
    if(flags.gps_position) {
        predicted.segment<3>(i) << in.position();
        i += 3;
    }

    /*
    GPS velocity prediction is just the previous estimate's NED velocity.
    */
    if(flags.gps_velocity) {
        predicted.segment<3>(i) << in.velocity();
        i += 3;
    }

    /*
    The true airspeed is predicted by subtracting the wind velocity from the
    ECEF velocity, transforming it into body frame and then taking the +X
    component.
    */
    if(flags.pitot_tas) {
        Vector3r temp_3d = attitude * (in.velocity() - in.wind_velocity());
        predicted[i] = temp_3d[0];
        i++;
    }

    /*
    The barometric altitude prediction is just the height from the state
    vector position.
    */
    if(flags.barometer_amsl) {
        predicted[i] = in.position()[2];
        i++;
    }

    return predicted;
}

MeasurementVector IOBoardModel::get_covariance() const {
    MeasurementVector::Index max_size = (MeasurementVector::Index)size();
    MeasurementVector::Index i = 0;
    if(flags.accelerometer) {
        max_size -= 3;
    }
    MeasurementVector out_covariance(max_size);
    Vector3r temp_3d;

    if(flags.accelerometer) {
        out_covariance.segment<3>(i) << covariance.segment<3>(0);
        i += 3;
    }

    if(flags.gyroscope) {
        out_covariance.segment<3>(i) << covariance.segment<3>(3);
        i += 3;
    }

    if(flags.magnetometer) {
        out_covariance.segment<3>(i) << covariance.segment<3>(6);
        i += 3;
    }

    if(flags.gps_position) {
        out_covariance.segment<3>(i) << covariance.segment<3>(9);
        i += 3;
    }

    if(flags.gps_velocity) {
        out_covariance.segment<3>(i) << covariance.segment<3>(12);
        i += 3;
    }

    if(flags.pitot_tas) {
        out_covariance[i] = covariance[15];
        i++;
    }

    if(flags.barometer_amsl) {
        out_covariance[i] = covariance[16];
        i++;
    }

    return out_covariance;
}


/*
The mean is calculated differently because we first have to normalise both the
gravitational acceleration and the magnetometer value, then add the
gravitational acceleration to the kinematic acceleration to get the expected
accelerometer measurement.
*/
MeasurementVector IOBoardModel::calculate_mean(
const MatrixXr &in,
const VectorXr &weights) {
    MeasurementVector::Index max_size = (MeasurementVector::Index)size();
    MeasurementVector::Index i = 0, j = 0;
    if(flags.accelerometer) {
        max_size -= 3;
    }
    MeasurementVector mean(max_size);

    /* Calculate the mean. */
    MeasurementVector initial_mean = in * weights;

    if(flags.accelerometer) {
        /*
        Normalise the gravitational acceleration, and add the gravitational
        and kinematic accelerations.
        */
        mean.segment<3>(i) <<
            (initial_mean.segment<3>(j+3).normalized() * G_ACCEL) +
            initial_mean.segment<3>(j);

        i += 3;
        j += 6;
    }

    if(flags.gyroscope) {
        mean.segment<3>(i) << initial_mean.segment<3>(j);
        i += 3;
        j += 3;
    }

    if(flags.magnetometer) {
        /* Normalise the magnetometer. */
        mean.segment<3>(i) << initial_mean.segment<3>(j).normalized() *
            magnetometer.field.norm();
        i += 3;
        j += 3;
    }

    if(flags.gps_position) {
        mean.segment<3>(i) << initial_mean.segment<3>(j);
        i += 3;
        j += 3;
    }

    if(flags.gps_velocity) {
        mean.segment<3>(i) << initial_mean.segment<3>(j);
        i += 3;
        j += 3;
    }

    if(flags.pitot_tas) {
        mean[i] = initial_mean[j];
        i++;
        j++;
    }

    if(flags.barometer_amsl) {
        mean[i] = initial_mean[j];
        i++;
        j++;
    }

    return mean;
}

/*
This is specialised because we have to add the gravitational acceleration to
the kinematic acceleration in the sigma points before subtracting the mean
from them.
*/
MatrixXr IOBoardModel::calculate_deltas(
const MatrixXr &in,
const MeasurementVector &mean) {
    MatrixXr deltas;
    MeasurementVector::Index max_size = (MeasurementVector::Index)in.rows();

    if(flags.accelerometer) {
        max_size -= 3;
        deltas = MatrixXr(max_size, in.cols());

        deltas.block(0, 0, 3, in.cols()) =
            in.block(0, 0, 3, in.cols()) +
            in.block(3, 0, 3, in.cols());

        if(max_size > 3) {
            deltas.block(3, 0, max_size-3, in.cols()) =
                in.block(6, 0, max_size-3, in.cols());
        }

        deltas.colwise() -= mean;
    } else {
        deltas = in.colwise() - mean;
    }

    return deltas;
}
