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

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "types.h"
#include "state.h"
#include "integrator.h"
#include "dynamics.h"
#include "interface.h"
#include "sensors.h"
#include "ukf.h"

static IOBoardModel model = IOBoardModel(
    Quaternionr(1, 0, 0, 0),
    Vector3r(0, 0, 0),
    Quaternionr(1, 0, 0, 0),
    Quaternionr(1, 0, 0, 0),
    Vector3r(1, 0, 0));
static UnscentedKalmanFilter ukf = UnscentedKalmanFilter(model);
static CentripetalModel centripetal_model = CentripetalModel();
static FixedWingFlightDynamicsModel fixed_wing_model = FixedWingFlightDynamicsModel();

void set_position(real_t lat, real_t lon, real_t alt) {
    State temp = ukf.get_state();
    temp.position() << lat, lon, alt;
    ukf.set_state(temp);
}

void set_velocity(real_t x, real_t y, real_t z) {
    State temp = ukf.get_state();
    temp.velocity() << x, y, z;
    ukf.set_state(temp);
}

void set_acceleration(real_t x, real_t y, real_t z) {
    State temp = ukf.get_state();
    temp.acceleration() << x, y, z;
    ukf.set_state(temp);
}

void set_attitude(real_t w, real_t x, real_t y, real_t z) {
    State temp = ukf.get_state();
    temp.attitude() << x, y, z, w;
    ukf.set_state(temp);
}

void set_angular_velocity(real_t x, real_t y, real_t z) {
    State temp = ukf.get_state();
    temp.angular_velocity() << x, y, z;
    ukf.set_state(temp);
}

void set_angular_acceleration(real_t x, real_t y, real_t z) {
    State temp = ukf.get_state();
    temp.angular_acceleration() << x, y, z;
    ukf.set_state(temp);
}

void set_wind_velocity(real_t x, real_t y, real_t z) {
    State temp = ukf.get_state();
    temp.wind_velocity() << x, y, z;
    ukf.set_state(temp);
}

void set_gyro_bias(real_t x, real_t y, real_t z) {
    State temp = ukf.get_state();
    temp.gyro_bias() << x, y, z;
    ukf.set_state(temp);
}

void get_state(struct state *in) {
    State current = ukf.get_state();
    in->position[0] = current.position()[0];
    in->position[1] = current.position()[1];
    in->position[2] = current.position()[2];
    in->velocity[0] = current.velocity()[0];
    in->velocity[1] = current.velocity()[1];
    in->velocity[2] = current.velocity()[2];
    in->acceleration[0] = current.acceleration()[0];
    in->acceleration[1] = current.acceleration()[1];
    in->acceleration[2] = current.acceleration()[2];
    in->attitude[0] = current.attitude()[3];
    in->attitude[1] = current.attitude()[0];
    in->attitude[2] = current.attitude()[1];
    in->attitude[3] = current.attitude()[2];
    in->angular_velocity[0] = current.angular_velocity()[0];
    in->angular_velocity[1] = current.angular_velocity()[1];
    in->angular_velocity[2] = current.angular_velocity()[2];
    in->angular_acceleration[0] = current.angular_acceleration()[0];
    in->angular_acceleration[1] = current.angular_acceleration()[1];
    in->angular_acceleration[2] = current.angular_acceleration()[2];
    in->wind_velocity[0] = current.wind_velocity()[0];
    in->wind_velocity[1] = current.wind_velocity()[1];
    in->wind_velocity[2] = current.wind_velocity()[2];
    in->gyro_bias[0] = current.gyro_bias()[0];
    in->gyro_bias[1] = current.gyro_bias()[1];
    in->gyro_bias[2] = current.gyro_bias()[2];
}

void get_state_covariance(real_t state_covariance[24*24]) {
    Eigen::Map<StateVectorCovariance> covariance_map(state_covariance);
    covariance_map = ukf.get_state_covariance();
}

void sensor_clear() {
    model.clear();
}

void sensor_set_accelerometer(real_t x, real_t y, real_t z) {
    model.set_accelerometer(Vector3r(x, y, z));
}

void sensor_set_gyroscope(real_t x, real_t y, real_t z) {
    model.set_gyroscope(Vector3r(x, y, z));
}

void sensor_set_magnetometer(real_t x, real_t y, real_t z) {
    model.set_magnetometer(Vector3r(x, y, z));
}

void sensor_set_gps_position(real_t lat, real_t lon, real_t alt) {
    model.set_gps_position(Vector3r(lat, lon, alt));
}

void sensor_set_gps_velocity(real_t x, real_t y, real_t z) {
    model.set_gps_velocity(Vector3r(x, y, z));
}

void sensor_set_pitot_tas(real_t tas) {
    model.set_pitot_tas(tas);
}

void sensor_set_barometer_amsl(real_t amsl) {
    model.set_barometer_amsl(amsl);
}

void sensor_set_covariance(real_t sensor_covariance[17]) {
    Eigen::Map< Eigen::Matrix<real_t, 17, 1> > covariance_map =
        Eigen::Map< Eigen::Matrix<real_t, 17, 1> >(sensor_covariance);
    MeasurementVector covariance = covariance_map;
    model.set_covariance(covariance);
}

void ukf_set_params(struct ioboard_params *in) {
    model = IOBoardModel(
        Quaternionr(
            in->accel_orientation[0],
            in->accel_orientation[1],
            in->accel_orientation[2],
            in->accel_orientation[3]),
        Vector3r(
            in->accel_offset[0],
            in->accel_offset[1],
            in->accel_offset[2]),
        Quaternionr(
            in->gyro_orientation[0],
            in->gyro_orientation[1],
            in->gyro_orientation[2],
            in->gyro_orientation[3]),
        Quaternionr(
            in->mag_orientation[0],
            in->mag_orientation[1],
            in->mag_orientation[2],
            in->mag_orientation[3]),
        Vector3r(
            in->mag_field[0],
            in->mag_field[1],
            in->mag_field[2]));
}

void ukf_set_field(real_t x, real_t y, real_t z) {
    model.set_magnetic_field(Vector3r(x, y, z));
}

void ukf_choose_dynamics(enum model_types t) {
    switch(t) {
        case MODEL_NONE:
            ukf.set_dynamics_model((DynamicsModel *)NULL);
            break;
        case MODEL_CENTRIPETAL:
            ukf.set_dynamics_model(&centripetal_model);
            break;
        case MODEL_FIXED_WING:
            ukf.set_dynamics_model(&fixed_wing_model);
            break;
    }
}

void ukf_iterate(float dt, real_t control_vector[4]) {
    ukf.iterate(dt, Eigen::Matrix<real_t, 4, 1>(control_vector));
}

void ukf_set_process_noise(real_t process_noise_covariance[24]) {
    Eigen::Map< Eigen::Matrix<real_t, 24, 1> > covariance_map =
        Eigen::Map< Eigen::Matrix<real_t, 24, 1> >(process_noise_covariance);
    ProcessCovariance covariance = covariance_map;
    ukf.set_process_noise(covariance);
}

void fixedwingdynamics_set_mass(real_t mass) {
    fixed_wing_model.set_mass(mass);
}

void fixedwingdynamics_set_inertia_tensor(real_t inertia_tensor[9]) {
    fixed_wing_model.set_inertia_tensor(
        Eigen::Matrix<real_t, 3, 3>(inertia_tensor));
}

void fixedwingdynamics_set_prop_coeffs(real_t in_prop_area, real_t in_prop_cve){
    fixed_wing_model.set_prop_coeffs(in_prop_area, in_prop_cve);
}

void fixedwingdynamics_set_drag_coeffs(real_t coeffs[5]) {
    fixed_wing_model.set_drag_coeffs(Eigen::Matrix<real_t, 5, 1>(coeffs));
}

void fixedwingdynamics_set_lift_coeffs(real_t coeffs[5]) {
    fixed_wing_model.set_lift_coeffs(Eigen::Matrix<real_t, 5, 1>(coeffs));
}

void fixedwingdynamics_set_side_coeffs(real_t coeffs[12]) {
    fixed_wing_model.set_side_coeffs(Eigen::Matrix<real_t, 12, 1>(coeffs));
}

void fixedwingdynamics_set_pitch_moment_coeffs(real_t coeffs[6]) {
    fixed_wing_model.set_pitch_moment_coeffs(
        Eigen::Matrix<real_t, 6, 1>(coeffs));
}

void fixedwingdynamics_set_roll_moment_coeffs(real_t coeffs[5]) {
    fixed_wing_model.set_roll_moment_coeffs(
        Eigen::Matrix<real_t, 5, 1>(coeffs));
}

void fixedwingdynamics_set_yaw_moment_coeffs(real_t coeffs[6]) {
    fixed_wing_model.set_yaw_moment_coeffs(
        Eigen::Matrix<real_t, 6, 1>(coeffs));
}
