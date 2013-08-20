import os
from ctypes import *

class IOBoardParams(Structure):
    _fields_ = [
        ("accel_orientation", c_double * 4),
        ("accel_offset", c_double * 3),
        ("gyro_orientation", c_double * 4),
        ("mag_orientation", c_double * 4),
        ("mag_field", c_double * 3),
    ]

class State(Structure):
    _fields_ = [
        ("position", c_double * 3),
        ("velocity", c_double * 3),
        ("acceleration", c_double * 3),
        ("attitude", c_double * 4),
        ("angular_velocity", c_double * 3),
        ("angular_acceleration", c_double * 3),
        ("wind_velocity", c_double * 3),
        ("gyro_bias", c_double * 3)
    ]

    def __repr__(self):
        fields = {
            "position": tuple(self.position),
            "velocity": tuple(self.velocity),
            "acceleration": tuple(self.acceleration),
            "attitude": tuple(self.attitude),
            "angular_velocity": tuple(self.angular_velocity),
            "angular_acceleration": tuple(self.angular_acceleration),
            "wind_velocity": tuple(self.wind_velocity),
            "gyro_bias": tuple(self.gyro_bias)
        }

        return str(fields)

cukf = cdll.LoadLibrary(
    os.path.join(os.path.dirname(__file__), "c", "libcukf.dylib"))

cukf.set_position.argtypes = [c_double, c_double, c_double]
cukf.set_position.restype = None

cukf.set_velocity.argtypes = [c_double, c_double, c_double]
cukf.set_velocity.restype = None

cukf.set_acceleration.argtypes = [c_double, c_double, c_double]
cukf.set_acceleration.restype = None

cukf.set_attitude.argtypes = [c_double, c_double, c_double, c_double]
cukf.set_attitude.restype = None

cukf.set_angular_velocity.argtypes = [c_double, c_double, c_double]
cukf.set_angular_velocity.restype = None

cukf.set_angular_acceleration.argtypes = [c_double, c_double, c_double]
cukf.set_angular_acceleration.restype = None

cukf.set_wind_velocity.argtypes = [c_double, c_double, c_double]
cukf.set_wind_velocity.restype = None

cukf.set_gyro_bias.argtypes = [c_double, c_double, c_double]
cukf.set_gyro_bias.restype = None

cukf.get_state.argtypes = [POINTER(State)]
cukf.get_state.restype = None

cukf.get_state_covariance.argtypes = [POINTER(c_double * (24**2))]
cukf.get_state_covariance.restype = None

cukf.sensor_clear.argtypes = []
cukf.sensor_clear.restype = None

cukf.sensor_set_covariance.argtypes = [POINTER(c_double * 17)]
cukf.sensor_set_covariance.restype = None

cukf.sensor_set_accelerometer.argtypes = [c_double, c_double, c_double]
cukf.sensor_set_accelerometer.restype = None

cukf.sensor_set_gyroscope.argtypes = [c_double, c_double, c_double]
cukf.sensor_set_gyroscope.restype = None

cukf.sensor_set_magnetometer.argtypes = [c_double, c_double, c_double]
cukf.sensor_set_magnetometer.restype = None

cukf.sensor_set_gps_position.argtypes = [c_double, c_double, c_double]
cukf.sensor_set_gps_position.restype = None

cukf.sensor_set_gps_velocity.argtypes = [c_double, c_double, c_double]
cukf.sensor_set_gps_velocity.restype = None

cukf.sensor_set_pitot_tas.argtypes = [c_double]
cukf.sensor_set_pitot_tas.restype = None

cukf.sensor_set_barometer_amsl.argtypes = [c_double]
cukf.sensor_set_barometer_amsl.restype = None

cukf.ukf_set_params.argtypes = [POINTER(IOBoardParams)]
cukf.ukf_set_params.restype = None

cukf.ukf_set_field.argtypes = [c_double, c_double, c_double]
cukf.ukf_set_field.restype = None

cukf.ukf_choose_dynamics.argtypes = [c_int]
cukf.ukf_choose_dynamics.restype = None

cukf.ukf_iterate.argtypes = [c_float, POINTER(c_double * 4)]
cukf.ukf_iterate.restype = None

cukf.ukf_set_process_noise.argtypes = [POINTER(c_double * 24)]
cukf.ukf_set_process_noise.restype = None

cukf.fixedwingdynamics_set_mass.argtypes = [c_double]
cukf.fixedwingdynamics_set_mass.restype = None

cukf.fixedwingdynamics_set_inertia_tensor.argtypes = [POINTER(c_double * 9)]
cukf.fixedwingdynamics_set_inertia_tensor.restype = None

cukf.fixedwingdynamics_set_prop_coeffs.argtypes = [c_double, c_double]
cukf.fixedwingdynamics_set_prop_coeffs.restype = None

cukf.fixedwingdynamics_set_drag_coeffs.argtypes = [POINTER(c_double * 5)]
cukf.fixedwingdynamics_set_drag_coeffs.restype = None

cukf.fixedwingdynamics_set_lift_coeffs.argtypes = [POINTER(c_double * 5)]
cukf.fixedwingdynamics_set_lift_coeffs.restype = None

cukf.fixedwingdynamics_set_side_coeffs.argtypes = [POINTER(c_double * 12)]
cukf.fixedwingdynamics_set_side_coeffs.restype = None

cukf.fixedwingdynamics_set_pitch_moment_coeffs.argtypes = [
    POINTER(c_double * 6)]
cukf.fixedwingdynamics_set_pitch_moment_coeffs.restype = None

cukf.fixedwingdynamics_set_roll_moment_coeffs.argtypes = [
    POINTER(c_double * 5)]
cukf.fixedwingdynamics_set_roll_moment_coeffs.restype = None

cukf.fixedwingdynamics_set_yaw_moment_coeffs.argtypes = [
    POINTER(c_double * 6)]
cukf.fixedwingdynamics_set_yaw_moment_coeffs.restype = None

current = State()
cukf.get_state(byref(current))
