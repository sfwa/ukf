#Copyright (C) 2013 Daniel Dyer
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

import os
from ctypes import *


# Taken from c/cukf.h
UKF_PRECISION_FLOAT = 0
UKF_PRECISION_DOUBLE = 1


state = None
covariance = None
sensor_errors = None


# Internal globals, set during init
_cukf = None
_REAL_T = None


class _SigmaPoint(object):
    def __init__(self, arr):
        self.attitude = (arr[0], arr[1], arr[2], arr[3])
        self.angular_velocity = (arr[4], arr[5], arr[6])
        self.acceleration = (arr[7], arr[8], arr[9])

    def __repr__(self):
        return str(self.__dict__)

# Internal classes, wrapping cukf structs directly
class _SensorParams(Structure):
    pass


class _State(Structure):
    def __repr__(self):
        fields = {
            "attitude": tuple(self.attitude),
            "angular_velocity": tuple(self.angular_velocity),
            "acceleration": tuple(self.acceleration)
        }
        return str(fields)

class _SensorErrors(Structure):
    def __repr__(self):
        field = {
            "accel_bias": tuple(self.accel_bias),
            "accel_scale": tuple(self.accel_scale),
            "gyro_bias": tuple(self.gyro_bias),
            "gyro_scale": tuple(self.gyro_scale),
            "mag_bias": tuple(self.mag_bias),
            "mag_scale": tuple(self.mag_scale),
        }
        return std(fields)

# Public interface
def iterate(dt):
    global _cukf, state, sensor_errors

    if not _cukf:
        raise RuntimeError("Please call ukf.init()")

    _cukf.ukf_set_state(state)
    _cukf.ukf_iterate(dt)
    _cukf.ukf_sensor_clear()
    _cukf.ukf_get_state(state)
    #_cukf.ukf_get_state_covariance(covariance)
    _cukf.ukf_get_sensor_errors(sensor_errors)


def set_sensors(accelerometer=None, gyroscope=None, magnetometer=None):
    if accelerometer is not None:
        _cukf.ukf_sensor_set_accelerometer(*accelerometer)
    if gyroscope is not None:
        _cukf.ukf_sensor_set_gyroscope(*gyroscope)
    if magnetometer is not None:
        _cukf.ukf_sensor_set_magnetometer(*magnetometer)


def configure_sensors(accelerometer_covariance=None,
        gyroscope_covariance=None, magnetometer_covariance=None):
    params = _SensorParams()

    if getattr(accelerometer_covariance, '__iter__', False):
        params.accel_covariance = accelerometer_covariance
    elif accelerometer_covariance is not None:
        params.accel_covariance = (accelerometer_covariance, ) * 3
    else:
        params.accel_covariance = (1.0, 1.0, 1.0)

    if getattr(gyroscope_covariance, '__iter__', False):
        params.gyro_covariance = gyroscope_covariance
    elif gyroscope_covariance is not None:
        params.gyro_covariance = (gyroscope_covariance, ) * 3
    else:
        params.gyro_covariance = (1.0, 1.0, 1.0)

    if getattr(magnetometer_covariance, '__iter__', False):
        params.mag_covariance = magnetometer_covariance
    elif magnetometer_covariance is not None:
        params.mag_covariance = (magnetometer_covariance, ) * 3
    else:
        params.mag_covariance = (1.0, 1.0, 1.0)

    _cukf.ukf_set_params(params)


def configure_process_noise(process_noise_covariance):
    _cukf.ukf_set_process_noise((_REAL_T * 9)(*process_noise_covariance))


def init():
    global _cukf, _REAL_T, state, sensor_errors

    lib = os.path.join(os.path.dirname(__file__), "libahrs.dylib")
    _cukf = cdll.LoadLibrary(lib)

    _cukf.ukf_init.argtypes = []
    _cukf.ukf_init.restype = None

    _cukf.ukf_config_get_precision.argtypes = []
    _cukf.ukf_config_get_precision.restype = c_long

    _cukf.ukf_config_get_state_dim.argtypes = []
    _cukf.ukf_config_get_state_dim.restype = c_long

    _cukf.ukf_config_get_measurement_dim.argtypes = []
    _cukf.ukf_config_get_measurement_dim.restype = c_long

    _PRECISION = _cukf.ukf_config_get_precision()
    _REAL_T = c_double if _PRECISION == UKF_PRECISION_DOUBLE else c_float
    _STATE_DIM = _cukf.ukf_config_get_state_dim()
    _MEASUREMENT_DIM = _cukf.ukf_config_get_measurement_dim()

    _SensorParams._fields_ = [
        ("accel_covariance", _REAL_T * 3),
        ("gyro_covariance", _REAL_T * 3),
        ("mag_covariance", _REAL_T * 3)
    ]

    _State._fields_ = [
        ("attitude", _REAL_T * 4),
        ("angular_velocity", _REAL_T * 3),
        ("acceleration", _REAL_T * 3)
    ]

    _SensorErrors._fields_ = [
        ("accel_bias", _REAL_T * 3),
        ("accel_scale", _REAL_T * 3),
        ("gyro_bias", _REAL_T * 3),
        ("gyro_scale", _REAL_T * 3),
        ("mag_bias", _REAL_T * 3),
        ("mag_scale", _REAL_T * 9)
    ]

    # Set up the function prototypes
    _cukf.ukf_set_acceleration.argtypes = [_REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_set_acceleration.restype = None

    _cukf.ukf_set_attitude.argtypes = [_REAL_T, _REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_set_attitude.restype = None

    _cukf.ukf_set_angular_velocity.argtypes = [_REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_set_angular_velocity.restype = None

    _cukf.ukf_get_state.argtypes = [POINTER(_State)]
    _cukf.ukf_get_state.restype = None

    _cukf.ukf_set_state.argtypes = [POINTER(_State)]
    _cukf.ukf_set_state.restype = None

    _cukf.ukf_get_state_covariance.argtypes = [
        POINTER(_REAL_T * (_STATE_DIM**2))]
    _cukf.ukf_get_state_covariance.restype = None

    _cukf.ukf_sensor_clear.argtypes = []
    _cukf.ukf_sensor_clear.restype = None

    _cukf.ukf_sensor_set_accelerometer.argtypes = [_REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_sensor_set_accelerometer.restype = None

    _cukf.ukf_sensor_set_gyroscope.argtypes = [_REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_sensor_set_gyroscope.restype = None

    _cukf.ukf_sensor_set_magnetometer.argtypes = [_REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_sensor_set_magnetometer.restype = None

    _cukf.ukf_set_params.argtypes = [POINTER(_SensorParams)]
    _cukf.ukf_set_params.restype = None

    _cukf.ukf_iterate.argtypes = [c_float]
    _cukf.ukf_iterate.restype = None

    _cukf.ukf_set_process_noise.argtypes = [POINTER(_REAL_T * _STATE_DIM)]
    _cukf.ukf_set_process_noise.restype = None

    _cukf.ukf_get_sensor_errors.argtypes = [POINTER(_SensorErrors)]
    _cukf.ukf_get_sensor_errors.restype = None

    # Initialize the library
    _cukf.ukf_init()

    # Set up the state
    state = _State()
    _cukf.ukf_get_state(state)

    # Set up the sensor errors
    sensor_errors = _SensorErrors()
    _cukf.ukf_get_sensor_errors(sensor_errors)