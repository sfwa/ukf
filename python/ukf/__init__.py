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
UKF_MODEL_NONE = 0
UKF_MODEL_CENTRIPETAL = 1
UKF_MODEL_CUSTOM = 2
UKF_MODEL_X8 = 3


state = None
covariance = None
model = UKF_MODEL_NONE
custom_model = None


# Internal globals, set during init
_cukf = None
_REAL_T = None
_CONTROL_DIM = None
_CUSTOM_MODEL_FUNC = None


class _SigmaPoint(object):
    def __init__(self, arr):
        self.position = (arr[0], arr[1], arr[2])
        self.velocity = (arr[3], arr[4], arr[5])
        self.acceleration = (arr[6], arr[7], arr[8])
        self.attitude = (arr[9], arr[10], arr[11], arr[12])
        self.angular_velocity = (arr[13], arr[14], arr[15])
        self.angular_acceleration = (arr[16], arr[17], arr[18])
        self.wind_velocity = (arr[19], arr[20], arr[21])
        self.gyro_bias = (arr[22], arr[23], arr[24])

    def __repr__(self):
        return str(self.__dict__)


# Wrapper around the custom model function
_custom_model_func_wrapper_cb = None # avoid garbage-collection
def _custom_model_wrapper(state, control, output):
    if not custom_model:
        return

    result = custom_model(_SigmaPoint(state.contents),
                          tuple(control.contents[0:3]))

    output.contents = (_REAL_T * 6)(*result[0:6])


# Internal classes, wrapping cukf structs directly
class _IOBoardParams(Structure):
    pass


class _State(Structure):
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


# Public interface
def iterate(dt, control=None):
    global _cukf, state, model

    if not _cukf:
        raise RuntimeError("Please call ukf.init()")

    if control is None:
        control = (0.0, ) * _CONTROL_DIM
    elif len(control) != _CONTROL_DIM:
        raise ValueError("Control vector must contain %d elements" %
                         _CONTROL_DIM)

    _cukf.ukf_choose_dynamics(model)
    if model == UKF_MODEL_CUSTOM:
        if not custom_model:
            raise RuntimeError(
                "Can't use ukf.model == UKF_MODEL_CUSTOM without a value " +
                "for ukf.custom_model"
            )

    _cukf.ukf_set_state(state)
    _cukf.ukf_iterate(dt, (_REAL_T * _CONTROL_DIM)(*control))
    _cukf.ukf_sensor_clear()
    _cukf.ukf_get_state(state)
    #_cukf.ukf_get_state_covariance(covariance)


def set_sensors(accelerometer=None, gyroscope=None, magnetometer=None,
        gps_position=None, gps_velocity=None, pitot_tas=None,
        barometer_amsl=None):
    if accelerometer is not None:
        _cukf.ukf_sensor_set_accelerometer(*accelerometer)
    if gyroscope is not None:
        _cukf.ukf_sensor_set_gyroscope(*gyroscope)
    if magnetometer is not None:
        _cukf.ukf_sensor_set_magnetometer(*magnetometer)
    if gps_position is not None:
        _cukf.ukf_sensor_set_gps_position(*gps_position)
    if gps_velocity is not None:
        _cukf.ukf_sensor_set_gps_velocity(*gps_velocity)
    if pitot_tas is not None:
        _cukf.ukf_sensor_set_pitot_tas(pitot_tas)
    if barometer_amsl is not None:
        _cukf.ukf_sensor_set_barometer_amsl(barometer_amsl)


def configure_sensors(accelerometer_offset=None,
        accelerometer_orientation=None, gyroscope_orientation=None,
        magnetometer_orientation=None, wmm_field=None,
        accelerometer_covariance=None, gyroscope_covariance=None,
        magnetometer_covariance=None, gps_position_covariance=None,
        gps_velocity_covariance=None, pitot_tas_covariance=None,
        barometer_amsl_covariance=None):
    params = _IOBoardParams()

    if accelerometer_offset is not None:
        params.accel_offset = accelerometer_offset
    else:
        params.accel_offset = (0.0, 0.0, 0.0)

    if accelerometer_orientation is not None:
        params.accel_orientation = accelerometer_orientation
    else:
        params.accel_orientation = (0.0, 0.0, 0.0, 1.0)

    if gyroscope_orientation is not None:
        params.gyro_orientation = gyroscope_orientation
    else:
        params.gyro_orientation = (0.0, 0.0, 0.0, 1.0)

    if magnetometer_orientation is not None:
        params.mag_orientation = magnetometer_orientation
    else:
        params.mag_orientation = (0.0, 0.0, 0.0, 1.0)

    if wmm_field is not None:
        params.mag_field = wmm_field
    else:
        params.mag_field = (1.0, 0.0, 0.0)

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

    if getattr(gps_position_covariance, '__iter__', False):
        params.gps_position_covariance = gps_position_covariance
    elif gps_position_covariance is not None:
        params.gps_position_covariance = (gps_position_covariance, ) * 3
    else:
        params.gps_position_covariance = (1.0, 1.0, 1.0)

    if getattr(gps_velocity_covariance, '__iter__', False):
        params.gps_velocity_covariance = gps_velocity_covariance
    elif gps_velocity_covariance is not None:
        params.gps_velocity_covariance = (gps_velocity_covariance, ) * 3
    else:
        params.gps_velocity_covariance = (1.0, 1.0, 1.0)

    if pitot_tas_covariance is not None:
        params.pitot_covariance = pitot_tas_covariance
    else:
        params.pitot_covariance = 1.0

    if barometer_amsl_covariance is not None:
        params.barometer_amsl_covariance = barometer_amsl_covariance
    else:
        params.barometer_amsl_covariance = 1.0

    _cukf.ukf_set_params(params)


def configure_process_noise(process_noise_covariance):
    _cukf.ukf_set_process_noise((_REAL_T * 24)(*process_noise_covariance))


def init(implementation="c"):
    global _cukf, _REAL_T, _CONTROL_DIM, _CUSTOM_MODEL_FUNC, state, \
           _custom_model_func_wrapper_cb


    # Load the requested library and determine configuration parameters
    if implementation == "c":
        lib = os.path.join(os.path.dirname(__file__), "c", "libcukf.dylib")
    elif implementation == "c66x":
        lib = os.path.join(os.path.dirname(__file__), "ccs-c66x",
                           "libc66ukf.dylib")
    else:
        raise NameError(
            "Unknown UKF implementation: %s (options are 'c', 'c66x')" %
            implementation)

    _cukf = cdll.LoadLibrary(lib)

    _cukf.ukf_init.argtypes = []
    _cukf.ukf_init.restype = None

    _cukf.ukf_config_get_precision.argtypes = []
    _cukf.ukf_config_get_precision.restype = c_long

    _cukf.ukf_config_get_state_dim.argtypes = []
    _cukf.ukf_config_get_state_dim.restype = c_long

    _cukf.ukf_config_get_control_dim.argtypes = []
    _cukf.ukf_config_get_control_dim.restype = c_long

    _cukf.ukf_config_get_measurement_dim.argtypes = []
    _cukf.ukf_config_get_measurement_dim.restype = c_long

    _PRECISION = _cukf.ukf_config_get_precision()
    _REAL_T = c_double if _PRECISION == UKF_PRECISION_DOUBLE else c_float
    _CONTROL_DIM = _cukf.ukf_config_get_control_dim()
    _STATE_DIM = _cukf.ukf_config_get_state_dim()
    _MEASUREMENT_DIM = _cukf.ukf_config_get_measurement_dim()

    _IOBoardParams._fields_ = [
        ("accel_orientation", _REAL_T * 4),
        ("accel_offset", _REAL_T * 3),
        ("gyro_orientation", _REAL_T * 4),
        ("mag_orientation", _REAL_T * 4),
        ("mag_field", _REAL_T * 3),

        ("accel_covariance", _REAL_T * 3),
        ("gyro_covariance", _REAL_T * 3),
        ("mag_covariance", _REAL_T * 3),
        ("gps_position_covariance", _REAL_T * 3),
        ("gps_velocity_covariance", _REAL_T * 3),
        ("pitot_covariance", _REAL_T),
        ("barometer_amsl_covariance", _REAL_T)
    ]

    _State._fields_ = [
        ("position", _REAL_T * 3),
        ("velocity", _REAL_T * 3),
        ("acceleration", _REAL_T * 3),
        ("attitude", _REAL_T * 4),
        ("angular_velocity", _REAL_T * 3),
        ("angular_acceleration", _REAL_T * 3),
        ("wind_velocity", _REAL_T * 3),
        ("gyro_bias", _REAL_T * 3)
    ]

    # Set up the function prototypes
    _cukf.ukf_set_position.argtypes = [_REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_set_position.restype = None

    _cukf.ukf_set_velocity.argtypes = [_REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_set_velocity.restype = None

    _cukf.ukf_set_acceleration.argtypes = [_REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_set_acceleration.restype = None

    _cukf.ukf_set_attitude.argtypes = [_REAL_T, _REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_set_attitude.restype = None

    _cukf.ukf_set_angular_velocity.argtypes = [_REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_set_angular_velocity.restype = None

    _cukf.ukf_set_angular_acceleration.argtypes = [_REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_set_angular_acceleration.restype = None

    _cukf.ukf_set_wind_velocity.argtypes = [_REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_set_wind_velocity.restype = None

    _cukf.ukf_set_gyro_bias.argtypes = [_REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_set_gyro_bias.restype = None

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

    _cukf.ukf_sensor_set_gps_position.argtypes = [_REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_sensor_set_gps_position.restype = None

    _cukf.ukf_sensor_set_gps_velocity.argtypes = [_REAL_T, _REAL_T, _REAL_T]
    _cukf.ukf_sensor_set_gps_velocity.restype = None

    _cukf.ukf_sensor_set_pitot_tas.argtypes = [_REAL_T]
    _cukf.ukf_sensor_set_pitot_tas.restype = None

    _cukf.ukf_sensor_set_barometer_amsl.argtypes = [_REAL_T]
    _cukf.ukf_sensor_set_barometer_amsl.restype = None

    _cukf.ukf_set_params.argtypes = [POINTER(_IOBoardParams)]
    _cukf.ukf_set_params.restype = None

    _cukf.ukf_choose_dynamics.argtypes = [c_int]
    _cukf.ukf_choose_dynamics.restype = None

    _CUSTOM_MODEL_FUNC = CFUNCTYPE(None,
                                   POINTER(_REAL_T * (_STATE_DIM + 1)),
                                   POINTER(_REAL_T * _CONTROL_DIM),
                                   POINTER(_REAL_T * 6))

    _cukf.ukf_set_custom_dynamics_model.argtypes = [_CUSTOM_MODEL_FUNC]
    _cukf.ukf_set_custom_dynamics_model.restype = None

    _cukf.ukf_iterate.argtypes = [c_float, POINTER(_REAL_T * _CONTROL_DIM)]
    _cukf.ukf_iterate.restype = None

    _cukf.ukf_set_process_noise.argtypes = [POINTER(_REAL_T * _STATE_DIM)]
    _cukf.ukf_set_process_noise.restype = None

    # Initialize the library
    _cukf.ukf_init()

    # Set the custom model callback
    _custom_model_func_wrapper_cb = _CUSTOM_MODEL_FUNC(_custom_model_wrapper)
    _cukf.ukf_set_custom_dynamics_model(_custom_model_func_wrapper_cb)

    # Set up the state
    state = _State()
    _cukf.ukf_get_state(state)
