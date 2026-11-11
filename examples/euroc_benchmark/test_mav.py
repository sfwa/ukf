"""Independent sensor-only checks of the MAV attitude estimator.

The flight generator integrates Newton's law in world coordinates. It supplies
only gyro and specific force to the executable: neither reference orientation
nor simulated velocity is an input to the estimator.
"""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation


ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = Path(os.environ.get(
    "EUROC_BENCHMARK_EXECUTABLE",
    ROOT / "build/examples/euroc_benchmark/euroc_benchmark",
)).resolve()
HEADER = "time,gx,gy,gz,ax,ay,az,mx,my,mz"
G = 9.80665


def aircraft_to_imu(thrust_axis):
    """A fixed sensor mounting, with aircraft +X along collective thrust."""
    x = np.asarray(thrust_axis, dtype=float).copy()
    x /= np.linalg.norm(x)
    y = np.array([0., 1., 0.])
    if abs(y @ x) > .9:
        y = np.array([0., 0., 1.])
    y -= (y @ x)*x
    y /= np.linalg.norm(y)
    return Rotation.from_matrix(np.column_stack((x, y, np.cross(x, y))))


def rotor_flight(time, thrust_axis, drag_rate, tilt_degrees=8.0, bias=None,
                 vibration_amplitude=1.2):
    """Smooth pitch excursion with a collective-thrust altitude controller.

    Aerodynamic drag opposes velocity in the rotor plane. Integrating world
    velocity avoids reusing the estimator's rotating-frame process equations.
    The collective controller preserves zero vertical velocity; it does not
    constrain either horizontal velocity or the orientation estimate. Powered
    flight adds zero-mean 70 Hz accelerometer mount vibration along thrust.
    This measurement perturbation does not change the mean-force trajectory.
    """
    time = np.asarray(time)
    amplitude = np.deg2rad(tilt_degrees)
    bias = np.zeros(3) if bias is None else np.asarray(bias)
    sensor_mount = aircraft_to_imu(thrust_axis)
    initial = Rotation.from_euler('z', 37, degrees=True) * Rotation.from_euler('y', -90, degrees=True)

    def smooth(t, start, stop):
        x = np.clip((t-start)/(stop-start), 0, 1)
        value = x*x*(3-2*x)
        rate = 6*x*(1-x)/(stop-start) if start < t < stop else 0
        return value, rate

    def aircraft_attitude(t):
        up, up_rate = smooth(t, 6, 8)
        down, down_rate = smooth(t, 15, 17)
        return initial * Rotation.from_euler('y', amplitude*(up-down)), amplitude*(up_rate-down_rate)

    def specific_force(t, velocity):
        orientation, _ = aircraft_attitude(t)
        body_velocity = orientation.inv().apply(velocity)
        drag = -drag_rate*body_velocity
        drag[0] = 0
        drag_world = orientation.apply(drag)
        thrust_world = orientation.apply([1., 0., 0.])
        collective = (G-drag_world[2])/thrust_world[2]
        return drag + [collective, 0., 0.]

    def derivative(t, velocity):
        orientation, _ = aircraft_attitude(t)
        return orientation.apply(specific_force(t, velocity)) - [0., 0., G]

    velocity = solve_ivp(derivative, (float(time[0]), float(time[-1])),
                         np.zeros(3), t_eval=time, rtol=1e-10, atol=1e-11,
                         max_step=.025).y.T
    rotations, gyro, acceleration = [], [], []
    for t, v in zip(time, velocity):
        orientation, rate = aircraft_attitude(float(t))
        rotations.append((orientation * sensor_mount.inv()).as_quat())
        gyro.append(sensor_mount.apply([0., rate, 0.])+bias)
        vibration = vibration_amplitude*np.sin(2*np.pi*70*t) if t >= 6 else 0
        acceleration.append(sensor_mount.apply(specific_force(float(t), v)+[vibration, 0., 0.]))
    samples = np.column_stack((time, gyro, acceleration, np.zeros((len(time), 3))))
    return samples, Rotation.from_quat(rotations), velocity


class MAVEstimatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        result = subprocess.run([str(EXECUTABLE), "--describe-filter", "mav-srukf"],
                                capture_output=True, text=True, check=True, timeout=30)
        cls.parameters = json.loads(result.stdout)
        cls.thrust_axis = np.asarray(cls.parameters["thrust_axis_imu"], dtype=float)
        cls.drag_rate = float(cls.parameters["drag_rate_s_inv"])
        cls.flight_time = np.arange(6001)/200
        cls.flight_samples, cls.flight_truth, cls.flight_velocity = rotor_flight(
            cls.flight_time, cls.thrust_axis, cls.drag_rate)

    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.run_number = 0

    def replay(self, samples, algorithm="mav-srukf"):
        self.run_number += 1
        stem = self.directory / str(self.run_number)
        source, output, diagnostic = (stem.with_suffix(ext) for ext in (".in.csv", ".out.csv", ".diag.csv"))
        np.savetxt(source, samples, delimiter=",", header=HEADER, comments="", fmt="%.17g")
        result = subprocess.run([
            str(EXECUTABLE), "--input", str(source), "--output", str(output),
            "--filter", algorithm, "--axes", "6", "--diagnostics", str(diagnostic),
        ], capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr)
        values = np.loadtxt(output, delimiter=",", skiprows=1, ndmin=2)
        self.assertEqual(values.shape, (len(samples), 5))
        np.testing.assert_allclose(values[:, 0], samples[:, 0], atol=1e-10, rtol=0)
        self.assertTrue(np.isfinite(values).all())
        np.testing.assert_allclose(np.linalg.norm(values[:, 1:], axis=1), 1, atol=1e-10, rtol=0)
        details = np.genfromtxt(diagnostic, delimiter=",", names=True, ndmin=1)
        self.assertEqual(len(details), len(samples))
        return Rotation.from_quat(values[:, 1:], scalar_first=True), details

    @staticmethod
    def bias(details):
        return np.column_stack([details[name] for name in ("bias_x", "bias_y", "bias_z")])

    @staticmethod
    def inclination(estimate, truth):
        vertical = np.tile([0., 0., 1.], (len(estimate), 1))
        cosine = np.sum(estimate.inv().apply(vertical)*truth.inv().apply(vertical), axis=1)
        return np.rad2deg(np.arccos(np.clip(cosine, -1, 1)))

    def test_stationary_arbitrary_mount_and_large_gyro_bias(self):
        time = np.arange(4001)/100
        initial = Rotation.from_euler('xyz', [20, -50, 37], degrees=True)
        truth = Rotation.from_quat(np.tile(initial.as_quat(), (len(time), 1)))
        bias = np.array([.01, -.015, -.08])
        samples = np.column_stack((time, np.tile(bias, (len(time), 1)),
                                   truth.inv().apply(np.tile([0., 0., G], (len(time), 1))),
                                   np.zeros((len(time), 3))))
        estimate, details = self.replay(samples)
        self.assertLess(self.inclination(estimate, truth)[time > 30].max(), .5)
        np.testing.assert_allclose(self.bias(details)[-1], bias, atol=.002, rtol=0)

    def test_rotor_force_model_rejects_false_gravity_during_translation(self):
        samples, truth, velocity = self.flight_samples, self.flight_truth, self.flight_velocity
        estimate, details = self.replay(samples)
        accelerating = (self.flight_time > 8) & (self.flight_time < 12)
        model_rmse = np.sqrt(np.mean(self.inclination(estimate, truth)[accelerating]**2))
        # Check a physical prediction, not a ranking on recorded data. During
        # acceleration, the specific-force direction differs from world up, so
        # an estimator which took the accelerometer for gravity would be wrong
        # by the angle between them.
        force_world = truth.apply(samples[:, 4:7])
        apparent = np.rad2deg(np.arctan2(np.linalg.norm(force_world[:, :2], axis=1), force_world[:, 2]))
        apparent_rmse = np.sqrt(np.mean(apparent[accelerating]**2))
        self.assertGreater(np.linalg.norm(velocity, axis=1).max(), 1)
        self.assertLess(abs(velocity[:, 2]).max(), 1e-7)
        self.assertFalse(details["rest"][accelerating].any(),
                         "Smooth accelerating flight must not be mistaken for rest")
        self.assertGreater(details["rotor_model_active"][accelerating].mean(), .9)
        self.assertGreater(apparent_rmse, 2)
        self.assertLess(model_rmse, .5*apparent_rmse)
        self.assertLess(model_rmse, 1.0)

    def test_handheld_rotation_without_motor_vibration_uses_gravity_fallback(self):
        time = np.arange(2001)/200
        initial = Rotation.from_euler('xyz', [20, -50, 37], degrees=True)
        omega = np.array([.3, -.2, .4])
        truth = initial*Rotation.from_rotvec(time[:, None]*omega)
        acceleration = truth.inv().apply(np.tile([0., 0., G], (len(time), 1)))
        samples = np.column_stack((time, np.tile(omega, (len(time), 1)), acceleration,
                                   np.zeros((len(time), 3))))
        estimate, details = self.replay(samples)
        self.assertFalse(details["rotor_model_active"].any())
        # The sensor is never at rest here, so the gyroscope bias starts with
        # its full five degrees per second of prior uncertainty; the tilt
        # wanders a little while the filter learns it from the motion.
        self.assertLess(self.inclination(estimate, truth).max(), .5)

    def test_prefix_causality_and_unused_magnetometer(self):
        samples = self.flight_samples.copy()
        full, _ = self.replay(samples)
        prefix, _ = self.replay(samples[:2001])
        np.testing.assert_allclose(full[:2001].as_quat(), prefix.as_quat(), atol=1e-12, rtol=0)
        samples[:, 7:10] = [123., -45., 67.]
        changed, _ = self.replay(samples)
        np.testing.assert_allclose(full.as_quat(), changed.as_quat(), atol=1e-12, rtol=0)

    def test_rotor_motion_observes_horizontal_bias_after_takeoff(self):
        samples = self.flight_samples.copy()
        bias_axis = aircraft_to_imu(self.thrust_axis).apply([0., 1., 0.])
        injected_bias = .01*bias_axis
        in_flight = self.flight_time >= 8
        # The bias appears after startup calibration, with no subsequent rest.
        # It is perpendicular to thrust: yaw bias at level hover is not
        # observable from this six-axis experiment and is not claimed here.
        samples[in_flight, 1:4] += injected_bias
        estimate, details = self.replay(samples)
        gravity_control, _ = self.replay(samples, "mav-srukf-gravity")
        self.assertFalse(details["rest"][in_flight].any())
        learned = self.bias(details)[-1]
        self.assertLess(np.linalg.norm(learned-injected_bias), .8*np.linalg.norm(injected_bias))
        late = self.flight_time > 20
        model_rmse = np.sqrt(np.mean(self.inclination(estimate, self.flight_truth)[late]**2))
        gravity_rmse = np.sqrt(np.mean(self.inclination(gravity_control, self.flight_truth)[late]**2))
        self.assertLess(model_rmse, gravity_rmse)

    def test_yaw_spin_uses_original_imu_frame(self):
        time = np.arange(2001)/100
        mount = aircraft_to_imu(self.thrust_axis)
        initial = Rotation.from_euler('y', -90, degrees=True) * mount.inv()
        omega = .4*self.thrust_axis
        truth = initial * Rotation.from_rotvec(time[:, None]*omega)
        acceleration = np.tile(G*self.thrust_axis, (len(time), 1))
        samples = np.column_stack((time, np.tile(omega, (len(time), 1)), acceleration,
                                   np.zeros((len(time), 3))))
        estimate, _ = self.replay(samples)
        aligned = (truth[0]*estimate[0].inv())*estimate
        error = np.rad2deg((truth.inv()*aligned).magnitude())
        self.assertLess(error.max(), .2)

    def test_zero_collective_freefall_coasts_on_gyro(self):
        time = np.arange(1001)/100
        falling = (time >= 5) & (time < 7)
        gyro = np.zeros((len(time), 3))
        omega = np.array([.3, -.4, .2])
        gyro[falling] = omega
        acceleration = np.tile(G*self.thrust_axis, (len(time), 1))
        acceleration[falling] = 0
        samples = np.column_stack((time, gyro, acceleration, np.zeros((len(time), 3))))
        estimate, details = self.replay(samples)
        start = np.flatnonzero(falling)[0]-1
        # The first interval uses the midpoint of the previous/current gyro.
        elapsed = time[falling]-time[start]-(time[start+1]-time[start])/2
        expected = estimate[start]*Rotation.from_rotvec(elapsed[:, None]*omega)
        error = np.rad2deg((expected.inv()*estimate[falling]).magnitude())
        self.assertFalse(details["rotor_model_active"][falling].any())
        # The rest detector needs a few samples to notice that the fall began.
        self.assertFalse(details["rest"][falling & (time >= 5.1)].any())
        # In free fall the accelerometer says nothing about which way is down,
        # and the velocity really does run away from zero. The attitude must
        # nevertheless follow the gyroscope, and the bias must hardly move.
        bias = self.bias(details)
        np.testing.assert_allclose(bias[falling], np.tile(bias[start], (falling.sum(), 1)),
                                   atol=2e-3, rtol=0)
        self.assertLess(error.max(), 1)

    def test_timestamp_origin_does_not_change_state(self):
        samples = self.flight_samples.copy()
        original, _ = self.replay(samples)
        samples[:, 0] += 10000
        shifted, _ = self.replay(samples)
        np.testing.assert_allclose(original.as_quat(), shifted.as_quat(), atol=1e-8, rtol=0)


if __name__ == "__main__":
    unittest.main()
