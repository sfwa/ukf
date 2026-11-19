#include <cmath>
#include <cstdio>

#include "ahrs.h"

namespace {
bool near(real_t a, real_t b) { return std::abs(a - b) < real_t(1e-4); }

bool same_vector(const real_t (&a)[3], const real_t (&b)[3]) {
  return near(a[0], b[0]) && near(a[1], b[1]) && near(a[2], b[2]);
}

ukf_innovation_t sample(bool gyro_first) {
  ukf_init();
  if (gyro_first) {
    ukf_sensor_set_gyroscope(0.1f, 0.2f, 0.3f);
  }
  ukf_sensor_set_accelerometer(1, 0, -9.8f);
  if (!gyro_first) {
    ukf_sensor_set_gyroscope(0.1f, 0.2f, 0.3f);
  }
  ukf_iterate(0.01f);
  ukf_innovation_t innovation;
  ukf_get_innovation(&innovation);
  return innovation;
}
}  // namespace

int main() {
  const real_t zero[3] = {};
  const auto first = sample(false);
  const auto reversed = sample(true);
  if (!same_vector(first.gyro, reversed.gyro) ||
      !same_vector(first.accel, reversed.accel) ||
      same_vector(first.accel, first.gyro) ||
      !same_vector(reversed.mag, zero)) {
    std::fputs("Incorrect AHRS innovation mapping\n", stderr);
    return 1;
  }

  ukf_state_t previous, current;
  ukf_get_state(&previous);
  ukf_sensor_clear();
  ukf_sensor_set_gyroscope(1, 2, 3);
  ukf_iterate(0.01f);
  ukf_get_state(&current);
  ukf_innovation_t innovation;
  ukf_get_innovation(&innovation);
  if (!same_vector(previous.acceleration, current.acceleration) ||
      !same_vector(innovation.accel, zero) ||
      !same_vector(innovation.mag, zero) ||
      same_vector(innovation.gyro, zero)) {
    std::fputs("Incorrect gyro-only AHRS update\n", stderr);
    return 1;
  }

  ukf_init();
  real_t noise[UKF_STATE_DIM] = {0, 0, 0, 4, 0, 0};
  ukf_set_process_noise(noise);
  noise[0] = -1;
  ukf_set_process_noise(noise);
  ukf_iterate(0.01f);
  real_t covariance[UKF_STATE_DIM];
  ukf_get_state_covariance_diagonal(covariance);
  ukf_get_innovation(&innovation);
  if (!near(covariance[3], 4.000001f) || !same_vector(innovation.accel, zero) ||
      !same_vector(innovation.gyro, zero) ||
      !same_vector(innovation.mag, zero)) {
    std::fputs("Incorrect prediction with semidefinite process noise\n",
               stderr);
    return 1;
  }
  return 0;
}
