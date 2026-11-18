#include <cmath>
#include <cstdio>

#include "cukf.h"

namespace {
int calls = 0;
bool correct_control = true;
const real_t input[UKF_CONTROL_DIM] = {0.1, 0.2, 0.3, 0.4};

void custom(const real_t*, const real_t* control, real_t* output) {
  ++calls;
  for (int i = 0; i < UKF_CONTROL_DIM; ++i) {
    correct_control = correct_control && control[i] == input[i];
  }
  for (int i = 0; i < 6; ++i) {
    output[i] = control[i % UKF_CONTROL_DIM];
  }
}

ukf_state_t model_step(ukf_model_t model, real_t throttle) {
  ukf_init();
  ukf_set_velocity(4, 0, 0);
  ukf_set_angular_velocity(0, 0, 1);
  ukf_choose_dynamics(model);
  real_t control[UKF_CONTROL_DIM] = {throttle, 0.5, 0.5, 0};
  ukf_iterate(0.01f, control);
  ukf_state_t state;
  ukf_get_state(&state);
  return state;
}
}  // namespace

int main() {
  ukf_init();
  ukf_set_custom_dynamics_model(custom);
  real_t control[UKF_CONTROL_DIM] = {input[0], input[1], input[2], input[3]};
  ukf_iterate(0.01f, control);
  ukf_state_t state;
  ukf_get_state(&state);
  if (calls != 2 * UKF_STATE_DIM + 1 || !correct_control ||
      std::abs(state.acceleration[2] - input[2]) > 1e-9 ||
      std::abs(state.angular_acceleration[0] - input[3]) > 1e-9) {
    std::fputs("Custom dynamics did not propagate the supplied controls\n",
               stderr);
    return 1;
  }
  ukf_choose_dynamics(UKF_MODEL_NONE);
  const int previous_calls = calls;
  ukf_iterate(0.01f, nullptr);
  if (calls != previous_calls) {
    std::fputs("Switching dynamics did not disable the callback\n", stderr);
    return 1;
  }
  const auto none = model_step(UKF_MODEL_NONE, 0);
  const auto centripetal = model_step(UKF_MODEL_CENTRIPETAL, 0);
  const auto idle = model_step(UKF_MODEL_X8, 0);
  const auto powered = model_step(UKF_MODEL_X8, 1);
  if (!(std::abs(none.acceleration[1]) < 1e-9 &&
        std::abs(centripetal.acceleration[1] - none.acceleration[1]) > 0.1 &&
        powered.acceleration[0] > idle.acceleration[0] + 1)) {
    std::fputs("Built-in dynamics did not affect the predicted acceleration\n",
               stderr);
    return 1;
  }
  return 0;
}
