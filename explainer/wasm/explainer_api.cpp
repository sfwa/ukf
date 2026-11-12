/*
The explainer's interface to the multirotor estimator: see explainer_api.h.

Two things happen here.

ex_run() runs the estimator over a recording exactly as the benchmark
program does (examples/euroc_benchmark/estimator.cpp), and keeps what the
page plots: the estimate and its standard deviations at every sample.

ex_explain() takes one filter step apart. MAV_Estimator hides its
intermediate values, and so does UKF::SquareRootCore, so this file rebuilds
them from the state and covariance before the step. It calls the same
public functions of the state and measurement vectors which the core calls,
in the same order. The result is then compared with the real step, so a
change to the estimator which this file does not follow shows up as a large
"check" value (and fails the explainer's tests).
*/

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/StdVector>
#include "mav_estimator.h"
#include "explainer_api.h"

#if defined(__EMSCRIPTEN__)
#include <emscripten/emscripten.h>
#define EX_EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EX_EXPORT
#endif

namespace {

constexpr int STATE_SIZE = static_cast<int>(MAV_StateVector::size());
constexpr int ERROR_SIZE = static_cast<int>(MAV_StateVector::covariance_size());
constexpr int SIGMA_COUNT = static_cast<int>(MAV_StateVector::num_sigma());
constexpr int MAX_ROWS = static_cast<int>(MAV_MeasurementVector::max_size());

const double NOT_A_NUMBER = std::numeric_limits<double>::quiet_NaN();

using StateCovariance = MAV_StateVector::CovarianceMatrix;
using StateSigmaPoints = MAV_StateVector::SigmaPointDistribution;
using StateSigmaDeltas = MAV_StateVector::SigmaPointDeltas;
using MeasurementSigmaPoints = MAV_MeasurementVector::SigmaPointDistribution<MAV_StateVector>;
using MeasurementSigmaDeltas = MAV_MeasurementVector::SigmaPointDeltas<MAV_StateVector>;
using MeasurementCovariance = MAV_MeasurementVector::CovarianceMatrix;
using CrossCovariance = Eigen::Matrix<real_t, ERROR_SIZE, Eigen::Dynamic, Eigen::ColMajor, ERROR_SIZE, MAX_ROWS>;

/* A named block of numbers in a result buffer. */
struct Field {
    const char* name;
    int rows;
    int cols;
};

int field_size(const Field& field) {
    return field.rows * field.cols;
}

template <std::size_t N>
int field_offset(const Field (&fields)[N], int index) {
    int offset = 0;
    for (int i = 0; i < index; i++) {
        offset += field_size(fields[i]);
    }
    return offset;
}

template <std::size_t N>
int field_total(const Field (&fields)[N]) {
    return field_offset(fields, static_cast<int>(N));
}

/* Model constants fixed in mav_estimator.h. */
struct Constant {
    const char* name;
    double value;
};

const Constant CONSTANTS[] = {
    {"gravity", G_ACCEL},
    {"drag_rate", MAV_DRAG_RATE},
    {"thrust_tilt", MAV_THRUST_TILT},
    {"thrust_axis_x", MAV_ThrustAxis[0]},
    {"thrust_axis_y", MAV_ThrustAxis[1]},
    {"thrust_axis_z", MAV_ThrustAxis[2]},
    {"lateral_axis_1_x", MAV_LateralAxis1[0]},
    {"lateral_axis_1_y", MAV_LateralAxis1[1]},
    {"lateral_axis_1_z", MAV_LateralAxis1[2]},
    {"lateral_axis_2_x", MAV_LateralAxis2[0]},
    {"lateral_axis_2_y", MAV_LateralAxis2[1]},
    {"lateral_axis_2_z", MAV_LateralAxis2[2]},
    {"state_size", STATE_SIZE},
    {"error_size", ERROR_SIZE},
    {"sigma_count", SIGMA_COUNT},
    {"sigma_spread", std::sqrt(UKF::Parameters::ScaledDimension<MAV_StateVector>)},
    {"weight_mean_centre", UKF::Parameters::Sigma_WM0<MAV_StateVector>},
    {"weight_mean_outer", UKF::Parameters::Sigma_WMI<MAV_StateVector>},
    {"weight_covariance_centre", UKF::Parameters::Sigma_WC0<MAV_StateVector>},
    {"weight_covariance_outer", UKF::Parameters::Sigma_WCI<MAV_StateVector>},
    {"max_rows", MAX_ROWS},
    {"checkpoint_interval", EX_CHECKPOINT_INTERVAL},
};

/* Every numeric member of MAV_Tuning. */
struct TuningField {
    const char* name;
    real_t& (*value)(MAV_Tuning&);
};

const TuningField TUNING_FIELDS[] = {
    {"initial_tilt", [](MAV_Tuning& t) -> real_t& { return t.initial_tilt; }},
    {"initial_velocity", [](MAV_Tuning& t) -> real_t& { return t.initial_velocity; }},
    {"initial_gyro_bias", [](MAV_Tuning& t) -> real_t& { return t.initial_gyro_bias; }},
    {"gyroscope_noise", [](MAV_Tuning& t) -> real_t& { return t.gyroscope_noise; }},
    {"gyroscope_scale_error", [](MAV_Tuning& t) -> real_t& { return t.gyroscope_scale_error; }},
    {"accelerometer_noise", [](MAV_Tuning& t) -> real_t& { return t.accelerometer_noise; }},
    {"rotor_model_noise", [](MAV_Tuning& t) -> real_t& { return t.rotor_model_noise; }},
    {"gyro_bias_walk", [](MAV_Tuning& t) -> real_t& { return t.gyro_bias_walk; }},
    {"rotor_drag_noise", [](MAV_Tuning& t) -> real_t& { return t.rotor_drag_noise; }},
    {"zero_velocity_noise", [](MAV_Tuning& t) -> real_t& { return t.zero_velocity_noise; }},
    {"rest_gyroscope_noise", [](MAV_Tuning& t) -> real_t& { return t.rest_gyroscope_noise; }},
    {"heading_sigma_limit", [](MAV_Tuning& t) -> real_t& { return t.heading_sigma_limit; }},
    {"rest_cutoff_frequency", [](MAV_Tuning& t) -> real_t& { return t.rest_cutoff_frequency; }},
    {"rest_averaging_time", [](MAV_Tuning& t) -> real_t& { return t.rest.averaging_time; }},
    {"rest_gyroscope_change", [](MAV_Tuning& t) -> real_t& { return t.rest.gyroscope_change; }},
    {"rest_accelerometer_change", [](MAV_Tuning& t) -> real_t& { return t.rest.accelerometer_change; }},
    {"rest_accelerometer_minimum", [](MAV_Tuning& t) -> real_t& { return t.rest.accelerometer_minimum; }},
    {"rest_accelerometer_maximum", [](MAV_Tuning& t) -> real_t& { return t.rest.accelerometer_maximum; }},
    {"rest_gyroscope_bias", [](MAV_Tuning& t) -> real_t& { return t.rest.gyroscope_bias; }},
    {"rest_minimum_time", [](MAV_Tuning& t) -> real_t& { return t.rest.minimum_time; }},
    {"rotors_cutoff_frequency", [](MAV_Tuning& t) -> real_t& { return t.rotors.cutoff_frequency; }},
    {"rotors_averaging_time", [](MAV_Tuning& t) -> real_t& { return t.rotors.averaging_time; }},
    {"rotors_on_vibration", [](MAV_Tuning& t) -> real_t& { return t.rotors.on_vibration; }},
    {"rotors_off_vibration", [](MAV_Tuning& t) -> real_t& { return t.rotors.off_vibration; }},
    {"minimum_thrust", [](MAV_Tuning& t) -> real_t& { return t.minimum_thrust; }},
    {"maximum_thrust", [](MAV_Tuning& t) -> real_t& { return t.maximum_thrust; }},
};

constexpr int TUNING_NUMBERS = static_cast<int>(sizeof(TUNING_FIELDS) / sizeof(TUNING_FIELDS[0]));

/* The one boolean member of MAV_Tuning follows the numbers. */
constexpr int TUNING_USE_ROTOR_MODEL = TUNING_NUMBERS;
constexpr int TUNING_COUNT = TUNING_NUMBERS + 1;

/*
One record of the trace. "sd" holds the standard deviation of each of the
nine error states: the attitude error about the world X, Y and Z axes (rad;
Z is the heading), then the velocity (m/s) and gyro bias (rad/s) errors.
*/
const Field TRACE_FIELDS[] = {
    {"attitude", 4, 1},
    {"velocity", 3, 1},
    {"gyro_bias", 3, 1},
    {"sd", ERROR_SIZE, 1},
    {"rotor_vibration", 1, 1},
    {"resting", 1, 1},
    {"rotor_model_active", 1, 1},
    {"rows", 1, 1},
    {"innovation", MAX_ROWS, 1},
    {"innovation_sd", MAX_ROWS, 1},
};

/*
One filter step taken apart, in the order in which the values are made.
Each sigma point is one column. A state is written as attitude (w, x, y,
z), velocity, gyro bias; a "delta" is the difference between a sigma point
and its centre, in the nine error coordinates.
*/
const Field STEP_FIELDS[] = {
    {"index", 1, 1},
    {"dt", 1, 1},
    {"resting", 1, 1},
    {"rotor_model_active", 1, 1},
    {"rows", 1, 1},
    {"check", 1, 1},
    {"gyroscope", 3, 1},
    {"accelerometer", 3, 1},
    {"process_sd", ERROR_SIZE, 1},
    {"measurement", MAX_ROWS, 1},
    {"measurement_sd", MAX_ROWS, 1},

    /* Before the step. */
    {"prior_state", STATE_SIZE, 1},
    {"prior_covariance", ERROR_SIZE, ERROR_SIZE},
    {"prior_sigma", STATE_SIZE, SIGMA_COUNT},
    {"prior_delta", ERROR_SIZE, SIGMA_COUNT},

    /* Prediction: every sigma point goes through the process model. */
    {"moved_sigma", STATE_SIZE, SIGMA_COUNT},
    {"moved_delta", ERROR_SIZE, SIGMA_COUNT},
    {"predicted_state", STATE_SIZE, 1},
    {"predicted_covariance", ERROR_SIZE, ERROR_SIZE},

    /* New sigma points, which include the process noise. */
    {"predicted_sigma", STATE_SIZE, SIGMA_COUNT},
    {"predicted_delta", ERROR_SIZE, SIGMA_COUNT},

    /* Correction: what each sigma point expects the sensors to read. */
    {"expected", MAX_ROWS, SIGMA_COUNT},
    {"expected_mean", MAX_ROWS, 1},
    {"innovation", MAX_ROWS, 1},
    {"innovation_covariance", MAX_ROWS, MAX_ROWS},
    {"cross_covariance", ERROR_SIZE, MAX_ROWS},
    {"gain", ERROR_SIZE, MAX_ROWS},
    {"correction", ERROR_SIZE, 1},
    {"corrected_state", STATE_SIZE, 1},
    {"corrected_covariance", ERROR_SIZE, ERROR_SIZE},

    /* What the real step produced, after the heading limit. */
    {"final_state", STATE_SIZE, 1},
    {"final_covariance", ERROR_SIZE, ERROR_SIZE},
};

constexpr int TRACE_FIELD_COUNT = static_cast<int>(sizeof(TRACE_FIELDS) / sizeof(TRACE_FIELDS[0]));
constexpr int STEP_FIELD_COUNT = static_cast<int>(sizeof(STEP_FIELDS) / sizeof(STEP_FIELDS[0]));

/* Everything the module remembers between calls. */
struct Session {
    std::vector<double> sample_values;
    std::vector<ImuSample> samples;
    std::vector<double> tuning_values;
    std::vector<double> trace;
    std::vector<double> step;
    std::vector<double> prediction;
    std::vector<MAV_Estimator, Eigen::aligned_allocator<MAV_Estimator>> checkpoints;
    MAV_Tuning tuning;
    int valid = 0;
};

Session& session() {
    static Session instance;
    return instance;
}

void load_tuning(const MAV_Tuning& tuning, std::vector<double>& values) {
    MAV_Tuning copy = tuning;
    values.resize(TUNING_COUNT);
    for (int i = 0; i < TUNING_NUMBERS; i++) {
        values[i] = TUNING_FIELDS[i].value(copy);
    }
    values[TUNING_USE_ROTOR_MODEL] = copy.use_rotor_model ? 1.0 : 0.0;
}

MAV_Tuning store_tuning(const std::vector<double>& values) {
    MAV_Tuning tuning;
    for (int i = 0; i < TUNING_NUMBERS; i++) {
        TUNING_FIELDS[i].value(tuning) = static_cast<real_t>(values[i]);
    }
    tuning.use_rotor_model = values[TUNING_USE_ROTOR_MODEL] != 0.0;
    return tuning;
}

/* Writes into one named field of a result buffer. */
template <std::size_t N>
double* field_pointer(std::vector<double>& buffer, std::size_t record, const Field (&fields)[N],
        const char* name) {
    int offset = 0;
    for (const Field& field : fields) {
        if (std::string(field.name) == name) {
            return buffer.data() + record + offset;
        }
        offset += field_size(field);
    }
    return nullptr;
}

/* A state as attitude (w, x, y, z), velocity, gyro bias. */
void write_state(double* out, const MAV_StateVector& state) {
    UKF::Quaternion attitude = state.get_field<Attitude>();
    UKF::Vector<3> velocity = state.get_field<Velocity>();
    UKF::Vector<3> gyro_bias = state.get_field<GyroBias>();
    out[0] = attitude.w();
    out[1] = attitude.x();
    out[2] = attitude.y();
    out[3] = attitude.z();
    for (int i = 0; i < 3; i++) {
        out[4 + i] = velocity[i];
        out[7 + i] = gyro_bias[i];
    }
}

void write_sigma_points(double* out, const StateSigmaPoints& points) {
    for (int i = 0; i < SIGMA_COUNT; i++) {
        write_state(out + i * STATE_SIZE, MAV_StateVector(points.col(i)));
    }
}

/* Any matrix, column by column. Unused rows of a padded block stay NaN. */
template <typename Derived>
void write_matrix(double* out, const Eigen::MatrixBase<Derived>& matrix, int padded_rows) {
    for (int column = 0; column < matrix.cols(); column++) {
        for (int row = 0; row < matrix.rows(); row++) {
            out[column * padded_rows + row] = matrix(row, column);
        }
    }
}

template <typename Derived>
void write_matrix(double* out, const Eigen::MatrixBase<Derived>& matrix) {
    write_matrix(out, matrix, static_cast<int>(matrix.rows()));
}

/* The measurements which MAV_Estimator::update() uses at this step. */
MAV_MeasurementVector step_measurement(const MAV_Estimator& estimator, const ImuSample& sample) {
    MAV_MeasurementVector measurement;
    if (estimator.rotor_model_active()) {
        measurement.set_field<RotorDrag>(mav_lateral(sample.accelerometer));
    } else {
        measurement.set_field<ZeroVelocity>(UKF::Vector<3>(0, 0, 0));
    }
    if (estimator.resting()) {
        measurement.set_field<RestGyroscope>(sample.gyroscope);
    }
    return measurement;
}

void write_trace_record(Session& s, int index, const MAV_Estimator& estimator) {
    const int stride = field_total(TRACE_FIELDS);
    const std::size_t record = static_cast<std::size_t>(index) * stride;
    const MAV_Filter& filter = estimator.ukf();

    UKF::Quaternion attitude = estimator.attitude();
    double* out = field_pointer(s.trace, record, TRACE_FIELDS, "attitude");
    out[0] = attitude.w();
    out[1] = attitude.x();
    out[2] = attitude.y();
    out[3] = attitude.z();
    write_matrix(field_pointer(s.trace, record, TRACE_FIELDS, "velocity"), estimator.velocity());
    write_matrix(field_pointer(s.trace, record, TRACE_FIELDS, "gyro_bias"), estimator.gyro_bias());

    /* The covariance is root * root^T, so each variance is a row's squared length. */
    out = field_pointer(s.trace, record, TRACE_FIELDS, "sd");
    for (int i = 0; i < ERROR_SIZE; i++) {
        out[i] = filter.root_covariance.row(i).norm();
    }

    *field_pointer(s.trace, record, TRACE_FIELDS, "rotor_vibration") = estimator.rotor_vibration();
    *field_pointer(s.trace, record, TRACE_FIELDS, "resting") = estimator.resting() ? 1.0 : 0.0;
    *field_pointer(s.trace, record, TRACE_FIELDS, "rotor_model_active") =
        estimator.rotor_model_active() ? 1.0 : 0.0;

    /* The first sample only initialises the filter, so it has no innovation. */
    const int rows = index > 0 ? static_cast<int>(filter.innovation.size()) : 0;
    *field_pointer(s.trace, record, TRACE_FIELDS, "rows") = rows;
    double* innovation = field_pointer(s.trace, record, TRACE_FIELDS, "innovation");
    double* innovation_sd = field_pointer(s.trace, record, TRACE_FIELDS, "innovation_sd");
    for (int i = 0; i < rows; i++) {
        innovation[i] = filter.innovation[i];
        innovation_sd[i] = filter.innovation_root_covariance.row(i).norm();
    }
}

}

extern "C" {

EX_EXPORT int ex_abi_version(void) {
    return EX_ABI_VERSION;
}

EX_EXPORT int ex_constant_count(void) {
    return static_cast<int>(sizeof(CONSTANTS) / sizeof(CONSTANTS[0]));
}

EX_EXPORT const char* ex_constant_name(int index) {
    return index >= 0 && index < ex_constant_count() ? CONSTANTS[index].name : nullptr;
}

EX_EXPORT double ex_constant(int index) {
    return index >= 0 && index < ex_constant_count() ? CONSTANTS[index].value : NOT_A_NUMBER;
}

EX_EXPORT double* ex_samples(int count) {
    Session& s = session();
    s.valid = 0;
    s.samples.clear();
    s.trace.clear();
    s.step.clear();
    s.prediction.clear();
    s.checkpoints.clear();
    if (count < 2) {
        s.sample_values.clear();
        return nullptr;
    }
    s.sample_values.assign(static_cast<std::size_t>(count) * EX_SAMPLE_STRIDE, 0.0);
    return s.sample_values.data();
}

EX_EXPORT int ex_sample_count(void) {
    return static_cast<int>(session().sample_values.size() / EX_SAMPLE_STRIDE);
}

EX_EXPORT int ex_tuning_count(void) {
    return TUNING_COUNT;
}

EX_EXPORT const char* ex_tuning_name(int index) {
    if (index >= 0 && index < TUNING_NUMBERS) {
        return TUNING_FIELDS[index].name;
    }
    return index == TUNING_USE_ROTOR_MODEL ? "use_rotor_model" : nullptr;
}

EX_EXPORT double* ex_tuning(void) {
    Session& s = session();
    if (s.tuning_values.empty()) {
        load_tuning(MAV_Tuning(), s.tuning_values);
    }
    return s.tuning_values.data();
}

EX_EXPORT void ex_tuning_reset(void) {
    load_tuning(MAV_Tuning(), session().tuning_values);
}

EX_EXPORT int ex_run(void) {
    Session& s = session();
    const int count = ex_sample_count();
    if (count < 2) {
        return -1;
    }

    /* Copy the samples into the form which the estimator takes. */
    s.samples.resize(count);
    for (int i = 0; i < count; i++) {
        const double* values = s.sample_values.data() + static_cast<std::size_t>(i) * EX_SAMPLE_STRIDE;
        s.samples[i].time = values[0];
        s.samples[i].gyroscope << values[1], values[2], values[3];
        s.samples[i].accelerometer << values[4], values[5], values[6];
        s.samples[i].magnetometer.setZero();
    }

    ex_tuning();
    s.tuning = store_tuning(s.tuning_values);
    s.trace.assign(static_cast<std::size_t>(count) * field_total(TRACE_FIELDS), NOT_A_NUMBER);
    s.checkpoints.clear();
    s.valid = 0;

    /* The same loop as the benchmark program. */
    MAV_Estimator estimator(s.tuning);
    if (!estimator.initialise(s.samples[0])) {
        return 0;
    }
    for (int i = 0; i < count; i++) {
        if (i > 0 && !estimator.update(s.samples[i - 1], s.samples[i])) {
            break;
        }
        write_trace_record(s, i, estimator);
        if (i % EX_CHECKPOINT_INTERVAL == 0) {
            s.checkpoints.push_back(estimator);
        }
        s.valid = i + 1;
    }

    return s.valid;
}

EX_EXPORT const double* ex_trace(void) {
    return session().trace.data();
}

EX_EXPORT int ex_trace_stride(void) {
    return field_total(TRACE_FIELDS);
}

EX_EXPORT int ex_trace_field_count(void) {
    return TRACE_FIELD_COUNT;
}

EX_EXPORT const char* ex_trace_field_name(int index) {
    return index >= 0 && index < TRACE_FIELD_COUNT ? TRACE_FIELDS[index].name : nullptr;
}

EX_EXPORT int ex_trace_field_offset(int index) {
    return index >= 0 && index < TRACE_FIELD_COUNT ? field_offset(TRACE_FIELDS, index) : -1;
}

EX_EXPORT int ex_trace_field_size(int index) {
    return index >= 0 && index < TRACE_FIELD_COUNT ? field_size(TRACE_FIELDS[index]) : -1;
}

EX_EXPORT int ex_explain(int index) {
    Session& s = session();
    if (index < 1 || index >= s.valid) {
        return 1;
    }
    s.step.assign(field_total(STEP_FIELDS), NOT_A_NUMBER);
    auto out = [&s](const char* name) { return field_pointer(s.step, 0, STEP_FIELDS, name); };

    /* Bring a copy of the estimator to the sample before this step. */
    const int checkpoint = (index - 1) / EX_CHECKPOINT_INTERVAL;
    MAV_Estimator estimator = s.checkpoints[checkpoint];
    for (int i = checkpoint * EX_CHECKPOINT_INTERVAL + 1; i < index; i++) {
        if (!estimator.update(s.samples[i - 1], s.samples[i])) {
            return 2;
        }
    }

    /* Take the real step, and keep the estimator from before it. */
    const MAV_Estimator before = estimator;
    const ImuSample& previous = s.samples[index - 1];
    const ImuSample& sample = s.samples[index];
    if (!estimator.update(previous, sample)) {
        return 2;
    }
    const MAV_Filter& prior = before.ukf();
    const MAV_Filter& result = estimator.ukf();

    /* The inputs of the step, as MAV_Estimator::update() forms them. */
    const real_t dt = sample.time - previous.time;
    MAV_Input input;
    input.gyroscope = (previous.gyroscope + sample.gyroscope) * 0.5;
    input.accelerometer = (previous.accelerometer + sample.accelerometer) * 0.5;
    input.powered = estimator.rotor_model_active();
    const MAV_MeasurementVector measurement = step_measurement(estimator, sample);
    const int rows = static_cast<int>(measurement.size());

    *out("index") = index;
    *out("dt") = dt;
    *out("resting") = estimator.resting() ? 1.0 : 0.0;
    *out("rotor_model_active") = estimator.rotor_model_active() ? 1.0 : 0.0;
    *out("rows") = rows;
    write_matrix(out("gyroscope"), input.gyroscope);
    write_matrix(out("accelerometer"), input.accelerometer);
    write_matrix(out("process_sd"), result.process_noise_root_covariance.diagonal());
    write_matrix(out("measurement"), measurement);

    /* 1. Sigma points around the estimate from before the step. */
    const real_t spread = std::sqrt(UKF::Parameters::ScaledDimension<MAV_StateVector>);
    const StateCovariance prior_covariance = prior.root_covariance * prior.root_covariance.transpose();
    const StateSigmaPoints prior_sigma = prior.state.calculate_sigma_point_distribution(
        prior.root_covariance * spread);
    const StateSigmaDeltas prior_delta = prior.state.calculate_sigma_point_deltas(prior_sigma);
    write_state(out("prior_state"), prior.state);
    write_matrix(out("prior_covariance"), prior_covariance);
    write_sigma_points(out("prior_sigma"), prior_sigma);
    write_matrix(out("prior_delta"), prior_delta);

    /* 2. Each one goes through the process model. */
    StateSigmaPoints moved_sigma;
    for (int i = 0; i < SIGMA_COUNT; i++) {
        moved_sigma.col(i) << MAV_StateVector(prior_sigma.col(i)).process_model<UKF::IntegratorRK4>(dt, input);
    }

    /* 3. Their weighted mean and spread, plus process noise, is the prediction. */
    const MAV_StateVector predicted_state = MAV_StateVector::calculate_sigma_point_mean(moved_sigma);
    const StateSigmaDeltas moved_delta = predicted_state.calculate_sigma_point_deltas(moved_sigma);
    StateCovariance predicted_covariance = MAV_StateVector::calculate_sigma_point_covariance(moved_delta);
    predicted_covariance += result.process_noise_root_covariance * result.process_noise_root_covariance.transpose();
    write_sigma_points(out("moved_sigma"), moved_sigma);
    write_matrix(out("moved_delta"), moved_delta);
    write_state(out("predicted_state"), predicted_state);
    write_matrix(out("predicted_covariance"), predicted_covariance);

    /* 4. New sigma points around the prediction. */
    const Eigen::LLT<StateCovariance> predicted_factor(predicted_covariance);
    if (predicted_factor.info() != Eigen::Success) {
        return 2;
    }
    const StateCovariance predicted_root = predicted_factor.matrixL();
    const StateSigmaPoints predicted_sigma = predicted_state.calculate_sigma_point_distribution(
        predicted_root * spread);
    const StateSigmaDeltas predicted_delta = predicted_state.calculate_sigma_point_deltas(predicted_sigma);
    write_sigma_points(out("predicted_sigma"), predicted_sigma);
    write_matrix(out("predicted_delta"), predicted_delta);

    /* 5. What each sigma point expects the sensors to read. */
    const MeasurementSigmaPoints expected =
        measurement.calculate_sigma_point_distribution<MAV_StateVector>(predicted_sigma);
    const MAV_MeasurementVector expected_mean =
        measurement.calculate_sigma_point_mean<MAV_StateVector>(expected);
    const MeasurementSigmaDeltas expected_delta =
        expected_mean.calculate_sigma_point_deltas<MAV_StateVector>(expected);
    const MAV_MeasurementVector innovation = expected_mean.calculate_innovation(measurement);
    const MeasurementCovariance measurement_root =
        measurement.calculate_measurement_root_covariance(result.measurement_root_covariance, expected_mean);
    MeasurementCovariance innovation_covariance =
        expected_mean.calculate_sigma_point_covariance<MAV_StateVector>(expected_delta);
    innovation_covariance += measurement_root * measurement_root.transpose();
    write_matrix(out("measurement_sd"), measurement_root.diagonal());
    write_matrix(out("expected"), expected, MAX_ROWS);
    write_matrix(out("expected_mean"), expected_mean);
    write_matrix(out("innovation"), innovation);
    write_matrix(out("innovation_covariance"), innovation_covariance, MAX_ROWS);

    /* 6. How the state error and the measurement error move together. */
    CrossCovariance cross_covariance = CrossCovariance::Zero(ERROR_SIZE, rows);
    for (int i = 0; i < SIGMA_COUNT; i++) {
        const real_t weight = i == 0 ? UKF::Parameters::Sigma_WC0<MAV_StateVector> :
            UKF::Parameters::Sigma_WCI<MAV_StateVector>;
        cross_covariance += weight * predicted_delta.col(i) * expected_delta.col(i).transpose();
    }

    /* 7. The gain, the correction and the corrected covariance. */
    const Eigen::LLT<MeasurementCovariance> innovation_factor(innovation_covariance);
    if (innovation_factor.info() != Eigen::Success) {
        return 2;
    }
    const CrossCovariance gain = innovation_factor.solve(cross_covariance.transpose()).transpose();
    const MAV_StateVector::StateVectorDelta correction = gain * innovation;
    MAV_StateVector corrected_state = predicted_state;
    corrected_state.apply_delta(correction);
    const StateCovariance corrected_covariance =
        predicted_covariance - gain * innovation_covariance * gain.transpose();
    write_matrix(out("cross_covariance"), cross_covariance);
    write_matrix(out("gain"), gain);
    write_matrix(out("correction"), correction);
    write_state(out("corrected_state"), corrected_state);
    write_matrix(out("corrected_covariance"), corrected_covariance);

    /* 8. The estimator's limit on the heading uncertainty. */
    StateCovariance limited_covariance = corrected_covariance;
    const real_t heading_sigma = std::sqrt(limited_covariance(2, 2));
    if (heading_sigma > s.tuning.heading_sigma_limit) {
        limited_covariance.row(2) *= s.tuning.heading_sigma_limit / heading_sigma;
        limited_covariance.col(2) *= s.tuning.heading_sigma_limit / heading_sigma;
    }

    /* Compare with the real step. */
    const StateCovariance final_covariance = result.root_covariance * result.root_covariance.transpose();
    write_state(out("final_state"), result.state);
    write_matrix(out("final_covariance"), final_covariance);
    const real_t state_difference = (corrected_state - result.state).cwiseAbs().maxCoeff();
    const real_t covariance_difference = (limited_covariance - final_covariance).cwiseAbs().maxCoeff();
    *out("check") = state_difference > covariance_difference ? state_difference : covariance_difference;

    return 0;
}

EX_EXPORT int ex_predict(int index, int steps) {
    Session& s = session();
    if (index < 0 || index >= s.valid || steps < 0) {
        return 0;
    }

    /* Bring a copy of the estimator to this sample. */
    const int checkpoint = index / EX_CHECKPOINT_INTERVAL;
    MAV_Estimator estimator = s.checkpoints[checkpoint];
    for (int i = checkpoint * EX_CHECKPOINT_INTERVAL + 1; i <= index; i++) {
        if (!estimator.update(s.samples[i - 1], s.samples[i])) {
            return 0;
        }
    }

    /* Its sigma points, as the filter itself would make them for the next step. */
    const MAV_Filter& filter = estimator.ukf();
    const real_t spread = std::sqrt(UKF::Parameters::ScaledDimension<MAV_StateVector>);
    StateSigmaPoints sigma = filter.state.calculate_sigma_point_distribution(filter.root_covariance * spread);

    const int last = index + steps < s.valid ? index + steps : s.valid - 1;
    const int snapshots = last - index + 1;
    const std::size_t snapshot_size = static_cast<std::size_t>(SIGMA_COUNT) * STATE_SIZE;
    s.prediction.assign(snapshot_size * snapshots, NOT_A_NUMBER);
    write_sigma_points(s.prediction.data(), sigma);

    const int stride = field_total(TRACE_FIELDS);
    for (int i = index + 1; i <= last; i++) {
        const ImuSample& previous = s.samples[i - 1];
        const ImuSample& sample = s.samples[i];
        MAV_Input input;
        input.gyroscope = (previous.gyroscope + sample.gyroscope) * 0.5;
        input.accelerometer = (previous.accelerometer + sample.accelerometer) * 0.5;
        input.powered = *field_pointer(s.trace, static_cast<std::size_t>(i) * stride, TRACE_FIELDS,
            "rotor_model_active") != 0.0;
        const real_t dt = sample.time - previous.time;
        for (int k = 0; k < SIGMA_COUNT; k++) {
            MAV_StateVector moved = MAV_StateVector(sigma.col(k)).process_model<UKF::IntegratorRK4>(dt, input);
            /* The filter renormalises the attitude at every step; do the same here. */
            moved.set_field<Attitude>(moved.get_field<Attitude>().normalized());
            sigma.col(k) << moved;
        }
        write_sigma_points(s.prediction.data() + snapshot_size * (i - index), sigma);
    }

    return snapshots;
}

EX_EXPORT const double* ex_prediction(void) {
    return session().prediction.data();
}

EX_EXPORT const double* ex_step(void) {
    return session().step.data();
}

EX_EXPORT int ex_step_size(void) {
    return field_total(STEP_FIELDS);
}

EX_EXPORT int ex_step_field_count(void) {
    return STEP_FIELD_COUNT;
}

EX_EXPORT const char* ex_step_field_name(int index) {
    return index >= 0 && index < STEP_FIELD_COUNT ? STEP_FIELDS[index].name : nullptr;
}

EX_EXPORT int ex_step_field_offset(int index) {
    return index >= 0 && index < STEP_FIELD_COUNT ? field_offset(STEP_FIELDS, index) : -1;
}

EX_EXPORT int ex_step_field_rows(int index) {
    return index >= 0 && index < STEP_FIELD_COUNT ? STEP_FIELDS[index].rows : -1;
}

EX_EXPORT int ex_step_field_cols(int index) {
    return index >= 0 && index < STEP_FIELD_COUNT ? STEP_FIELDS[index].cols : -1;
}

}
