/*
C interface between the interactive explainer and the multirotor estimator
of the EuRoC benchmark (examples/euroc_benchmark/mav_estimator.h).

The explainer page runs the estimator in the browser. This file is the whole
surface between the two: the page writes recorded IMU samples and tuning
values into the module's memory, asks for a run, and reads the results back.
The estimator itself is not changed or copied; explainer_api.cpp calls
MAV_Estimator::initialise() and MAV_Estimator::update(), as the benchmark
program does.

All values are doubles. Matrices are stored column by column. Quaternions
are written scalar first (w, x, y, z), as in the benchmark's output files.

The two result buffers describe their own layout (see ex_trace_field_* and
ex_step_field_*), so the page looks fields up by name and never hard-codes
an offset.
*/

#ifndef EXPLAINER_API_H
#define EXPLAINER_API_H

#ifdef __cplusplus
extern "C" {
#endif

/* Changes whenever the meaning of this interface changes. */
#define EX_ABI_VERSION (1)

/* Numbers per IMU sample: time, gyroscope x y z, accelerometer x y z. */
#define EX_SAMPLE_STRIDE (7)

/* The estimator is saved every this many samples, so that ex_explain() can
restart close to any step instead of at the beginning of the recording. */
#define EX_CHECKPOINT_INTERVAL (200)

int ex_abi_version(void);

/*
Model constants which are fixed in mav_estimator.h, so that the page shows
the true values: see ex_constant_name() for the list.
*/
int ex_constant_count(void);
const char* ex_constant_name(int index);
double ex_constant(int index);

/*
The recording. ex_samples() makes room for "count" samples and returns the
buffer to fill, EX_SAMPLE_STRIDE numbers per sample. Returns null if there
is not enough memory. Calling it discards earlier results.
*/
double* ex_samples(int count);
int ex_sample_count(void);

/*
Tuning: every numeric member of MAV_Tuning, by name. ex_tuning() returns
the values, which the page may change before a run. A boolean member is 0
or 1. ex_tuning_reset() restores the benchmark's values.
*/
int ex_tuning_count(void);
const char* ex_tuning_name(int index);
double* ex_tuning(void);
void ex_tuning_reset(void);

/*
Run the estimator over the whole recording with the current tuning. Returns
the number of samples which have an estimate: the sample count if every
step succeeded, less if a filter step failed (the first sample without an
estimate is the one whose step failed), or -1 if there is nothing to run.
*/
int ex_run(void);

/*
The result of the last run: one record per sample. Records after a failed
step are filled with NaN.
*/
const double* ex_trace(void);
int ex_trace_stride(void);
int ex_trace_field_count(void);
const char* ex_trace_field_name(int index);
int ex_trace_field_offset(int index);
int ex_trace_field_size(int index);

/*
One filter step taken apart. ex_explain(index) repeats the step which moved
the estimate from sample index-1 to sample index, and fills the step buffer
with every intermediate value: the sigma points before and after the
process model, the predicted state, the expected measurements, the
innovation, the gain and the corrected state.

The real step is run as well, and the two results are compared. The
largest difference is stored in the "check" field; it should be at the
level of rounding error.

Returns 0 on success, 1 if the index is out of range or there is no run,
and 2 if the step fails.
*/
int ex_explain(int index);
const double* ex_step(void);
int ex_step_size(void);
int ex_step_field_count(void);
const char* ex_step_field_name(int index);
int ex_step_field_offset(int index);
int ex_step_field_rows(int index);
int ex_step_field_cols(int index);

/*
Looking ahead. ex_predict(index, steps) takes the estimate at sample
"index", makes its sigma points, and moves each of them through the process
model for the next "steps" samples of the recording. No measurement is used
and no process noise is added, so this shows what the process model alone
does to each sigma point. Each step uses the IMU samples of the recording
and the rotor model decision of the last run.

The prediction buffer then holds one snapshot per sample, starting with the
sigma points at "index" itself. A snapshot is the sigma points one after
another, each written as a state (attitude w x y z, velocity, gyro bias).

Returns the number of snapshots (steps + 1, or fewer at the end of the
recording), or 0 if the index is out of range or there is no run.
*/
int ex_predict(int index, int steps);
const double* ex_prediction(void);

#ifdef __cplusplus
}
#endif

#endif
