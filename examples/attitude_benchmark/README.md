# Recorded-data attitude benchmark

This example compares covariance UKF and square-root UKF with the same attitude
model against Madgwick, Mahony, and online VQF on BROAD and Sassari. Both six-axis
(gyro + accelerometer) and nine-axis (+ magnetometer) configurations are run.
See [the comparison and methodology](../../docs/attitude-benchmark.md).

From the repository root, using Python 3.12 and a C++17 compiler:

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r examples/attitude_benchmark/requirements.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target attitude_benchmark --parallel 2
python -m unittest discover -s examples/attitude_benchmark -p 'test_*.py' -v
python examples/attitude_benchmark/run.py --smoke
python examples/attitude_benchmark/run.py
```

The first run downloads datasets into `.benchmark-cache` (~215 MB for the full
suite). Every file is checked against its SHA-256 and size in
`datasets.lock.json`, including on cache hits. Do not commit dataset binaries.
Network access is needed for initial dependency and data downloads; a verified
populated dataset cache supports offline replay.

Full runs use all 39 BROAD trials and all 18 Sassari speed/sensor combinations.
Smoke runs use the complete BROAD 01 and Sassari slow/XS1 recordings; they are
integration checks, not representative accuracy rankings. Select subsets with
`--datasets broad`, `--axes 6`, or `--algorithms ukf srukf vqf`. Use
`--cache`, `--output`, and `--executable` to change paths.

The runner writes `benchmark-results/comparison.md`, `results.json`, and
`per_trial.csv`. JSON contains source and executable hashes, package versions,
per-trial RMSE/P95, sample counts, alignment windows, and diagnostic timing.
Failed filters, malformed datasets, and checksum failures stop the run instead
of disappearing from the average. Sources must not change during a run.

The C++ executable is independently usable:

```sh
build/examples/attitude_benchmark/attitude_benchmark \
  --input sensors.csv --output attitude.csv --filter srukf --axes 9
```

Input header: `time,gx,gy,gz,ax,ay,az,mx,my,mz`. Time is seconds; gyro is rad/s;
acceleration is specific force in m/s², including gravity; magnetic values may
use any consistent scale. Supply all ten columns even in six-axis mode (zeros
are acceptable for the unused magnetometer). Samples must have finite values,
strictly increasing timestamps, and nonzero used direction vectors.

Output header: `time,qw,qx,qy,qz`. Quaternions rotate body to a right-handed
z-up world; nine-axis initialization puts horizontal magnetic north on +X.
Six-axis world yaw is arbitrary. The first output is sensor-only initialization;
each subsequent output incorporates the sample at that timestamp. Timing JSON
is written to stderr; CSV parsing/writing is outside the timed update loop.

Data attribution and upstream licenses are recorded in `datasets.lock.json`.
BROAD is CC BY 4.0; the Sassari repository provides an MIT license. Reference
implementations remain installed dependencies rather than vendored source:
AHRS 0.4.0 (MIT) and VQF 2.0.1 (MIT).
