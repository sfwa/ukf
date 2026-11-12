# Interactive explainer

A long-form article that explains the unscented Kalman filter with
interactive figures, from a single noisy sensor up to the multirotor
estimator of the [EuRoC benchmark](../docs/euroc-benchmark.md). The
estimator runs in the reader's browser, on recorded flights.

The page is a static site: HTML, CSS, JavaScript modules, one WebAssembly
file and a few data files. It has no framework, no bundler and no server
side, and it loads nothing from other sites.

## What runs in the browser

The filter in the page is `MAV_Estimator` from
[`examples/euroc_benchmark/mav_estimator.h`](../examples/euroc_benchmark/mav_estimator.h),
compiled to WebAssembly. It is not copied or rewritten.

- [`wasm/explainer_api.h`](wasm/explainer_api.h) is the whole interface
  between the page and the estimator. The page gives it IMU samples and
  tuning values, and reads back the estimate and its uncertainty at every
  sample.
- [`wasm/explainer_api.cpp`](wasm/explainer_api.cpp) runs the estimator with
  the same loop as the benchmark program. It can also take one filter step
  apart for the figures, using only the public functions of the library, and
  it checks the result against the real step.
- The estimator gained one read-only accessor, `MAV_Estimator::ukf()`, so that
  the page can show the covariance. Nothing else in the benchmark changed.

The input is also the benchmark's. `tools/pack_euroc.py` reads the pinned
EuRoC files through the benchmark's own loader and stores the IMU samples as
the sensor's integer counts. The page decodes them to the same doubles that
the benchmark reads, bit for bit.

## Build and view

You need CMake, Python 3 with numpy, and
[Emscripten](https://emscripten.org/docs/getting_started/downloads.html):

```sh
git clone https://github.com/emscripten-core/emsdk.git
emsdk/emsdk install 7.0.0
emsdk/emsdk activate 7.0.0
source emsdk/emsdk_env.sh
```

Then, from the repository root:

```sh
python explainer/build.py
python -m http.server --directory explainer/_site
```

Open <http://localhost:8000>. The first build downloads about 40 MB of EuRoC
files into `.benchmark-cache/euroc`, the same cache that the benchmark uses.

Opening `index.html` straight from the disk does not work. Browsers do not
load WebAssembly or JavaScript modules from `file://` addresses.

While you work on the page, rebuild only what changed:

```sh
python explainer/build.py --skip-wasm --skip-data --dev
```

`--dev` also copies two test pages, `dev/concepts.html` and `dev/real.html`,
which show the figures without the article around them.

## Tests

The tests need the site built as above, and Node.js 22 or newer.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUKF_BUILD_EXPLAINER=ON
cmake --build build --target euroc_benchmark explainer_check
EXPLAINER_REQUIRE_ALL=1 python -m unittest discover -s explainer/test -v
npm install --no-save playwright@1.63.0
npx playwright install chromium
node explainer/test/smoke.cjs
```

A Python test is skipped if what it needs has not been built.
`EXPLAINER_REQUIRE_ALL=1` turns a skip into a failure, as in CI.

The Python tests hold the page to its claim that it runs the benchmark's
estimator on the benchmark's inputs:

- the packed recordings decode to the numbers the benchmark reads;
- the interface, built natively, writes the same attitude file as
  `euroc_benchmark`, byte for byte;
- the WebAssembly module, driven by the page's own JavaScript in Node, agrees
  with the native program to rounding error;
- a step taken apart for the figures agrees with the real step;
- the page's tilt and heading errors are the benchmark's.

`smoke.cjs` opens the site in Chromium (it needs the `playwright` package).
It checks that every figure draws, that nothing reports an error, and that
every number written in the article is the number the estimator gives. If
the estimator changes, this test says which sentences to update.

## Publishing

[`.github/workflows/explainer.yml`](../.github/workflows/explainer.yml)
builds and tests the site on pull requests, and publishes it with GitHub
Pages on every push to `master`.

One setting is needed, once: in the repository's settings, under Pages, set
the source to "GitHub Actions". The site then appears at
`https://<owner>.github.io/<repository>/`. The page uses only relative
addresses, so it works under any path.

## The recordings

The page ships five EuRoC recordings by default; the first one in the list
is the one the article is written around. Choose others with

```sh
python explainer/build.py --sequences V1_02_medium MH_01_easy V2_03_difficult
```

Two figures use `MH_01_easy` by name, so keep it in the list.

The recordings are not in this repository. The build downloads them and
checks them against `examples/euroc_benchmark/datasets.lock.json`. They are
part of the [EuRoC MAV dataset](https://doi.org/10.1177/0278364915620033),
which ETH Zurich publishes as
[In Copyright, Non-Commercial Use Permitted](https://doi.org/10.3929/ethz-b-000690084).
A published copy of this site contains the IMU samples and the reference
trajectory of the chosen recordings. The article names the dataset and its
terms. The library's licence does not apply to the dataset.

## Layout

```
explainer/
  build.py              builds the site into explainer/_site
  CMakeLists.txt        the WebAssembly module, and its native check
  wasm/                 the C interface to the estimator
  tools/pack_euroc.py   packs the recordings; describes the file format
  web/                  the page
    index.html          the article
    css/style.css       styles and the colours of both themes
    js/core/            figures, controls, plots, 3D drawing
    js/data/            recordings, the estimator, the benchmark's metrics
    js/figures/         one file for each figure
    dev/                test pages for the figures
  test/                 the tests described above
```

## Writing

The article is written in plain language on purpose: short sentences, the
active voice, everyday words, and one name for each thing. A term is
explained where it first appears. Please keep to this when you change the
text.

A number in the text that depends on the estimator or on a recording is
written as `<span class="live" data-live="name">value</span>`. The page
replaces the value with the real one when it loads
(`web/js/figures/facts.js`), and the browser test fails if the two differ.

In every figure, blue is what the filter believes, orange is what a sensor
reads, and green is what a model predicts. Everything else is neutral.
