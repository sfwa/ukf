#!/usr/bin/env python3
"""Build the explainer site: one folder of static files, ready to serve.

    python explainer/build.py

does three things:

  1. compiles the estimator to WebAssembly (needs Emscripten's emcmake on
     the PATH, and CMake),
  2. packs the EuRoC recordings (needs numpy; downloads the pinned files on
     the first run, about 30 MB),
  3. copies the page (explainer/web) next to them.

The result is explainer/_site. Serve it with any static web server, for
example:

    python -m http.server --directory explainer/_site

and open http://localhost:8000. Opening index.html straight from the disk
does not work: browsers do not load WebAssembly or JavaScript modules from
file:// addresses.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
WASM_NAME = "ukf_explainer.wasm"


def run(command, **options):
    print("+", " ".join(str(part) for part in command), flush=True)
    subprocess.run(command, check=True, **options)


def build_wasm(build_directory: Path, output: Path) -> None:
    """Compile the estimator and its interface (explainer/wasm) with Emscripten."""
    if shutil.which("emcmake") is None:
        raise SystemExit("emcmake is not on the PATH. Install Emscripten and activate it "
                         "(see explainer/README.md), or pass --skip-wasm to keep an existing module.")
    run(["emcmake", "cmake", "-S", ROOT, "-B", build_directory,
         "-DCMAKE_BUILD_TYPE=Release",
         "-DUKF_BUILD_EXPLAINER=ON",
         "-DUKF_BUILD_TESTS=OFF", "-DUKF_BUILD_BENCHMARKS=OFF", "-DUKF_BUILD_EXAMPLES=OFF"])
    run(["cmake", "--build", build_directory, "--target", "ukf_explainer", "--parallel", "2"])
    shutil.copyfile(build_directory / "explainer" / WASM_NAME, output / WASM_NAME)


def pack_data(output: Path, cache: Path, sequences: list[str] | None) -> None:
    """Pack the recordings with explainer/tools/pack_euroc.py."""
    command = [sys.executable, HERE / "tools" / "pack_euroc.py", "--output", output / "data", "--cache", cache]
    if sequences:
        command += ["--sequences", *sequences]
    run(command)


def copy_page(output: Path, include_dev: bool) -> None:
    """Copy explainer/web. The test pages under dev/ are left out unless asked for."""
    def ignore(directory, names):
        if include_dev:
            return []
        return [name for name in names if name == "dev" or name.startswith("dev-")]
    shutil.copytree(HERE / "web", output, dirs_exist_ok=True, ignore=ignore)
    # GitHub Pages must serve the files as they are.
    (output / ".nojekyll").write_text("")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--output", type=Path, default=HERE / "_site", help="where to put the site (default: explainer/_site)")
    parser.add_argument("--wasm-build", type=Path, default=ROOT / "build-wasm", help="CMake build folder for the WebAssembly module")
    parser.add_argument("--cache", type=Path, default=ROOT / ".benchmark-cache" / "euroc",
                        help="where verified EuRoC files are kept (default: the benchmark's cache)")
    parser.add_argument("--sequences", nargs="+", help="recordings to include; the first is the one the page opens with")
    parser.add_argument("--skip-wasm", action="store_true", help="keep the module already in the output folder")
    parser.add_argument("--skip-data", action="store_true", help="keep the recordings already in the output folder")
    parser.add_argument("--dev", action="store_true", help="also copy the figure test pages (explainer/web/dev)")
    args = parser.parse_args()

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    copy_page(output, args.dev)
    if not args.skip_wasm:
        build_wasm(args.wasm_build.resolve(), output)
    if not args.skip_data:
        pack_data(output, args.cache.resolve(), args.sequences)

    missing = [name for name in (WASM_NAME, "data/recordings.json") if not (output / name).is_file()]
    if missing:
        raise SystemExit(f"The site is not complete: {', '.join(missing)} missing from {output}")
    print(f"The site is in {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
