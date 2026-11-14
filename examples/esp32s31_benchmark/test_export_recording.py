"""Exporter checks: timing fidelity, bounded input, and C++ representation."""

import hashlib
import importlib.util
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("export_recording", HERE / "export_recording.py")
EXPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXPORT)


class ExportRecordingTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.source = self.directory / "sensors.csv"
        self.output = self.directory / "recording.h"

    def write(self, times, sensor="0,0,0,0,0,9.80665,0.6,0,-0.8"):
        self.source.write_text(",".join(EXPORT.HEADER) + "\n" +
                               "".join(f"{time},{sensor}\n" for time in times))

    def test_rebases_epoch_without_retiming_irregular_samples(self):
        self.write(["2000000000.000000000", "2000000000.002000000",
                    "2000000000.003500000", "2000000000.007000000"])
        rows, intervals, digest = EXPORT.load_recording(self.source, start=1, count=3)
        self.assertEqual([row[0] for row in rows], [0.0, 0.0015, 0.005])
        self.assertAlmostEqual(intervals[1], 0.0035)
        self.assertEqual(digest, hashlib.sha256(self.source.read_bytes()).hexdigest())
        EXPORT.export_recording(self.source, self.output, start=1, count=3)
        header = self.output.read_text()
        self.assertIn("source_start_row = 1;", header)
        self.assertIn("dt_min_seconds", header)
        self.assertIn("dt_mean_seconds", header)
        self.assertIn("no resampling or retiming", header)

    def test_does_not_silently_truncate(self):
        self.write(range(EXPORT.MAX_SAMPLES + 1))
        with self.assertRaisesRegex(ValueError, "use --start/--count"):
            EXPORT.load_recording(self.source)
        rows, _, _ = EXPORT.load_recording(self.source, count=EXPORT.MAX_SAMPLES)
        self.assertEqual(len(rows), EXPORT.MAX_SAMPLES)

    def test_rejects_invalid_data(self):
        cases = [
            ("time,gx\n0,0\n1,0\n", "header"),
            (",".join(EXPORT.HEADER) + "\n0,0\n", "ten finite"),
            (",".join(EXPORT.HEADER) + "\n0,0,0,0,0,0,1,1,1,nan\n", "ten finite"),
        ]
        for text, message in cases:
            with self.subTest(text=text):
                self.source.write_text(text)
                with self.assertRaisesRegex(ValueError, message):
                    EXPORT.load_recording(self.source)
        for times in ([0, 0], [1, 0]):
            self.write(times)
            with self.assertRaisesRegex(ValueError, "Time must increase"):
                EXPORT.load_recording(self.source)
        self.write([0, 1], sensor="1e39,0,0,0,0,1,1,0,0")
        with self.assertRaisesRegex(ValueError, "finite float"):
            EXPORT.load_recording(self.source)
        self.write([0, "1e-50"])
        with self.assertRaisesRegex(ValueError, "positive in single precision"):
            EXPORT.load_recording(self.source)

    def test_rejects_short_and_out_of_range_selection(self):
        self.write([0, 1, 2])
        for kwargs in ({"start": -1}, {"count": 1}, {"count": EXPORT.MAX_SAMPLES + 1},
                       {"start": 2}, {"count": 4}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                EXPORT.load_recording(self.source, **kwargs)

    def test_generated_header_is_valid_cpp_with_unusual_filename(self):
        compiler = shutil.which("c++")
        if compiler is None:
            self.skipTest("No C++ compiler")
        self.source = self.directory / 'sensors "é😀\\u00e9\n.csv'
        self.write([0, 0.002])
        EXPORT.export_recording(self.source, self.output)
        program = self.directory / "check.cpp"
        program.write_text('#include "recording.h"\n#include <cstdio>\n'
                           'int main() { using namespace AttitudeBenchmark::Recording;\n'
                           'static_assert(sizeof(samples[0].gyroscope[0]) == 4);\n'
                           'static_assert(samples[1].time == 0.002);\n'
                           'std::printf("%s", name); }\n')
        executable = self.directory / "check"
        build = subprocess.run([compiler, "-std=c++17", str(program), "-o", str(executable)],
                               capture_output=True, text=True)
        self.assertEqual(build.returncode, 0, build.stderr)
        run = subprocess.run([str(executable)], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0)
        self.assertEqual(run.stdout, self.source.name)

    def test_cli_budget_note_and_no_input_overwrite(self):
        self.write([0, 0.01])
        import sys
        result = subprocess.run([sys.executable, str(HERE / "export_recording.py"),
                                 str(self.source), str(self.output)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("does not imply", result.stdout)
        original = self.source.read_bytes()
        result = subprocess.run([sys.executable, str(HERE / "export_recording.py"),
                                 str(self.source), str(self.source)], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.source.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
