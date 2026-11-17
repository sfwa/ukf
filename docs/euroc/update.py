#!/usr/bin/env python3
"""Record a complete, consistently configured EuRoC replay in documentation."""
import argparse
import csv
import datetime
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ALGORITHMS = {'ukf', 'srukf', 'mav-srukf-gravity', 'mav-srukf', 'madgwick', 'mahony', 'vqf', 'vqf-mav'}


def validate(report):
    manifest = json.loads((ROOT / 'examples/euroc_benchmark/datasets.lock.json').read_text())
    sequences = set(manifest['sequences'])
    if (report['scope'] != 'full EuRoC dataset' or report['split'] != 'all'
            or report['sequence_patterns'] or report['axes'] != [6]
            or set(report['algorithms']) != ALGORITHMS
            or set(report['selected_sequences']) != sequences):
        raise ValueError(f'Only a complete, six-axis, {len(ALGORITHMS)}-algorithm EuRoC replay can replace the snapshot')
    rows = report['per_sequence']
    wanted = {(name, algorithm) for name in sequences for algorithm in ALGORITHMS}
    actual = {(row['sequence'], row['algorithm']) for row in rows}
    if actual != wanted or len(rows) != len(wanted) or any(
            row['dataset'] != 'euroc' or row['axes'] != 6 for row in rows):
        raise ValueError('Missing, duplicate or unexpected sequence/algorithm results')
    if not report['source_sha256'] or set(report['binary_sha256']) != {'attitude_benchmark', 'euroc_benchmark'}:
        raise ValueError('Source and both native executable hashes are required')
    return rows


def table(summary):
    lines = ['| Algorithm | Sequences | Inclination | Heading | Total |',
             '|---|---:|---:|---:|---:|']
    for row in summary:
        lines.append(f"| {row['algorithm']} | {row['sequences']} | "
                     f"{row['inclination_rmse_deg']:.3f} | {row['heading_rmse_deg']:.3f} | "
                     f"{row['total_rmse_deg']:.3f} |")
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    parser.add_argument('--build-description', required=True)
    args = parser.parse_args()
    report = json.loads(args.report.read_text())
    rows = validate(report)
    destination = Path(__file__).resolve().parent
    metadata = {key: value for key, value in report.items() if key != 'per_sequence'}
    metadata['build_description'] = args.build_description
    (destination / 'provenance.json').write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n')
    with (destination / 'per-sequence.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    lines = ['# Recorded EuRoC MAV attitude comparison', '',
             f'Recorded {datetime.date.today().isoformat()}. **{len(rows)} evaluations**: all eleven sequences, '
             f'{len(ALGORITHMS)} algorithms/variants, six-axis IMU input only.', '',
             'Values are arithmetic means of per-sequence RMSEs in degrees. Heading and total error use '
             'one fixed initial yaw alignment; absolute heading is unobservable. The primary process-model '
             'target is inclination accuracy.', '',
             '## Primary validation: externally constrained Vicon attitude', '',
             'V1_03_difficult and V2_01_easy, V2_02_medium, V2_03_difficult were not used for tuning; the '
             'methodology records two structural choices made after validation results had been seen. '
             'The reference combines Vicon poses and IMU in an offline batch estimate; it is not wholly '
             'independent of the input IMU.', '']
    lines += table(report['summary_primary_validation'])
    lines += ['', '## Primary Vicon subset, including development', '',
              'The five Vicon sequences excluding V1_01_easy. V1_02_medium is a development recording.', '']
    lines += table(report['summary_primary'])
    lines += ['', '## All eleven recordings (descriptive)', '',
              'Includes five Machine Hall recordings with no independent optical orientation observations '
              'and V1_01_easy with its documented canonical-reference orientation issue.', '']
    lines += table(report['summary'])
    lines += ['', '## Per-sequence model comparison', '',
              'Inclination RMSE in degrees. The gravity control is the MAV estimator with its rotor model '
              'switched off; the default and adapted VQF configurations remain separately visible.', '',
              '| Sequence | Split | Gravity control | MAV SR-UKF | VQF | Adapted VQF |',
              '|---|---|---:|---:|---:|---:|']
    for sequence in report['selected_sequences']:
        selected = {row['algorithm']: row for row in rows if row['sequence'] == sequence}
        lines.append(f"| {sequence} | {selected['mav-srukf']['split']} | "
                     f"{selected['mav-srukf-gravity']['inclination_rmse_deg']:.3f} | "
                     f"{selected['mav-srukf']['inclination_rmse_deg']:.3f} | "
                     f"{selected['vqf']['inclination_rmse_deg']:.3f} | "
                     f"{selected['vqf-mav']['inclination_rmse_deg']:.3f} |")
    lines += ['', 'Development used MH_01_easy and V1_02_medium. All estimators process complete sensor '
              'recordings; no ground-truth position, velocity or bias enters an estimator, and no '
              'reference-derived flight mask selects favorable intervals. Configuration changes and '
              'model limitations are described in the methodology.', '',
              '[Methodology](../euroc-benchmark.md) · [Per-sequence metrics](per-sequence.csv) · '
              '[Configuration and provenance](provenance.json)', '',
              args.build_description+' The recorded source and executable hashes identify the replayed '
              'implementation. Timing is diagnostic, not a cross-language speed ranking.', '']
    (destination / 'comparison.md').write_text('\n'.join(lines))


if __name__ == '__main__':
    main()
