#!/usr/bin/env python3
"""Compare native transport outputs and reconstruct moments from saved histories."""
import argparse
import csv
import json
import math
from pathlib import Path


def rows(path):
    return [[float(value.replace('D', 'E')) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]


def check(a, b, label, *, atol=1.e-20, rtol=1.e-8):
    if not math.isclose(a, b, abs_tol=atol, rel_tol=rtol):
        raise ValueError(f'{label}: {a} != {b}')


def curve(path):
    result = rows(path)
    if result[-1] != [0., 0.]:
        raise ValueError(f'{path.name}: missing terminal sentinel')
    return result[:-1]


def reconstruct(history):
    initial, samples, previous = {}, {}, {}
    for _, marker, step, t, s in history:
        marker, step = int(marker), int(step)
        if step != previous.get(marker, -1) + 1:
            raise ValueError('history is not a retained contiguous prefix')
        previous[marker] = step
        if step == 0:
            initial[marker] = s
        else:
            samples.setdefault(step, []).append((t, (s-initial[marker])**2))
    if len(initial) != 16:
        raise ValueError('the fixed public fixture requires sixteen markers')
    result = []
    for step, values in sorted(samples.items()):
        n = len(values)
        second = math.fsum(q for _, q in values)/n
        fourth = math.fsum(q*q for _, q in values)/n
        error = math.sqrt(max(0., fourth-second*second)/n)
        t = values[0][0]
        for sample_time, _ in values:
            check(sample_time, t, 'history time', atol=1.e-18, rtol=1.e-12)
        result.append((step, t, n, second, error))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scalar', type=Path, required=True)
    parser.add_argument('--cuda', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('transport.csv'))
    args = parser.parse_args()
    cpu_history = rows(args.scalar/'particle_histories.dat')
    gpu_history = rows(args.cuda/'particle_histories.dat')
    if len(cpu_history) != len(gpu_history):
        raise ValueError('CPU/CUDA history lengths differ')
    for cpu, gpu in zip(cpu_history, gpu_history):
        if cpu[:3] != gpu[:3]:
            raise ValueError('CPU/CUDA history identifiers differ')
        for a, b in zip(cpu[3:], gpu[3:]):
            check(a, b, 'CPU/CUDA history', atol=1.e-12)
    for directory in (args.scalar, args.cuda):
        if rows(directory/'n_lost_particles.dat') != [[9.]]:
            raise ValueError('the fixed fixture requires nine physical losses')
    cpu_events = rows(args.scalar/'lost_particle_events.dat')
    gpu_events = rows(args.cuda/'lost_particle_events.dat')
    if len(cpu_events) != 9 or len(gpu_events) != 9:
        raise ValueError('incorrect loss-event count')
    for cpu, gpu in zip(cpu_events, gpu_events):
        if cpu[1:4] != gpu[1:4]:
            raise ValueError('CPU/CUDA lost marker, reason or step differs')
        for a, b in zip(cpu, gpu):
            check(a, b, 'CPU/CUDA loss event', atol=1.e-12)
    reference = reconstruct(gpu_history)
    curves = [curve(directory/name) for directory in (args.scalar, args.cuda)
              for name in ('mean_square_displacement.dat', 'mean_square_uncertainty.dat')]
    if any(len(values) != len(reference) for values in curves):
        raise ValueError('reported curves do not cover retained measurement steps')
    output = []
    max_cpu, max_cuda = 0., 0.
    for expected, cpu_mean, cpu_error, gpu_mean, gpu_error in zip(reference, *curves):
        step, t, n, mean, error = expected
        for measured in (cpu_mean, cpu_error, gpu_mean, gpu_error):
            check(measured[0], t, 'reported time', atol=1.e-18, rtol=1.e-12)
        for measured, target in ((cpu_mean[1], mean), (gpu_mean[1], mean),
                                 (cpu_error[1], error), (gpu_error[1], error)):
            check(measured, target, 'independent retained-history estimator')
        max_cpu = max(max_cpu, abs(cpu_mean[1]-mean), abs(cpu_error[1]-error))
        max_cuda = max(max_cuda, abs(gpu_mean[1]-mean), abs(gpu_error[1]-error))
        output.append((step, t, n, cpu_mean[1], gpu_mean[1], mean,
                       cpu_error[1], gpu_error[1], error))
    with args.output.open('w', newline='') as destination:
        writer = csv.writer(destination)
        writer.writerow(('step', 'time_seconds', 'contributors', 'scalar_second',
                         'cuda_second', 'history_second', 'scalar_standard_error',
                         'cuda_standard_error', 'history_standard_error'))
        writer.writerows(output)
    print(json.dumps(dict(measurements=len(output), physical_losses=9,
                          maximum_scalar_reference_error=max_cpu,
                          maximum_cuda_reference_error=max_cuda)))


if __name__ == '__main__':
    main()
