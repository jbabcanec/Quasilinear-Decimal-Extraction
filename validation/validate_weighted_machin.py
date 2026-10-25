"""Validate arithmetic, complete cores, decimal intervals, retries, and CLI limits.

Build the native checker sources first; see validation/README.md. Existing
reports are never overwritten. All checks run outside performance samples.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--executable', type=Path, default=ROOT / 'extract_weighted_machin.exe')
    parser.add_argument('--mpfr', type=Path, default=ROOT / 'tmp/weighted_machin_final/mpfr_window_interval.exe')
    args = parser.parse_args()
    env = os.environ.copy()
    if os.name == 'nt':
        env['PATH'] = r'C:\msys64\ucrt64\bin' + os.pathsep + env.get('PATH', '')
    arithmetic = ROOT / 'tmp/weighted_machin_arithmetic_checks_final.exe'
    cores = ROOT / 'tmp/weighted_machin_position_core_checks.exe'
    files = [args.executable.resolve(), args.mpfr.resolve(), arithmetic, cores,
             ROOT / 'decimal_extract_weighted_machin.cpp', Path(__file__).resolve(),
             ROOT / 'validation/weighted_machin_arithmetic_checks.cpp',
             ROOT / 'validation/weighted_machin_position_core_checks.cpp']
    hashes = {str(p): digest(p) for p in files}
    report = {'status': 'running', 'started_utc': datetime.now(timezone.utc).isoformat(),
              'file_sha256_before': hashes, 'native_checks': [], 'windows': [],
              'guard_retries': [], 'block_overrides': [], 'rejected_arguments': []}
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2)

    def run(exe, *arguments, success=True):
        result = subprocess.run([str(exe), *map(str, arguments)], cwd=ROOT,
                                env=env, capture_output=True, text=True, timeout=300)
        if success and result.returncode:
            raise RuntimeError(result.stderr)
        if not success and (result.returncode == 0 or result.stdout.strip()):
            raise AssertionError(f'Invalid arguments accepted: {arguments}')
        return result

    references = {}

    def reference(position):
        if position not in references:
            record = json.loads(run(args.mpfr, position).stdout)
            assert record['method'] == 'MPFR directed-rounding enclosure'
            references[position] = record
        return references[position]['digits']

    try:
        for executable, arguments in [(args.executable, ['--selftest']),
                                      (arithmetic, []), (cores, [])]:
            result = run(executable, *arguments)
            records = [json.loads(line) for line in result.stdout.splitlines()]
            assert all('PASS' in record.values() for record in records)
            report['native_checks'].append({'executable': str(executable), 'results': records})
        for position in [*range(1, 41), 100, 997, 10000]:
            digits = run(args.executable, position).stdout.strip()
            assert digits == reference(position), (position, digits)
            report['windows'].append({'N': position, 'digits': digits, 'MPFR_interval_match': True})
        for guard in [1, 8, 16]:
            result = run(args.executable, 997, '--timing', f'--guard={guard}')
            stats = json.loads(result.stderr)
            assert result.stdout.strip() == reference(997) and stats['passes'] > 1
            report['guard_retries'].append({'initial_guard': guard, 'stats': stats})
        for block in [1, 7, 128, 1024]:
            result = run(args.executable, 997, f'--block={block}')
            assert result.stdout.strip() == reference(997)
            report['block_overrides'].append({'N': 997, 'block': block, 'MPFR_interval_match': True})
        for arguments in [['0'], ['-1'], ['10', '--guard=0'], ['10', '--block=0'],
                          ['10', '--block=-1'], ['10', '--euler'],
                          ['--core-check', '500000001']]:
            result = run(args.executable, *arguments, success=False)
            report['rejected_arguments'].append({'arguments': arguments, 'stderr': result.stderr.strip()})
        after = {str(p): digest(p) for p in files}
        assert after == hashes, 'A checked source or executable changed during validation'
        report.update(status='PASS', file_sha256_after=after, references=references,
                      completed_utc=datetime.now(timezone.utc).isoformat())
    except Exception as error:
        report.update(status='FAIL', error=repr(error))
        raise
    finally:
        args.output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'status': report['status'], 'windows': len(report['windows']),
                      'guard_retry_checks': len(report['guard_retries']),
                      'output': str(args.output.resolve())}))


if __name__ == '__main__':
    main()
