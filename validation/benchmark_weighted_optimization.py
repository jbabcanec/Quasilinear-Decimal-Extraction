"""Compare native extractors in fresh processes and independently verify each window.

Example (from the project root)::

    python validation/benchmark_weighted_optimization.py \
        --engine original=extract_weighted.exe \
        --engine weighted_fast=extract_weighted_optimized.exe \
        --engine euler=extract_v5.exe --positions 10000 100000 1000000 \
        --repeats 3 --output validation/weighted_optimization_run1.json \
        --check-baseline-core extract_weighted.exe --core-engine weighted_fast

Each engine must accept POSITION --timing, print a 12-digit window to stdout,
and emit a JSON timing object to stderr. Core checks require --core POSITION.
Repeat --engine-arg NAME=ARG to append engine-specific benchmark options,
for example --engine-arg weighted_fast=--block=15625. Core checks do not use them.
Only the Python standard library is used. Existing reports are never overwritten.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
PHASES = {"pre", "word", "core", "vphase", "total", "factor", "coefficients",
          "real", "rational", "assembly"}
CORE_POSITIONS = (100, 997, 10000, 100000)
DEFAULT_SOURCES = (
    "decimal_extract_weighted.cpp", "decimal_extract_v3.cpp",
    "decimal_extract_v5.cpp", "decimal_extract_shared_residues.cpp",
    "validation/mpfr_window_interval.cpp",
    "validation/benchmark_weighted_optimization.py",
)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_path(value):
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"File does not exist: {path}")
    return path


def parse_engines(values):
    engines = {}
    for value in values:
        name, separator, filename = value.partition("=")
        if not separator or not name.strip() or not filename:
            raise ValueError(f"Expected --engine NAME=PATH, got {value!r}")
        if name in engines:
            raise ValueError(f"Duplicate engine name: {name}")
        engines[name] = checked_path(filename)
    return engines


def parse_engine_arguments(values, engines):
    arguments = {name: [] for name in engines}
    for value in values:
        name, separator, argument = value.partition("=")
        if not separator or not argument:
            raise ValueError(f"Expected --engine-arg NAME=ARG, got {value!r}")
        if name not in arguments:
            raise ValueError(f"Unknown engine in --engine-arg: {name}")
        arguments[name].append(argument)
    return arguments


def call(executable, arguments, env, timeout):
    command = [str(executable), *map(str, arguments)]
    start = time.perf_counter()
    try:
        result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True,
                                text=True, timeout=timeout, check=False)
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(f"Timed out after {timeout:g}s: {command}") from error
    wall = time.perf_counter() - start
    if result.returncode:
        raise RuntimeError(f"Exit {result.returncode}: {command}\n"
                           f"stdout: {result.stdout[-2000:]}\n"
                           f"stderr: {result.stderr[-4000:]}")
    return result.stdout.strip(), result.stderr.strip(), wall


def timing_json(stderr):
    # Accept one complete object or a JSON line alongside diagnostic text.
    candidates = [stderr, *reversed(stderr.splitlines())]
    for candidate in candidates:
        try:
            stats = json.loads(candidate)
        except json.JSONDecodeError:
            continue
        if isinstance(stats, dict) and "total" in stats:
            if not numeric_time(stats["total"]):
                raise ValueError("Engine timing 'total' must be finite and nonnegative")
            return stats
    raise ValueError(f"No JSON timing object with 'total' on stderr: {stderr[-4000:]}")


def numeric_time(value):
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value) and value >= 0)


def timed_run(name, executable, position, repeat, order, env, timeout, engine_arguments):
    digits, stderr, wall = call(executable, [position, "--timing", *engine_arguments], env, timeout)
    if not re.fullmatch(r"[0-9]{12}", digits):
        raise ValueError(f"{name} at {position}: expected 12 digits, got {digits!r}")
    stats = timing_json(stderr)
    if "N" in stats and stats["N"] != position:
        raise ValueError(f"{name}: timing position {stats['N']} != {position}")
    return {"engine": name, "N": position, "repeat": repeat, "order": order,
            "digits": digits, "stats": stats, "stderr": stderr,
            "engine_arguments": engine_arguments,
            "wall_seconds": wall}


def interval_reference(executable, position, env, timeout):
    stdout, stderr, wall = call(executable, [position], env, timeout)
    reference = json.loads(stdout)
    if (reference.get("N") != position
            or not re.fullmatch(r"[0-9]{12}", reference.get("digits", ""))
            or reference.get("method") != "MPFR directed-rounding enclosure"):
        raise ValueError(f"Unexpected independent interval reference: {reference!r}")
    reference["verification_wall_seconds"] = wall
    if stderr:
        reference["stderr"] = stderr
    return reference


def exact_core(executable, position, env, timeout):
    stdout, _, _ = call(executable, ["--core", position], env, timeout)
    core = json.loads(stdout)
    precision = core.get("bits")
    value = core.get("V_hex", "")
    if (not isinstance(precision, int) or isinstance(precision, bool)
            or precision < 1 or not re.fullmatch(r"[0-9a-fA-F]+", value)):
        raise ValueError(f"Invalid --core response from {executable}")
    residue = int(value, 16)
    if residue.bit_length() > precision:
        raise ValueError(f"Unreduced --core residue from {executable}")
    return precision, residue


def summarize(samples, reference):
    stats = [sample["stats"] for sample in samples]
    phase_keys = set.intersection(*(set(record) for record in stats))
    phase_keys = sorted(key for key in phase_keys
                        if (key in PHASES or isinstance(stats[0][key], float))
                        and all(numeric_time(record[key]) for record in stats))
    compute = [record["total"] for record in stats]
    walls = [sample["wall_seconds"] for sample in samples]
    return {
        "engine": samples[0]["engine"], "N": samples[0]["N"],
        "repetitions": len(samples), "digits": reference["digits"],
        "compute_median": statistics.median(compute),
        "compute_min": min(compute), "compute_max": max(compute),
        "wall_median": statistics.median(walls),
        "wall_min": min(walls), "wall_max": max(walls),
        "phase_medians": {key: statistics.median(record[key] for record in stats)
                          for key in phase_keys},
        "MPFR_interval_match": True,
    }


def benchmark(args):
    engines = parse_engines(args.engine)
    engine_arguments = parse_engine_arguments(args.engine_arg, engines)
    verifier = checked_path(args.mpfr_exe)
    baseline = checked_path(args.check_baseline_core) if args.check_baseline_core else None
    sources = ([checked_path(path) for paths in args.source for path in paths]
               if args.source else [ROOT / path for path in DEFAULT_SOURCES
                                    if (ROOT / path).is_file()])
    sources = sorted(set([*sources, Path(__file__).resolve()]), key=str)
    core_engines = args.core_engine or [name for name in engines if "weighted" in name.lower()]
    if baseline and not core_engines:
        raise ValueError("Specify --core-engine NAME for the exact-core comparison")
    if any(name not in engines for name in core_engines):
        raise ValueError("Each --core-engine must name an --engine entry")
    if args.core_engine and not baseline:
        raise ValueError("--core-engine requires --check-baseline-core")

    output = args.output.expanduser().resolve()
    env = os.environ.copy()
    runtime_dirs = [str(Path(path).expanduser().resolve()) for path in args.runtime_dir]
    if runtime_dirs:
        env["PATH"] = os.pathsep.join([*runtime_dirs, env.get("PATH", "")])
    tracked = set([*engines.values(), verifier, *sources])
    if baseline:
        tracked.add(baseline)
    hashes = {str(path): sha256(path) for path in sorted(tracked, key=str)}
    report = {
        "schema_version": 1, "started_utc": utc_now(), "verified": False,
        "status": "running", "output": str(output),
        "method": "Sequential fresh processes; rotate engine order each repetition and "
                  "position; no discarded warm-ups. Engine JSON times exclude process "
                  "startup; wall times include it. Exact core comparisons and independent "
                  "MPFR directed-rounding checks run outside benchmark samples. Per-engine "
                  "options follow --timing and are not passed to exact core checks.",
        "command": sys.argv, "python": sys.version,
        "platform": platform.platform(), "processor": platform.processor(),
        "runtime_directories_prepended": runtime_dirs,
        "engines": {name: {"executable": str(path), "sha256": hashes[str(path)],
                           "arguments": engine_arguments[name]}
                    for name, path in engines.items()},
        "reference_executable": str(verifier),
        "core_baseline": str(baseline) if baseline else None,
        "source_sha256": {str(path): hashes[str(path)] for path in sources},
        "file_sha256_before": hashes, "positions": args.positions,
        "repeats": args.repeats, "timeout_seconds": args.timeout,
        "reference_timeout_seconds": args.reference_timeout,
        "runs": [], "summary": [], "references": {}, "core_checks": [],
    }

    # Exclusive creation protects historical reports, including the default name.
    with output.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")

    def save():
        output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    try:
        if baseline:
            for position in CORE_POSITIONS:
                wanted = exact_core(baseline, position, env, args.timeout)
                for name in core_engines:
                    actual = exact_core(engines[name], position, env, args.timeout)
                    check = {"engine": name, "N": position, "bits": actual[0],
                             "exact_match": actual == wanted,
                             "residue_hex_sha256": hashlib.sha256(
                                 format(actual[1], "x").encode("ascii")).hexdigest()}
                    report["core_checks"].append(check)
                    save()
                    print(json.dumps({"core_check": check}), flush=True)
                    if not check["exact_match"]:
                        raise ValueError(f"Exact core mismatch: {name} at {position}")

        names = list(engines)
        for position_index, position in enumerate(args.positions):
            samples = {name: [] for name in names}
            for repeat in range(args.repeats):
                offset = (position_index + repeat) % len(names)
                ordering = names[offset:] + names[:offset]
                for order, name in enumerate(ordering):
                    sample = timed_run(name, engines[name], position, repeat, order,
                                       env, args.timeout, engine_arguments[name])
                    samples[name].append(sample)
                    report["runs"].append(sample)
                    save()
                    print(json.dumps(sample), flush=True)
            reference = interval_reference(verifier, position, env, args.reference_timeout)
            report["references"][str(position)] = reference
            save()
            for name in names:
                if any(sample["digits"] != reference["digits"] for sample in samples[name]):
                    raise ValueError(f"MPFR interval mismatch: {name} at {position}; "
                                     f"expected {reference['digits']}")
            summaries = [summarize(samples[name], reference) for name in names]
            for summary in summaries:
                summary["speedup_reference_engine"] = names[0]
                for kind in ("compute", "wall"):
                    duration = summary[f"{kind}_median"]
                    summary[f"{kind}_speedup_vs_first_engine"] = (
                        summaries[0][f"{kind}_median"] / duration if duration else None)
                report["summary"].append(summary)
            save()
            print(json.dumps({"verified_position": position, "summary": summaries}), flush=True)

        hashes_after = {str(path): sha256(path) for path in sorted(tracked, key=str)}
        report["file_sha256_after"] = hashes_after
        report["files_changed_during_run"] = [path for path in hashes if hashes[path] != hashes_after[path]]
        if report["files_changed_during_run"]:
            raise ValueError("Sources or executables changed during the benchmark: "
                             + ", ".join(report["files_changed_during_run"]))
        report.update(verified=True, status="complete", completed_utc=utc_now())
        save()
        print(json.dumps({"verified": True, "output": str(output)}), flush=True)
    except Exception as error:
        report.update(status="failed", failed_utc=utc_now(), error=f"{type(error).__name__}: {error}")
        save()
        raise
    except KeyboardInterrupt:
        report.update(status="interrupted", failed_utc=utc_now(), error="KeyboardInterrupt")
        save()
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--engine", action="append", required=True, metavar="NAME=PATH")
    parser.add_argument("--engine-arg", action="append", default=[], metavar="NAME=ARG",
                        help="Append a per-engine benchmark argument after --timing (repeatable); not used by core checks")
    parser.add_argument("--positions", nargs="+", type=int, default=[10000, 100000, 1000000])
    parser.add_argument("--repeats", type=int, default=3)
    default_name = "weighted_optimization_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + ".json"
    parser.add_argument("--output", type=Path, default=ROOT / "validation" / default_name)
    parser.add_argument("--timeout", type=float, default=600, help="Per extractor/core process timeout (seconds)")
    parser.add_argument("--reference-timeout", type=float, default=600, help="Per MPFR verifier timeout (seconds)")
    parser.add_argument("--mpfr-exe", type=Path, default=ROOT / "validation" / "mpfr_window_interval.exe")
    parser.add_argument("--source", type=Path, nargs="+", action="append", default=[],
                        help="Source files to hash (repeatable); defaults to standard main sources")
    parser.add_argument("--check-baseline-core", type=Path, metavar="PATH",
                        help="Exact --core baseline; checks positions 100, 997, 10000, 100000")
    parser.add_argument("--core-engine", action="append", default=[], metavar="NAME",
                        help="Engine to compare with core baseline (repeatable); default: names containing weighted")
    default_runtime = [r"C:\msys64\ucrt64\bin"] if os.name == "nt" and Path(r"C:\msys64\ucrt64\bin").is_dir() else []
    parser.add_argument("--runtime-dir", action="append", default=default_runtime, metavar="PATH",
                        help="Prepend a DLL/shared-library executable search directory to PATH")
    args = parser.parse_args()
    if args.repeats < 1 or not all(1 <= position <= 10000000 for position in args.positions):
        parser.error("repeats must be positive and positions must be in the MPFR verifier range 1..10000000")
    if len(set(args.positions)) != len(args.positions):
        parser.error("positions must be unique")
    if not all(math.isfinite(value) and value > 0 for value in (args.timeout, args.reference_timeout)):
        parser.error("timeouts must be finite and positive")
    try:
        benchmark(args)
    except (OSError, ValueError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
