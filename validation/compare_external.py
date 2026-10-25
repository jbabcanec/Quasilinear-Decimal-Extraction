"""Sequential, same-machine extraction comparisons; no reference in timed runs.

Build (including a hash-checked download of Bellard's source):
    python -B validation/compare_external.py --build
Check adapters:
    python -B validation/compare_external.py --check
Benchmark:
    python -B validation/compare_external.py --bench 10000 100000 --repeats 3
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess
import sys
import time
import urllib.request
import gmpy2
import psutil

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "validation" / "baselines"
DEFAULT_OUT = ROOT / "validation" / "external_comparison_2026-09-07.json"
sys.path.insert(0, str(ROOT))
from validation.compare_cpp_euler import mpfr_window

ENV = os.environ.copy()
if os.name == "nt":
    ENV["PATH"] = r"C:\msys64\ucrt64\bin" + os.pathsep + ENV["PATH"]
ENGINES = ("bellard", "gourdon", "v5", "weighted", "mpfr_prefix")

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def build():
    provenance = json.loads((BASE / "provenance.json").read_text(encoding="utf-8"))
    bellard = BASE / "bellard_pi1.c"
    if not bellard.exists():
        raw = urllib.request.urlopen(provenance["bellard"]["url"], timeout=30).read()
        if hashlib.sha256(raw).hexdigest() != provenance["bellard"]["sha256"]:
            raise RuntimeError("Bellard source hash differs; inspect it before running.")
        bellard.write_bytes(raw)
    assert digest(bellard) == provenance["bellard"]["sha256"]
    gcc = shutil.which("gcc", path=ENV["PATH"])
    gxx = shutil.which("g++", path=ENV["PATH"])
    if not gcc or not gxx:
        raise RuntimeError("gcc and g++ must be on PATH.")
    common = ["-O3", "-march=native"]
    commands = [
        [gcc, *common, "-std=c17", "-DHAS_LONG_LONG", str(bellard), "-lm",
         "-o", str(BASE / "bellard_pi1.exe")],
        [gxx, *common, "-std=c++17", str(BASE / "gourdon_window.cpp"),
         "-o", str(BASE / "gourdon_window.exe")],
        [gxx, *common, "-std=c++17", str(BASE / "mpfr_prefix.cpp"), "-lmpfr", "-lgmp",
         "-o", str(BASE / "mpfr_prefix.exe")],
        [gxx, *common, "-std=c++17", str(ROOT / "decimal_extract_v5.cpp"),
         "-lgmpxx", "-lgmp", "-o", str(ROOT / "extract_v5.exe")],
    ]
    for command in commands:
        subprocess.run(command, cwd=ROOT, env=ENV, check=True, capture_output=True)
    result = {"commands": commands,
              "gcc": subprocess.check_output([gcc, "--version"], env=ENV, text=True).splitlines()[0],
              "g++": subprocess.check_output([gxx, "--version"], env=ENV, text=True).splitlines()[0]}
    (BASE / "build.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result), flush=True)

def run(engine, position):
    if engine == "weighted":
        command = [sys.executable, "-B", str(ROOT / "validation/weighted_block_crt.py"),
                   "--pos", str(position)]
    elif engine == "v5":
        command = [str(ROOT / "extract_v5.exe"), str(position), "--timing"]
    else:
        exe = {"bellard": "bellard_pi1.exe", "gourdon": "gourdon_window.exe",
               "mpfr_prefix": "mpfr_prefix.exe"}[engine]
        command = [str(BASE / exe), str(position)]
    start = time.perf_counter()
    proc = subprocess.Popen(command, cwd=ROOT, env=ENV, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    monitor = psutil.Process(proc.pid)
    peak = 0
    while proc.poll() is None:
        try:
            mem = monitor.memory_info()
            peak = max(peak, getattr(mem, "peak_wset", mem.rss))
        except psutil.NoSuchProcess:
            break
        time.sleep(0.001)
    stdout, stderr = proc.communicate()
    wall = time.perf_counter() - start
    if proc.returncode:
        raise RuntimeError(f"{command}: {stderr}")
    result = {"engine": engine, "N": position, "wall_seconds": wall,
              "peak_working_set_bytes": peak, "command": command}
    if engine == "weighted":
        data = json.loads(stdout)
        digits = data["digits"]
        result["engine_stats"] = data["stats"]
    elif engine == "bellard":
        match = re.fullmatch(r"Decimal digits of pi at position (\d+): (\d{9})\s*", stdout)
        assert match and int(match[1]) == position, stdout
        digits = match[2]
    else:
        digits = stdout.strip()
        if engine == "v5":
            result["engine_stats"] = json.loads(stderr)
    width = 12 if engine in ("weighted", "v5") else 9
    assert re.fullmatch(rf"\d{{{width}}}", digits), (engine, stdout)
    result.update(digits=digits, computed_width=width)
    return result

def check():
    checks = []
    for position in (100, 137, 1000, 2718, 9999):
        want = mpfr_window(position)
        for engine in ("bellard", "gourdon", "mpfr_prefix"):
            record = run(engine, position)
            assert record["digits"] == want[:record["computed_width"]], (record, want)
            checks.append({"N": position, "engine": engine, "digits": record["digits"]})
    report = {"verified": True, "reference": "MPFR, explicit floor, two guard precisions",
              "checks": checks}
    (BASE / "adapter_checks.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report), flush=True)

def benchmark(positions, repeats, output):
    sources = ["decimal_extract_v5.cpp", "validation/weighted_block_crt.py",
               "validation/harmonic_prime_crt.py", "decimal_extract_v5.py",
               "decimal_bbp_extract.py", "validation/compare_cpp_euler.py",
               "validation/compare_external.py", "validation/baselines/bellard_pi1.c",
               "validation/baselines/gourdon_window.cpp", "validation/baselines/mpfr_prefix.cpp"]
    report = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "method": "Fresh processes, sequential runs, rotating engine order. Wall time includes startup and output; peak process working set sampled every 1 ms. Reference calculations excluded. No warm-up runs.",
        "comparison_width": 9,
        "precision": {"bellard": "Original double accumulator, n+20 truncation parameter, nine digits",
                      "gourdon": "64-bit fractional accumulator, default n0=12, nine digits",
                      "v5": "G=192, twelve digits (all retained and verified)",
                      "weighted": "G=192, twelve digits (all retained and verified)",
                      "mpfr_prefix": "Prefix through N+8, 192 additional bits, nine digits"},
        "hardware": {"processor": "AMD Ryzen 7 5700", "physical_cores": psutil.cpu_count(logical=False),
                     "logical_processors": psutil.cpu_count(), "physical_memory_bytes": psutil.virtual_memory().total},
        "os": platform.platform(), "python": sys.version,
        "gmpy2": gmpy2.version(), "gmp": gmpy2.mp_version(), "mpfr": gmpy2.mpfr_version(),
        "compiler": json.loads((BASE / "build.json").read_text(encoding="utf-8")),
        "provenance": json.loads((BASE / "provenance.json").read_text(encoding="utf-8")),
        "source_sha256": {name: digest(ROOT / name) for name in sources},
        "runs": [], "references": {}, "summary": [], "verified": False,
    }
    def save():
        output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for position in positions:
        if not 11 <= position <= 1000000:
            raise ValueError("Benchmark position outside audited range 11..1000000.")
        for repeat in range(repeats):
            offset = repeat % len(ENGINES)
            order = ENGINES[offset:] + ENGINES[:offset]
            for engine in order:
                record = run(engine, position)
                record["repeat"] = repeat
                report["runs"].append(record)
                save()
                print(json.dumps(record), flush=True)
        want = mpfr_window(position)
        report["references"][str(position)] = want
        for engine in ENGINES:
            samples = [r for r in report["runs"] if r["engine"] == engine and r["N"] == position]
            assert all(r["digits"] == want[:r["computed_width"]] for r in samples), (samples, want)
            times = [r["wall_seconds"] for r in samples]
            report["summary"].append({
                "engine": engine, "N": position, "repetitions": len(samples),
                "computed_width": samples[0]["computed_width"],
                "wall_median": statistics.median(times), "wall_min": min(times), "wall_max": max(times),
                "peak_bytes_max": max(r["peak_working_set_bytes"] for r in samples),
            })
        save()
        print(json.dumps({"N": position, "verification": "PASS", "summary": report["summary"][-5:]}), flush=True)
    report["verified"] = True
    report["completed_utc"] = datetime.now(timezone.utc).isoformat()
    save()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--bench", type=int, nargs="+")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    if args.build: build()
    if args.check: check()
    if args.bench: benchmark(args.bench, args.repeats, args.output)
