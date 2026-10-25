"""Measure the fixed-guard weighted-block prototype in fresh processes.

This is a performance experiment, not a proof of the complexity exponent.
Run from any directory: python -B validation/benchmark_weighted_blocks.py
"""
from pathlib import Path
import hashlib
import json
import platform
import statistics
import subprocess
import sys
import time
import gmpy2
import psutil

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation" / "weighted_block_performance_2026-09-06.json"
sys.path.insert(0, str(ROOT))
from validation.compare_cpp_euler import mpfr_window

def run(position):
    start = time.perf_counter()
    proc = subprocess.Popen(
        [sys.executable, "-B", str(ROOT / "validation/weighted_block_crt.py"),
         "--pos", str(position)], cwd=ROOT,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    monitor = psutil.Process(proc.pid)
    peak = 0
    while proc.poll() is None:
        try:
            memory = monitor.memory_info()
            peak = max(peak, getattr(memory, "peak_wset", memory.rss))
        except psutil.NoSuchProcess:
            break
        time.sleep(0.01)
    stdout, stderr = proc.communicate()
    wall = time.perf_counter() - start
    if proc.returncode:
        raise RuntimeError(stderr)
    result = json.loads(stdout)
    result.update(N=position, wall_seconds=wall, peak_working_set_bytes=peak)
    return result

report = {
    "date": "2026-09-06",
    "engine": "Python weighted-block CRT, fixed guard G=192, width=12",
    "method": "Three sequential fresh processes per position; computation time from the engine; wall time includes startup; peak working set sampled every 10 ms.",
    "hardware": {
        "processor": "AMD Ryzen 7 5700",
        "physical_cores": 8, "logical_processors": 16,
        "physical_memory_bytes": psutil.virtual_memory().total,
    },
    "os": platform.platform(), "python": sys.version,
    "gmpy2": gmpy2.version(), "gmp": gmpy2.mp_version(),
    "mpfr": gmpy2.mpfr_version(),
    "source_sha256": {
        name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        for name in ["validation/weighted_block_crt.py",
                     "validation/harmonic_prime_crt.py",
                     "decimal_extract_v5.py", "decimal_bbp_extract.py"]
    },
    "runs": [], "summary": [], "references": {},
}
for position in [10_000, 100_000, 1_000_000]:
    samples = []
    for repeat in range(3):
        result = run(position)
        result["repeat"] = repeat
        samples.append(result)
        report["runs"].append(result)
        OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result), flush=True)
    want = mpfr_window(position)
    assert all(sample["digits"] == want for sample in samples)
    report["references"][str(position)] = {
        "digits": want, "method": "MPFR const_pi with explicit floor at two guard precisions; reference time excluded"
    }
    compute = [s["stats"]["elapsed_seconds"] for s in samples]
    report["summary"].append({
        "N": position, "repetitions": 3,
        "computation_seconds_median": statistics.median(compute),
        "computation_seconds_min": min(compute),
        "computation_seconds_max": max(compute),
        "wall_seconds_median": statistics.median(s["wall_seconds"] for s in samples),
        "peak_working_set_bytes_max": max(s["peak_working_set_bytes"] for s in samples),
        "blocks": samples[0]["stats"]["blocks"],
        "forests": samples[0]["stats"]["forests"],
        "digits": want,
    })
    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
report["verified"] = True
OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"verified": True, "summary": report["summary"]}), flush=True)
