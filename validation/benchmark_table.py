#!/usr/bin/env python3
"""benchmark_table.py -- one consistent measurement of the paper's timing table.
(2026-07-31)

Re-measures every cell in one session on one machine, rather than mixing runs
taken on different days under different load. Reports wall time, peak working
set, and the digit window, so the two engines also cross-check each other.

Peak memory matters here because it is what decides how far each engine can
go: the Python engine materializes the band as a list of tuples, so its
ceiling arrives well before the compiled engine's.

Usage:  python validation/benchmark_table.py [maxN_python] [maxN_cpp]
"""

import subprocess
import sys
import time

import psutil

ROOT = r"C:\Users\josep\dropbox\babcanec works\mathematics\bbp algos"
# the child times only the computation, so interpreter and gmpy2 import
# startup (~0.15 s) does not inflate the small-N cells
PY_SNIPPET = (
    "import sys, time; sys.path.insert(0, r'%s');"
    "import decimal_bbp_extract as D;"
    "t=time.perf_counter(); d=D.extract_v3(%d, window=12)[0];"
    "print('%%.3f %%s' %% (time.perf_counter()-t, d))"
)


def run(cmd, cwd=ROOT):
    """Run cmd, returning (wall_seconds, peak_working_set_bytes, stdout)."""
    t0 = time.perf_counter()
    p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True)
    peak = 0
    try:
        ps = psutil.Process(p.pid)
        while p.poll() is None:
            try:
                mi = ps.memory_info()
                peak = max(peak, getattr(mi, "peak_wset", mi.rss))
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break
            time.sleep(0.02)
    except psutil.NoSuchProcess:
        pass
    out, _ = p.communicate()
    wall = time.perf_counter() - t0
    toks = (out or "").strip().split()
    # a python child prints "<self-timed seconds> <digits>"; the exe prints digits
    if len(toks) == 2:
        return float(toks[0]), peak, [toks[1]]
    return wall, peak, toks[-1:]


def human(nbytes):
    for unit, div in (("GB", 1 << 30), ("MB", 1 << 20), ("KB", 1 << 10)):
        if nbytes >= div:
            return "%.1f %s" % (nbytes / div, unit)
    return "%d B" % nbytes


def main():
    max_py = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 7
    max_cpp = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10 ** 8

    print("%-10s | %-12s %-10s %-14s | %-12s %-10s %-14s | agree"
          % ("N", "Python", "peak", "digits",
             "C++", "peak", "digits"))
    print("-" * 104)
    for e in range(4, 9):
        N = 10 ** e
        row = ["10^%d" % e]

        if N <= max_py:
            t, m, d = run([sys.executable, "-c", PY_SNIPPET % (ROOT, N)])
            dpy = d[0] if d else "(none)"
            row += ["%.2f s" % t, human(m), dpy]
        else:
            dpy = None
            row += ["not run", "-", "-"]

        if N <= max_cpp:
            t, m, d = run([r".\extract_v3.exe", str(N)])
            dcc = d[0] if d else "(none)"
            row += ["%.2f s" % t, human(m), dcc]
        else:
            dcc = None
            row += ["not run", "-", "-"]

        agree = "yes" if (dpy and dcc and dpy == dcc) else (
            "-" if not (dpy and dcc) else "MISMATCH")
        print("%-10s | %-12s %-10s %-14s | %-12s %-10s %-14s | %s"
              % tuple(row + [agree]))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
