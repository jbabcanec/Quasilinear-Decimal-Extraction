"""Fresh native weighted-block performance samples, with separate MPFR checks."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time
import psutil

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from validation.compare_cpp_euler import mpfr_window
env=os.environ.copy()
if os.name=="nt":
    env["PATH"]=r"C:\msys64\ucrt64\bin"+os.pathsep+env["PATH"]

def run(N):
    start=time.perf_counter()
    p=subprocess.Popen([str(ROOT/"extract_weighted.exe"),str(N),"--timing"],env=env,
                       stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    process=psutil.Process(p.pid); peak=0
    while p.poll() is None:
        try:
            mem=process.memory_info(); peak=max(peak,getattr(mem,"peak_wset",mem.rss))
        except psutil.NoSuchProcess: break
        time.sleep(.001)
    stdout,stderr=p.communicate()
    wall=time.perf_counter()-start
    if p.returncode: raise RuntimeError(stderr)
    stats=json.loads(stderr)
    return {"N":N,"digits":stdout.strip(),"stats":stats,
            "wall_seconds":wall,"peak_working_set_bytes":peak}

def bench(positions,repeats,output):
    report={"started_utc":datetime.now(timezone.utc).isoformat(),
            "engine":"Native C++ weighted-block CRT, G=192, width=12",
            "method":"Sequential fresh processes; no discarded warm-ups; 1 ms peak working-set sampling. Engine times exclude process startup. MPFR reference outside timed runs.",
            "processor":"AMD Ryzen 7 5700","os":platform.platform(),
            "compiler":"GCC 16.1.0, -O3 -march=native -std=c++17, GMP 6.3.0",
            "source_sha256":{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                             for name in ("decimal_extract_weighted.cpp","validation/benchmark_weighted_cpp.py")},
            "runs":[],"summary":[],"references":{},"verified":False}
    def save(): output.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    for N in positions:
        samples=[]
        for rep in range(repeats):
            record=run(N); record["repeat"]=rep
            samples.append(record); report["runs"].append(record); save()
            print(json.dumps(record),flush=True)
        want=mpfr_window(N); report["references"][str(N)]=want
        assert all(s["digits"]==want for s in samples),(N,samples,want)
        report["summary"].append({
            "N":N,"repetitions":repeats,"digits":want,
            "compute_median":statistics.median(s["stats"]["total"] for s in samples),
            "compute_min":min(s["stats"]["total"] for s in samples),
            "compute_max":max(s["stats"]["total"] for s in samples),
            "wall_median":statistics.median(s["wall_seconds"] for s in samples),
            "wall_min":min(s["wall_seconds"] for s in samples),
            "wall_max":max(s["wall_seconds"] for s in samples),
            "peak_bytes_max":max(s["peak_working_set_bytes"] for s in samples),
            "phase_medians":{key:statistics.median(s["stats"][key] for s in samples)
                             for key in ("word","core","factor","coefficients","real","rational","assembly")},
        })
        save()
    report["verified"]=True; report["completed_utc"]=datetime.now(timezone.utc).isoformat(); save()
    print(json.dumps({"verified":True,"summary":report["summary"]}),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--positions",nargs="+",type=int,default=[10000,100000,1000000])
    ap.add_argument("--repeats",type=int,default=3)
    ap.add_argument("--output",type=Path,default=ROOT/"validation/weighted_cpp_performance_2026-09-07.json")
    args=ap.parse_args()
    bench(args.positions,args.repeats,args.output)
