"""Independent correctness checks for the native weighted-block port."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import time
from gmpy2 import mpz

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from decimal_bbp_extract import cutoff
from decimal_extract_v5 import build_v_euler
from validation.compare_cpp_euler import mpfr_window

env=os.environ.copy()
if os.name=="nt":
    env["PATH"]=r"C:\msys64\ucrt64\bin"+os.pathsep+env["PATH"]
exe=ROOT/"extract_weighted.exe"
report={"source_sha256":hashlib.sha256((ROOT/"decimal_extract_weighted.cpp").read_bytes()).hexdigest(),
        "native_selftest":None,"core_checks":[],"digit_checks":[],"verified":False}
def call(*args):
    return subprocess.check_output([str(exe),*map(str,args)],env=env,text=True)
report["native_selftest"]=json.loads(call("--selftest"))
assert report["native_selftest"]["selftest"]=="PASS"
for N in (100,997,10000,100000,1000000):
    native=json.loads(call("--core",N))
    X=mpz(5)**(N-1)
    ref,m,_,_=build_v_euler(N,cutoff(N,192,"pi"),X.bit_length()+200)
    assert native["bits"]==m and mpz(native["V_hex"],16)==ref,N
    report["core_checks"].append({"N":N,"bits":m,"exact_match":True,
                                  "reference":"Separate 2-adic Euler construction, Python/GMP"})
    print(json.dumps(report["core_checks"][-1]),flush=True)
for N in [*range(1,41),100,997,10000,100000,1000000]:
    got=call(N).strip()
    want=mpfr_window(N)
    assert got==want,(N,got,want)
    report["digit_checks"].append({"N":N,"digits":got,"MPFR_match":True})
report["verified"]=True
(ROOT/"validation/weighted_cpp_checks_2026-09-07.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"verified":True,"native_selftest":report["native_selftest"],
                  "exact_large_cores":len(report["core_checks"]),"decimal_windows":len(report["digit_checks"])}),flush=True)
