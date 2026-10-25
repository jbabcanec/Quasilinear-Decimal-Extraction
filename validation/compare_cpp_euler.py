"""Same-machine, sequential C++ v3/v5 comparison with independent MPFR checks."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time

import gmpy2 as g
import psutil

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from decimal_bbp_extract import pi_reference


def run(engine,N):
    env=os.environ.copy()
    env['PATH']=r'C:\msys64\ucrt64\bin'+os.pathsep+env['PATH']
    env['V3_TIMING']='1'
    exe='extract_v3_control.exe' if engine=='v3' else 'extract_v5.exe'
    cmd=[str(ROOT/exe),str(N)]
    if engine=='v5': cmd.append('--timing')
    start=time.perf_counter()
    p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,env=env)
    proc=psutil.Process(p.pid)
    peak=0
    while p.poll() is None:
        try:
            mem=proc.memory_info()
            peak=max(peak,getattr(mem,'peak_wset',mem.rss))
        except psutil.NoSuchProcess:
            break
        time.sleep(0.01)
    out,err=p.communicate()
    wall=time.perf_counter()-start
    if p.returncode:
        raise RuntimeError(f'{cmd}: {err}')
    digits=out.strip()
    assert re.fullmatch(r'\d{12}',digits),(cmd,out,err)
    if engine=='v5':
        stats=json.loads(err.strip())
    else:
        fields=dict(re.findall(r'\[\s*(\w+)\s+([0-9.]+)s\]',err))
        stats=dict(pre=float(fields['X']),word=float(fields['word']),vphase=float(fields['V']))
        stats['total']=sum(stats.values())
    stats.update(engine=engine,N=N,digits=digits,wall_seconds=wall,peak_working_set_bytes=peak)
    return stats


def mpfr_window(N,window=12):
    # For verification ONLY: this computes a prefix independently with MPFR.
    # The extractor itself never calls MPFR or this reference function.
    shift=N+window-1
    scale=g.mpz(10)**shift
    answers=[]
    for extra in (128,192):
        with g.context(precision=scale.bit_length()+extra):
            # mpz(mpfr) follows the context's rounding mode; use floor
            # explicitly because digit extraction truncates, never rounds.
            answers.append(str(g.floor(g.const_pi()*scale)%(10**window)).zfill(window))
    assert answers[0]==answers[1]
    return answers[0]


def selftest():
    assert mpfr_window(4,1)=='5'  # truncation of 3.141592..., not rounding to 6
    positions=list(range(1,61))+[100,997,2500,10000]
    ref=pi_reference(max(positions)+20)
    for N in positions:
        new=run('v5',N)['digits']
        assert new==ref[N-1:N+11],(N,new,ref[N-1:N+11])
    print(json.dumps(dict(cpp_digit_selftest='PASS',positions=len(positions),
                         reference='integer Chudnovsky')),flush=True)


def benchmark(sizes,repeats,outfile):
    report=dict(method='Sequential runs, alternating engine order; process peak working set sampled every 10 ms.',
                python=sys.version,gmp=g.mp_version(),mpfr=g.mpfr_version(),
                compiler_flags='-O3 -march=native -std=c++17',
                source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                               for name in ('decimal_extract_v3.cpp','decimal_extract_v5.cpp')},
                runs=[],references={})
    for N in sizes:
        for rep in range(repeats):
            order=['v3','v5'] if rep%2==0 else ['v5','v3']
            for name in order:
                record=run(name,N)
                record['repeat']=rep
                report['runs'].append(record)
                print(json.dumps(record),flush=True)
                if outfile:
                    Path(outfile).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        want=mpfr_window(N)
        report['references'][str(N)]=dict(digits=want,method='MPFR const_pi at two guard precisions')
        for record in report['runs']:
            if record['N']==N:
                assert record['digits']==want,(N,record,want)
        print(json.dumps(dict(N=N,independent_MPFR='PASS',digits=want)),flush=True)
        if outfile:
            Path(outfile).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')


def verify_report(path):
    report=json.loads(Path(path).read_text(encoding='utf-8'))
    for N in sorted({r['N'] for r in report['runs']}):
        want=mpfr_window(N)
        for record in report['runs']:
            if record['N']==N:
                assert record['digits']==want,(N,record,want)
        report['references'][str(N)]=dict(digits=want,
            method='MPFR const_pi with explicit floor, at two guard precisions')
        print(json.dumps(dict(N=N,independent_MPFR='PASS',digits=want)),flush=True)
    report['verified']=True
    Path(path).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--selftest',action='store_true')
    ap.add_argument('--bench',nargs='+',type=int)
    ap.add_argument('--repeats',type=int,default=2)
    ap.add_argument('--output')
    ap.add_argument('--verify-report')
    args=ap.parse_args()
    if args.selftest: selftest()
    if args.bench: benchmark(args.bench,args.repeats,args.output)
    if args.verify_report: verify_report(args.verify_report)
