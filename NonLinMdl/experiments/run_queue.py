"""Run qLPV-MPC experiments from a queue file in parallel MATLAB batch processes.

Queue file: one experiment per line, tab separated:
    name<TAB>plant<TAB>wind<TAB>MATLAB struct expression for the parameters
Lines starting with # are ignored. Experiments with an existing result file are skipped.

Usage: python run_queue.py queue.txt [max_parallel]
"""
import os
import subprocess
import sys
import time

MATLAB = r'C:\Program Files\MATLAB\R2023a\bin\matlab.exe'
EXP_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    queue_file = sys.argv[1]
    max_parallel = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    jobs = []
    for line in open(queue_file, encoding='utf-8'):
        line = line.rstrip('\n')
        if not line.strip() or line.startswith('#'):
            continue
        name, plant, wind, params = line.split('\t')
        if os.path.isfile(os.path.join(EXP_DIR, 'results', name + '.mat')):
            print('skip (exists):', name, flush=True)
            continue
        jobs.append((name, plant, wind, params))
    os.makedirs(os.path.join(EXP_DIR, 'logs'), exist_ok=True)
    running = []
    while jobs or running:
        while jobs and len(running) < max_parallel:
            name, plant, wind, params = jobs.pop(0)
            cmd = "addpath('%s'); run_qlmpc_experiment('%s','%s','%s',%s)" % (EXP_DIR, name, plant, wind, params)
            log = os.path.join(EXP_DIR, 'logs', name + '.txt')
            p = subprocess.Popen([MATLAB, '-batch', cmd, '-logfile', log],
                                 stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            running.append((name, p, time.time()))
            print('start', name, flush=True)
        time.sleep(10)
        for item in running[:]:
            name, p, t0 = item
            if p.poll() is not None:
                ok = os.path.isfile(os.path.join(EXP_DIR, 'results', name + '.mat'))
                print('done %s (%.1f min, %s)' % (name, (time.time() - t0) / 60, 'result' if ok else 'NO RESULT'), flush=True)
                running.remove(item)
    print('queue finished', flush=True)


if __name__ == '__main__':
    main()
