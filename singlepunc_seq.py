# Build a sequence of single spinning punctures

import subprocess
import numpy as np
import os
from multiprocessing import Pool


steps=10

TEMPLATE = """\
# singlepunc_seq.py
npoints_A=60
npoints_B=60
npoints_phi=40
par_b=1.e-5
par_m_plus=1.
par_P_plus1=0.
par_P_plus2=0.
par_P_plus3=0.
par_S_plus1=0.
par_S_plus2=0.
par_S_plus3={spin}
par_m_minus=1.e-12
par_P_minus1=0.
par_P_minus2=0.
par_P_minus3=0.
par_S_minus1=0.
par_S_minus2=0.
par_S_minus3=0.
do_solution_file_output=1
do_bam_file_output=1
verbose=1
#
"""

def run(spin):
    """
    Run TwoPuncture.c
    """
    dir="SinglePuncture"
    parfile = "singlepunc_S{0:.3f}.par".format(spin)
    os.system("mkdir -p {}".format(dir))
    open("{}/".format(dir) + parfile, "w").write(TEMPLATE.format(spin=spin))
    cmd = "../TwoPuncturesRun.x " + parfile
    cwd=dir+"/"
    return subprocess.call(cmd, cwd=cwd, shell=True)

if __name__ == '__main__':

    spins = np.unique(np.concatenate((np.linspace(0.3,0.4,5,endpoint=False)[1::], np.linspace(0.4,0.5,5,endpoint=False)[1::],np.linspace(0.5,0.6,10,endpoint=False)[1::],np.linspace(0.6,0.7,10,endpoint=False)[1::])))
    spins = np.concatenate((np.array([0.59]), np.linspace(0.6,0.7,10,endpoint=False)[1::]))
    spins = np.array([0.02,0.04,0.06,0.08,0.12,0.14,0.16,0.18,0.22,0.24,0.26,0.28])
    pool = Pool(processes=4)
    pool.map(run, spins)
    #for s in spins: run(s)
