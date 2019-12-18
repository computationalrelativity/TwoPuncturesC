# build a sequence of single spining punctures

import subprocess
import numpy as np
import os

steps=10

TEMPLATE = """\
# singlepunc_seq.py
npoints_A=60
npoints_B=60
npoints_phi=40
par_b=1.e-10
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
    parfile = "singlepunc_S{0:.3f}.par".format(s)
    os.system("mkdir -p {}".format(dir))
    open("{}/".format(dir) + parfile, "w").write(TEMPLATE.format(spin=spin))
    cmd = "../TwoPuncturesRun.x " + parfile
    cwd=dir+"/"
    return subprocess.call(cmd, cwd=cwd, shell=True)

if __name__ == '__main__':

    spins = np.linspace(0.,0.9, steps)
    for s in spins: run(s)
