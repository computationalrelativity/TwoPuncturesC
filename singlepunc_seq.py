# build a sequence of single spining punctures

import subprocess
import numpy as np

steps=10

TEMPLATE = """\
# singlepunc_seq.py
bar_b=1
par_m_plus=1.0
par_S_plus3={spin}
#
"""

def run(spin):
    """
    Run TwoPuncture.c
    """
    parfile = "singlepunc_S{0:.3f}.par".format(s)
    open(parfile, "w").write(TEMPLATE.format(spin=spin)
    cmd = "TwoPuncturesRun.x " + parfile
    return subprocess.call(cmd, shell=True)

if __name__ == '__main__':

    spins = np.linespace(0.,0.9, steps)
    for s in spins: run(s)
