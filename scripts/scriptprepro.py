import os
import sys
import re
from multiprocessing import Pool, Process, JoinableQueue, Pipe
import time
import numpy as np
import subprocess
import shutil
import time
from pathlib import Path
##Constant
# Number of process to parallelization
NUM_PROCESS = 2
# LIST of nuclides to make data
nuclidlist = []
# Type of cross-section lib, allowed only "tendl" and "endfb". By default "tendl"
_MOD = "tendl"
# PREPRO KIT
PREPROEXEC = ["endf2c", "linear", "recent", "sigma1", "activate",
              "legend", "fixup", "dictin", "groupie", "MT.DAT"]
# INPUTS , dict key - filename : value - file pattern

INPUTS = {"LINEAR.INP" : """          0          0 1.000Ee-30          0
ENDFB.OUT
LINEAR.OUT
     0 0  0  999999999
                        (BLANK CARD TERMINATES MAT REQUEST RANGES)
 0.0000E-00 1.0000E-03
                        (BLANK CARD TERMINATES FILE 3 ERROR LAW)
 0.0000E-00 1.0000E-05
 1.0000E+00 1.0000E-05
 2.0000E+00 1.0000E-04
 2.0000E+07 1.0000E-04
                        (BLANK CARD TERMINATES FILE 3 ERROR LAW)
====(this line and below is not read)===================================
za092235
LINEAR.OUT
"""
, "RECENT.INP" :   """           0 1.0000E-30          1          1          1          1
LINEAR.OUT
RECENT.OUT
          1       9999
                        (BLANK CARD TERMINATES MAT REQUEST RANGES)
 0.0000E-00 1.0000E-02
                        (BLANK CARD TERMINATES FILE 3 ERROR LAW)
 0.0000E-00 5.0000E-04
 1.0000E+00 5.0000E-04
 2.0000E+00 5.0000E-03
 2.0000E+07 5.0000E-03
                        (BLANK CARD TERMINATES FILE 3 ERROR LAW)
"""
, "SIGMA1.INP" : """          0          2 3.0000E+02 1.0000E-30          1
RECENT.OUT
SIGMA1.OUT
          1       9999
                       (BLANK CARD TERMINATES MAT RANGE REQUESTS)
 0.0000E+ 0 1.0000E-02
                       (BLANK CARD TERMINATES FILE 3 ERROR LAW)
 0.0000E+ 0 5.0000E-04
 1.0000E+ 0 5.0000E-04
 2.0000E+ 0 5.0000E-03
 2.0000E+ 7 5.0000E-03
                       (BLANK CARD TERMINATES FILE 3 ERROR LAW)
 0.0000E+ 0 1.0000E-04
                       (BLANK CARD TERMINATES FILE 3 ERROR LAW)



====(this line and below is not read)==============================
          0          2 1.2000E+06 1.0000E-10          1
tape.135.RECENT
tape.135.SIGMA1
"""
, "ACTIVATE.INP" : """SIGMA1.OUT
ACTIVATE.OUT """
, "LEGEND.INP" : """ 1.0000E-02      20000          2          1          2          0
ACTIVATE.OUT
LEGEND.OUT
     0 0  0  999999999"""
, "FIXUP.INP" : """10002111110001
LEGEND.OUT
FIXUP.OUT

---(this line and below is not read)-----------------------------------
Standard Input is Below
12202111110001"""
, "DICTIN.INP" :  """FIXUP.OUT
PREPRO19.OUT"""
, "GROUPIE.INP" : """          0        -9         2           -2 1.0000E-03           0
PREPRO19.OUT
GROUPIE.OUT
          1          1        1           1          1
Groupie Test Run
     1 1  1  999999999"""
, "verify" :   """#!/bin/sh
./endf2c
./linear
./recent
./sigma1
./activate
./legend
./fixup
./dictin
./groupie
"""}

class freenumber:
    """Generate a current free number for a new process

    Parametres:
    -----------
    numbers : int
       number of the process

    Mathods:
    --------
    free(int index)
        get free the index process
    busy(int index)
        get busy the index process
    get_free()
        return : int free number"""
    def __init__(self, numbers):
        self._status = [True for i in range(numbers)]
        self._index = 0
    def _findfree(self):
        for i, v in enumerate(self._status):
            if (v):
                return i
        return -1
    def free(self, index):
        self._status[index] = True
        self._index = self._findfree()
        pass
    def busy(self, index):
        self._status[index] = False
        self._index = self._findfree()
        pass
    def get_free(self):
        return self._index

def preproc(preprodir, execdir, numprocess=1):
    """Make a directories to store multigroup data

    Paramteres:
    -----------
    preprodir : str
       path to prepro executive
    execdir : str
       path to execution directory
    numprocess : int
       number of the process"""
    if not os.path.isdir(os.path.join(execdir, "proceed")):
        os.mkdir(os.path.join(execdir, "proceed"))
    for n in range(numprocess):
        if not os.path.isdir(os.path.join(execdir, "proceed", str(n + 1))):
           os.mkdir(os.path.join(execdir, "proceed", str(n + 1)))
        for pe in PREPROEXEC:
            shutil.copyfile(os.path.join(preprodir, pe),
                            os.path.join(execdir, "proceed", str(n + 1), pe))
        for i, v in INPUTS.items():
            with open(os.path.join(preprodir, i), 'w') as f:
                f.write(v)
            shutil.copyfile(os.path.join(preprodir, i),
                            os.path.join(execdir, "proceed", str(n + 1), i))

def wrapp_func(pdir, nuclidname):
    """Start the execution

    Parametres:
    pdir : str
        path to nuclidname library directory
    nuclidname : str
        name of nuclide to proceed"""
    os.chdir(pdir)
    DEVNULL = os.open(os.devnull, os.O_WRONLY)
    ps = subprocess.call([os.path.join(pdir, "verify")], stdout = DEVNULL)
    time.sleep(0.1)
    try:
        r = next(Path.cwd().glob('*.UNSHIELD.LST'))
        os.rename(r.name, "{}.tab".format(nuclidname))
    except(StopIteration, UnboundLocalError):
        print("Error file {} not proceeding".format(nuclidname))

def parse_endfblib(libdir):
    """Parse ENDF/B library
    Parametres:
    -----------
    libdir : str
       directory with ENDFB file structure"""
    filepaths = []
    nuclidnames = []
    endf_dir = Path(libdir)
    neutron_files = tuple((endf_dir / "neutrons").glob("*endf"))
    for n in neutron_files:
        filepaths.append(n.absolute())
        nuclidnames.append(n.name.split('_')[1] +
         re.split("^0*", n.name.split('_')[2][:-5])[-1])
    return nuclidnames, filepaths

def parse_tendlib(libdir):
    """Parse TENDL library
    Parametres:
    -----------
    libdir : str
       directory with TENDL file structure"""
    filepaths = []
    nuclidnames = []
    for dirs, subdir, files in os.walk(libdir):
        for f in files:
            name_isotope = f.split('.')[0][2:]
            numpart_ = int(re.search("\d+", name_isotope).group())
            name_isotope = re.split("\d+", name_isotope)[0] + str(numpart_)
            nuclidnames.append(name_isotope)
            filepaths.append(os.path.join(dirs, f))
    print(nuclidnames)
    print(filepaths)
    return nuclidnames, filepaths

def exec(nuclidnames, filepaths, execdir, nuclidlist):
    """The main execution loop

    Parametres:
    -----------
    nuclidnames : iterable of str
       name of nuclides to cross-section processing
    filepaths : iterable of str/path
       name of file path to nuclude pointed in nuclidnames"""
    workers = [None]* NUM_PROCESS
    curfree = 0
    gen = freenumber(NUM_PROCESS)
    for (n, f) in zip(nuclidnames, filepaths):
        if (n in nuclidlist):
            if (gen.get_free() == -1):
                cycle = True
                while cycle:
                    for w in workers:
                        if not w.is_alive():
                            cycle = False
                            curfree = workers.index(w)
                            workers[curfree] = None
                            gen.free(curfree)
                    time.sleep(0.5)
            if (gen.get_free() != -1):
                curfree = gen.get_free()
                gen.busy(curfree)
                shutil.copyfile(f, os.path.join(execdir,
                 str(curfree + 1), "ENDFB.IN"))
                worker=Process(target = wrapp_func, args=(os.path.join(execdir,
                str(curfree + 1)) ,n))
                worker.start()
                print('isotope with name' + n + ' added')
                workers[curfree] = worker
    while (len(workers) > 0):
        if (workers == [None]* NUM_PROCESS):
            print("ALL NONE")
            return
        for w in workers:
            if (w is not None):
                if not w.is_alive():
                    workers.pop(workers.index(w))
                else:
                    time.sleep(5.)

if __name__ == '__main__':
    preprodir = sys.argv[-3]
    libdir = sys.argv[-2]
    targetdir = sys.argv[-1]
    execdir = os.path.join(targetdir, "proceed")
    if (_MOD == "endfb"):
        nucls, paths = parse_endfblib(libdir)
    else:
        nucls, paths = parse_tendlib(libdir)
    if (len(nuclidlist) == 0):
        nuclidlist = nucls
    preproc(preprodir,
            targetdir, NUM_PROCESS)
    exec(nucls, paths, execdir, nuclidlist)

