import os
import subprocess
import ctypes
import time
import matplotlib.pyplot as plt
from pybench import Benchmark,parser

def compile(bench,compiler,opt_level,basename):

    cmd_opt = 'polycc --tile --parallel %s.cpp'%basename

    name_pluto = basename+".cpp.pluto.c"
    cmd_comp = switchcompiler(compiler,opt_level,name_pluto)

    try:
        subprocess.check_call(cmd_opt,shell=True)
        subprocess.check_call(cmd_comp,shell=True)
        return name_pluto
    except subprocess.CalledProcessError as e:
        print "Compilation error: ", e
        raise Exception("Failed to compile ")

# produce command to run
def switchcompiler(compiler,opt_level,name):
# -fPIC -shared
    basename = name
    cc = 'g++'
    cc += ' -fopenmp'
    cc += ' -O%d'% opt_level
    cc += " -o %s_g++" % basename
    cc += " %s" % basename
    cc += " test.h"

    clang = 'clang++'
    clang += ' -fopenmp'
    clang += ' -o %s_clang' % basename
    clang += ' %s'% basename
    clang += ' -O%d' % opt_level

    return {
        'g++':cc,
        'clang':clang,
    }[compiler]

def runlib(bench,basename,compiler):
    try: 
        with bench.timed_region("running time with %s" % compiler):
            perf = "" # "perf stat -e task-clock,cycles,instructions,cache-references,cache-misses "
            log = "" #"-o " + "log_default " + "--append "
            cmd = perf + log + os.getcwd() + "/" + basename + '_'+ compiler
            subprocess.check_call(cmd,shell=True)

    except OSError as e:
        print "running error:", e


class myBench(Benchmark):

    warmups = 0
    repeats = 1
    method = 'benchmarking'
    benchmark = 'myBench'

    compilers = ['g++']
    sizes = ['4 4 1000','4 1000','4 32','4 4 32']
    params =  [('opt_level', range(3,4)),('compiler',compilers),('tile',sizes)]
    #when adding new options , add in switchcompiler, pollyoptions, 
 
    basename = 'eigenwave3d'
    
    def benchmarking(self, opt_level=0,times=3 , compiler='',tile = '32 32 32'):
        # #set(4)
        self.series['times'] = times

        f = open('tile.sizes','w')
        f.write(tile)
        f.close()

        #if not self.compiled[compiler]:
        pluto_name = compile(self,compiler,opt_level,self.basename)

        runlib(self,pluto_name,compiler)
        

# os.environ["OMP_NUM_THREADS"] = "4"
# os.environ["KMP_AFFINITY"]="granularity=thread,scatter"
myBench().main()