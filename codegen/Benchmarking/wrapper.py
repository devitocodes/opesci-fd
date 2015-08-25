import os
import subprocess
import ctypes
import time
import matplotlib.pyplot as plt

from pybench import Benchmark,parser

def reset():

    os.environ["OMP_NUM_THREADS"] =""
    os.environ["KMP_AFFINITY"]=""


def set(aff):
    reset()
    if(aff>=1):
        os.environ["OMP_NUM_THREADS"] = str(aff)
        os.environ["KMP_AFFINITY"]="verbose,granularity=thread,scatter"
    


def compile(bench,compiler,opt_level,basename):

    
    cmd = switchcompiler(compiler,opt_level,basename)
    try:
        #with bench.timed_region("compile time with %s" % compiler):
        print cmd
        subprocess.check_call(cmd,shell=True)
    except subprocess.CalledProcessError as e:
        print "Compilation error: ", e
        raise Exception("Failed to compile ")

# produce command to run
def switchcompiler(compiler,opt_level,name):
# -fPIC -shared
    basename = name
    cc = 'g++'
    cc += ' -fopenmp -O%d'% opt_level
    cc += " -o %s_g++" % basename
    cc += " %s.cpp" % basename
    cc += " test.h"

    clang = 'clang++'
    clang += ' -fopenmp'
    clang += ' -o %s_clang' % basename
    clang += ' %s.cpp'% basename
    clang += ' -O%d' % opt_level

    polly = 'clang++'
    polly += ' -fopenmp'
    polly += ' -Xclang -load -Xclang LLVMPolly.so'
    polly += ' -O3 -mllvm'

    # pluto = polly

    polly += ' -polly' # -mllvm --polly-tile-sizes=8,8,8'

    pollywithcmd = polly

    polly += ' -o %s_polly' % basename
    polly += ' %s.cpp' % basename

    extracmd = pollyoption(compiler)    # extra options for polly

    pollywithcmd += ' -mllvm -polly'
    pollywithcmd += extracmd    # pollywithcmd += ' -fPIC -shared'
    pollywithcmd += ' -o %s_%s' % (basename, compiler)
    pollywithcmd += ' %s.cpp' % basename

    # pluto += ' --tile'
    # pluto += ' -fPIC -shared -o %s_pluto.so %s.cpp' % (basename,basename)

    return {
        'g++':cc,
        'clang':clang,
        'polly':polly,
        'pollyvector':pollywithcmd,
        'pollyparallel':pollywithcmd,
        'pollyboth':pollywithcmd,
        'pollynotiling':pollywithcmd,
        'pollynoaliasing':pollywithcmd,
        'pollypocc':pollywithcmd,
        'pollyshow':pollywithcmd,
        'pollyfunc':pollywithcmd,
        'pollyexport':pollywithcmd,
        'pollymain':pollywithcmd
        # 'pluto':pluto           ##### careful with this, pluto already in polly
    }[compiler]


def runlib(bench,basename,compiler):
    try: 
        with bench.timed_region("running time with %s" % compiler):
            # my_lib = ctypes.cdll.LoadLibrary(os.getcwd() + "/" + basename + '_'+ compiler)
            # my_lib.main()
            # #my_lib.func()
            # ctypes.cdll.LoadLibrary('libdl.so').dlclose(my_lib._handle)
            perf = "perf stat -e task-clock,cycles,instructions,cache-references,cache-misses "
            log = "-o " + "log_default " + "--append "
            cmd = perf + log + os.getcwd() + "/" + basename + '_'+ compiler
            subprocess.check_call(cmd,shell=True)

    except OSError as e:
        print "running error:", e

    #with bench.timed_region('running time with %s  ' % compiler):
    

def pollyoption(x):
    if(x=='pollyvector'):
        return '-vectorizer=stripmine'
    elif(x== 'pollyparallel'):
        return '-parallel -lgomp'
    elif(x== 'pollyboth' ):
        return '-only-func=critical -mllvm -polly-no-tiling' #-mllvm -polly-ignore-aliasing 
    elif(x== 'pollynotiling'):
        return '-no-tiling'
    elif(x== 'pollynoaliasing'):
        return '-ignore-aliasing'
    elif(x== 'pollypocc'):
        return '-optimizer=none'    #doesn't work properly because didn't compile with scoplib
    elif(x== 'pollyshow'):
        return '-show'          # doesn't work somehow
    elif(x== 'pollyfunc'):
        return '-only-func=critical'    
    elif(x== 'pollymain'):
        return '-only-func=noncrit'
    elif(x== 'pollyexport'):
        return '-export'
    else :
        return ''
    

class myBench(Benchmark):

    warmups = 0
    repeats = 1
    method = 'benchmarking'
    benchmark = 'myBench'

    # compilers = ['g++','clang','polly','pollyvector','pollyparallel','pollynotiling','pollyboth','pollynoaliasing',
        # 'pollymain','pollyfunc']
   
    # compilers = ['pollyvector','pollyparallel','polly','clang','pollynotiling','pollyfunc' ]
    # compilers = ['g++' , 'clang','polly','pollynotiling','pollyfunc','pollyparallel','pollynoaliasing','pollyvector' ]
    compilers = ['clang', 'polly', 'pollynotiling']#,'pollynotiling','pollyparallel','pollyvector']

    params =  [('opt_level', range(3,4)),('compiler',compilers)]
    #when adding new options , add in switchcompiler, pollyoptions, 
 
    basename = 'test3d_new'
    # only compile once 
    compiled =dict ((comp,False) for comp in compilers)

    #params automatically change which parameter of the function,
    #and runs, gives total time
    #This saves user writting loops inside the function 

    def benchmarking(self, opt_level=0,times=3 , compiler=''):
        # #set(4)
        self.series['times'] = times
        if not self.compiled[compiler]:
            compile(self,compiler,opt_level,self.basename)

        self.compiled[compiler] = True
        runlib(self,self.basename,compiler)
        

# os.environ["OMP_NUM_THREADS"] = "4"
# os.environ["KMP_AFFINITY"]="granularity=thread,scatter"
myBench().main()

# remainder : python wrapper.py -b  -s    -l -- opt_level=3
#                               run save        gives argument to the method 
