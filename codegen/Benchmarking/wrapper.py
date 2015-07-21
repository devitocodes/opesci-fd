import os
import subprocess
import ctypes
import time
import matplotlib.pyplot as plt

from pybench import Benchmark,parser


def compile(bench,compiler,opt_level,basename):

    
    cmd = switchcompiler(compiler,opt_level,basename)
    try:
        #with bench.timed_region("compile time with %s" % compiler):
        subprocess.check_call(cmd,shell=True)
    except subprocess.CalledProcessError as e:
        print "Compilation error: ", e
        raise Exception("Failed to compile ")

def switchcompiler(compiler,opt_level,name):
    basename = name
    cc = 'g++'
    cc += ' -fPIC -shared -O%d'% opt_level
    cc += " -o %s_g++.so" % basename
    cc += " %s.cpp" % basename
    cc += " test.h"

    clang = '~/./polly/llvm_build/bin/clang++'
    clang += ' -fPIC -shared'
    clang += ' -o %s_clang.so' % basename
    clang += ' %s.cpp'% basename
    clang += ' -O%d' % opt_level

    polly = '~/./polly/llvm_build/bin/clang++'
    polly += ' -Xclang -load -Xclang LLVMPolly.so'
    polly += ' -O3 -mllvm'

    # pluto = polly

    polly += ' -polly'

    pollywithcmd = polly

    polly += ' -fPIC -shared'
    polly += ' -o %s_polly.so' % basename
    polly += ' %s.cpp' % basename

    extracmd = pollyoption(compiler)    # extra options for polly

    pollywithcmd += extracmd        
    pollywithcmd += ' -fPIC -shared'
    pollywithcmd += ' -o %s_%s.so' % (basename, compiler)
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
        'pollyexport':pollywithcmd
        # 'pluto':pluto           ##### careful with this, pluto already in polly
    }[compiler]


def runlib(bench,basename,compiler):
    try: 
        with bench.timed_region("running time with %s" % compiler):
            my_lib = ctypes.cdll.LoadLibrary(os.getcwd() + "/" + basename + '_'+ compiler + '.so')
            my_lib.main()
            #my_lib.func()
            ctypes.cdll.LoadLibrary('libdl.so').dlclose(my_lib._handle)

    except OSError as e:
        print "running error:", e

    #with bench.timed_region('running time with %s  ' % compiler):
    

def pollyoption(x):
    if(x=='pollyvector'):
        return '-vectorizer=stripmine'
    elif(x== 'pollyparallel'):
        return '-parallel -lgomp'
    elif(x== 'pollyboth' ):
        return '-vectorizer=stripmine -mllvm -polly-parallel -mllvm -polly-ignore-aliasing -lgomp ' #-parallel -lgomp
    elif(x== 'pollynotiling'):
        return '-no-tiling'
    elif(x== 'pollynoaliasing'):
        return '-ignore-aliasing'
    elif(x== 'pollypocc'):
        return '-optimizer=none'    #doesn't work properly because didn't compile with scoplib
    elif(x== 'pollyshow'):
        return '-show'          # doesn't work somehow
    elif(x== 'pollyfunc'):
        return '-only-func=name'    # have to add a functionname here
    elif(x== 'pollyexport'):
        return '-export'
    else :
        return ''
    

class myBench(Benchmark):

    warmups = 0
    repeats = 3
    method = 'benchmarking'
    benchmark = 'myBench'

    # compilers = ['g++','clang','polly','pollyvector','pollyparallel','pollynotiling','pollyboth','pollynoaliasing']
    # compilers = ['pollyparallel','pollynotiling','pollynoaliasing']
    # compilers = ['pollyboth']
    compilers = ['g++']
    params =  [('opt_level', range(2,3)),('compiler',compilers)]
    # opt_level affect clang++ and g++ 

    #when adding new options , add in switchcompiler, pollyoptions, 


    basename = 'test'

    # ('extra',['-parallel -lgomp ','-vectorizer=stripmine','-vectorizer=stripmine-parallel -lgomp'])
    compiled =dict ((comp,False) for comp in compilers)

    #params automatically change which parameter of the function,
    #and runs, gives total time
    #This saves user writting loops inside the function 

    def benchmarking(self, opt_level=0,times=3 , compiler=''):

        
        self.series['times'] = times
        if not self.compiled[compiler]:
            compile(self,compiler,opt_level,self.basename)

        self.compiled[compiler] = True
        runlib(self,self.basename,compiler)

myBench().main()

# remainder : python wrapper.py -b  -s    -l -- opt_level=3
#                               run save        gives argument to the method 
