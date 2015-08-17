import os
import subprocess
from collections import defaultdict

from pybench import Benchmark


def switchcompiler(compiler, basename, opt_level):
    """ Produce compiler command to run """
    cc = 'g++'
    cc += ' -fopenmp -O%d' % opt_level
    cc += " -o %s_g++" % basename
    cc += " %s.cpp" % basename

    clang = 'clang++'
    clang += ' -fopenmp'
    clang += ' -o %s_clang' % basename
    clang += ' %s.cpp'% basename
    clang += ' -O%d' % opt_level

    return {
        'g++' : cc,
        'clang' : clang,
    }[compiler]

class PropagatorBench(Benchmark):
    """
    Benchmarking tool for FD Propagator code

    Execute benchmark runs with:
    python prop_bench.py -b -s -l -- basename=<basename> compiler=<compiler> opt_level=<opt_level>
    """
    warmups = 0
    repeats = 3

    method = 'propagator'
    benchmark = 'Propagator'

    compiled = defaultdict(lambda: False)

    def compile(self, compiler, basename, opt_level):
        if self.compiled[compiler]:
            return

        cmd = switchcompiler(compiler, basename, opt_level)
        try:
            print "Compiling:", cmd
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            print "Compilation error: ", e
            raise Exception("Failed to compile ")

        self.compiled[compiler] = True

    def runlib(self, basename, compiler):
        try:
            cmd = os.getcwd() + "/" + basename + '_'+ compiler
            print "Running:", cmd
            with self.timed_region("total"):
                # This is still rather crude since we're timing the OS
                # and Python overheads of invoking the execution command
                subprocess.check_call(cmd,shell=True)
        except OSError as e:
            print "running error:", e

    def propagator(self, basename='test3d', compiler='g++', opt_level=3):
        self.series['compiler'] = compiler
        self.series['basename'] = basename
        self.series['opt_level'] = opt_level

        self.compile(compiler, basename, opt_level)
        self.runlib(basename, compiler)


PropagatorBench().main()
