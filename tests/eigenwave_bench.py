import os
import subprocess
from collections import defaultdict
import platform

from pybench import Benchmark


def switchcompiler(compiler, basename, opt_level):
    """ Produce compiler command to run """
    cc = 'g++'
    cc += ' -fopenmp -lpapi -std=c++11 -mavx -DOPESCI_HAVE_PAPI -O%d' % opt_level
    cc += " -o %s_g++" % basename
    cc += " opesciProfiling.cpp %s.cpp 2> %s.log" % (basename, basename)

    clang = 'clang++'
    clang += ' -fopenmp -DOPESCI_HAVE_PAPI'
    clang += ' -o %s_clang' % basename
    clang += ' %s.cpp 2> %s.log'% (basename, basename)
    clang += ' -O%d' % opt_level

    intel = 'icpc'
    intel += ' -openmp -lpapi -std=c++11 -mavx -DOPESCI_HAVE_PAPI'
    intel += ' -o %s_icpc' % basename
    intel += ' opesciProfiling.cpp %s.cpp 2> %s.log'% (basename, basename)
    intel += ' -O%d' % opt_level

    return {
        'g++' : cc,
        'clang' : clang,
        'icpc' : intel,
    }[compiler]

class PropagatorBench(Benchmark):
    """
    Benchmarking tool for FD Propagator code

    Execute benchmark runs with:
    python prop_bench.py -b -s -l -- basename=<basename> compiler=<compiler> opt_level=<opt_level>
    """
    warmups = 0
    repeats = 1

    method = 'propagator'
    benchmark = 'Propagator'
    stream_filename = 'stream_file'
    trial_filename = ''

    compiled = defaultdict(lambda: False)
    is_outfile_opened = False
    outfile = ''

    def __del__(self):
	if self.is_outfile_opened:
 	   self.outfile.close()

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

        #self.compiled[compiler] = True

    def runlib(self, basename, compiler):
        try:
            cmd = os.getcwd() + "/" + basename + '_'+ compiler
            print "Running:", cmd
            with self.timed_region("total"):
                # This is still rather crude since we're timing the OS
                # and Python overheads of invoking the execution command
                #subprocess.check_call(cmd,shell=True)
		output = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).stdout.read()
		self.save_file(output)			
        except OSError as e:
            print "running error:", e

    def propagator(self, basename='eigenwave', compiler='g++', opt_level=3,
                   nthreads=1, affinity='close'):
        self.series['compiler'] = compiler
        self.series['basename'] = basename
        self.series['opt_level'] = opt_level
        self.series['nthreads'] = nthreads
        self.series['affinity'] = affinity

        # Parallel thread settings
        os.environ["OMP_NUM_THREADS"] = str(nthreads)
        if affinity in ['close', 'spread']:
            os.environ["OMP_PROC_BIND"] = affinity
        elif affinity in ['compact', 'scatter']:
            os.environ["KMP_AFFINITY"] = "granularity=thread,%s" % affinity
        else:
            print """ERROR: Only the following affinity settings are supported:
 * OMP_PROC_BIND: 'close', 'spread'
 * KMP_AFFINITY: 'compact', 'scatter'"""
            raise ValueError("Unknown thread affinity setting: %s")
	
	if not self.is_outfile_opened:		
	   self.trial_filename = os.getcwd() + "/results/" + basename + '_' + compiler + '_O' + str(opt_level) + '_' + str(nthreads) + '_' + affinity + '.dat'
	   self.is_outfile_opened = True
	   self.outfile = open(self.trial_filename,'w')
        
	self.compile(compiler, basename, opt_level)
        self.runlib(basename, compiler)
	self.save_file('rpeak' + ' ' + str(self.get_theoretical_peak()) + '\n')
	self.compileStream(compiler,self.stream_filename, opt_level)
	self.runStream(compiler)
	

    def save_file(self, string):
	#outfile = open(self.trial_filename,'a')
	self.outfile.write(string)
	#outfile.close()

    def compileStream(self, compiler, basename, opt_level):
	""" Compile stream to get the memory bandwith. """
	if self.compiled[compiler]:
           return

        cmd = switchcompiler(compiler, self.stream_filename, opt_level)
        try:
            print "Compiling:", cmd
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            print "Compilation error: ", e
            raise Exception("Failed to compile Stream benchmark! ")

	self.compiled[compiler] = True

    def runStream(self,compiler):
	""" Run the stream binary. """
	try:
            print "Running: streamfile"
            cmd = os.getcwd() + '/' + self.stream_filename + '_'+ compiler
	    output = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).stdout.read()
	    self.save_file(output)
        except subprocess.CalledProcessError as e:
            print "Running error: ", e
            raise Exception("Failed to run stream! ")


    def calculate_theoretical_peak(self,cores,frequency,flops_per_cycle,fma_instruction):
        """
        Generic Calculator the processor theoretical peak.
        :param cores: number total of cores
        :param frequency: processor's frequency in Ghz
        :param flops_per_cycle: total number of flops computed in one cycle
        :param fma_instruction: Extension that enable add and mult at the same cycle
        """   
        return cores*frequency*flops_per_cycle*fma_instruction
   

    def get_theoretical_peak(self):
        """
        Get characteristics from the processor in order to calculate theoretical peak.
        """
        cores = 0
        frequency = 0
        ##sockets = 0
        sp_flops_per_cycle = 0
        fma_instruction=0
        mhz_to_ghz=1000
        if platform.system() == "Linux":
            cmd_cores = "cat /proc/cpuinfo | grep processor | wc -l"
            cmd_frequency = "lscpu | grep -i mhz | awk -F : '{print $2}'"
            cmd_threads = "lscpu | grep Thread | awk -F : '{print $2}'"
            ##cmd_sockets = "lscpu | grep Socket | awk -F : '{print $2}'"
            cmd_has_avx = "cat /proc/cpuinfo | grep avx"
	    cmd_has_fma = "cat /proc/cpuinfo | grep 'fma\|fma4'"
           
            cores_total = subprocess.Popen(cmd_cores, stdout=subprocess.PIPE, shell=True).stdout.read()
            frequency = subprocess.Popen(cmd_frequency, stdout=subprocess.PIPE, shell=True).stdout.read()
            threads = subprocess.Popen(cmd_threads, stdout=subprocess.PIPE, shell=True).stdout.read()
            ##sockets = subprocess.Popen(cmd_sockets, stdout=subprocess.PIPE, shell=True).stdout.read()
            is_avx = subprocess.Popen(cmd_has_avx, stdout=subprocess.PIPE, shell=True).stdout.read()
	    is_fma = subprocess.Popen(cmd_has_fma, stdout=subprocess.PIPE, shell=True).stdout.read()       

            cores = float(cores_total)/float(threads)
            frequency = float(frequency)/mhz_to_ghz
            if is_avx != "":
                sp_flops_per_cycle = 8
            else:
                sp_flops_per_cycle = 4
	    if is_fma != "":
	       fma_instruction=2
	    else:
	       fma_instruction=1           
            #DEBUG
            #print cores,frequency,dp_flops_per_cycle,sp_flops_per_cycle
        else:
            raise ValueError("Platform not supported!")
        return self.calculate_theoretical_peak(cores,frequency,sp_flops_per_cycle,fma_instruction)

PropagatorBench().main()
