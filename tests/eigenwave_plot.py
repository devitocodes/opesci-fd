from pybench import Benchmark, parser

class PropagatorPlot(Benchmark):
    figsize = (6, 4)
    profileregions = ['total']

    def plot_compiler_comparison(self, opt_level):
        groups = ['compiler']
        opt_str = ['-O%s' % o for o in opt_level]
        for region in self.profileregions:
            self.plot(figsize=self.figsize, format='pdf', figname='PropCC_%s'%region,
                      xaxis='opt_level', xvals=opt_level, xticklabels=opt_str,
                      xlabel='Compiler configuration', groups=groups, regions=[region],
                      kinds='bar', title='Performance: %s'%region, legend={'loc': 'best'})

    def plot_parallel_scaling(self, nthreads):
        groups = ['compiler', 'opt_level', 'affinity']
        for region in self.profileregions:
            self.plot(figsize=self.figsize, format='pdf', figname='PropOMP_%s'%region,
                      xaxis='nthreads', xticklabels=nthreads, xlabel='Number of threads',
                      regions=[region], groups=groups, xmax=nthreads[-1], trendline='Perfect speedup',
                      kinds='loglog', title='Performance: %s'%region, legend={'loc': 'best'})


def switchStreamCompiler(compiler):
    """ Produce compiler command to run """
    cc = 'gcc'
    cc += ' -fopenmp'
    cc += ' -o streamfile'
    cc += ' stream_file.c'
    cc += ' -O3 -mavx'

    clang = 'clang'
    clang += ' -fopenmp'
    clang += ' -o streamfile'
    clang += ' stream_file.cpp'
    clang += ' -O3 -mavx'

    intel = 'icc'
    intel += ' -openmp'
    intel += ' -o streamfile'
    intel += ' stream_file.c'
    intel += ' -O3 -xAVX'

    return {
        'gcc' : cc,
        'clang' : clang,
        'icc' : intel,
    }[compiler]

def compileStream(compiler):
	""" Compile stream to get the memory bandwith. """
        cmd = switchStreamCompiler(compiler)
        try:
            print "Compiling:", cmd
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            print "Compilation error: ", e
            raise Exception("Failed to compile! ")


def runStream():
	""" Run the stream binary. """
	try:
            print "Running: streamfile"
            subprocess.check_call('./streamfile', shell=True)
        except subprocess.CalledProcessError as e:
            print "Running error: ", e
            raise Exception("Failed to run stream! ")


def calculate_theoretical_peak(cores,frequency,flops_per_cycle,fma_instruction):
    """
        Generic Calculator the processor theoretical peak.
    :param cores: number total of cores
    :param frequency: processor's frequency in Ghz
    :param flops_per_cycle: total number of flops computed in one cycle
    :param fma_instruction: Extension that enable add and mult at the same cycle
        """   
    return cores*frequency*flops_per_cycle*fma_instruction
   

def get_theoretical_peak(is_peak_double):
    """
        Get characteristics from the processor in order to calculate theoretical peak.
    :param is_peak_double: True or False for theoretical peak in double or single precision
        """
    cores = 0
    frequency = 0
    ##sockets = 0
    dp_flops_per_cycle = 0
    sp_flops_per_cycle = 0
    fma_instruction=2
    mhz_to_ghz=1000
    if platform.system() == "Linux":
        cmd_cores = "cat /proc/cpuinfo | grep processor | wc -l"
        cmd_frequency = "lscpu | grep -i mhz | awk -F : '{print $2}'"
        cmd_threads = "lscpu | grep Thread | awk -F : '{print $2}'"
        ##cmd_sockets = "lscpu | grep Socket | awk -F : '{print $2}'"
        cmd_is_avx = "cat /proc/cpuinfo | grep avx"
       
        cores_total = subprocess.Popen(cmd_cores, stdout=subprocess.PIPE, shell=True).stdout.read()
        frequency = subprocess.Popen(cmd_frequency, stdout=subprocess.PIPE, shell=True).stdout.read()
        threads = subprocess.Popen(cmd_threads, stdout=subprocess.PIPE, shell=True).stdout.read()
        ##sockets = subprocess.Popen(cmd_sockets, stdout=subprocess.PIPE, shell=True).stdout.read()
        is_avx = subprocess.Popen(cmd_is_avx, stdout=subprocess.PIPE, shell=True).stdout.read()
       
        cores = float(cores_total)/float(threads)
        frequency = float(frequency)/mhz_to_ghz
        if is_avx != "":
            dp_flops_per_cycle = 4
        else:
            dp_flops_per_cycle = 2
        sp_flops_per_cycle = dp_flops_per_cycle*2
        #DEBUG
        #print cores,frequency,dp_flops_per_cycle,sp_flops_per_cycle
    else:
        raise ValueError("Platform not supported!")
    if is_peak_double:
        return calculate_theoretical_peak(cores,frequency,dp_flops_per_cycle,fma_instruction)
    else:
        return calculate_theoretical_peak(cores,frequency,sp_flops_per_cycle,fma_instruction)


if __name__ == '__main__':
    p = parser(description="Performance benchmark plotter for 3D Propagator test.")
    p.add_argument('--basename', type=str, nargs='+',
                   help='Name of test used for benchmarking')
    p.add_argument('--compiler', type=str, nargs='+',
                   help='Compilers to compare')
    p.add_argument('--opt_level', type=int, nargs='+',
                   help='Compiler optimisation levels to compare')
    p.add_argument('--affinity', type=str, nargs='+',
                   help='Parallel thread affintiy setting used')
    args = p.parse_args()
    basename = args.basename or ['test3d']
    compiler = args.compiler or ['g++', 'clang', 'icpc']
    opt_level = args.opt_level or [2, 3]
    nthreads = args.parallel or [1]
    affinity = args.affinity or ['close']

    b = PropagatorPlot(benchmark='Propagator-Performance',
                      resultsdir=args.resultsdir, plotdir=args.plotdir)
    b.combine_series([('basename', basename), ('compiler', compiler),
                      ('opt_level', opt_level), ('nthreads', nthreads),
                      ('affinity', affinity)], filename='Propagator')

    if len(nthreads) > 1:
        b.plot_parallel_scaling(nthreads)
    else:
        b.plot_compiler_comparison(opt_level)
