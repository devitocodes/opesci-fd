from pybench import Benchmark, parser
import subprocess
from jinja2 import Template

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

def roofline(basename, compiler, opt_level, nthreads, affinity, arithmetic_intensity, file_path, precision):
    """
    Roofline plot.
    :param basename: basename from the benchmark
    :param compiler: compiler used on the benchmark
    :param opt_level: optimization level used on the benchmark
    :param nthreads: number of threads used on the benchmark
    :param affinity: affinity used on the benchmark 
    :param arithmetic_intensity: arithmetic intensity of the code
    :param file_path: path to result file of the benchmark
    :param precision: double or single
    """
    rpeak = []
    bw = []
    flops = []
    basenameCount = len(basename)
    for i in range(0, basenameCount):	  
        filename = basename[i] + '_' + compiler[0] + '_O' + str(opt_level[0]) + '_' + str(nthreads[0])  + '_' + affinity[0] #+ '.dat'
        #print filename	
    
        if file_path != '':
           print "Working with result file:" + file_path[0] + '/' + filename
           filename = file_path[0] + '/' + filename
        else:
	    filename = 'results/' + filename	 
        try:
           subprocess.check_call('ls ' + filename + '.dat', shell=True)
           rpeak_ = "grep rpeak " + filename + '.dat' + " | awk  '{sum +=$2}; END {print sum/NR}'"
           rpeak.append(float(subprocess.Popen(rpeak_ , stdout=subprocess.PIPE, shell=True).stdout.read()))
           if precision[0] == 'double':
              rpeak[i] = float(rpeak[i])/2
           bw_ = "grep Triad " + filename + '.dat' + " | awk  '{sum +=$2}; END {print sum/NR/1000}'"
           bw.append(float(subprocess.Popen(bw_, stdout=subprocess.PIPE, shell=True).stdout.read()))
           if bw == '':
              raise Exception("Failed! No Bandwidth output!")	
           flops_ = "grep MFlops/s " + filename + '.dat' + " | awk  '{sum +=$4}; END {print sum/NR/10^3}'"
           flops.append(float(subprocess.Popen(flops_ , stdout=subprocess.PIPE, shell=True).stdout.read()))
           if flops == '':
              raise Exception("Failed! No Mflops output!")	
        except subprocess.CalledProcessError as e:
              print "Result file error: ", e
              raise Exception("Failed to find result file!")
        #print rpeak, bw, flops, arithmetic_intensity
    
    if basenameCount > 1:
       filename = 'results/compareCodes'	
    bwTotal = sum(bw)/len(bw)
    rpeakTotal = sum(rpeak)/len(rpeak)
    generatePlotDatFile('result.dat', rpeakTotal/bwTotal, bwTotal, rpeakTotal, arithmetic_intensity)
    generateGnuplotScript(filename + '.pdf', basename, precision[0],float(rpeakTotal)*2,'result.dat')
    plotRoofline()

def generatePlotDatFile(filename, machineAI, machineBW, machineTopGFps, codeAI):
    """
    Generate a dat file with data to generate the roofline plot.
    :param filename: name to result file
    :param machineAI: arithmetic intensity of the code 
    :param machineTopGFps: machine theoretical performance
    :param machineBW: machine bandwidth performance
    :param codeAI: arithmetic intensity of the codes
    """  
    # Use these default values, that are in range with Xeon computers, the values are chosen to give a smooth roofline plot. 
    xAxis = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32]
	
    # Insert the machine AI in the xAxis values.
    insertIndex = 0
    for x in xAxis:
	if x > machineAI:
	   insertIndex = xAxis.index(x)
	   break
	
    if insertIndex != 0:
       xAxis.insert(insertIndex, machineAI)

    generateDatFile(filename, xAxis, machineAI, machineBW, machineTopGFps, codeAI)


def generateDatFile(filename, xAxis, machineAI, machineBW, machineTopGFps, codeAI):
    """
    Generate a dat file with data to generate the roofline plot.
    :param filename: name to result file
    :param xAxis: plot x range
    :param machineAI: arithmetic intensity of the code 
    :param machineTopGFps: machine theoretical performance
    :param machineBW: machine bandwidth performance
    :param codeAI: arithmetic intensity of the codes
    """
    
    # The y axis to be filled.
    yAxis = []	

    for x in xAxis:
	if x < machineAI:
 	   yAxis.append(x * machineBW)
	else:
           yAxis.append(machineTopGFps)
		
    writeAxisValuesToFile(filename, xAxis, yAxis, "#AI", "#Theo", machineTopGFps, machineBW, codeAI)

def writeAxisValuesToFile(filename, xAxis, yAxis, xLabel, yLabel, machineTopGFps, machineBW, codeAI):
    """
    Has the logic to write the axis to the dat file.
    :param filename: name to result file
    :param xAxis: plot x range
    :param yAxis: plot y range
    :param xLabel: file label 
    :param yLabel: file label
    :param machineTopGFps: machine theoretical performance
    :param machineBW: machine bandwidth performance
    :param codeAI: arithmetic intensity of the codes
    """
    try:
	aiValue = []
	codeAIBW = []	
	f = open(filename, 'w')
	f.write(xLabel.ljust(11) + yLabel.ljust(11))
		
	count = len(xAxis)
	ailength = len(codeAI)
	for i in range(0, ailength):
	   f.write((' #AI_code' + str(i)).ljust(11) + (' #AIcode' + str(i) + '*bw').ljust(11))
	f.write('\n');
	for i in range(0, ailength):
	   aiValue.append(('%g' % float(codeAI[i])).ljust(11))
	   codeAIBW.append(float(codeAI[i])*machineBW)
	   if codeAIBW[i] > machineTopGFps:
              codeAIBW[i] = machineTopGFps
        for i in range(0, count):
	    xValue = '%g' % xAxis[i]
	    # Align the columns!
	    xValue = xValue.ljust(11)
	    #f.write(xValue + ('%g \n' % yAxis[i]))
	    yValue = ('%g ' % yAxis[i]).ljust(11)
	    f.write(xValue + yValue)
#	    print "aa"
	    if i != 0:
	       for i in range(0, ailength):
                 f.write(aiValue[i].ljust(12) + ('%g ' % codeAIBW[i]).ljust(11))		  
	    else:
	       for i in range(0, ailength):
	         f.write(aiValue[i].ljust(12) + ('%g ' % 0).ljust(11))
	    f.write('\n')

	f.close()
    except Exception as e:
	   print e
	   print "Error writing to file"

def generateGnuplotScript(scriptName, basename, precision_, rpeak_, data_):
    """
    Generates the gnuplot script based on the template found in the roofline folder.
    :param scriptName: name of the output file
    :param basename: list with name of the codes
    :param precision_: double or single
    :param rpeak_: machine theoretical performance
    :param data_: filename of the result file
    """
    fileContent = ""
    try:
	f = open('roofline/roofline_plot.in', 'r')

	fileContent = f.read()

	f.close()
    except Exception as e:
	print e
	print "Error reading the template file"
	return

    plotTemplate = Template(fileContent,trim_blocks=True,keep_trailing_newline=True)
    gnuplotScript = plotTemplate.render({'plotFilename': scriptName,'precision': precision_, 'rpeak': rpeak_, 'data' : data_, 'basename' : basename})

    try:
	f = open('roofline_plot.plt', 'w')
	f.write(gnuplotScript)

	f.close()
    except Exception as e:
	print e
	print "Error writing to file"

def plotRoofline():
    """
    Generates the pdf file.
    """
    try:
        print "Ploting roofline... \n"
        subprocess.check_call("gnuplot roofline_plot.plt", shell=True)
    except OSError as e:
        print "running error:", e


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
    p.add_argument('--ai', type=str, nargs='+',
                   help='Code AI')
    p.add_argument('--precision', type=str, nargs='+',
                   help='Single or double precision values')
    p.add_argument('--file_path', type=str, nargs='+',
                   help='Path to benchmark result file')
    p.add_argument('mode', choices=('roofline', 'normal'),
                   nargs='?', default='normal', help='Mode choice')
    args = p.parse_args()
    basename = args.basename or ['eigenwave']
    compiler = args.compiler or ['g++', 'clang', 'icpc']
    opt_level = args.opt_level or [2, 3]
    nthreads = args.parallel or [1]
    affinity = args.affinity or ['close']
    arithmetic_intensity = args.ai or [0]
    file_path = args.file_path or ''
    precision = args.precision or ['single']	

    if args.mode == 'roofline':
       if len(basename) != len(arithmetic_intensity):
	  raise Exception("Failed! Number of basename and ai must be the same!")
       else:	 
	  roofline(basename, compiler, opt_level, nthreads, affinity, arithmetic_intensity, file_path, precision)
    else:
        b = PropagatorPlot(benchmark='Propagator-Performance',
                      resultsdir=args.resultsdir, plotdir=args.plotdir)
        b.combine_series([('basename', basename), ('compiler', compiler),
                      ('opt_level', opt_level), ('nthreads', nthreads),
                      ('affinity', affinity)], filename='Propagator')
        if len(nthreads) > 1:
           b.plot_parallel_scaling(nthreads)
        else:
           b.plot_compiler_comparison(opt_level)
