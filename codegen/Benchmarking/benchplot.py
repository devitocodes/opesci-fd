from pybench import parser
from wrapper import myBench
import pylab as plt
import affinity as aff
import os



class myPlot(myBench):
    figsize = (6, 4)

if __name__ == '__main__':
    p = parser(description="Performance benchmark for optlevel compile time and running time.")
    p.add_argument('-tm', '--times', type=int, nargs='+',
            help='times to plot')
    p.add_argument('-comp','--compiler', type=str,nargs='+',
            help = 'compiler to plot')

    args = p.parse_args()
   # dim = args.dim
    tms = args.times
    compiler = args.compiler


    regions = ['compile time', 'running time']
    labels = {(False, False): 'Implicit', (False, True): 'Explicit',
              (True, True): 'Explicit, zero-tracking'}

    b = myPlot(benchmark='myBench',resultsdir=args.resultsdir, plotdir=args.plotdir)

    # print tms
    # WARNING: have to use -tm 3 to make sure it combines with series 3....
    b.combine_series([('times', tms)], filename='myBench')

    compilers_ = [x for x in b.result['params'] if x[0]=='compiler'][0][1]
    compilers_check = ['running time with '+x for x in compilers_ ]
    # get the right parameter , because in the results the order is sometimes read in a wrong way.

    

    times = []
    # print b.result
    for data in b.result['timings'].values():
        # print data.values()
        for x in data.keys():
            if x in compilers_check:
                times.append(data[x])

    # checking it is the right results, by checking the key name, 


    labels_ = b.result['timings'].keys()
    labels = [x[0] for x in labels_]

    fig = plt.figure();
    ax  = fig.add_axes((.1,.4,.8,.5))

    detail ={'clang':'using clang++ with opt_level given in wrapper , no other options enabled.'
                , 'pollyvector':'polly with vectoer=stripemine '
                , 'pollyparallel':'polly with parallel'
                , 'polly':'polly default'
                , 'pollynotiling':'notiling'
                , 'pollyboth':'all polly'
                , 'g++':'g++ default with opt_level in wrapper'
                , 'pollynoaliasing':'polly with noaliasing option'
                , 'pollyfunc': 'only optimise specific func'
                }


    text = ''
    for l in labels:
        text += l + ': ' + (detail[l]) + '\n'


    print '+++'
    # print aff.get_process_affinity_mask(5)
    # print os.environ['GOMP_CPU_AFFINITY=\"0 1 2 3 \"']
    print os.getpid()
    print '---'
    # plot bar first then add ticks and labels , otherwise no extra space 
    plt.bar(range(0,len(times)),times,width = 0.5, align = 'center')
    plt.xticks(range(0,len(times)), labels, rotation = 17)
    plt.ylabel('running time ')
    fig.text(.1,.0, text)
    plt.savefig('my_fig%s'%b.basename)







    
    # compiler_str = ['%s-haha' % d for d in compilers_check]
    # #for region in regions:
    # region = regions[0]
    # b.plot(figsize=b.figsize, format='pdf',figname='myBench',
    #     xaxis='opt_level', xvals=labels,regions=[region],xticklabels=compiler_str,
    #     xlabel='x discretisation', labels = labels,
    #     kinds='bar', title="")

        #  groups=groups, labels=labels, legend={'loc': 'best'},
        # xticklabels=degree_str, 






    # print b.result



