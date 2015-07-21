from pybench import parser
from wrapper import myBench
import pylab as plt



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

    plt.bar(range(0,len(times)),times,align = 'center')
    plt.xticks(range(0,len(times)), labels)
    plt.ylabel('running time ')
    plt.savefig('my_fig')
    
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



