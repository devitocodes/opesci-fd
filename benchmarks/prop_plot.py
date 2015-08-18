from pybench import Benchmark, parser

class PropagatorPlot(Benchmark):
    figsize = (6, 4)


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

    regions = ['total']
    if len(nthreads) > 1:
        groups = ['compiler', 'opt_level', 'affinity']
        for region in regions:
            b.plot(figsize=b.figsize, format='pdf', figname='Prop3D_%s'%region,
                   xaxis='nthreads', xlabel='Number of threads', xticklabels=nthreads,
                   regions=[region], groups=groups,
                   kinds='semilogx', title='Performance: %s'%region, legend={'loc': 'best'})
    else:
        groups = ['compiler']
        opt_str = ['-O%s' % o for o in opt_level]
        for region in regions:
            b.plot(figsize=b.figsize, format='pdf', figname='Prop3D_%s'%region,
                   xaxis='opt_level', xvals=opt_level, xticklabels=opt_str,
                   xlabel='Compiler configuration', groups=groups, regions=[region],
                   kinds='bar', title='Performance: %s'%region, legend={'loc': 'best'})
