from pybench import Benchmark, parser

class Eigenwave3DPlot(Benchmark):
    benchmark = 'Eigenwave3D'
    figsize = (6, 4)
    profileregions = ['execute']

    def plot_compiler_comparison(self, vectorise):
        groups = ['compiler']
        opt_str = ['#pragma %s' % o for o in vectorise]
        for region in self.profileregions:
            self.plot(figsize=self.figsize, format='pdf', figname='Eigen3D_CC_%s'%region,
                      xaxis='vectorise', xvals=vectorise, xticklabels=opt_str,
                      xlabel='Compiler configuration', groups=groups, regions=[region],
                      kinds='bar', title='Performance: %s'%region, legend={'loc': 'best'})

    def plot_parallel_scaling(self, nthreads):
        groups = ['compiler', 'affinity']
        for region in self.profileregions:
            self.plot(figsize=self.figsize, format='pdf', figname='Eigen3D_OMP_%s'%region,
                      xaxis='nthreads', xticklabels=nthreads, xlabel='Number of threads',
                      regions=[region], groups=groups, xmax=nthreads[-1], trendline='Perfect speedup',
                      kinds='loglog', title='Performance: %s'%region, legend={'loc': 'best'})


if __name__ == '__main__':
    p = parser(description="Benchmark plotter for Eigenwave3D test.")
    p.add_argument('--basename', type=str, nargs='+',
                   help='Name of test used for benchmarking')
    p.add_argument('--compiler', type=str, nargs='+',
                   help='Compilers to compare')
    p.add_argument('--affinity', type=str, nargs='+',
                   help='Parallel thread affintiy setting used')
    p.add_argument('--vectorise', type=str, nargs='+',
                   help='Vectorisation pragma used, can be None')
    args = p.parse_args()
    basename = args.basename or ['eigenwave3d']
    compiler = args.compiler or ['g++', 'icpc']
    nthreads = args.parallel or [1]
    affinity = args.affinity or ['close']
    vectorise = args.vectorise or ['ivdep', 'simd', 'none']

    b = Eigenwave3DPlot(resultsdir=args.resultsdir, plotdir=args.plotdir)
    b.combine_series([('basename', basename), ('compiler', compiler),
                      ('nthreads', nthreads), ('affinity', affinity),
                      ('vectorise', vectorise)])
    if len(nthreads) > 1:
        b.plot_parallel_scaling(nthreads)
    else:
        b.plot_compiler_comparison(vectorise)
