from os import path, environ
from collections import defaultdict

from eigenwave3d import eigenwave3d
from pybench import Benchmark


class Eigenwave3DBench(Benchmark):
    """
    Benchmarking tool for Opesci Eigenwave3D tests

    Execute benchmark runs with:
    python eigenwave3d_bench.py -b -s -l -- basename=<basename> compiler=<compiler>
    """
    warmups = 0
    repeats = 3

    method = 'opesci_run'
    benchmark = 'Eigenwave3D'

    profileregions = ['execute']

    _compiled = defaultdict(lambda: False)

    def __init__(self, domain_size=(1.0, 1.0, 1.0),
                 grid_size=(100, 100, 100),
                 dt = 0.002, tmax = 1.0):
        super(Eigenwave3DBench, self).__init__()
        self.grid = eigenwave3d(domain_size, grid_size, dt, tmax)

    def opesci_run(self, basename='eigenwave3d', compiler='g++',
                   nthreads=1, affinity='close'):
        self.series['basename'] = basename
        self.series['compiler'] = compiler
        self.series['nthreads'] = nthreads
        self.series['affinity'] = affinity

        # Parallel thread settings
        environ["OMP_NUM_THREADS"] = str(nthreads)
        if affinity in ['close', 'spread']:
            environ["OMP_PROC_BIND"] = affinity
        elif affinity in ['compact', 'scatter']:
            environ["KMP_AFFINITY"] = "granularity=thread,%s" % affinity
        else:
            print """ERROR: Only the following affinity settings are supported:
 * OMP_PROC_BIND: 'close', 'spread'
 * KMP_AFFINITY: 'compact', 'scatter'"""
            raise ValueError("Unknown thread affinity setting: %s")

        # Generate and compile the test model
        testdir = path.join(path.dirname(__file__), "src")
        filename = path.join(testdir, '%s_%s.cpp' % (basename, compiler))
        if not self._compiled[compiler]:
            self.grid.compile(filename, compiler=compiler, shared=True)
            self._compiled[compiler] = True

        # Timed model execution
        with self.timed_region("execute"):
            self.grid.execute(filename)

Eigenwave3DBench().main()
