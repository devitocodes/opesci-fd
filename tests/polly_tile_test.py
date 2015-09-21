from pybench import Benchmark
from eigenwave3d import default
from os import path
import subprocess

_test_dir = path.join(path.dirname(__file__), "src")


class polly_tile_bench(Benchmark):
    """
    tile_size has to be a string contains the tile sizes,
    e.g. "4 8", "4 4 1000", "nopolly", "default"
    example of using this script with bash:
    declare -a arr=("2 4 8" "4 4 16")
    for i in "${arr[@]}"
    do
        python tests/polly_tile_test.py -b -s -l -- nthreads=4 tile_size="$i"
    done
    python tests/pluto_plot.py -ts "${arr[@]}"
    """

    warmups = 0
    repeats = 1
    method = 'benchmarking'
    benchmark = 'polly_tile_bench'

    def test(self, compiler=None, execute=False, nthreads=4,
             output=False, profiling=False, papi_events=[],
             tile_size=''):
        print 'tile in test', tile_size
        grid = default(compiler=compiler, nthreads=nthreads,
                       output=output, profiling=profiling, execute=False,
                       papi_events=papi_events, tile=tile_size)
        with self.timed_region("tiling size: %s" % tile_size):
            cmd = str(grid)
            print 'cmd', cmd
            subprocess.check_call(cmd)

    def benchmarking(self, compiler='polly', nthreads=4,
                     output=False, profiling=False, papi_events=[],
                     polly=True, tile_size=''):
        print "benchmarking with tile size: ", tile_size
        self.series['tile_size'] = tile_size
        tile_size = tile_size
        compiler = compiler
        if(tile_size == "default"):
            tile_size = None
        elif(tile_size == "nopolly"):
            compiler = 'clang'

        print 'tile in benchmarking:', tile_size
        self.test(compiler=compiler, nthreads=nthreads, output=output, profiling=profiling,
                  papi_events=papi_events, tile_size=tile_size)

polly_tile_bench().main()
