from pybench import Benchmark
from eigenwave3d import default
from os import path
import subprocess

_test_dir = path.join(path.dirname(__file__), "src")


class pluto_tile_bench(Benchmark):
    """
    tile_size has to be a string contains the tile sizes,
    e.g. "4 8", "4 4 1000", "nopluto", "default"

    example of using this script with bash:
    declare -a arr=("2 4 8" "4 4 16")
    for i in "${arr[@]}"
    do
        python tests/pluto_tile_test.py -b -s -l -- nthreads=4 tile_size="$i"
    done

    python tests/pluto_plot.py -ts "${arr[@]}"
    """

    warmups = 0
    repeats = 1
    method = 'benchmarking'
    benchmark = 'pluto_tile_bench'

    def test(self, compiler=None, execute=False, nthreads=1,
             output=False, profiling=False, papi_events=[],
             pluto=False, tile_size=''):
        grid = default(compiler=compiler, nthreads=nthreads, output=output, profiling=profiling,
                       papi_events=papi_events, pluto=pluto, tile=tile_size, execute=False)
        with self.timed_region("tiling size: %s" % tile_size):
            cmd = str(grid)
            print cmd
            subprocess.check_call(cmd)

    def benchmarking(self, compiler='gnu', nthreads=1,
                     output=False, profiling=False, papi_events=[],
                     pluto=True, tile_size=''):
        print "benchmarking with tile size: ", tile_size
        self.series['tile_size'] = tile_size
        tile_size = tile_size
        pluto = pluto
        if(tile_size == "default"):
            tile_size = ''
        elif(tile_size == "nopluto"):
            pluto = False

        self.test(compiler=compiler, nthreads=nthreads, output=output, profiling=profiling,
                  papi_events=papi_events, pluto=pluto, tile_size=tile_size)

pluto_tile_bench().main()
