from pybench import Benchmark
from eigenwave3d import default
from os import path

_test_dir = path.join(path.dirname(__file__), "src")


class pluto_tile_bench(Benchmark):
    """
    tile_size has to be a string contains the tile sizes,
    e.g. "4 8", "4 4 1000", "default"

    example of using this script with bash:
    declare -a arr=("2 4 8" "4 4 16")
    for i in "${arr[@]}"
    do
        python tests/pluto_tile_test.py -b -s -l -- nthreads=4 tile_size="$i"
    done

    declare -a ser="";
    for i in "${arr[@]}"
    do
        ser+="\"$i\" "
    done
    python tests/pluto_plot.py -ts $ser;
    """

    warmups = 0
    repeats = 1
    method = 'benchmarking'
    benchmark = 'pluto_tile_bench'

    def benchmarking(self, compiler='gnu', nthreads=1,
                     output=False, profiling=False, papi_events=[],
                     pluto=False, tile_size=''):
        print "benchmarking with tile size: ", tile_size
        self.series['tile_size'] = tile_size
        default(compiler=compiler, nthreads=nthreads, output=output, profiling=profiling,
                papi_events=papi_events, pluto=pluto, tile=tile_size)

pluto_tile_bench().main()
