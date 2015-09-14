from pybench import Benchmark
from eigenwave3d import eigenwave3d
from os import path

_test_dir = path.join(path.dirname(__file__), "src")


def generate_pluto_grid():
    domain_size = (1.0, 1.0, 1.0)
    grid_size = (100, 100, 100)
    dt = 0.002
    tmax = 1.0
    filename = path.join(_test_dir, 'eigenwave3d.cpp')
    grid = eigenwave3d(domain_size, grid_size, dt, tmax,
                       o_converge=True, omp=True, simd=False,
                       ivdep=True, filename=filename, pluto=True)
    return grid


def default(bench, compiler='gnu', nthreads=1,):
    """Eigenwave test case on a unit cube grid (100 x 100 x 100)
    """
    domain_size = (1.0, 1.0, 1.0)
    grid_size = (100, 100, 100)
    dt = 0.002
    tmax = 1.0
    filename = path.join(_test_dir, 'eigenwave3d.cpp')

    grid = eigenwave3d(domain_size, grid_size, dt, tmax,
                       o_converge=True, omp=True, simd=False,
                       ivdep=True, filename=filename)
    grid.compile(filename, compiler=compiler, shared=False)
    with bench.timed_region("default "):
        grid.execute(filename, compiler=compiler, nthreads=nthreads)
    grid.convergence()


def pluto(bench, grid, compiler='gnu', nthreads=1, tile='32 32 32'):
    if tile:
        f = open('tile.sizes', 'w')
        f.write(tile)
        f.close()

    filename = path.join(_test_dir, 'eigenwave3d.cpp')
    grid.compile(filename, compiler=compiler, shared=False)
    filename_p = grid.pluto_op(filename)
    grid.src_file = filename_p
    if compiler in ['clang', 'clang++']:
        # an ugly fix, but pluto will always attack <omp.h> to the first line
        # which would fail clang
        with open(filename_p, 'r') as fin:
            data = fin.read().splitlines(True)
        with open(filename_p, 'w') as fout:
            fout.writelines(data[1:])

    grid.compile(filename_p, compiler=compiler, shared=False)

    with bench.timed_region("pluto with tile size%s" % tile):
        grid.execute(filename_p, compiler=compiler, nthreads=nthreads)
    grid.convergence()


class pluto_tile_bench(Benchmark):
    """
    tile_size has to be a string contains the tile sizes,
    e.g. "4 8", "4 4 1000", "nopluto"

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

    sizes = ['4 4 1000']
    #params = [('tile_size', sizes)]
    grid = generate_pluto_grid()

    def benchmarking(self, compiler='gnu', nthreads=1, tile_size='32 32 32'):
        print "benchmarking: ", nthreads, compiler
        self.series['tile_size'] = tile_size
        if tile_size == 'nopluto':
            default(self, compiler=compiler, nthreads=nthreads)
        else:
            pluto(self, self.grid, compiler=compiler, nthreads=nthreads, tile=tile_size)

pluto_tile_bench().main()
