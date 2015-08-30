import pytest
from os import path
from eigenwave3d import eigenwave3d


@pytest.fixture
def testdir():
    return path.join(path.abspath(path.dirname(__file__)), "src")


@pytest.fixture
def domain_size():
    return (1.0, 1.0, 1.0)


@pytest.fixture
def grid_size():
    return (100, 100, 100)


@pytest.fixture
def dt():
    return 0.002


@pytest.fixture
def tmax():
    return 1.0


@pytest.fixture
def grid(domain_size, grid_size, dt, tmax):
    return eigenwave3d(domain_size, grid_size, dt, tmax)


@pytest.mark.parametrize(('ivdep', 'simd', 'output'),
                         [(True, False, False),
                          (False, True, False),
                          (False, False, True)])
def test_parallel_single(grid, testdir, ivdep, simd, output):
    filename = path.join(testdir, 'eigenwave3d.cpp')
    grid.set_switches(omp=True, double=False, output_vts=output,
                      simd=simd, ivdep=ivdep)
    grid.compile(filename, compiler='g++')
    grid.execute(filename, compiler='g++')
    grid.convergence()


@pytest.mark.parametrize(('ivdep', 'simd', 'output'),
                         [(True, False, False),
                          (False, True, False),
                          (False, False, True)])
def test_parallel_double(grid, testdir, ivdep, simd, output):
    filename = path.join(testdir, 'eigenwave3d.cpp')
    grid.set_switches(omp=True, double=True, output_vts=output,
                      simd=simd, ivdep=ivdep)
    grid.compile(filename, compiler='g++')
    grid.execute(filename, compiler='g++')
    grid.convergence()
