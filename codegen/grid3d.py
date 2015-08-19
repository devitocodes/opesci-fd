"""
main entry script to create c++ source code for pre-defined test cases
- calls run_test in grid_test.py
- travis.yml script runs the default test case
"""

from grid_test import *
import sys


def default():
    """
    default test case
    eigen wave in unit cube, 100 x 100 x 100 grids
    output L2 norm
    """
    domain_size = (1.0, 1.0, 1.0)
    grid_size = (100, 100, 100)
    dt = 0.002
    tmax = 1.0
    run_test(domain_size, grid_size, dt, tmax, o_step=False, o_converge=True,
             omp=True, simd=False, ivdep=True, io=False,
             filename='src/tests/test3d.cpp')


def default_vtk():
    """
    default test case
    eigen wave in unit cube, 100 x 100 x 100 grids
    output L2 norm
    """
    domain_size = (1.0, 1.0, 1.0)
    grid_size = (100, 100, 100)
    dt = 0.002
    tmax = 1.0
    run_test(domain_size, grid_size, dt, tmax, o_step=True, o_converge=True,
             omp=True, simd=False, ivdep=True, io=True,
             filename='src/tests/test3d_vtk.cpp')


def read_data():
    """
    test case for reading data (media parameters) from input file
    eigen wave in unit cube, 200 x 200 x 200 grids
    """
    domain_size = (1.0, 1.0, 1.0)
    grid_size = (195, 195, 195)
    dt = 0.002
    tmax = 1.0
    run_test(domain_size, grid_size, dt, tmax, o_step=True, o_converge=False,
             omp=True, simd=False, ivdep=True, io=True, read=True,
             filename='src/tests/test3d_read.cpp',
             rho_file='RHOhomogx200',
             vp_file='VPhomogx200',
             vs_file='VShomogx200')


def cx1():
    """
    test case for comparison between pragma simd and pragma ivdep on cx1
    """
    domain_size = (1.0, 1.0, 1.0)
    grid_size = (200, 200, 200)
    dt = 0.001
    tmax = 5.0
    run_test(domain_size, grid_size, dt, tmax, o_step=False, o_converge=False,
             omp=True, simd=False, ivdep=True, io=False,
             filename='src/tests/test3d_ivdep.cpp')
    run_test(domain_size, grid_size, dt, tmax, o_step=False, o_converge=False,
             omp=True, simd=True, ivdep=False, io=False,
             filename='src/tests/test3d_simd.cpp')


def converge_test():
    """
    - test case for convergence analysis of (2,4) scheme
    - eigen wave in unit cube
    - start with spacing 0.1, spacing halves for each test
    - dt reduce by factor of 4 for each test
    """
    domain_size = (1.0, 1.0, 1.0)
    s = 10
    c = 0.4*s
    dt = c/(s**2)
    tmax = 5.0
    run_test(domain_size, (s, s, s), dt, tmax, o_step=False, o_converge=True,
             omp=True, simd=False, ivdep=True, io=False,
             filename='tmp/test3d_'+str(s)+'.cpp')

    s = s*2
    dt = c/(s**2)
    run_test(domain_size, (s, s, s), dt, tmax, o_step=False, o_converge=True,
             omp=True, simd=False, ivdep=True, io=False,
             filename='tmp/test3d_'+str(s)+'.cpp')

    s = s*2
    dt = c/(s**2)
    run_test(domain_size, (s, s, s), dt, tmax, o_step=False, o_converge=True,
             omp=True, simd=False, ivdep=True, io=False,
             filename='tmp/test3d_'+str(s)+'.cpp')

    s = s*2
    dt = c/(s**2)
    run_test(domain_size, (s, s, s), dt, tmax, o_step=False, o_converge=True,
             omp=True, simd=False, ivdep=True, io=False,
             filename='tmp/test3d_'+str(s)+'.cpp')


def main():
    if len(sys.argv) > 2:
        print 'too many arguments!'
        return

    if len(sys.argv) == 1:
        default()
        return

    if sys.argv[1] == 'cx1':
        cx1()
        return

    if sys.argv[1] == 'converge':
        converge_test()
        return

    if sys.argv[1] == 'read':
        read_data()
        return

    if sys.argv[1] == 'vtk':
        default_vtk()
        return

if __name__ == "__main__":
    main()
