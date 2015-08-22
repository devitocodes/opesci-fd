from opesci import *
from sympy import symbols, Eq, sqrt, pi, cos, sin, Float
from os import path
import sys

_test_dir = path.join(path.dirname(__file__), "src")


def run_test(domain_size, grid_size, dt, tmax, o_step=False, o_converge=True,
             omp=True, simd=False, ivdep=True, io=False, double=False,
             filename='test.cpp', read=False, expand=True, eval_const=True,
             rho_file='', vp_file='', vs_file=''):
    """
    create 3D eigen waves and run FD simulation
    :param domain_size: define size of domain
    e.g. (1.0, 1.0, 1.0) for unit cube
    :param spacing: define grid spacing, e.g. (0.1, 0.1, 0.1)
    :param dt: define time step
    :param tmax: define simulation time
    :param o_step: switch for code to output at each time step
    (such as save field as vtk file), default False (no output)
    :param o_converge: switch for code to calculate L2 norms
    for convergence test, default True (output)
    :param omp: swtich for inserting #pragma omp for before outter loop
    default True (use omp)
    :param simd: switch for inserting #pragma simd before inner loop
    default False (do not use simd)
    :param ivddp: switch for inserting #praga ivdep before inner loop
    default True (use ivdep)
    :param io: switch for include io header files
    default False (not include vtk header files)
    :param double: switch for using double as real number variables
    default False (use float)
    :param expand: expand kernel fully (no factorisation), default True
    :param eval_const: evaluate all constants in kernel in generated code default True
    :param filename: output source code file name as string
    :param read: switch for reading meida parameters from input files
    default False (not reading)
    :param rho_file: file name for input file of rho (density)
    :param vp_file: file name for input file of Vp (primary velocity)
    :param vs_file: file name for input file of Vs (secondary velocity)
    """

    print 'domain size: ' + str(domain_size)
    print 'grid size: ' + str(grid_size)
    print 'dt: ' + str(dt)
    print 'tmax: ' + str(tmax)

    # declare fields
    Txx = SField('Txx')
    Txx.set(dimension=3, direction=(1, 1))
    Tyy = SField('Tyy')
    Tyy.set(dimension=3, direction=(2, 2))
    Tzz = SField('Tzz')
    Tzz.set(dimension=3, direction=(3, 3))
    Txy = SField('Txy')
    Txy.set(dimension=3, direction=(1, 2))
    Tyz = SField('Tyz')
    Tyz.set(dimension=3, direction=(2, 3))
    Txz = SField('Txz')
    Txz.set(dimension=3, direction=(1, 3))
    U = VField('U')
    U.set(dimension=3, direction=1)
    V = VField('V')
    V.set(dimension=3, direction=2)
    W = VField('W')
    W.set(dimension=3, direction=3)

    grid = StaggeredGrid(dimension=3)

    # set switches
    grid.set_omp(omp)
    grid.set_simd(simd)
    grid.set_ivdep(ivdep)
    grid.set_io(io)
    grid.set_double(double)
    grid.set_expand(expand)
    grid.set_eval_const(eval_const)

    grid.set_stress_fields([Txx, Tyy, Tzz, Txy, Tyz, Txz])
    grid.set_velocity_fields([U, V, W])
    grid.set_domain_size(domain_size)
    grid.set_grid_size(grid_size)
    grid.set_time_step(dt, tmax)

    # define parameters
    rho, beta, lam, mu = symbols('rho beta lambda mu')
    t, x, y, z = symbols('t x y z')
    grid.set_index([x, y, z])

    if read:
        grid.set_media_params(read=True, rho_file=rho_file,
                              vp_file=vp_file, vs_file=vs_file)
    else:
        grid.set_media_params(read=False, rho=1.0, vp=1.0, vs=0.5)

    print 'require dt < ' + str(grid.get_time_step_limit())

    # define eigen waves
    Omega = pi*sqrt(2*mu*beta)
    A = sqrt(2*mu/beta)
    U_func = cos(pi*x)*(sin(pi*y)-sin(pi*z))*cos(Omega*t)
    V_func = cos(pi*y)*(sin(pi*z)-sin(pi*x))*cos(Omega*t)
    W_func = cos(pi*z)*(sin(pi*x)-sin(pi*y))*cos(Omega*t)
    Txx_func = -A*sin(pi*x)*(sin(pi*y)-sin(pi*z))*sin(Omega*t)
    Tyy_func = -A*sin(pi*y)*(sin(pi*z)-sin(pi*x))*sin(Omega*t)
    Tzz_func = -A*sin(pi*z)*(sin(pi*x)-sin(pi*y))*sin(Omega*t)
    Txy_func = Float(0)
    Tyz_func = Float(0)
    Txz_func = Float(0)

    U.set_analytic_solution(U_func)
    V.set_analytic_solution(V_func)
    W.set_analytic_solution(W_func)
    Txx.set_analytic_solution(Txx_func)
    Tyy.set_analytic_solution(Tyy_func)
    Tzz.set_analytic_solution(Tzz_func)
    Txy.set_analytic_solution(Txy_func)
    Tyz.set_analytic_solution(Tyz_func)
    Txz.set_analytic_solution(Txz_func)

    grid.calc_derivatives()

    # PDEs: momentum equations
    eq1 = Eq(U.d[0][1], beta*(Txx.d[1][2] + Txy.d[2][2] + Txz.d[3][2]))
    eq2 = Eq(V.d[0][1], beta*(Txy.d[1][2] + Tyy.d[2][2] + Tyz.d[3][2]))
    eq3 = Eq(W.d[0][1], beta*(Txz.d[1][2] + Tyz.d[2][2] + Tzz.d[3][2]))
    # PDEs: stress-strain equations
    eq4 = Eq(Txx.d[0][1], (lam + 2*mu)*U.d[1][2] + lam*(V.d[2][2]+W.d[3][2]))
    eq5 = Eq(Tyy.d[0][1], (lam + 2*mu)*V.d[2][2] + lam*(U.d[1][2]+W.d[3][2]))
    eq6 = Eq(Tzz.d[0][1], (lam + 2*mu)*W.d[3][2] + lam*(U.d[1][2]+V.d[2][2]))
    eq7 = Eq(Txy.d[0][1], mu*(U.d[2][2] + V.d[1][2]))
    eq8 = Eq(Tyz.d[0][1], mu*(V.d[3][2] + W.d[2][2]))
    eq9 = Eq(Txz.d[0][1], mu*(U.d[3][2] + W.d[1][2]))

    grid.solve_fd([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9])

    Txx_expr = (lam + 2*mu)*U.d[1][1] + lam*(V.d[2][1]+W.d[3][1])
    Txx.set_dt(Txx_expr)
    Tyy_expr = (lam + 2*mu)*V.d[2][1] + lam*(U.d[1][1]+W.d[3][1])
    Tyy.set_dt(Tyy_expr)
    Tzz_expr = (lam + 2*mu)*W.d[3][1] + lam*(U.d[1][1]+V.d[2][1])
    Tzz.set_dt(Tzz_expr)
    Txy_expr = mu*(U.d[2][1] + V.d[1][1])
    Txy.set_dt(Txy_expr)
    Tyz_expr = mu*(V.d[3][1] + W.d[2][1])
    Tyz.set_dt(Tyz_expr)
    Txz_expr = mu*(U.d[3][1] + W.d[1][1])
    Txz.set_dt(Txz_expr)

    grid.set_free_surface_boundary(dimension=1, side=0)
    grid.set_free_surface_boundary(dimension=1, side=1)
    grid.set_free_surface_boundary(dimension=2, side=0)
    grid.set_free_surface_boundary(dimension=2, side=1)
    grid.set_free_surface_boundary(dimension=3, side=0)
    grid.set_free_surface_boundary(dimension=3, side=1)

    # Generate code and write to output file
    grid.generate(filename, o_step, o_converge)


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
             filename=path.join(_test_dir, 'test_grid3d.cpp'))


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
             filename=path.join(_test_dir, 'test_grid3d_vtk.cpp'))


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
             filename=path.join(_test_dir, 'test_grid3d_read.cpp'),
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
             filename=path.join(_test_dir, 'test_grid3d_ivdep.cpp'))
    run_test(domain_size, grid_size, dt, tmax, o_step=False, o_converge=False,
             omp=True, simd=True, ivdep=False, io=False,
             filename=path.join(_test_dir, 'test_grid3d_simd.cpp'))


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
