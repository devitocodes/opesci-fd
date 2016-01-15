from opesci import *
from os import path
from argparse import ArgumentParser, RawTextHelpFormatter


_test_dir = path.join(path.dirname(__file__), "src")


def eigenwave3d(domain_size, grid_size, dt, tmax, output_vts=False, o_converge=True,
                accuracy_order=[1, 2, 2, 2],
                omp=True, simd=False, ivdep=True, double=False, pluto=False,
                filename='test.cpp', read=False, expand=True, eval_const=True,
                rho_file='', vp_file='', vs_file='', fission=False):
    """
    create 3D eigen waves and run FD simulation
    :param domain_size: define size of domain
    e.g. (1.0, 1.0, 1.0) for unit cube
    :param grid_size: define grid size, e.g. (10, 10, 10)
    :param dt: define time step
    :param tmax: define simulation time
    :param output_vts: switch for code to output at each time step
    (such as save field as vtk file), default False (no output)
    :param o_converge: switch for code to calculate L2 norms
    for convergence test, default True (output)
    :param omp: swtich for inserting #pragma omp for before outter loop
    default True (use omp)
    :param simd: switch for inserting #pragma simd before inner loop
    default False (do not use simd)
    :param ivddp: switch for inserting #praga ivdep before inner loop
    default True (use ivdep)
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
    :param pluto: switch for inserting #pragma scop and #pragma endscop for
    pluto optimisation
    :param fission: switch for doing loop fission optimisation
    """

    print 'domain size: ' + str(domain_size)
    print 'grid size: ' + str(grid_size)
    print 'approximation order: ' + str(accuracy_order)
    print 'dt: ' + str(dt)
    print 'tmax: ' + str(tmax)

    # Declare fields
    Txx = SField('Txx', dimension=3, direction=(1, 1))
    Tyy = SField('Tyy', dimension=3, direction=(2, 2))
    Tzz = SField('Tzz', dimension=3, direction=(3, 3))
    Txy = SField('Txy', dimension=3, direction=(1, 2))
    Tyz = SField('Tyz', dimension=3, direction=(2, 3))
    Txz = SField('Txz', dimension=3, direction=(1, 3))
    U = VField('U', dimension=3, direction=1)
    V = VField('V', dimension=3, direction=2)
    W = VField('W', dimension=3, direction=3)

    grid = StaggeredGrid(dimension=3, domain_size=domain_size,
                         grid_size=grid_size,
                         stress_fields=[Txx, Tyy, Tzz, Txy, Tyz, Txz],
                         velocity_fields=[U, V, W], pluto=pluto, fission=fission)
    grid.set_time_step(dt, tmax)

    grid.set_switches(omp=omp, simd=simd, ivdep=ivdep, double=double,
                      expand=expand, eval_const=eval_const,
                      output_vts=output_vts, converge=o_converge)

    # define parameters
    rho, beta, lam, mu = symbols('rho beta lambda mu')
    t, x, y, z = symbols('_t x y z')
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

    grid.set_order(accuracy_order)
    grid.calc_derivatives()

    # PDEs: momentum equations
    eq1 = Eq(U.d[0][1], beta*(Txx.d[1][1] + Txy.d[2][1] + Txz.d[3][1]))
    eq2 = Eq(V.d[0][1], beta*(Txy.d[1][1] + Tyy.d[2][1] + Tyz.d[3][1]))
    eq3 = Eq(W.d[0][1], beta*(Txz.d[1][1] + Tyz.d[2][1] + Tzz.d[3][1]))
    # PDEs: stress-strain equations
    eq4 = Eq(Txx.d[0][1], (lam + 2*mu)*U.d[1][1] + lam*(V.d[2][1]+W.d[3][1]))
    eq5 = Eq(Tyy.d[0][1], (lam + 2*mu)*V.d[2][1] + lam*(U.d[1][1]+W.d[3][1]))
    eq6 = Eq(Tzz.d[0][1], (lam + 2*mu)*W.d[3][1] + lam*(U.d[1][1]+V.d[2][1]))
    eq7 = Eq(Txy.d[0][1], mu*(U.d[2][1] + V.d[1][1]))
    eq8 = Eq(Tyz.d[0][1], mu*(V.d[3][1] + W.d[2][1]))
    eq9 = Eq(Txz.d[0][1], mu*(U.d[3][1] + W.d[1][1]))

    grid.solve_fd([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9])

    grid.set_free_surface_boundary(dimension=1, side=0)
    grid.set_free_surface_boundary(dimension=1, side=1)
    grid.set_free_surface_boundary(dimension=2, side=0)
    grid.set_free_surface_boundary(dimension=2, side=1)
    grid.set_free_surface_boundary(dimension=3, side=0)
    grid.set_free_surface_boundary(dimension=3, side=1)

    # print 'kernel Weighted AI: ' + '%.2f' % grid.get_overall_kernel_ai()[1]
    print 'stress kernel AI'
    print '%.2f, %.2f (weighted), %d ADD, %d MUL, %d LOAD, %d STORE' % grid.get_stress_kernel_ai()
    print 'velocity kernel AI'
    print '%.2f, %.2f (weighted), %d ADD, %d MUL, %d LOAD, %d STORE' % grid.get_velocity_kernel_ai()
    # print 'stress bc AI'
    # print grid.get_stress_bc_ai()
    # print 'velocity bc AI'
    # print grid.get_velocity_bc_ai()
    print 'overall algorithm AI'
    print '%.2f, %.2f (weighted)' % grid.get_overall_kernel_ai()

    return grid


def default(compiler=None, execute=False, nthreads=1,
            accuracy_order=[2, 4, 4, 4],
            output=False, profiling=False, papi_events=[], pluto=False, tile=' ', fission=False):
    """
    Eigenwave test case on a unit cube grid (100 x 100 x 100)
    """
    if pluto:
        f = open('tile.sizes', 'w')
        f.write(str(tile))
        f.close()
    domain_size = (1.0, 1.0, 1.0)
    grid_size = (100, 100, 100)
    dt = 0.002
    tmax = 1.0
    filename = path.join(_test_dir, 'eigenwave3d.cpp')
    grid = eigenwave3d(domain_size, grid_size, dt, tmax,
                       accuracy_order=accuracy_order,
                       o_converge=True, omp=True, simd=False,
                       ivdep=True, filename=filename, pluto=pluto, fission=fission)
    grid.set_switches(output_vts=output, profiling=profiling)
    grid.set_papi_events(papi_events)
    out = None
    if pluto:
        grid.generate(filename)
        filename_p = grid.pluto_op(filename)
        if compiler in ['clang', 'clang++']:
            # an ugly fix, but pluto will always attack <omp.h> to the first line
            # which would fail clang
            with open(filename_p, 'r') as fin:
                data = fin.read().splitlines(True)
            with open(filename_p, 'w') as fout:
                fout.writelines(data[1:])
        grid.src_file = filename_p
        if compiler is None:
            share = True
        else:
            share = False
            out = grid.compile(filename_p, compiler=compiler, shared=share)
    else:
        filename_p = filename
        if compiler is None:
            grid.generate(filename_p)
        else:
            out = grid.compile(filename_p, compiler=compiler, shared=False)

    if execute:
        # Test Python-based execution for the base test
        grid.execute(filename_p, compiler=compiler, nthreads=nthreads)
        grid.convergence()
    return out


def read_data(compiler=None, execute=False, nthreads=1,
              accuracy_order=[2, 4, 4, 4],
              output=False, profiling=False, papi_events=[], fission=False):

    """Test for model intialisation from input file

    Computes eigenwave on a unit cube grid (200 x 200 x 200)
    """
    domain_size = (1.0, 1.0, 1.0)
    grid_size = (195, 195, 195)
    dt = 0.002
    tmax = 1.0
    filename = path.join(_test_dir, 'eigenwave3d_read.cpp')
    grid = eigenwave3d(domain_size, grid_size, dt, tmax, read=True,
                       accuracy_order=accuracy_order,
                       o_converge=False, omp=True, simd=False, ivdep=True,
                       filename=filename, rho_file='RHOhomogx200',
                       vp_file='VPhomogx200', vs_file='VShomogx200', fission=fission)
    grid.set_switches(output_vts=output, profiling=profiling)
    grid.set_papi_events(papi_events)
    if compiler is None:
        grid.generate(filename)
    else:
        grid.compile(filename, compiler=compiler, shared=False)
    if execute:
        # Test Python-based execution for the base test
        grid.execute(filename, compiler=compiler, nthreads=nthreads)
        grid.convergence()


def cx1():
    """
    test case for comparison between pragma simd and pragma ivdep on cx1
    """
    domain_size = (1.0, 1.0, 1.0)
    grid_size = (200, 200, 200)
    dt = 0.001
    tmax = 5.0
    eigenwave3d(domain_size, grid_size, dt, tmax, output_vts=False, o_converge=False,
                omp=True, simd=False, ivdep=True,
                filename=path.join(_test_dir, 'eigenwave3d_ivdep.cpp'))
    eigenwave3d(domain_size, grid_size, dt, tmax, output_vts=False, o_converge=False,
                omp=True, simd=True, ivdep=False,
                filename=path.join(_test_dir, 'eigenwave3d_simd.cpp'))


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
    eigenwave3d(domain_size, (s, s, s), dt, tmax, output_vts=False, o_converge=True,
                omp=True, simd=False, ivdep=True,
                filename='tmp/test3d_'+str(s)+'.cpp')

    s = s*2
    dt = c/(s**2)
    eigenwave3d(domain_size, (s, s, s), dt, tmax, output_vts=False, o_converge=True,
                omp=True, simd=False, ivdep=True,
                filename='tmp/test3d_'+str(s)+'.cpp')

    s = s*2
    dt = c/(s**2)
    eigenwave3d(domain_size, (s, s, s), dt, tmax, output_vts=False, o_converge=True,
                omp=True, simd=False, ivdep=True,
                filename='tmp/test3d_'+str(s)+'.cpp')

    s = s*2
    dt = c/(s**2)
    eigenwave3d(domain_size, (s, s, s), dt, tmax, output_vts=False, o_converge=True,
                omp=True, simd=False, ivdep=True,
                filename='tmp/test3d_'+str(s)+'.cpp')


def main():
    ModeHelp = """Avalable testing modes:
default:   Eigenwave test case on a unit cube grid (100 x 100 x 100)

read:      Test for model intialisation from input file; computes
           eigenwave on a unit cube grid (200 x 200 x 200)

converge:  Convergence test of the (2,4) scheme, which is 2nd order
           in time and 4th order in space. The test halves spacing
           starting from 0.1 and reduces dt by a factor of 4 for
           each step
"""
    p = ArgumentParser(description="Standalone testing script for the Eigenwave3D example",
                       formatter_class=RawTextHelpFormatter)
    p.add_argument('mode', choices=('default', 'read', 'converge', 'cx1'),
                   nargs='?', default='default', help=ModeHelp)
    p.add_argument('-so', '--spatial_order', default=4, type=int, dest='so',
                   help='order of the spatial discretisation to use for code generation * 2, eg. order=4 to use 4th order approximation in x,y,z')
    p.add_argument('-c', '--compiler', default=None,
                   help='C++ Compiler to use for model compilation, eg. g++ or icpc')
    p.add_argument('-x', '--execute', action='store_true', default=False,
                   help='Dynamically execute the generated model')
    p.add_argument('-n', '--nthreads', type=int, default=1,
                   help='Number of threads for dynamic execution')
    p.add_argument('-o', '--output', action='store_true', default=False,
                   help='Activate solution output in .vts format')
    p.add_argument('-p', '--profiling', action='store_true', default=False,
                   help='Activate performance profiling from PAPI')
    p.add_argument('--papi-events', dest='papi_events', nargs='+', default=[],
                   help='Specific PAPI events to measure')
    p.add_argument('--tile', default=None,
                   help="tile-size for pluto optimisation e.g. --tile '4 4 32'")
    p.add_argument('--pluto', action='store_true', default=False,
                   help="Apply pluto optimisation ")
    p.add_argument('--fission', action='store_true', default=False,
                   help="Apply loop fission optimisation")

    args = p.parse_args()
    print "Eigenwave3D example (mode=%s)" % args.mode

    if args.mode == 'default':
        default(compiler=args.compiler, execute=args.execute,
                nthreads=args.nthreads, output=args.output,
                accuracy_order=[2, args.so, args.so, args.so],
                profiling=args.profiling, papi_events=args.papi_events,
                pluto=args.pluto, tile=args.tile, fission=args.fission)
    elif args.mode == 'read':
        read_data(compiler=args.compiler, execute=args.execute,
                  nthreads=args.nthreads, output=args.output,
                  accuracy_order=[2, args.so, args.so, args.so],
                  profiling=args.profiling, papi_events=args.papi_events, fission=args.fission)
    elif args.mode == 'converge':
        converge_test()
    elif args.mode == 'cx1':
        cx1()

if __name__ == "__main__":
    main()
