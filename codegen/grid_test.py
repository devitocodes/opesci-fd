from sympy import symbols, Eq, sqrt, pi, cos, sin, Float
from mako.lookup import TemplateLookup
from mako.runtime import Context
from StringIO import StringIO
from grid import *

def run_test(domain_size, grid_size, dt, tmax, o_step=False, o_converge=True, omp=True, simd=False, ivdep=True, io=False, double=False, filename='test.cpp'):
    """
    create 3D eigen waves and run FD simulation
    :param domain_size: define size of domain, e.g. (1.0, 1.0, 1.0) for unit cube
    :param spacing: define grid spacing, e.g. (0.1, 0.1, 0.1)
    :param dt: define time step
    :param tmax: define simulation time
    :param o_step: switch for code to output at each time step (such as save field as vtk file), default False (no output)
    :param o_converge: switch for code to calculate L2 norms for convergence test, default True (output)
    :param omp: swtich for inserting #pragma omp for before outter loop, default True (use omp)
    :param simd: switch for inserting #pragma simd before inner loop, default False (do not use simd)
    :param ivddp: switch for inserting #praga ivdep before inner loop, default True (use ivdep)
    :param io: switch for include io header files, default False (not include vtk header files)
    :param double: switch for using double as real number variables, default False (use float)
    :param filename: output source code file name as string
    """

    print 'domain size: ' + str(domain_size)
    print 'gird size: ' + str(grid_size)
    print 'dt: ' + str(dt)
    print 'tmax: ' + str(tmax)

    # declare fields
    Txx = SField('Txx'); Txx.set(dimension=3, direction=(1, 1))
    Tyy = SField('Tyy'); Tyy.set(dimension=3, direction=(2, 2))
    Tzz = SField('Tzz'); Tzz.set(dimension=3, direction=(3, 3))
    Txy = SField('Txy'); Txy.set(dimension=3, direction=(1, 2))
    Tyz = SField('Tyz'); Tyz.set(dimension=3, direction=(2, 3))
    Txz = SField('Txz'); Txz.set(dimension=3, direction=(1, 3))
    U = VField('U'); U.set(dimension=3, direction=1)
    V = VField('V'); V.set(dimension=3, direction=2)
    W = VField('W'); W.set(dimension=3, direction=3)

    grid = StaggeredGrid(dimension=3)

    # set switches
    grid.set_omp(omp)
    grid.set_simd(simd)
    grid.set_ivdep(ivdep)
    grid.set_io(io)
    grid.set_double(double)

    grid.set_stress_fields([Txx, Tyy, Tzz, Txy, Tyz, Txz])
    grid.set_velocity_fields([U, V, W])
    grid.set_domain_size(domain_size)
    grid.set_grid_size(grid_size)
    grid.set_time_step(dt, tmax)

    # define parameters
    rho, beta, lam, mu = symbols('rho beta lambda mu')
    t, x, y, z = symbols('t x y z')
    grid.set_index([x, y, z])

    grid.set_media_params(read=False, rho=1.0, vp=1.0, vs=0.5, rho_file='', vp_file='', vs_file='')

    print 'require dt < ' + str(grid.get_time_step_limit())

    # define eigen waves
    Omega = pi*sqrt(2*mu/rho)
    A = sqrt(2*rho*mu)
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

    Txx_expr = (lam + 2*mu)*U.d[1][1] + lam*(V.d[2][1]+W.d[3][1]); Txx.set_dt(Txx_expr)
    Tyy_expr = (lam + 2*mu)*V.d[2][1] + lam*(U.d[1][1]+W.d[3][1]); Tyy.set_dt(Tyy_expr)
    Tzz_expr = (lam + 2*mu)*W.d[3][1] + lam*(U.d[1][1]+V.d[2][1]); Tzz.set_dt(Tzz_expr)
    Txy_expr = mu*(U.d[2][1] + V.d[1][1]); Txy.set_dt(Txy_expr)
    Tyz_expr = mu*(V.d[3][1] + W.d[2][1]); Tyz.set_dt(Tyz_expr)
    Txz_expr = mu*(U.d[3][1] + W.d[1][1]); Txz.set_dt(Txz_expr)

    grid.set_free_surface_boundary(dimension=1, side=0); grid.set_free_surface_boundary(dimension=1, side=1)
    grid.set_free_surface_boundary(dimension=2, side=0); grid.set_free_surface_boundary(dimension=2, side=1)
    grid.set_free_surface_boundary(dimension=3, side=0); grid.set_free_surface_boundary(dimension=3, side=1)

    if o_step:
        output_step = grid.output_step()
    else:
        output_step = ''

    if o_converge:
        output_final = grid.converge_test()
    else:
        output_final = ''

    # write to template file
    mylookup = TemplateLookup(directories=['templates/staggered', 'templates/'])
    mytemplate = mylookup.get_template('staggered3d_tmpl.cpp')
    buf = StringIO()
    dict1 = {'io': io, 'time_stepping': grid.time_stepping(),
             'define_constants': grid.define_variables(),
             'declare_fields': grid.declare_fields(),
             'initialise': grid.initialise(),
             'initialise_bc': grid.initialise_boundary(),
             'stress_loop': grid.stress_loop(),
             'velocity_loop': grid.velocity_loop(),
             'stress_bc': grid.stress_bc(),
             'velocity_bc': grid.velocity_bc(),
             'output_step': output_step,
             'output_final': output_final}
    ctx = Context(buf, **dict1)
    mytemplate.render_context(ctx)
    code = buf.getvalue()

    # generate compilable C++ source code
    f = open(filename, 'w')
    f.write(code)
    f.close()

    print filename + ' generated!'
