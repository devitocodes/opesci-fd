from sympy import Symbol, factorial, Matrix, Rational, Indexed, IndexedBase
from sympy import solve, Eq, expand
from sympy.mpmath.libmp import prec_to_dps
from sympy.printing.ccode import CCodePrinter

import sympy.mpmath.libmp as mlib
from mako.lookup import TemplateLookup
from mako.runtime import Context
from StringIO import StringIO


hf = Rational(1, 2)  # 1/2


class MyCPrinter(CCodePrinter):
    """
    MyCPrinter extends sympy.printing.ccode.CCodePrinter
    """
    def __init__(self, settings={}):
        CCodePrinter.__init__(self, settings)

    def _print_Indexed(self, expr):
        """
        override method in CCodePrinter
        Print field as C style multidimensional array
        e.g. U[t,x,y,z] -> U[t][x][y][z]
        """
        output = self._print(expr.base.label) \
            + ''.join(['[' + self._print(x) + ']' for x in expr.indices])
        return output

    def _print_Rational(self, expr):
        """
        override method in CCodePrinter
        print fractional number as float/float
        (default was long double/long double)
        """
        p, q = int(expr.p), int(expr.q)
        return '%d.0F/%d.0F' % (p, q)  # float precision by default

    def _print_Mod(self, expr):
        """
        override method in CCodePrinter
        print mod using % operator in C++
        """
        args = map(ccode, expr.args)
        args = ['('+x+')' for x in args]
        result = '%'.join(args)
        return result

    def _print_Float(self, expr):
        """
        override method in StrPrinter
        always printing floating point numbers in scientific notation
        """
        prec = expr._prec
        if prec < 5:
            dps = 0
        else:
            dps = prec_to_dps(expr._prec)
        if self._settings["full_prec"] is True:
            strip = False
        elif self._settings["full_prec"] is False:
            strip = True
        elif self._settings["full_prec"] == "auto":
            strip = self._print_level > 1
        rv = mlib.to_str(expr._mpf_, dps, strip_zeros=strip, max_fixed=-2, min_fixed=2)
        if rv.startswith('-.0'):
            rv = '-0.' + rv[3:]
        elif rv.startswith('.0'):
            rv = '0.' + rv[2:]
        return rv


def ccode(expr, **settings):
    """
    generate C++ code from expression expr
    calling MyCPrinter class
    """
    return MyCPrinter(settings).doprint(expr, None)


def ccode_eq(eq, **settings):
    """
    genereate C++ assignment from equation eq
    assigning RHS to LHS
    """
    return MyCPrinter(settings).doprint(eq.lhs, None) \
        + ' = ' + MyCPrinter(settings).doprint(eq.rhs, None)


def render(tmpl, dict1):
    """
    render Mako template using dictionary dict1
    returns the result as string
    """
    buf = StringIO()
    ctx = Context(buf, **dict1)
    tmpl.render_context(ctx)
    return buf.getvalue()


def IndexedBases(s):
    """
    declare multiple IndexedBase objects
    :param s: string of names seperated by white space
    returns IndxedBase objects as tuple
    """
    l = s.split()
    bases = [IndexedBase(x) for x in l]
    return tuple(bases)


def tc(dx, n):
    """
    return coefficient of power n term in Taylor series expansion
    :param n: power
    :param dx: distance from expansion reference
    """
    return (dx**n)/factorial(n)


def Taylor(dx, n):
    """
    create Matrix of Taylor Coefficients M, such that M * D = R
    where D is list of derivatives at x [f, f', f'' ..]
    R is list of values at neighbour grid point [.. f(x-dx), f(x), f(x+dx) ..]
    :param dx: spacing between grid points
    :param n: length of D and R (i.e. number of derivatives to calculate)
    returns Maxtrix object M
    """
    l = []
    for i in range(-n+1, n):
        ll = [tc((i*dx), j) for j in range(2*n-1)]
        l.append(ll)
    return Matrix(l)


def Taylor_half(dx, n):
    """
    similar function as Taylor() for staggered grids
    create Matrix of Taylor Coefficients M, such that M * D = R
    where D is list of derivatives at x [f, f', f'' ..]
    R is list of values at neighbour grid point
    e.g. R = [.. f(x-3/2*dx), f(x-dx/2), f(x), f(x+dx/2), f(x+3/2*dx) ..]
    :param dx: spacing between grid points
    :param n: length of D and R (i.e. number of derivatives to calculate)
    returns Maxtrix object M
    """
    l = []
    for i in range(-n*2+1, n*2, 2):
        ll = [tc((i*dx/2), j) for j in range(2*n)]
        l.append(ll)
    return Matrix(l)


def Deriv(U, i, k, d, n):
    """
    calculate the FD approximation for nth derivative
    the function works by using M * D = R, as specified in Taylor()
    then inverting M, so that D = M_inverse * R
    :param U: the field
    :param i: list of indices (Symbol objects) of field U
    e.g. possibley [t,x,y,z] for 3D
    :param k: determine the derivative is with respect to which dimension
    e.g. dU/dt if k=0, dU/dx if k=0, assuming i=[t,x,y,z]
    :param d: spacing of the grid, e.g. dx, dt
    :param n: accuracy of approximation
    e.g. n=2 for 2nd-order FD approximation
    returns D=list of expressions for FD approxiation
    D[1] for first derivative, D[2] for 2nd derivative etc
    raises NotImplementedError exception if dimension is more than 4
    (i.e. 1 time and 3 space dimensions)
    """
    M = Taylor(d, n)
    s = [0]*len(i)
    s[k] = 1  # switch on kth dimension
    # generate matrix of RHS, i.e. [ ... U[x-1], U[x], U[x+1] ... ]
    if len(i) == 1:
        RX = Matrix([U[i[0]+s[0]*x] for x in range(-n+1, n)])
    elif len(i) == 2:
        RX = Matrix([U[i[0]+s[0]*x, i[1]+s[1]*x] for x in range(-n+1, n)])
    elif len(i) == 3:
        RX = Matrix([U[i[0]+s[0]*x,
                    i[1]+s[1]*x,
                    i[2]+s[2]*x] for x in range(-n+1, n)])
    elif len(i) == 4:
        RX = Matrix([U[i[0]+s[0]*x,
                    i[1]+s[1]*x,
                    i[2]+s[2]*x,
                    i[3]+s[3]*x] for x in range(-n+1, n)])
    else:
        raise NotImplementedError(">4 dimensions, need to fix")

    return M.inv() * RX


def Deriv_half(U, i, k, d, n):
    """
    similar function as Deriv() for staggered grids
    calculate the FD approximation for nth derivative
    the function works by using M * D = R, as specified in Taylor()
    then inverting M, so that D = M_inverse * R
    :param U: the field
    :param i: list of indices (Symbol objects) of field U
    e.g. possibley [t,x,y,z] for 3D
    :param k: determine the derivative is with respect to which dimension
    e.g. dU/dt if k=0, dU/dx if k=0, assuming i=[t,x,y,z]
    :param d: spacing of the grid, e.g. dx, dt
    :param n: accuracy of approximation
    e.g. n=1 for 2nd-order FD approximation, n=2 for 4th-order FD approximation
    returns D=list of expressions for FD approxiation
    D[1] for first derivative, D[2] for 2nd derivative etc
    raises NotImplementedError exception if dimension is more than 4
    (1 time and 3 space dimensions)
    """
    M = Taylor_half(d, n)
    s = [0]*len(i)
    s[k] = 1  # switch on kth dimension
    # generate matrix of RHS, i.e. [ ... U[x-1], U[x], U[x+1] ... ]
    if len(i) == 1:
        RX = Matrix([U[i[0]+s[0]*x*hf] for x in range(-n*2+1, n*2, 2)])
    elif len(i) == 2:
        RX = Matrix([U[i[0]+s[0]*x*hf,
                    i[1]+s[1]*x*hf] for x in range(-n*2+1, n*2, 2)])
    elif len(i) == 3:
        RX = Matrix([U[i[0]+s[0]*x*hf,
                    i[1]+s[1]*x*hf,
                    i[2]+s[2]*x*hf] for x in range(-n*2+1, n*2, 2)])
    elif len(i) == 4:
        RX = Matrix([U[i[0]+s[0]*x*hf,
                    i[1]+s[1]*x*hf,
                    i[2]+s[2]*x*hf,
                    i[3]+s[3]*x*hf] for x in range(-n*2+1, n*2, 2)])
    else:
        raise NotImplementedError(">4 dimensions, need to fix")

    result = M.inv() * RX
    return result


def is_half(expr):
    """
    test if constants (non symbols) in an expression is integer+1/2 or integer
    (indices on staggered grid are either integer or integer+1/2)
    e.g. returns True for expression x+1/2, returns False for expression x+y+1
    used to test for index in staggered grid
    """
    d = {x: 0 for x in expr.free_symbols}
    return not expr.subs(d).is_Integer


def shift_grid(expr):
    """
    shift all field indices in the input expression to whole number
    if it is staggered (i.e. has index that is in the middle of grids)
    the result expression can then be converted to C code
    """
    if expr.is_Symbol:
        return expr
    if expr.is_Number:
        return expr
    if isinstance(expr, Indexed):
        b = expr.base
        idx = [x-hf if is_half(x) else x for x in list(expr.indices)]
        t = Indexed(b, *idx)
        return t
    args = tuple([shift_grid(arg) for arg in expr.args])  # recursive call
    expr2 = expr.func(*args)
    return expr2


def shift_index(expr, k, s):
    """
    shift the k-th index of all fields in input expression by s
    return the shifted expression
    :param expr: input expression
    k: the index number to be shifted
    s: the shifted amount
    e.g. k=1, s=1, U[x,y,z] -> U[x,y+1,z]
    """
    if expr.is_Symbol:
        return expr
    if expr.is_Number:
        return expr
    if isinstance(expr, Indexed):
        b = expr.base
        idx = list(expr.indices)
        idx[k] += s
        t = Indexed(b, *idx)
        return t
    # recursive call
    args = tuple([shift_index(arg, k, s) for arg in expr.args])
    expr2 = expr.func(*args)
    return expr2


class Field(IndexedBase):
    """
    - Class to represent fields on staggered grid
    - Extends sympy IndexedBase class, can be indexed with [] operator
    - Parent class to VFeild (velocity) and SField (stress)
    - Holds relevant information such as dimension, staggered-ness,
    expression for FD approximation of derivatives,
    code to calculate boundary cells
    """

    def set(self, dimension, staggered):
        self.dimension = dimension
        self.staggered = staggered
        self.lookup = TemplateLookup(directories=['templates/staggered/'])
        # list of list to store derivative expressions
        self.d = [[None]*4 for x in range(dimension+1)]
        # list of list to store boundary ghost cell code
        self.bc = [[None]*2 for x in range(dimension+1)]

    def set_analytic_solution(self, function):
        """
        set analytical function (exact solution) of the field
        used to compare against numerical solution
        param function: expression of the analytical function
        """
        self.sol = function

    def calc_derivative(self, l, k, d, n):
        """
        assign list self.d with FD approximations of 1st derivatives
        such that self.d[k][n] = FD approximation of 1st derivative,
        in kth index, of nth order accuracy
        input param description same as Deriv_half()
        """
        self.d[k][n] = Deriv_half(self, l, k, d, n)[1]

    def align(self, expr):
        """
        - shift the indices of fields in input expression
        according to the staggered-ness of this field
        - used to convert relative offset reference between fields
        to absolute reference (prepare to be converted to array)
        - return the modified expression
        """
        if expr.is_Symbol or expr.is_Number:
            return expr
        if isinstance(expr, Indexed):
            b = expr.base
            if not (isinstance(b, VField) or isinstance(b, SField)):
                return expr
            # align indices if input field staggered different from this field
            idx = []
            for k in range(len(expr.indices)):
                # if both are staggered or unstaggered in direction k
                # index is unchanged
                if self.staggered[k] == b.staggered[k]:
                    idx += [expr.indices[k]]
                # if this field is staggered but target field is unstaggered
                # in direction k, shift by +1/2
                elif self.staggered[k]:
                    idx += [expr.indices[k]+hf]
                # if this field is unstaggered but target field is staggered
                # in direction k, shift by -1/2
                else:
                    idx += [expr.indices[k]-hf]
            tmp = Indexed(b, *idx)
            return tmp
        # recursive call for all arguments of expr
        args = tuple([self.align(arg) for arg in expr.args])
        result = expr.func(*args)
        return result

    def set_fd_kernel(self, fd):
        """
        set the updating kernel of this field
        e.g. the expression to calculate U[t+1,x,y,z]
        store the kernel to self.fd
        store the aligned expression to self.fd_align
        """
        self.fd = fd
        tmp = self.align(fd)
        self.fd_align = tmp

    def vtk_save_field(self):
        """
        generate code to output this field with vtk
        uses Mako template
        returns the generated code as string
        """
        tmpl = self.lookup.get_template('save_field.txt')
        result = ''
        dict1 = {'filename': ccode(self.label)+'_', 'field': ccode(self.label)}
        result = render(tmpl, dict1)

        return result

    def set_dt(self, dt):
        """
        set the expression of first time derivative of the field, e.g. dU/dt
        used to calculate ghost cells for free-surface boundary condition
        """
        self.dt = dt

    def set_media_param(self, media_param):
        """
        - set the media parameters beta, lambda, mu
        - these parameters should be used for the updating kernel
        - only needed when read data from file
        """
        self.media_param = media_param


class VField(Field):
    """
    Class to represent velocity field on staggered grid
    subclass of Field
    """

    def set(self, dimension, direction):
        """
        - set number of dimensions and direction of the velocity field
        - work out the staggered-ness according to the direction
        - a velocity field is only staggered in the spatial index
        same as its direction
        - a velocity field is always staggered in time index
        i.e. in 3D, field Vx will have staggered = [True, True, False, False]
        :param dimension: number of dimensions, e.g. 3
        :param direction: the direction of the field, e.g. 1
        """
        self.direction = direction
        staggered = [False] * (dimension+1)
        staggered[0] = True
        staggered[direction] = True
        Field.set(self, dimension, staggered)

    def associate_stress_fields(self, sfields):
        """
        link this velocity field to a list of stress field
        e.g. Vx will be associated with Txx, Txy, Txz
        :param sfields: list of associated stress field
        """
        self.sfields = sfields

    def set_free_surface(self, indices, d, b, side, read=False):
        """
        - set free surface boundary condition to boundary d, at index b
        :param indices: list of indices, e.g. [t,x,y,z] for 3D
        :param d: direction of the boundary surface normal
        :param b: location of the boundary (index)
        :param side: lower boundary (0) or upper boundary (1)
        - e.g. set_free_surface([t,x,y,z],1,2,0)
        set y-z plane at x=2 to be lower free surface
        - ghost cells are calculated using reflection of stress fields
        - store the code to populate ghost cells to self.bc
        """
        # use this stress field to solve for ghost cell expression
        field = self.sfields[d]
        idx = list(indices)
        if not field.staggered[d]:
            # if not staggered, solve ghost cell using T'[b]=0
            eq = Eq(field.dt)
            shift = hf
            t = b - hf
        else:
            # if staggered, solve ghost cell using T'[b-1/2]=T[b+1/2]
            eq = Eq(field.dt.subs(indices[d], indices[d]-hf),
                    field.dt.subs(indices[d], indices[d]+hf))
            shift = 1
            t = b

        idx[d] -= ((-1)**side)*shift
        lhs = self[idx]
        rhs = solve(eq, lhs)[0]
        lhs = lhs.subs(indices[d], t)
        rhs = self.align(rhs.subs(indices[d], t))
        # if read data from file, replace media parameters with array
        # replace particular index with boundary
        if read:
            rhs = rhs.subs(self.media_param)
            rhs = rhs.subs(indices[d], b)

        self.bc[d][side] = ccode(lhs) + ' = ' + ccode(rhs) + ';\n'


class SField(Field):
    """
    Class to represent stress fields on staggered grid
    subclass of Field
    """

    def set(self, dimension, direction):
        """
        - set number of dimensions and direction of the stress field
        - work out the staggered-ness according to the direction
        - compression stress fields are not staggered
        - sheer stress fields are staggered in the surface
        normal and force direction
        - stress fields are not staggered in time index
        i.e. in 3D, field Txx has staggered = [False, False, False, False]
        Txy has staggered = [False, True, True, False]
        :param dimension: number of dimensions, e.g. 3
        :param direction: the direction of the field, e.g. (1,1) for Txx
        """
        self.direction = direction
        staggered = [False] * (dimension+1)
        if direction[0] == direction[1]:
            # compression stress, not staggered
            Field.set(self, dimension, staggered)
        else:
            # sheer stress, staggered
            for i in range(len(direction)):
                staggered[direction[i]] = True
            Field.set(self, dimension, staggered)

    def set_free_surface(self, indices, d, b, side, read=False):
        """
        set free surface boundary condition to boundary d, at index b
        :param indices: list of indices, e.g. [t,x,y,z] for 3D
        :param d: direction of the boundary surface normal
        :param b: location of the boundary (index)
        side: lower boundary (0) or upper boundary (1)
        e.g. set_free_surface([t,x,y,z],1,2,0)
        set y-z plane at x=2 to be lower free surface
        ghost cells are calculated using reflection of stress fields
        store the code to populate ghost cells to self.bc
        """
        if d not in self.direction:
            # e.g. x boundary for field Tyy is not required
            self.bc[d][side] = ''
            return

        idx = list(indices)  # ghost cell
        idx2 = list(indices)  # cell inside domain

        if not self.staggered[d]:
            # if not staggered, assign T[d]=0, assign T[d-1]=-T[d+1]
            idx[d] = b
            idx2[d] = b
            eq1 = Eq(self[idx])
        else:
            # if staggered, assign T[d-1/2]=T[d+1/2], assign T[d-3/2]=T[d+3/2]
            idx[d] = b - (1-side)
            idx2[d] = idx[d] + (-1)**side
            eq1 = Eq(self[idx], -self[idx2])

        idx[d] -= (-1)**side
        idx2[d] += (-1)**side
        eq2 = Eq(self[idx], -self[idx2])
        self.bc[d][side] = ccode_eq(eq1) + ';\n' + ccode_eq(eq2) + ';\n'


class Variable(Symbol):
    """
    wrapper class of Symbol to store extra information
    """
    def __new__(cls, name, *args):
        return Symbol.__new__(cls, name)

    def __init__(self, name, value=0, type='int', constant=False):
        """
        :param type: type string of the variable to be used in generated code
        :param value: associated value of the variable
        :param constant: if the variable is a constant
        """
        self.type = type
        self.constant = constant
        self.value = value


class StaggeredGrid:
    """
    - Class to represent velocity-stress method on staggered grid
    - calculates the computation kernel
    - generate compilable C++ code

    steps (refer to run_test() in grid_test.py):
    - create grid, set switches and :parameters
    (e.g. time step, time length, grid spacing, order of accuracy)
    - create velocity fields (VField objects) and stress fields
    (SField objects), link fields to grid
    - assign analytical solution to the fields (initial condition)
    - populate FD derivative approximations with calc_derivatives()
    - input the VS PDEs using FD derivative approximations
    - solve_fd()
    - set boundary conditions
    - use functions such as grid.stress_loop() to generate code fragments
    - insert the fragments into Mako templates

    switches can be set with set_switch_name methods:
    omp: insert #pragma omp for before outer loops, default True
    ivdep: insert #praga ivdep before inner loop, default True
    simd: insert #pragma simd before inner loop, default False
    double: use float (False) or double (True) for real numbers, default False
    io: include header files for io (e.g. vtk support), default False
    read: whether to read media parameters from input file, default False
    expand: expand kernel fully (no factorisation), default True
    eval_const: evaluate all constants in kernel in generated code default True
    """

    def __init__(self, dimension):
        self.dimension = dimension
        self.lookup = TemplateLookup(directories=['templates/staggered/'])

        # switches
        self.omp = True  # switch for inserting #pragma omp for
        self.ivdep = True  # switch for inserting #pragma ivdep to inner loop
        self.simd = False  # switch for inserting #pragma simd to inner loop
        self.double = False  # use float (False) or double (True)
        self.io = False  # include io header files
        self.expand = True  # expand kernel fully (no factorisation)
        self.eval_const = True  # evaluate all constants in kernel
        self.real_t = 'double' if self.double else 'float'

        # number of ghost cells for boundary
        self.margin = Variable('margin', 2, 'int', True)
        self.size = [1.0] * dimension  # default domain size
        # grid size symbols, dim1, dim2, ...
        self.dim = [Variable('dim'+str(k+1),
                             10+1+self.margin.value*2, 'int', True)
                    for k in range(self.dimension)]
        # spacing symbols, dx1, dx2, ...
        self.spacing = [Variable('dx'+str(k+1),
                        int(self.size[k] /
                            (self.dim[k].value-1-self.margin.value*2)),
                        self.real_t, True) for k in range(self.dimension)]
        # indices symbols, x1, x2 ...
        self.index = [Symbol('x'+str(k+1)) for k in range(dimension)]
        # default 2nd order in time, 4th order in space, i.e. (2,4) scheme
        default_accuracy = [1] + [2]*self.dimension
        self.set_accuracy(default_accuracy)

        self.t = Symbol('t')
        self.dt = Variable('dt', 0.01, self.real_t, True)
        self.ntsteps = Variable('ntsteps', 100, 'int', True)
        # user defined variables
        # use dictionary because of potential duplication
        self.defined_variable = {}

        self._update_spacing()

    def set_accuracy(self, accuracy):
        """
        - set the accuracy of the scheme
        - create the t variables for time stepping
        e.g. t0, t1 for 2nd order scheme, t0, t1, t2, t3 for 4th order scheme
        :param accuracy: list of time accuracy followed by spatial accuracy
        e.g. [1,2,2,2] for (2,4) scheme
        """
        self.order = accuracy
        # periodicity for time stepping
        self.tp = Variable('tp', self.order[0]*2, 'int', True)
        # add time variables for time stepping
        self.time = []
        for k in range(self.order[0]*2):
            name = 't' + str(k)
            v = Variable(name, 0, 'int', False)
            self.time.append(v)

    def set_io(self, io):
        """
        set io swtich
        """
        assert io is True or io is False
        self.io = io

    def set_omp(self, omp):
        """
        set omp swtich
        """
        assert omp is True or omp is False
        self.omp = omp

    def set_ivdep(self, ivdep):
        """
        set ivdep swtich
        """
        assert ivdep is True or ivdep is False
        self.ivdep = ivdep

    def set_simd(self, simd):
        """
        set simd swtich
        """
        assert simd is True or simd is False
        self.simd = simd

    def set_expand(self, expand):
        """
        set expand swtich
        """
        assert expand is True or expand is False
        self.expand = expand

    def set_eval_const(self, eval_const):
        """
        set eval_const swtich
        """
        assert eval_const is True or eval_const is False
        self.eval_const = eval_const

    def set_double(self, double):
        """
        set double swtich
        redefine the spacing varialbes with new type
        """
        assert double is True or double is False
        self.double = double
        self.real_t = 'double' if self.double else 'float'
        # spacing symbols, dx1, dx2, ...
        self.spacing = [Variable('dx' + str(k+1), 0.1, self.real_t, True)
                        for k in range(self.dimension)]

    def _update_spacing(self):
        """
        - update spacing variables (dx1, dx2 etc) with new values
        - called when domain size or grid size changes
        - set vec_size variable value to be total number of grid points needed
        used to define field as std::vector<>
        """
        self.spacing = [Variable('dx'+str(k+1),
                        self.size[k]/(self.dim[k].value-1-self.margin.value*2),
                        self.real_t, True) for k in range(self.dimension)]

        expr = 2*self.order[0]
        for d in self.dim:
            expr *= d
        self.vec_size = Variable('vec_size', expr, 'int', True)

    def set_domain_size(self, size):
        """
        set the (physical) domain size
        :param size: domain size as tuple, e.g. (1.0, 1.0, 1.0) for unit cube
        """
        self.size = size
        self._update_spacing()

    def set_grid_size(self, size):
        """
        set grid size (number of grids in each dimension)
        update the spacing variables with new value (domain size / grid size)
        :param size: grid spacing as tuple, e.g. (100, 100, 100)
        """
        self.dim = [Variable('dim'+str(k+1), size[k]+1+2*self.margin.value,
                    'int', True) for k in range(self.dimension)]

        self._update_spacing()

    def set_index(self, indices):
        """
        set indices of the grid
        :param indices: list of symbols as indices, e.g. [t,x,y,z]
        """
        self.index = indices

    def set_variable(self, var, value=0, type='int', constant=False):
        """
        set user defined variables, update defined_variable list
        :param var: name of variable
        value: variable value
        type: variable type as string
        constant: whether variable is constant
        """
        if isinstance(var, Symbol):
            var = var.name
        self.defined_variable[var] = Variable(var, value, type, constant)

    def get_time_step_limit(self):
        """
        return the maximum time step size for stability
        """
        if self.read:
            return 'physical parameters to be read from file,'\
                + ' please compile and run the executable'
        l = self.defined_variable['lambda'].value
        m = self.defined_variable['mu'].value
        r = self.defined_variable['rho'].value
        Vp = ((l + 2*m)/r)**0.5
        h = min([sp.value for sp in self.spacing])
        if self.order[1] == 1:
            return h/Vp/(3**0.5)
        elif self.order[1] == 2:
            return 0.495*h/Vp
        else:
            return 'not implemented yet'

    def set_time_step(self, dt, tmax):
        """
        set the time step size
        :param dt: time step size
        tmax: maximum simulation time
        """
        self.dt.value = dt
        self.ntsteps.value = int(tmax/dt)

    def set_stress_fields(self, sfields):
        """
        assign stress fields
        :param sfields: list of stress fields
        """
        num = self.dimension+self.dimension*(self.dimension-1)/2
        if not len(sfields) == num:
            raise Exception('wrong number of stress fields: '
                            + str(num) + ' fields required.')
        self.sfields = sfields

    def set_velocity_fields(self, vfields):
        """
        assign velocity fields
        :param vfields: list of velocity fields
        """
        num = self.dimension
        if not len(vfields) == num:
            raise Exception('wrong number of velocity fields: '
                            + str(num) + ' fields required.')
        self.vfields = vfields

    def calc_derivatives(self):
        """
        calculate the FD approximation of derivatives
        save the calculated expressions in field.d lists
        """
        l = [self.dt] + self.spacing
        for field in self.sfields+self.vfields:
            for k in range(self.dimension+1):
                # loop through dimensions
                h = Symbol(l[k].name)
                for o in range(1, self.order[k]+1):
                    # loop through order of derivatives
                    field.calc_derivative([self.t]+self.index, k, h, o)

    def get_all_variables(self):
        """
        return list of all variables defined
        """
        variables = self.dim + self.spacing + self.time\
            + [self.tp, self.dt, self.margin, self.ntsteps, self.vec_size]\
            + self.defined_variable.values()
        return variables

    def solve_fd(self, eqs):
        """
        - from the input PDEs, solve for the time stepping compuation kernel
        of all stress and velocity fields of the grid
        :param eqs: PDEs
        Note: sympy doesn't support solving simultaneous equations
        containing Indexed objects
        need to pass one equation to eqch field to solve
        """
        t = self.t
        t1 = t+hf+(self.order[0]-1)  # the most advanced time index
        index = [t1] + self.index
        self.eqs = eqs
        # populate substitution dictionary if evluating constants in kernel
        if self.eval_const:
            dict1 = {}
            variables = self.get_all_variables()
            for v in variables:
                if v.constant:
                    dict1[Symbol(v.name)] = v.value

        for field, eq in zip(self.vfields+self.sfields, eqs):
            # want the LHS of express to be at time t+1
            kernel = solve(eq, field[index])[0]
            kernel = kernel.subs({t: t+1-hf-(self.order[0]-1)})

            if self.expand:
                # expanding kernel (de-factorisation)
                kernel = expand(kernel)

            if self.eval_const:
                # substitute constants with values
                kernel = kernel.subs(dict1)

            if self.read:
                # replace with effective media parameters
                # if reading data from file (i.e. parameters not constant)
                kernel = kernel.subs(field.media_param)
            field.set_fd_kernel(kernel)

    def associate_fields(self):
        """
        - match up stress fields with velocity fields
        e.g. Vx will be matched with Txx, Txy, Txz
        - raise exception if stress fields cannot be matched
        with velocity fields (error with direction setting of fields)
        - associated stress fields are used to solve for boundary cells
        of velocity fields
        """
        for v in self.vfields:
            fields = [None]
            for d in range(1, self.dimension+1):
                lookfor = tuple(sorted([v.direction, d]))
                fields.append([f for f in self.sfields
                               if f.direction == lookfor][0])
            if not len(fields) == self.dimension+1:
                print len(fields)
                raise Exception('error in field directions')
            v.associate_stress_fields(fields)

    def set_free_surface_boundary(self, dimension, side):
        """
        set free surface boundary condition
        calls set_free_surface() method of SField and VField
        :param dimension: the normal direction to identify the surface
        e.g. dimension=1 for y-z place
        :param side: the side of the surface
        side=0 for bottom surface, side=1 for top surface
        """
        self.associate_fields()
        index = [self.t] + self.index
        for field in self.sfields+self.vfields:
            if side == 0:
                field.set_free_surface(index, dimension,
                                       self.margin.value, side, self.read)
            else:
                field.set_free_surface(index, dimension,
                                       self.dim[dimension-1]
                                       - self.margin.value - 1,
                                       side, self.read)

    def set_media_params(self, read=False, rho=1.0, vp=1.0, vs=0.5,
                         rho_file='', vp_file='', vs_file=''):
        """
        set media phyiscal parameters rho (density), Vp (primary velocity),
        Vs (secondary velocity) values of the grid
        or declare filename of the input data file
        :param read: whether to read the data from file or not
        default to False (not reading data from file)
        :param rho: the value of rho if read=False (constant across the grid)
        :param vp: the value of Vp if read=False (constant across the grid)
        :param vs: the value of Vs if read=False (constant across the grid)
        :param rho_file: name of input data file for rho, in binary form
        :param vp_file: name of input data file for Vp, in binary form
        :param vs_file: name of input data file for Vs, in binary form
        """
        self.read = read
        if self.read:
            self.rho_file = rho_file
            self.vp_file = vp_file
            self.vs_file = vs_file
            # input data
            self.rho = Field('rho')  # field to store rho
            self.vp = Field('vp')  # field to store vp
            self.vs = Field('vs')  # field to store vs
            # calculated from input data
            # list of fields to store beta (buoyancy and effective buoyancy)
            self.beta = [Field('beta')]\
                + [Field('beta'+str(k+1)) for k in range(self.dimension)]
            self.lam = Field('lambda')  # field to store lambda
            # list of fields to store mu
            # (shear modulus and effective shear modulus)
            self.mu = [Field('mu')] + [Field('mu12'), Field('mu13'),
                                       Field('mu23')]
            self.assign_media_param()

        else:
            self.set_variable('rho', rho, 'float', True)
            self.set_variable('beta', 1.0/rho, 'float', True)
            self.set_variable('lambda', rho*(vp**2-2*vs**2), 'float', True)
            self.set_variable('mu', rho*(vs**2), 'float', True)

    def assign_media_param(self):
        """
        - assign media parameters (beta, lambda, mu)
        to the stress fields and velocity fields
        - decide which media parameters should be used for the updating kernel
        - create dictionary used to map the parameters
        """
        for f in self.vfields + self.sfields:
            if isinstance(f, VField):
                if f.staggered[1]:
                    # Vx
                    f.set_media_param({Symbol('beta'):
                                       self.beta[1][self.index],
                                       Symbol('lambda'):
                                       self.lam[self.index],
                                       Symbol('mu'):
                                       self.mu[0][self.index]})
                elif f.staggered[2]:
                    # Vy
                    f.set_media_param({Symbol('beta'):
                                       self.beta[2][self.index],
                                       Symbol('lambda'):
                                       self.lam[self.index],
                                       Symbol('mu'):
                                       self.mu[0][self.index]})
                else:
                    # Vz
                    f.set_media_param({Symbol('beta'):
                                       self.beta[3][self.index],
                                       Symbol('lambda'):
                                       self.lam[self.index],
                                       Symbol('mu'):
                                       self.mu[0][self.index]})
            else:
                if f.staggered[1] and f.staggered[2]:
                    # Txy
                    f.set_media_param({Symbol('beta'):
                                       self.beta[0][self.index],
                                       Symbol('lambda'):
                                       self.lam[self.index],
                                       Symbol('mu'):
                                       self.mu[1][self.index]})
                elif f.staggered[1] and f.staggered[3]:
                    # Txz
                    f.set_media_param({Symbol('beta'):
                                       self.beta[0][self.index],
                                       Symbol('lambda'):
                                       self.lam[self.index],
                                       Symbol('mu'):
                                       self.mu[2][self.index]})
                elif f.staggered[2] and f.staggered[3]:
                    # Tyz
                    f.set_media_param({Symbol('beta'):
                                       self.beta[0][self.index],
                                       Symbol('lambda'):
                                       self.lam[self.index],
                                       Symbol('mu'):
                                       self.mu[3][self.index]})
                else:
                    # Txx, Tyy, Tzz
                    f.set_media_param({Symbol('beta'):
                                       self.beta[0][self.index],
                                       Symbol('lambda'):
                                       self.lam[self.index],
                                       Symbol('mu'):
                                       self.mu[0][self.index]})

    # ------------------- sub-routines for output -------------------- #

    def define_variables(self):
        """
        - generate code for declaring variables
        - return the generated code as string
        """
        result = ''
        variables = self.get_all_variables()
        for v in variables:
            line = ''
            if v.constant:
                line += 'const '
            line += v.type + ' ' + v.name + ' = ' + str(v.value) + ';\n'
            result += line
        return result

    def declare_fields(self):
        """
        - generate code for delcaring fields
        - the generated code first declare fields as std::vector
        of size=vec_size, then cast to multidimensional array
        - return the generated code as string
        """
        result = ''
        arr = ''  # = [dim1][dim2][dim3]...
        for d in self.dim:
            arr += '[' + d.name + ']'
        for field in self.sfields + self.vfields:
            vec = '_' + ccode(field.label) + '_vec'
            result += 'std::vector<' + self.real_t + '> ' + vec \
                + '(' + self.vec_size.name + ');\n'
            result += self.real_t + ' (*' + ccode(field.label) + ')' + arr \
                + '= (' + self.real_t + ' (*)' + arr + ') ' + vec \
                + '.data();\n'

        if self.read:
            # add code to read data
            result += self.read_data()

        return result

    def read_data(self):
        """
        - generate code for reading data (rho, Vp, Vs) from input files
        - calculate effective media parameters beta, lambda, mu from the data
        """
        result = ''
        if self.read:
            arr = ''  # =[dim2][dim3]...
            for d in self.dim[1:]:
                arr += '[' + d.name + ']'
            vsize = 1
            for d in self.dim:
                vsize *= d.value
            # declare fields to read physical parameters from file
            # always use float not double
            loop = [self.rho, self.vp, self.vs] + self.beta + [self.lam] \
                + self.mu
            for field in loop:
                vec = '_' + ccode(field.label) + '_vec'
                result += 'std::vector<float> ' + vec + '(' \
                    + str(vsize) + ');\n'
                result += 'float (*' + ccode(field.label) + ')' + arr \
                    + '= (float (*)' + arr + ') ' + vec + '.data();\n'
            # read from file
            result += 'opesci_read_simple_binary("' + self.rho_file + '",_' \
                + ccode(self.rho.label) + '_vec);\n'
            result += 'opesci_read_simple_binary("' + self.vp_file + '",_' \
                + ccode(self.vp.label) + '_vec);\n'
            result += 'opesci_read_simple_binary("' + self.vs_file + '",_' \
                + ccode(self.vs.label) + '_vec);\n'
            # calculated effective media parameter
            idx = self.index
            # make copies of index
            idx100 = list(idx)
            idx010 = list(idx)
            idx001 = list(idx)
            idx110 = list(idx)
            idx101 = list(idx)
            idx011 = list(idx)
            # shift the indices to obtain idx100=[x+1,y,z] etc
            idx100[0] += 1
            idx010[1] += 1
            idx001[2] += 1
            idx110[0] += 1
            idx110[1] += 1
            idx101[0] += 1
            idx101[2] += 1
            idx011[1] += 1
            idx011[2] += 1
            # beta
            kernel = ccode(self.beta[0][idx]) + '=' + ccode(1.0/self.rho[idx])
            result += self.simple_loop(kernel)
            # beta1 (effective bouyancy in x direction)
            kernel = ccode(self.beta[1][idx]) + '=' \
                + ccode((self.beta[0][idx] + self.beta[0][idx100])/2.0)
            result += self.simple_loop(kernel)
            # beta2 (effective bouyancy in y direction)
            kernel = ccode(self.beta[2][idx]) + '=' \
                + ccode((self.beta[0][idx] + self.beta[0][idx010])/2.0)
            result += self.simple_loop(kernel)
            # beta3 (effective bouyancy in z direction)
            kernel = ccode(self.beta[3][idx]) + '=' + \
                ccode((self.beta[0][idx] + self.beta[0][idx001])/2.0)
            result += self.simple_loop(kernel)
            # lambda
            kernel = ccode(self.lam[idx]) + '=' + \
                ccode(self.rho[idx]*(self.vp[idx]**2-2*self.vs[idx]**2))
            result += self.simple_loop(kernel)
            # mu
            kernel = ccode(self.mu[0][idx]) + '=' \
                + ccode(self.rho[idx]*(self.vs[idx]**2))
            result += self.simple_loop(kernel)
            # mu12 (effective shear modulus for shear stress sigma_xy)
            kernel = ccode(self.mu[1][idx]) + '=' \
                + ccode(1.0/(0.25*(1.0/self.mu[0][idx]+1.0/self.mu[0][idx100]
                        + 1.0/self.mu[0][idx010]+1.0/self.mu[0][idx110])))
            result += self.simple_loop(kernel)
            # mu13 (effective shear modulus for shear stress sigma_xz)
            kernel = ccode(self.mu[2][idx]) + '=' \
                + ccode(1.0/(0.25*(1.0/self.mu[0][idx]+1.0/self.mu[0][idx100]
                        + 1.0/self.mu[0][idx001]+1.0/self.mu[0][idx101])))
            result += self.simple_loop(kernel)
            # mu23 (effective shear modulus for shear stress sigma_yz)
            kernel = ccode(self.mu[3][idx]) + '=' \
                + ccode(1.0/(0.25*(1.0/self.mu[0][idx]+1.0/self.mu[0][idx010]
                        + 1.0/self.mu[0][idx001]+1.0/self.mu[0][idx011])))
            result += self.simple_loop(kernel)
        return result

    def simple_loop(self, kernel):
        """
        - helper function to generate simple nested loop over the entire domain
        (not including ghost cells) with kernel at the inner loop
        - variables defined in self.index are used as loop variables
        """
        result = ''
        tmpl = self.lookup.get_template('generic_loop.txt')
        m = self.margin.value
        for d in range(self.dimension-1, -1, -1):
            i = self.index[d]
            i0 = m
            i1 = ccode(self.dim[d]-m)
            if d == self.dimension-1:
                # inner loop
                result += kernel + ';\n'
            dict1 = {'i': i, 'i0': i0, 'i1': i1, 'body': result}
            result = render(tmpl, dict1)
        return result

    def initialise(self):
        """
        generate code for initialisation of the fields
        - substitute starting time to the analytical function of the fields
        - substitute field coordinates calculated from array indices
        to the analytical function of the fields
        - generate inner loop by inserting kernel into Mako template
        - recursive insertion to generate nested loop
        return generated code as string
        """

        tmpl = self.lookup.get_template('generic_loop.txt')
        result = ''
        m = self.margin.value
        loop = [Symbol('_'+x.name) for x in self.index]  # symbols for loop

        for field in self.sfields+self.vfields:
            body = ''
            if self.omp:
                result += '#pragma omp for\n'
            # populate xvalue, yvalue zvalue code
            for d in range(self.dimension-1, -1, -1):
                i = loop[d]
                i0 = m
                if field.staggered[d+1]:
                    i1 = ccode(self.dim[d]-m-1)
                    expr = self.spacing[d]*(loop[d] - self.margin.value + 0.5)
                else:
                    i1 = ccode(self.dim[d]-m)
                    expr = self.spacing[d]*(loop[d] - self.margin.value)
                pre = self.real_t + ' ' + self.index[d].name + '= ' \
                    + ccode(expr) + ';\n'
                post = ''
                if d == self.dimension-1:
                    # inner loop
                    # first time step
                    t0 = self.dt.value/2 if field.staggered[0] else 0
                    sol = field.sol.subs(self.t, t0)
                    if self.read:
                        sol = sol.subs(field.media_param)
                        for idx in self.index:
                            sol = sol.subs(idx, '_'+idx.name)
                    body = ccode(field[[0]+loop]) + '=' \
                        + ccode(sol) + ';\n'
                body = pre + body + post
                dict1 = {'i': i, 'i0': i0, 'i1': i1, 'body': body}
                body = render(tmpl, dict1)

            result += body
        return result

    def initialise_boundary(self):
        """
        - generate code for initialisation of boundary ghost cells
        - generate generic boundary cell code
        replace array indices [t] with [0]
        - return gerneated code as string
        """
        result = self.stress_bc().replace('[t1]', '[0]')
        result += self.velocity_bc().replace('[t1]', '[0]')
        return result

    def stress_loop(self):
        """
        generate code for stress field update loop
        - loop through stress fields to generate code of computation kernel
        - generate inner loop by inserting kernel into Mako template
        - recursive insertion to generate nested loop
        return generated code as string
        """
        tmpl = self.lookup.get_template('generic_loop.txt')
        m = self.margin.value
        body = ''
        for d in range(self.dimension-1, -1, -1):
            i = self.index[d]
            i0 = m
            i1 = ccode(self.dim[d]-m)
            if d == self.dimension-1:
                # inner loop
                idx = [self.time[1]] + self.index
                for field in self.sfields:
                    body += ccode(field[idx]) + '=' \
                        + ccode(field.fd_align.xreplace({self.t+1:
                                                        self.time[1],
                                                        self.t:
                                                        self.time[0]})) \
                        + ';\n'
            dict1 = {'i': i, 'i0': i0, 'i1': i1, 'body': body}
            body = render(tmpl, dict1)
            if self.ivdep and d == self.dimension-1:
                    body = '#pragma ivdep\n' + body
            if self.simd and d == self.dimension-1:
                    body = '#pragma simd\n' + body

        if self.omp:
            body = '#pragma omp for\n' + body

        return body

    def velocity_loop(self):
        """
        generate code for velocity field update loop
        - loop through velocity fields to generate code of computation kernel
        - generate inner loop by inserting kernel into Mako template
        - recursive insertion to generate nested loop
        return generated code as string
        """
        tmpl = self.lookup.get_template('generic_loop.txt')
        m = self.margin.value
        body = ''
        for d in range(self.dimension-1, -1, -1):
            i = self.index[d]
            i0 = m
            i1 = ccode(self.dim[d]-m)
            if d == self.dimension-1:
                # inner loop
                idx = [self.time[1]] + self.index
                for field in self.vfields:
                    body += ccode(field[idx]) + '=' \
                        + ccode(field.fd_align.xreplace({self.t+1:
                                                        self.time[1],
                                                        self.t:
                                                        self.time[0]})) \
                        + ';\n'
            dict1 = {'i': i, 'i0': i0, 'i1': i1, 'body': body}
            body = render(tmpl, dict1)
            if self.ivdep and d == self.dimension-1:
                    body = '#pragma ivdep\n' + body
            if self.simd and d == self.dimension-1:
                    body = '#pragma simd\n' + body

        if self.omp:
            body = '#pragma omp for\n' + body

        return body

    def stress_bc(self):
        """
        generate code for updating stress field boundary ghost cells
        - generate inner loop by inserting boundary code (saved in field.bc)
        - recursive insertion to generate nested loop
        - loop through all stress fields and sides
        return generated code as string
        """
        tmpl = self.lookup.get_template('generic_loop.txt')
        result = ''
        for field in self.sfields:
            for d in range(self.dimension):
                for side in range(2):
                    # skip if this boundary calculation is not needed
                    if field.bc[d+1][side] == '':
                        continue
                    if self.omp:
                        result += '#pragma omp for\n'
                    body = ''
                    for d2 in range(self.dimension-1, -1, -1):
                        # loop through other dimensions
                        if not d2 == d:
                            i = self.index[d2]
                            i0 = 0
                            i1 = self.dim[d2]
                            if body == '':
                                # inner loop, populate ghost cell calculation
                                body = field.bc[d+1][side]
                                dict1 = {'i': i, 'i0': i0,
                                         'i1': i1, 'body': body}
                                body = render(tmpl, dict1).replace('[t]',
                                                                   '[t1]')
                                if self.ivdep:
                                    body = '#pragma ivdep\n' + body
                                if self.simd:
                                    body = '#pragma simd\n' + body
                            else:
                                dict1 = {'i': i, 'i0': i0,
                                         'i1': i1, 'body': body}
                                body = render(tmpl, dict1).replace('[t]',
                                                                   '[t1]')

                    result += body

        return result

    def velocity_bc(self):
        """
        generate code for updating stress field boundary ghost cells
        - generate inner loop by inserting boundary code (saved in field.bc)
        - recursive insertion to generate nested loop
        - loop through all velocity fields and sides
        return generated code as string
        """
        tmpl = self.lookup.get_template('generic_loop.txt')
        result = ''
        for d in range(self.dimension):
            # update the staggered field first
            # because other fields depends on it
            sequence = [f for f in self.vfields if f.staggered[d+1]] \
                + [f for f in self.vfields if not f.staggered[d+1]]
            for field in sequence:
                for side in range(2):
                    # skip if this boundary calculation is not needed
                    if field.bc[d+1][side] == '':
                        continue
                    if self.omp:
                        result += '#pragma omp for\n'
                    body = ''
                    for d2 in range(self.dimension-1, -1, -1):
                        # loop through other dimensions
                        if not d2 == d:
                            i = self.index[d2]
                            i0 = 1
                            i1 = self.dim[d2]-1
                            if body == '':
                                # inner loop, populate ghost cell calculation
                                body = field.bc[d+1][side]
                                dict1 = {'i': i, 'i0': i0,
                                         'i1': i1, 'body': body}
                                body = render(tmpl, dict1).replace('[t]',
                                                                   '[t1]')
                                if self.ivdep:
                                    body = '#pragma ivdep\n' + body
                                if self.simd:
                                    body = '#pragma simd\n' + body
                            else:
                                dict1 = {'i': i, 'i0': i0,
                                         'i1': i1, 'body': body}
                                body = render(tmpl, dict1).replace('[t]',
                                                                   '[t1]')

                    result += body

        return result

    def time_stepping(self):
        """
        generate time index variable for time stepping
        e.g. for 2nd order time-accurate scheme, varibales are t0, t1
        for 4th order time-accurate scheme, variables are t0, t1, t2, t3
        the variables are used to address the field arrays
        e.g. in 2nd order scheme, U[t1] will be updated using U[t0]
        the variables are calculated by taking mod with time periodicity
        return generated code as string
        """

        result = ''
        tmpl = self.lookup.get_template('time_stepping.txt')
        _ti = Symbol('_ti')
        body = ''

        for i in range(len(self.time)):
            lhs = self.time[i].name
            if i == 0:
                rhs = ccode(_ti % self.tp)
            else:
                rhs = ccode((self.time[i-1]+1) % self.tp)
            body += lhs + ' = ' + rhs + ';\n'

        dict1 = {'body': body}
        result = render(tmpl, dict1)
        return result

    def output_step(self):
        """
        - generate code for output at each time step
        - typically output selected fields in vtk format
        - return generated code as string
        """
        result = ''
        result += self.vfields[0].vtk_save_field()
        return result

    def converge_test(self):
        """
        - generate code for convergence test
        - convergence test implemented by calculating L2 norm
        of the simulation against analytical solution
        - L2 norm of each field is calculated and output with printf()
        - return generated code as string
        """
        tmpl = self.lookup.get_template('generic_loop.txt')
        result = ''
        m = self.margin.value
        ti = self.ntsteps.value % 2  # last updated grid
        loop = [Symbol('_'+x.name) for x in self.index]  # symbols for loop

        for i in range(len(self.spacing)):
            result += 'printf("' + str(self.spacing[i].value) + '\\n");\n'

        for field in self.sfields+self.vfields:
            body = ''
            l2 = ccode(field.label)+'_l2'
            idx = [ti] + loop
            result += self.real_t + ' ' + l2 + ' = 0.0;\n'
            # populate xvalue, yvalue zvalue code
            for d in range(self.dimension-1, -1, -1):
                i = loop[d]
                i0 = m
                if field.staggered[d+1]:
                    i1 = ccode(self.dim[d]-m-1)
                    expr = self.spacing[d]*(loop[d] - self.margin.value + 0.5)
                else:
                    i1 = ccode(self.dim[d]-m)
                    expr = self.spacing[d]*(loop[d] - self.margin.value)
                pre = self.real_t + ' ' + self.index[d].name + '= ' \
                    + ccode(expr) + ';\n'
                if d == self.dimension-1:
                    # inner loop
                    tn = self.dt.value*self.ntsteps.value \
                        if not field.staggered[0] \
                        else self.dt.value*self.ntsteps.value \
                        + self.dt.value/2.0
                    body = l2 + '+=' \
                        + ccode((field[idx] -
                                (field.sol.subs(self.t, tn)))**2.0) + ';\n'
                body = pre + body
                dict1 = {'i': i, 'i0': i0, 'i1': i1, 'body': body}
                body = render(tmpl, dict1)

            result += body
            volume = 1.0
            for i in range(len(self.spacing)):
                volume *= self.spacing[i].value
            l2_value = 'pow(' + l2 + '*' + ccode(volume) + ', 0.5)'
            result += 'printf("' + l2 + '\\t%.10f\\n", ' + l2_value + ');\n'

        return result
