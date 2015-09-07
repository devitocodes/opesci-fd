from compilation import get_package_dir
from sympy import factorial, Matrix, Rational, Indexed, IndexedBase
from sympy import solve, Eq
from codeprinter import ccode, ccode_eq, render
from mako.lookup import TemplateLookup
from os import path

hf = Rational(1, 2)  # 1/2

__all__ = ['SField', 'VField']


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


def Deriv_half(U, i, k, d, n, p=False):
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
        if p:

            RX = Matrix([U[i[0]+s[0]*x*hf,#polly
                        i[1]+s[1]*x*hf,
                        i[2]+s[2]*x*hf,
                        i[3]+s[3]*x*hf] for x in range(-n*2+1, n*2, 2)])
   #         print RX
        else:
            RX = Matrix([U[i[0]+s[0]*x*hf,
                        i[1]+s[1]*x*hf,
                        i[2]+s[2]*x*hf,
                        i[3]+s[3]*x*hf] for x in range(-n*2+1, n*2, 2)])
  #          print RX
    else:
        raise NotImplementedError(">4 dimensions, need to fix")

    #print M.inv()
    result = M.inv() * RX
    #print result
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


def get_ops_expr(expr, arrays):
    """
    - get number of different operations in expression expr
    - types of operations are ADD (inc -) and MUL (inc /)
    - arrays (IndexedBase objects) in expr that are not in list arrays
    are added to the list
    - return tuple (#ADD, #MUL, list of unique names of fields)
    """
    # add array to list arrays if it is not in it
    if isinstance(expr, Indexed):
        base = str(expr.base.label)
        if base not in arrays:
            arrays += [base]
        return (0, 0, arrays)

    mul = 0
    add = 0
    if expr.is_Mul or expr.is_Add:
        args = expr.args
        # increment MUL or ADD by # arguments less 1
        # sympy multiplication and addition can have multiple arguments
        if expr.is_Mul:
            mul += len(args)-1
        else:
            add += len(args)-1
        arrays2 = arrays
        # recursive call of all arguments
        for expr2 in args:
            add2, mul2, arrays2 = get_ops_expr(expr2, arrays2)
            add += add2  # accumulate ADD
            mul += mul2  # acculate MUL

        return (add, mul, arrays2)
    # return zero and unchanged array if execution gets here
    return (0, 0, arrays)


class Field(IndexedBase):
    """
    - Class to represent fields on staggered grid
    - Extends sympy IndexedBase class, can be indexed with [] operator
    - Parent class to VFeild (velocity) and SField (stress)
    - Holds relevant information such as dimension, staggered-ness,
    expression for FD approximation of derivatives,
    code to calculate boundary cells
    """

    def __new__(typ, name, **kwargs):
        obj = IndexedBase.__new__(typ, name)
        return obj

    def __init__(self, *args, **kwargs):
        template_dir = path.join(get_package_dir(), "templates")
        staggered_dir = path.join(get_package_dir(), "templates/staggered")
        self.lookup = TemplateLookup(directories=[template_dir, staggered_dir])
        super(Field, self).__init__()

        # Pass additional arguments to self.set()
        if len(kwargs) > 0:
            self.set(**kwargs)

    def set(self, dimension, staggered):
        self.dimension = dimension
        self.staggered = staggered

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

    def calc_derivative(self, l, k, d, n, p=False):
        """
        assign list self.d with FD approximations of 1st derivatives
        such that self.d[k][n] = FD approximation of 1st derivative,
        in kth index, of nth order accuracy
        input param description same as Deriv_half()
        """
        self.d[k][n] = Deriv_half(self, l, k, d, n, p)[1]

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
       # print ', '.join("%s: %s" % item for item in vars(field).items())
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

        # print rhs

        # print lhs
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
