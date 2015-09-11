from sympy import Expr, Symbol, Indexed, Rational, IndexedBase, factorial, Matrix


hf = Rational(1, 2)  # 1/2


def get_all_objects(expr, typ):
    """
    helper function to return list of all objects in expr of type=typ
    """
    if isinstance(expr, typ):
        return [expr]
    if not isinstance(expr, Expr):
        return []
    else:
        args = expr.args
        result = []
        for arg in args:
            result += get_all_objects(arg, typ)
        return result


def variable_to_symbol(variables):
    """
    helper function to convert list of variables to symbols
    """
    return [Symbol(variable.name) for variable in variables]


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


def Taylor_generic(dx, left, right):
    """
    generic version of Taylor()
    expand neighbour grid point located from -left to +right
    """
    l = []
    n = left + right + 1
    for i in range(-left, right+1):
        ll = [tc((i*dx), j) for j in range(n)]
        l.append(ll)
    return Matrix(l)


def Deriv(U, index, k, d, n):
    """
    calculate the FD approximation for nth derivative
    the function works by using M * D = R, as specified in Taylor()
    then inverting M, so that D = M_inverse * R
    :param U: the field
    :param index: list of indices (Symbol objects) of field U
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
    s = [0]*len(index)
    s[k] = 1  # switch on kth dimension
    # generate matrix of RHS, i.e. [ ... U[x-1], U[x], U[x+1] ... ]
    if len(index) == 1:
        RX = Matrix([U[index[0]+s[0]*x] for x in range(-n+1, n)])
    elif len(index) == 2:
        RX = Matrix([U[index[0]+s[0]*x, index[1]+s[1]*x] for x in range(-n+1, n)])
    elif len(index) == 3:
        RX = Matrix([U[index[0]+s[0]*x,
                    index[1]+s[1]*x,
                    index[2]+s[2]*x] for x in range(-n+1, n)])
    elif len(index) == 4:
        RX = Matrix([U[index[0]+s[0]*x,
                    index[1]+s[1]*x,
                    index[2]+s[2]*x,
                    index[3]+s[3]*x] for x in range(-n+1, n)])
    else:
        raise NotImplementedError(">4 dimensions, need to fix")

    return M.inv() * RX


def Deriv_generic(U, index, dimension, delta, left, right):
    """
    generic version of Derive()
    """
    M = Taylor_generic(delta, left, right)
    mask = [0]*len(index)
    mask[dimension] = 1  # switch on kth dimension
    # generate matrix of RHS, i.e. [ ... U[x-1], U[x], U[x+1] ... ]
    index2 = list(index)
    ll = []
    for i in range(-left, right+1):
        index2[dimension] += i
        ll.append(U[index2])
    RX = Matrix(ll)

    return M.inv() * RX


def Deriv_half(U, index, dimension, delta, order):
    """
    similar function as Deriv() for staggered grids
    calculate the FD approximation for nth derivative
    the function works by using M * D = R, as specified in Taylor()
    then inverting M, so that D = M_inverse * R
    :param U: the field
    :param index: list of indices (Symbol objects) of field U
    e.g. possibley [t,x,y,z] for 3D
    :param k: determine the derivative is with respect to which dimension
    e.g. dU/dt if k=0, dU/dx if k=0, assuming i=[t,x,y,z]
    :param delta: spacing of the grid, e.g. dx, dt
    :param order: order of accuracy of approximation/2
    e.g. order=1 for 2nd-order FD approximation, order=2 for 4th-order FD approximation
    returns D=list of expressions for FD approxiation
    D[1] for first derivative, D[2] for 2nd derivative etc
    raises NotImplementedError exception if dimension is more than 4
    (1 time and 3 space dimensions)
    """
    M = Taylor_half(delta, order)
    mask = [0]*len(index)
    mask[dimension] = 1  # switch on kth dimension
    # generate matrix of RHS, i.e. [ ... U[x-1], U[x], U[x+1] ... ]
    if len(index) == 1:
        RX = Matrix([U[index[0]+mask[0]*x*hf] for x in range(-order*2+1, order*2, 2)])
    elif len(index) == 2:
        RX = Matrix([U[index[0]+mask[0]*x*hf,
                    index[1]+mask[1]*x*hf] for x in range(-order*2+1, order*2, 2)])
    elif len(index) == 3:
        RX = Matrix([U[index[0]+mask[0]*x*hf,
                    index[1]+mask[1]*x*hf,
                    index[2]+mask[2]*x*hf] for x in range(-order*2+1, order*2, 2)])
    elif len(index) == 4:
        RX = Matrix([U[index[0]+mask[0]*x*hf,
                    index[1]+mask[1]*x*hf,
                    index[2]+mask[2]*x*hf,
                    index[3]+mask[3]*x*hf] for x in range(-order*2+1, order*2, 2)])
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
