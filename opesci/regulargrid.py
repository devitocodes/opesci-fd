from grid import Grid
from variable import Variable
from fields import Media
from codeprinter import ccode, ccode_eq
from derivative import DDerivative
from util import *
from compilation import get_package_dir

from sympy import Symbol, Rational, solve, expand, Eq, symbols
from mako.lookup import TemplateLookup
import mmap
import cgen_wrapper as cgen
from os import path
from __builtin__ import str
from opesci.fields import RegularField
from argparse import ArgumentError

__all__ = ['RegularGrid']

class RegularGrid(Grid):
    template_base = 'regular3d_tmpl.cpp'
    
    template_keys = ['pluto', 'io', 'profiling', 'numevents_papi',
                     'time_stepping', 'define_constants', 'declare_fields',
                     'define_fields', 'store_fields', 'load_fields',
                     'initialise','define_profiling', 'define_papi_events', 'sum_papi_events', 'primary_loop', 'free_memory']

    _papi_events = []
    
    _switches = ['omp', 'ivdep', 'simd', 'double', 'expand', 'eval_const',
                 'output_vts', 'converge', 'profiling', 'pluto', 'fission']
    _params = ['c', 'v']
    def __init__(self, dimension, index=None, fields=None, double=False, profiling=False, pluto=False, fission=False, omp=True, ivdep=True, simd=False, io=False,
                 expand=True, eval_const=True, grid_size=(10, 10, 10), domain_size=None):
        super(RegularGrid, self).__init__()
        
        template_dir = path.join(get_package_dir(), "templates")
        regular_dir = path.join(get_package_dir(), "templates/regular")
        self.lookup = TemplateLookup(directories=[template_dir, regular_dir])
        
        self.dimension = dimension
        
        self.double = double
        self.real_t = 'double' if self.double else 'float'
        self.dt = Variable('dt', 0.01, self.real_t, True)
        
        
        
        
        
        
        self.ntsteps = Variable('ntsteps', 100, 'int', True)
        self.alignment = mmap.PAGESIZE  # default alignment for malloc
        # number of ghost cells for boundary
        self.margin = Variable('margin', 2, 'int', True)
        
        # grid size symbols, dim1, dim2, ...
        self.dim = [Variable('dim'+str(k+1),
                             10+1+self.margin.value*2, 'int', True)
                    for k in range(self.dimension)]
        
        self.size = [1.0] * dimension  # default domain size
        
        # default 2nd order in time, 4th order in space, i.e. (2,4) scheme
        default_order = [2] + [4]*self.dimension
        self.t = Symbol('_t')
        self.grid_size = grid_size
        self.max_derivative_order = 1
        if fields!=None:
            self.fields = fields
        self.set_order(default_order)
        self.set_grid_size(grid_size)
        
        # spacing symbols, dx1, dx2, ...
        self.spacing = [Variable('dx'+str(k+1),
                        int(self.size[k] /
                            (self.dim[k].value-1-self.margin.value*2)),
                        self.real_t, True) for k in range(self.dimension)]
        
        
        self.set_field_spacing()
        
        self.set_index(index)
        # user defined variables
        # use dictionary because of potential duplication
        self.defined_variable = {}
        self.defined_variable
        
        # Switches
        self.pluto = pluto
        self.omp = omp
        self.ivdep = ivdep
        self.simd = simd
        
        self.expand = expand
        self.eval_const = eval_const
        
        
        
        self.profiling = profiling
        self.fission = fission
        if domain_size:
            self.set_domain_size(domain_size)
        self.read = False
    
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
    
    def calc_derivatives(self, max_order=1):
        """
        populate field.d lists with Derivative objects
        """
        self.max_derivative_order = max_order
        self.set_order(self.order)
        for field in self.fields:
            field.populate_derivatives(max_order=max_order)  # VS scheme only requare 1st derivatives
    
    def set_order(self, order):
        """
        - set the order of approximation of the scheme
        - create the t variables for time stepping
        e.g. t0, t1 for 2nd order scheme, t0, t1, t2, t3 for 4th order scheme
        :param order: list of time order followed by spatial orders
        e.g. [2,4,4,4] for (2,4) scheme
        """
        for x in order:
            if x % 2 == 1:
                pass
                #raise ValueError(str(x) + ' is not a valid order (require even integer)')
        self.order = order
        num_time_vars = max(self.order[0], self.max_derivative_order+1)
        # periodicity for time stepping
        self.tp = Variable('tp', num_time_vars, 'int', True)
        # add time variables for time stepping: t0, t1 ...
        self.time = []
        
        for k in range(num_time_vars):
            name = str(self.t) + str(k)
            v = Variable(name, k, 'int', False)
            self.time.append(v)
        self.margin.value = self.order[1]/2
        self.set_grid_size(self.grid_size)
        self.update_field_order()
        self._update_spacing()
        
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
        self.grid_size = size
        self.dim = [Variable('dim'+str(k+1), size[k]+1+2*self.margin.value,
                    'int', True) for k in range(self.dimension)]

        self._update_spacing()
        
    def update_field_order(self):
        """
        update the order of acuracy of the fields
        """
        if hasattr(self, 'fields'):
            for field in self.fields:
                field.set_order(self.order)
    
    def set_time_step(self, dt, tmax):
        """
        set the time step size
        :param dt: time step size
        tmax: maximum simulation time
        """
        self.dt.value = dt
        self.ntsteps.value = int(tmax/dt)
        
    def set_switches(self, **kwargs):
        for switch, value in kwargs.items():
            if switch not in self._switches:
                raise KeyError("Unsupported switch: ", switch)
            if not isinstance(value, bool):
                raise ValueError("Only boolean values allowed for switches")

            setattr(self, switch, value)
            if switch == 'double':
                self.real_t = 'double' if self.double else 'float'
                self._update_spacing()
    
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

        expr = self.order[0]
        for d in self.dim:
            expr *= d
        self.vec_size = Variable('vec_size', expr, 'int', True)
    
    def set_index(self, index):
        """
        set indices of the grid
        :param indices: list of symbols as indices, e.g. [x,y,z]
        """
        if index is None:
            # default to indices symbols x1, x2 ...
            self.index = [Symbol('x'+str(k+1)) for k in range(self.dimension)]
        else:
            self.index = index

        if hasattr(self, 'fields'):
            for field in self.fields:
                # update the fields with new indices
                field.set_indices([self.t]+self.index)
    
    def set_params(self, **kwargs):
        for param, value in kwargs.items():
            if param not in self._params:
                raise KeyError("Unsupported parameter: ", param)
            setattr(self, param, value)
            self.set_variable(param, value, self.real_t, True)
    
    def set_field_spacing(self):
        for field in self.fields:
            field.set_spacing(variable_to_symbol([self.dt]+self.spacing))
    
    def solve_fd(self, equations):
        """
        - from the input PDEs, solve for the time stepping compuation kernel
        of all stress and velocity fields of the grid
        :param eqs: PDEs
        Note: sympy doesn't support solving simultaneous equations
        containing Indexed objects
        need to pass one equation to eqch field to solve
        """
        self.eq = []
        
        if not len(self.fields)==len(equations):
            raise KeyError("Number of equations must be the same as number of fields. Number of fields is ", len(self.fields))
        
        # save the equation
        for field, eq in zip(self.fields, equations):
            field.set_dt(eq.rhs)
            self.eq.append(eq)

        # replace derivatives in equations with FD approximations
        eqs = []
        for eq in self.eq:
            derivatives = get_all_objects(eq, DDerivative)
            for deriv in derivatives:
                eq = eq.subs(deriv, deriv.fd[deriv.max_accuracy])
            eqs += [eq]

        t = self.t
        t1 = t+hf+(self.order[0]/2-1)  # the most advanced time index
        t1_regular = t+1+(self.order[0]/2-1)
        index = [t1] + self.index
        index_regular = [t1_regular] + self.index
        
        simplify = True if max(self.order[1:]) <= 4 else False

        for field, eq in zip(self.fields, eqs):
            # want the LHS of express to be at time t+1
            if isinstance(field, RegularField):
                kernel = solve(eq, field[index_regular], simplify=simplify)[0] 
            else:
                kernel = solve(eq, field[index], simplify=simplify)[0]
            kernel = kernel.subs({t: t-(self.order[0]/2-1)})

            field.set_fd_kernel(kernel)
            
    def get_all_variables(self):
        """
        return list of all variables defined
        """
        variables = self.dim + self.spacing + self.time\
            + [self.tp, self.dt, self.margin, self.ntsteps]\
            + self.defined_variable.values()
        return variables
    
    def create_const_dict(self):
        """
        create the dictionary of all constants with their values, store in self.const_dict
        this is used to pre-evaluate all constants in the kernel (and bc)
        """
        dict1 = {}
        variables = self.get_all_variables()
        for v in variables:
            if v.constant:
                dict1[Symbol(v.name)] = v.value
        self.const_dict = dict1
    
    def get_kernel_ai(self, fields=None):
        """
        - get the arithmetic intensity of stress kernel
        - get the number of different operations of the stress field kernel
        - types of operations are ADD (inc -), MUL (inc /), LOAD, STORE
        - #LOAD = number of unique fields in the kernel
        - return tuple (#ADD, #MUL, #LOAD, #STORE)
        - arithmetic intensity AI = (ADD+MUL)/[(LOAD+STORE)*word size]
        - weighted AI, AI_w = (ADD+MUL)/(2*Max(ADD,MUL)) * AI
        """
        
        if fields==None:
            fields = self.fields
        
        if self.eval_const:
            self.create_const_dict()
        store = 0
        add = 0
        mul = 0
        arrays = []  # to store name of arrays loaded
        for field in fields:
            store += 1  # increment STORE by 1 (assignment)
            expr = self.transform_kernel(field)
            add2, mul2, arrays2 = get_ops_expr(expr, arrays)
            add += add2  # accumulate # ADD
            mul += mul2  # accumulate # MUL
            arrays = arrays2  # replace with new list of field names

        # 8 byte if double, 4 if float used
        # media parameter fields are always float, need to amend this
        word_size = 8 if self.double else 4
        load = len(arrays)
        ai = float(add+mul)/(load+store)/word_size
        ai_w = ai*(add+mul)/max(add, mul)/2.0

        return (ai, ai_w, add, mul, load, store)
    
    def transform_kernel(self, field):
        """
        transform the kernel of field based on setting of grid, such as expand, eval_const, read
        """
        kernel_new = field.kernel_aligned
        if self.expand:
            # expanding kernel (de-factorisation)
            kernel_new = expand(kernel_new)

        if self.eval_const:
            # substitute constants with values
            kernel_new = kernel_new.subs(self.const_dict)

        return kernel_new
    
    # ------------------- sub-routines for PAPI profiling ------------ #

    def set_papi_events(self, events=[]):
        self._papi_events = events

    @property
    def define_profiling(self):
        """Code fragment that defines global PAPI counters and events"""
        code = [cgen.Initializer(cgen.Value('float', 'g_%s' % v), 0.0) for v in ['rtime', 'ptime', 'mflops']]
        code += [cgen.Initializer(cgen.Value('long long', 'g_%s' % e), 0) for e in self._papi_events]
        return str(cgen.Module(code))

    @property
    def numevents_papi(self):
        return len(self._papi_events)

    @property
    def define_papi_events(self):
        """Code fragment that starts PAPI counters for specified events"""
        code = []
        code.append(cgen.Initializer(cgen.Value('int', 'numevents'), self.numevents_papi))
        code.append(cgen.ArrayOf(cgen.Value('int', 'events'), self.numevents_papi))
        code.append(cgen.ArrayOf(cgen.Value('long long', 'counters'), self.numevents_papi))
        code += [cgen.Statement('opesci_papi_name2event("%s", &(events[%d]))' % (e, i)) for i, e in enumerate(self._papi_events)]
        return str(cgen.Module(code))

    @property
    def sum_papi_events(self):
        """Code fragment that reads PAPI counters for specified events"""
        code = [cgen.Statement('profiling->g_%s += counters[%d]' % (e, i)) for i, e in enumerate(self._papi_events)]
        return str(cgen.Module(code))
    
    @property
    def io(self):
        """Flag whether to include I/O headers"""
        return False
    
    # ------------------- sub-routines for output -------------------- #

    @property
    def define_constants(self):
        """
        - generate code for declaring variables
        - return the generated code as string
        """

        result = []
        variables = self.get_all_variables()
        for v in variables:
            if v.constant:
                line = cgen.Initializer(cgen.Const(cgen.Value(v.type, v.name)), v.value)
            else:
                line = cgen.Initializer(cgen.Value(v.type, v.name), v.value)
            result.append(line)

        return str(cgen.Module(result))
    
    @property
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

        _ti = Symbol('_ti')
        body = []

        for i in range(len(self.time)):
            lhs = self.time[i].name
            if i == 0:
                rhs = ccode(_ti % self.tp)
            else:
                rhs = ccode((self.time[i-1]+1) % self.tp)
            body.append(cgen.Assign(lhs, rhs))

        body = cgen.Block(body)
        body = cgen.Module([cgen.Pragma('omp single'), body])
        return str(body)
    
    @property
    def define_fields(self):
        """Code fragment that defines field arrays"""
        result = []
        for f in self.fields:
            var = cgen.Pointer(cgen.Value(self.real_t, ccode(f.label)))
            result.append(var)

        return str(cgen.Module(result))

    @property
    def store_fields(self):
        """Code fragment that stores field arrays to 'grid' struct"""
        result = []
        for f in self.fields:
            assignment = cgen.Assign('grid->%s' % ccode(f.label), '(%s*) %s' % (self.real_t, ccode(f.label)))  # There must be a better way of doing this. This hardly seems better than string manipulation
            result.append(assignment)

        return str(cgen.Module(result))

    @property
    def load_fields(self):
        """Code fragment that loads field arrays from 'grid' struct"""
        idxs = ''.join(['[%d]' % d.value for d in self.dim])
        result = []
        for f in self.fields:
            back_assign = cgen.Initializer(cgen.Value(self.real_t, "(*%s)%s" % (ccode(f.label), idxs)), '(%s (*)%s) grid->%s' % (self.real_t, idxs, ccode(f.label)))  # Another hackish attempt.
            result.append(back_assign)

        return str(cgen.Module(result))
    
    
    @property
    def declare_fields(self):
        return self.declare_fields_raw(True)
    
    def declare_fields_raw(self, as_string=True):
        """
        - generate code for declaring fields
        - the generated code first declare fields as std::vector
        of size=vec_size, then cast to multidimensional array
        - return the generated code as string
        """
        result = []
        arr = ''  # = [dim1][dim2][dim3]...
        for d in self.dim:
            arr += '[' + d.name + ']'
        vsize = 1
        for d in self.dim:
            vsize *= d.value
        
        vsize *= len(self.time)
        statements = []
        for field in self.fields:
            vec = "_%s_vec" % ccode(field.label)
            vec_value = cgen.Pointer(cgen.Value(self.real_t, vec))
            # alloc aligned memory (on windows and linux)
            statements.append(vec_value)
            ifdef = cgen.IfDef('_MSC_VER', [cgen.Assign(vec, '(%s*) _aligned_malloc(%s*sizeof(%s), %s)' % (self.real_t, str(vsize), self.real_t, str(self.alignment)))],
                               [cgen.Statement('posix_memalign((void **)(&%s), %d, %d*sizeof(%s))' % (vec, self.alignment, vsize, self.real_t))])
            statements.append(ifdef)
            # cast pointer to multidimensional array
            cast_pointer = cgen.Initializer(cgen.Value(self.real_t, "(*%s)%s" % (ccode(field.label), arr)), '(%s (*)%s) %s' % (self.real_t, arr, vec))
            statements.append(cast_pointer)
        result += statements
        
        if as_string:
            return str(cgen.Module(result))
        else:
            return cgen.Module(result)
    
    @property
    def initialise(self):
        loop = [Symbol('_'+x.name) for x in self.index]  # symbols for loop

        statements = []
        for field in self.fields:
            body = []
            if self.omp:
                statements.append(cgen.Pragma('omp for schedule(static,1)'))
            # populate xvalue, yvalue zvalue code
            for d in range(self.dimension-1, -1, -1):
                i = loop[d]
                i0 = 0
                
                i1 = ccode(self.dim[d])
                expr = self.spacing[d]*(loop[d])
                pre = []

                post = []
                if d == self.dimension-1:
                    # inner loop
                    # first time step
                    t0 = 0
                    sol = field.sol.subs(self.t, t0)
                    if self.read:
                        sol = sol.subs(self.media_dict)
                        sol = self.resolve_media_params(sol)
                    for idx in self.index:
                        sol = sol.subs(idx, '_'+idx.name)
                    body = [cgen.Assign(ccode(field[[0]+loop]), ccode(sol))]
                body = pre + body + post
                body = [cgen.For(cgen.InlineInitializer(cgen.Value('int', i), i0), cgen.Line('%s<%s' % (i, i1)), cgen.Line('++%s' % i), cgen.Block(body))]

            statements.append(body[0])
            statements += self.generate_second_initialisation()
        return str(cgen.Module(statements))
    
    def generate_second_initialisation(self):
        loop = [Symbol('_'+x.name) for x in self.index]  # symbols for loop
        m = self.margin.value
        statements = []
        v=symbols("v")
        for field in self.fields:
            body = []
            if self.omp:
                statements.append(cgen.Pragma('omp for schedule(static,1)'))
            # populate xvalue, yvalue zvalue code
            for d in range(self.dimension-1, -1, -1):
                i = loop[d]
                i0 = m
                
                i1 = ccode(self.dim[d]-m)
                expr = self.spacing[d]*(loop[d])
                pre = []

                post = []
                if d == self.dimension-1:
                    # inner loop
                    # first time step
                    t0 = 0
                    kernel = self.transform_kernel(field)
                    for arg in kernel.args:
                        if str(arg).startswith("-") and str(self.t - 1) in str(arg):
                            kernel = kernel.subs({arg: 0}, simultaneous=True)
                            arg = 2*v*self.dt
                            kernel = 0.5*(kernel + arg)
                    
                    kernel = kernel.subs({self.t:self.time[0]})
                    for idx in self.index:
                        kernel = kernel.subs(idx, '_'+idx.name)
                    
                    
                    
                    body = [cgen.Assign(ccode(field[[self.time[1]]+loop]), ccode(kernel))]
                body = pre + body + post
                body = [cgen.For(cgen.InlineInitializer(cgen.Value('int', i), i0), cgen.Line('%s<%s' % (i, i1)), cgen.Line('++%s' % i), cgen.Block(body))]

            statements.append(body[0])
        return statements
    def generate_loop(self, fields):
        """
        The functions to generate stress loops and velocity loops are identical,
        save for a single parameter. Moved the common code to this function to reduce repetition of code.
        """
        if self.eval_const:
            self.create_const_dict()
        m = self.margin.value
        body = []
        for d in range(self.dimension-1, -1, -1):
            i = self.index[d]
            i0 = m
            i1 = ccode(self.dim[d]-m)
            if d == self.dimension-1:
                # inner loop
                if not self.fission:
                    body = self.simple_kernel(fields, [d, i, i0, i1])
                else:
                    body = self.fission_kernel(fields, [d, i, i0, i1])
            if not d == self.dimension-1:
                body = [cgen.For(cgen.InlineInitializer(cgen.Value('int', i), i0), cgen.Line('%s<%s' % (i, i1)), cgen.Line('++%s' % i), cgen.Block(body))]

        if not self.pluto and self.omp:
            body.insert(0, cgen.Pragma('omp for schedule(static,1)'))
        return str(cgen.Module(body))
    
    @property
    def primary_loop(self):
        return self.generate_loop(self.fields)
    
    def kernel_sympy(self, field):
        kernel = self.transform_kernel(field)
        kernel = kernel.xreplace({self.t+1: self.time[2], self.t: self.time[1], self.t-1:self.time[0]})
        return kernel
    
    def simple_kernel(self, grid_field, indexes):
        """
        Generate the inner loop with all fields from stress or velocity
        :param grid_field: stress or velocity field array
        :param indexes: array with dimension, dimension var, initial margin, final margin
        - iterate through fields and replace mako template
        - return inner loop code as string
        """
        body = []
        idx = [self.time[len(self.time)-1]] + self.index
        # This loop create the most inner loop with all fields
        for field in grid_field:
            body.append(cgen.Assign(ccode(field[idx]), ccode(self.kernel_sympy(field))))
        body = [cgen.For(cgen.InlineInitializer(cgen.Value('int', indexes[1]), indexes[2]), cgen.Line('%s<%s' % (indexes[1], indexes[3])), cgen.Line('++%s' % indexes[1]), cgen.Block(body))]
        if not self.pluto and self.ivdep and indexes[0] == self.dimension-1:
            body.insert(0, self.compiler._ivdep)
        if not self.pluto and self.simd and indexes[0] == self.dimension-1:
            body.insert(0, cgen.Pragma('simd'))

        return body
    
    @property
    def free_memory(self):
        """
       - generate code for free allocated memory
       - return the generated code as string
       """
        statements = []
        for field in self.fields:
            # alloc aligned memory (on windows and linux)
            ifdef = cgen.IfDef('_MSC_VER', [cgen.Statement('_aligned_free(grid->%s)' % (ccode(field.label)))],
                               [cgen.Statement('free(grid->%s)' % (ccode(field.label)))])
            statements.append(ifdef)

        return str(cgen.Module(statements))


