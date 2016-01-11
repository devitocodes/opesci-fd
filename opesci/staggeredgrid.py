from grid import Grid
from variable import Variable
from fields import Media
from codeprinter import ccode, ccode_eq
from derivative import DDerivative
from util import *
from compilation import get_package_dir

from sympy import Symbol, Rational, solve, expand, Eq
from mako.lookup import TemplateLookup
import mmap
import cgen_wrapper as cgen
from os import path
from __builtin__ import str

__all__ = ['StaggeredGrid']

hf = Rational(1, 2)  # 1/2


class StaggeredGrid(Grid):
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

    Switches can also be set with set_switches(switch1=va1, switch2=val2, ...).
    Supported switches are:
    * omp: insert #pragma omp for before outer loops, default True
    * ivdep: insert #praga ivdep before inner loop, default True
    * simd: insert #pragma simd before inner loop, default False
    * double: use float (False) or double (True) for real numbers, default False
    * io: include header files for io (e.g. vtk support), default False
    * read: whether to read media parameters from input file, default False
    * expand: expand kernel fully (no factorisation), default True
    * eval_const: evaluate all constants in kernel in generated code default True
    * output_vts: Output solution fields at every timestep
    * converge: Generate code for computing analutical solution and L2 norms
    * profiling: Generate code for gathering profiling information via PAPI
    * pluto: Define scop for pluto optimisation
    * fission: Define loop fission optimisation
    """
    template_base = 'staggered3d_tmpl.cpp'

    template_keys = ['pluto', 'io', 'profiling', 'numevents_papi',
                     'time_stepping', 'define_constants', 'declare_fields',
                     'define_fields', 'store_fields', 'load_fields',
                     'initialise', 'initialise_bc', 'stress_loop',
                     'velocity_loop', 'stress_bc', 'velocity_bc', 'output_step',
                     'define_convergence', 'converge_test', 'print_convergence',
                     'define_profiling', 'define_papi_events', 'sum_papi_events']

    _switches = ['omp', 'ivdep', 'simd', 'double', 'expand', 'eval_const',
                 'output_vts', 'converge', 'profiling', 'pluto', 'fission']

    _papi_events = []

    def __init__(self, dimension, index=None, domain_size=None, grid_size=None,
                 time_step=None, stress_fields=None, velocity_fields=None,
                 omp=True, ivdep=True, simd=False, double=False, io=False,
                 expand=True, eval_const=True, output_vts=False,
                 converge=False, profiling=False, pluto=False, fission=False):

        self.dimension = dimension

        template_dir = path.join(get_package_dir(), "templates")
        staggered_dir = path.join(get_package_dir(), "templates/staggered")
        self.lookup = TemplateLookup(directories=[template_dir, staggered_dir])

        # List of associated fields
        self.sfields = []
        self.vfields = []

        # Switches
        self.pluto = pluto
        self.omp = omp
        self.ivdep = ivdep
        self.simd = simd
        self.double = double
        self.expand = expand
        self.eval_const = eval_const
        self.real_t = 'double' if self.double else 'float'
        self.output_vts = output_vts
        self.converge = converge
        self.profiling = profiling
        self.fission = fission

        # number of ghost cells for boundary
        self.margin = Variable('margin', 2, 'int', True)
        self.grid_size = (10, 10, 10)
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

        self.t = Symbol('_t')
        self.set_index(index)

        # default 2nd order in time, 4th order in space, i.e. (2,4) scheme
        default_order = [2] + [4]*self.dimension
        self.set_order(default_order)

        self.dt = Variable('dt', 0.01, self.real_t, True)
        self.ntsteps = Variable('ntsteps', 100, 'int', True)
        self.alignment = mmap.PAGESIZE  # default alignment for malloc

        # Optional further grid settings
        if stress_fields:
            self.set_stress_fields(stress_fields)
        if velocity_fields:
            self.set_velocity_fields(velocity_fields)
        if domain_size:
            self.set_domain_size(domain_size)
        if grid_size:
            self.set_grid_size(grid_size)

        # user defined variables
        # use dictionary because of potential duplication
        self.defined_variable = {}
        self.update_field_order()
        self._update_spacing()

    @property
    def fields(self):
        return self.sfields + self.vfields

    @property
    def io(self):
        """Flag whether to include I/O headers"""
        return self.read or self.output_vts

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
                raise ValueError(str(x) + ' is not a valid order (require even integer)')
        self.order = order
        # periodicity for time stepping
        self.tp = Variable('tp', self.order[0], 'int', True)
        # add time variables for time stepping: t0, t1 ...
        self.time = []
        for k in range(self.order[0]):
            name = str(self.t) + str(k)
            v = Variable(name, 0, 'int', False)
            self.time.append(v)
        self.margin.value = self.order[1]/2
        self.set_grid_size(self.grid_size)
        self.update_field_order()

    def update_field_order(self):
        """
        update the order of acuracy of the fields
        """
        if hasattr(self, 'sfields') and hasattr(self, 'vfields'):
            for field in self.sfields + self.vfields:
                field.set_order(self.order)

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

        if hasattr(self, 'sfields') and hasattr(self, 'vfields'):
            for field in self.sfields + self.vfields:
                # update the fields with new indices
                field.set_indices([self.t]+self.index)

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

    def set_alignment(self, alignment):
        """
        set alignment size to be used for malloc alignment
        :param alignment: new alignment size in bytes
        """
        self.alignment = alignment

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
        if self.order[1] == 2:
            return h/Vp/(3**0.5)
        elif self.order[1] == 4:
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
        for field in self.sfields:
            field.set_spacing(variable_to_symbol([self.dt]+self.spacing))

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
        for field in self.vfields:
            field.set_spacing(variable_to_symbol([self.dt]+self.spacing))

    def calc_derivatives(self):
        """
        populate field.d lists with Derivative objects
        """
        for field in self.sfields+self.vfields:
            field.populate_derivatives(max_order=1)  # VS scheme only requare 1st derivatives

    def get_all_variables(self):
        """
        return list of all variables defined
        """
        variables = self.dim + self.spacing + self.time\
            + [self.tp, self.dt, self.margin, self.ntsteps, self.vec_size]\
            + self.defined_variable.values()
        return variables

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
        # save the equation
        for field, eq in zip(self.vfields+self.sfields, equations):
            if self.read:
                eq = eq.subs(self.media_dict)
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
        index = [t1] + self.index

        simplify = True if max(self.order[1:]) <= 4 else False

        for field, eq in zip(self.vfields+self.sfields, eqs):
            # want the LHS of express to be at time t+1
            kernel = solve(eq, field[index], simplify=simplify)[0]
            kernel = kernel.subs({t: t+hf-(self.order[0]/2-1)})

            field.set_fd_kernel(kernel)

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

    def transform_bc(self, field, dimension, side):
        """
        transform the bc[dimension][side] item of field based on setting of grid, such as expand, eval_const, read
        """
        bc_new_list = []
        for bc in field.bc[dimension][side]:
            bc_new = bc
            if self.expand:
                # expanding kernel (de-factorisation)
                bc_new = expand(bc_new)

            if self.eval_const:
                # substitute constants with values
                bc_new = bc_new.subs(self.const_dict)

            bc_new_list.append(bc_new)
        return bc_new_list

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

        # all stress fields
        fields = [None] + [f for f in self.sfields if f.direction[0] == f.direction[1]]
        for s in fields[1:]:
            s.associate_stress_fields(fields)

    def set_free_surface_boundary(self, dimension, side):
        """
        set free surface boundary condition
        calls set_free_surface() method of SField and VField
        :param dimension: the normal direction to identify the surface
        e.g. dimension=1 for y-z place
        :param side: the side of the surface
        side=0 for bottom surface, side=1 for top surface
        """
        algo = 'robertsson'
        if self.order[dimension] == 4:
            # using different algorithm for free surface for 4th order
            algo = 'levander'
        self.associate_fields()
        for field in self.sfields+self.vfields:
            if side == 0:
                field.set_free_surface(dimension, self.margin.value, side, algo=algo)
            else:
                field.set_free_surface(dimension, self.dim[dimension-1]-self.margin.value-1, side, algo=algo)

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
            # field to store rho
            self.rho = Media('rho', dimension=3, staggered=[False, False, False], index=self.index)
            # field to store vp
            self.vp = Media('vp', dimension=3, staggered=[False, False, False], index=self.index)
            # field to store vs
            self.vs = Media('vs', dimension=3, staggered=[False, False, False], index=self.index)

            # calculated from input data
            # list of fields to store beta (buoyancy and effective buoyancy)
            self.beta = [Media('beta', dimension=3, staggered=[False, False, False], index=self.index),
                         Media('beta1', dimension=3, staggered=[False, False, False], index=self.index),
                         Media('beta2', dimension=3, staggered=[False, False, False], index=self.index),
                         Media('beta3', dimension=3, staggered=[False, False, False], index=self.index)]
            # field to store lambda
            self.lam = Media('lambda', dimension=3, staggered=[False, False, False], index=self.index)
            # list of fields to store mu
            # (shear modulus and effective shear modulus)
            self.mu = [Media('mu', dimension=3, staggered=[False, False, False], index=self.index),
                       Media('mu12', dimension=3, staggered=[False, False, False], index=self.index),
                       Media('mu13', dimension=3, staggered=[False, False, False], index=self.index),
                       Media('mu23', dimension=3, staggered=[False, False, False], index=self.index)]
            self.create_media_dict()

        else:
            self.set_variable('rho', rho, 'float', True)
            self.set_variable('beta', 1.0/rho, 'float', True)
            self.set_variable('lambda', rho*(vp**2-2*vs**2), 'float', True)
            self.set_variable('mu', rho*(vs**2), 'float', True)

    def create_media_dict(self):
        """
        create the dictionary of media parameters (for reading data)
        the parameter symbols in PDEs are substituted with Media objects
        """
        self.media_dict = {Symbol('beta'): self.beta[0][self.index],
                           Symbol('lambda'): self.lam[self.index],
                           Symbol('mu'): self.mu[0][self.index]}

    def resolve_individual_media_param(self, expr):
        """
        map individual media parameter to the effective media parameter
        """
        idx = list(expr.indices)
        if str(expr.base.label) == 'beta':
            if is_half(idx[0]):
                # replace with beta1
                idx[0] -= hf
                return self.beta[1][idx]
            elif is_half(idx[1]):
                # replace with beta2
                idx[1] -= hf
                return self.beta[2][idx]
            elif is_half(idx[2]):
                # replace with beta1
                idx[2] -= hf
                return self.beta[3][idx]
            else:
                return expr
        if str(expr.base.label) == 'mu':
            if is_half(idx[0]):
                idx[0] -= hf
                if is_half(idx[1]):
                    # replace with mu12
                    idx[1] -= hf
                    return self.mu[1][idx]
                else:
                    # replace with m13
                    idx[2] -= hf
                    return self.mu[2][idx]
            elif is_half(idx[1]):
                # replace with mu23
                idx[1] -= hf
                idx[2] -= hf
                return self.mu[3][idx]
            else:
                return expr
        return expr

    def resolve_media_params(self, expr):
        """
        resolve the media parameters in a kernel by replacing with effective media parameters
        e.g. mu[x,y+1/2,z+1/2] replaced with mu23[x,y,z]
        """
        if expr.is_Equality:
            rhs = self.resolve_media_params(expr.rhs)
            return Eq(expr.lhs, rhs)
        if expr.is_Symbol or expr.is_Number:
            return expr
        if isinstance(expr, Indexed):
            b = expr.base
            if not (isinstance(b, Media)):
                return expr
            else:
                return self.resolve_individual_media_param(expr)

        args = tuple([self.resolve_media_params(arg) for arg in expr.args])
        result = expr.func(*args)
        return result

    def get_velocity_kernel_ai(self):
        """
        - get the arithmetic intensity of velocity kernel
        - get the number of different operations of the velocity field kernel
        - types of operations are ADD (inc -), MUL (inc /), LOAD, STORE
        - #LOAD = number of unique fields in the kernel
        - return tuple (AI, AI_w, #ADD, #MUL, #LOAD, #STORE)
        - arithmetic intensity AI = (ADD+MUL)/[(LOAD+STORE)*word size]
        - weighted AI, AI_w = (ADD+MUL)/(2*Max(ADD,MUL)) * AI
        """
        if self.eval_const:
            self.create_const_dict()
        store = 0
        add = 0
        mul = 0
        arrays = []  # to store name of arrays loaded
        for field in self.vfields:
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

    def get_stress_kernel_ai(self):
        """
        - get the arithmetic intensity of stress kernel
        - get the number of different operations of the stress field kernel
        - types of operations are ADD (inc -), MUL (inc /), LOAD, STORE
        - #LOAD = number of unique fields in the kernel
        - return tuple (#ADD, #MUL, #LOAD, #STORE)
        - arithmetic intensity AI = (ADD+MUL)/[(LOAD+STORE)*word size]
        - weighted AI, AI_w = (ADD+MUL)/(2*Max(ADD,MUL)) * AI
        """
        if self.eval_const:
            self.create_const_dict()
        store = 0
        add = 0
        mul = 0
        arrays = []  # to store name of arrays loaded
        for field in self.sfields:
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

    def get_velocity_bc_ai(self):
        """
        - get the arithmetic intensity of velocity boundary conditions
        - get the number of different operations of the velocity boundary condition computation
        - return list of AIs and # ops for each boundary
        - types of operations are ADD (inc -), MUL (inc /), LOAD, STORE
        - #LOAD = number of unique fields in the kernel
        - return tuple (#ADD, #MUL, #LOAD, #STORE)
        - arithmetic intensity AI = (ADD+MUL)/[(LOAD+STORE)*word size]
        - weighted AI, AI_w = (ADD+MUL)/(2*Max(ADD,MUL)) * AI
        """
        result = []
        if self.eval_const:
            self.create_const_dict()

        # 8 byte if double, 4 if float used
        # media parameter fields are always float, need to amend this
        word_size = 8 if self.double else 4

        for field in self.vfields:
            for dimension in range(1, 4):
                for side in range(2):
                    store = 0
                    add = 0
                    mul = 0
                    load = 0
                    arrays = []  # to store name of arrays loaded
                    bc_list = self.transform_bc(field, dimension, side)
                    for bc in bc_list:
                        store += 1  # increment STORE by 1 (assignment)
                        add2, mul2, arrays2 = get_ops_expr(bc.rhs, arrays)
                        add += add2  # accumulate # ADD
                        mul += mul2  # accumulate # MUL
                        arrays = arrays2  # replace with new list of field names
                    load = len(arrays)
                    if (store == 0):
                        ai = 0
                        weight = 0  # weight of AI in overall AI calculation
                    else:
                        ai = float(add+mul)/(load+store)/word_size
                        weight = 1.0/(self.dim[dimension-1].value-self.margin.value*2)
                        if (add == 0 and mul == 0):
                            ai_w = ai
                        else:
                            ai_w = ai*(add+mul)/max(add, mul)/2.0
                    result.append({'weight': weight, 'ai': ai, 'ai_w': ai_w, 'add': add, 'mul': mul, 'load': load, 'store': store})
        return result

    def get_stress_bc_ai(self):
        """
        - get the arithmetic intensity of stress boundary conditions
        - get the number of different operations of the stress boundary condition computation
        - return list of AIs and # ops for each boundary
        - types of operations are ADD (inc -), MUL (inc /), LOAD, STORE
        - #LOAD = number of unique fields in the kernel
        - return tuple (#ADD, #MUL, #LOAD, #STORE)
        - arithmetic intensity AI = (ADD+MUL)/[(LOAD+STORE)*word size]
        - weighted AI, AI_w = (ADD+MUL)/(2*Max(ADD,MUL)) * AI
        """
        result = []
        if self.eval_const:
            self.create_const_dict()

        # 8 byte if double, 4 if float used
        # media parameter fields are always float, need to amend this
        word_size = 8 if self.double else 4

        for field in self.sfields:
            for dimension in range(1, 4):
                for side in range(2):
                    store = 0
                    add = 0
                    mul = 0
                    load = 0
                    arrays = []  # to store name of arrays loaded
                    bc_list = self.transform_bc(field, dimension, side)
                    for bc in bc_list:
                        store += 1  # increment STORE by 1 (assignment)
                        add2, mul2, arrays2 = get_ops_expr(bc.rhs, arrays)
                        add += add2  # accumulate # ADD
                        mul += mul2  # accumulate # MUL
                        arrays = arrays2  # replace with new list of field names
                    load = len(arrays)
                    if (store == 0):
                        ai = 0
                        weight = 0  # weight of AI in overall AI calculation
                    else:
                        ai = float(add+mul)/(load+store)/word_size
                        weight = 1.0/(self.dim[dimension-1].value-self.margin.value*2)
                        if (add == 0 and mul == 0):
                            ai_w = ai
                        else:
                            ai_w = ai*(add+mul)/max(add, mul)/2.0
                    result.append({'weight': weight, 'ai': ai, 'ai_w': ai_w, 'add': add, 'mul': mul, 'load': load, 'store': store})
        return result

    def get_overall_kernel_ai(self):
        """
        - get the overall arithmetic intensity of the kernel (velocity and stress)
        - arithmetic intensity AI = (ADD+MUL)/[(LOAD+STORE)*word size]
        - weighted AI, AI_w = (ADD+MUL)/(2*Max(ADD,MUL)) * AI
        """
        # get the AI of kernels and boundary conditions
        velocity_ai = self.get_velocity_kernel_ai()
        stress_ai = self.get_stress_kernel_ai()
        velocity_bc_ai = self.get_velocity_bc_ai()
        stress_bc_ai = self.get_stress_bc_ai()
        total_weight = 2.0  # velocity kernel and stress kernel
        overall_ai = velocity_ai[0] + stress_ai[0]
        overall_ai_w = velocity_ai[1] + stress_ai[1]

        for ai in velocity_bc_ai+stress_bc_ai:
            total_weight += ai['weight']
            overall_ai += ai['ai']*ai['weight']
            overall_ai_w += ai['ai_w']*ai['weight']

        overall_ai /= total_weight
        overall_ai_w /= total_weight

        # calculate adjustment due to ghost cells
        ghost_adj = 1.0
        for d in self.dim[1:]:
            ghost_adj *= 1 - float(self.margin.value)/d.value
        overall_ai *= ghost_adj
        overall_ai_w *= ghost_adj
        return overall_ai, overall_ai_w

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
        vsize *= self.order[0]
        statements = []
        for field in self.sfields + self.vfields:
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
        if self.read:
            # add code to read data
            result += self.read_data()
        return str(cgen.Module(result))

    def read_data(self):
        """
        - generate code for reading data (rho, Vp, Vs) from input files
        - calculate effective media parameters beta, lambda, mu from the data
        """
        statements = []
        if self.read:
            arr = ''  # =[dim2][dim3]...
            for d in self.dim[1:]:
                arr += '[' + d.name + ']'
            vsize = 1
            for d in self.dim:
                vsize *= d.value
            # declare fields to read physical parameters from file
            # always use float not double
            loop = [self.rho, self.vp, self.vs] + self.beta + [self.lam] + self.mu
            for field in loop:
                vec = "_%s_vec" % ccode(field.label)
                vec_value = cgen.Pointer(cgen.Value(self.real_t, vec))
                # alloc aligned memory (on windows and linux)
                statements.append(vec_value)
                ifdef = cgen.IfDef('_MSC_VER', [cgen.Assign(vec, '(%s*) _aligned_malloc(%d * sizeof(%s), %d)' % (self.real_t, vsize, self.real_t, self.alignment))],
                                   [cgen.Statement('posix_memalign((void **)(&%s), %d, %d*sizeof(%s))' % (vec, self.alignment, vsize, self.real_t))])
                statements.append(ifdef)
                cast_pointer = cgen.Initializer(cgen.Value(self.real_t, "(*%s)%s" % (ccode(field.label), arr)), '(%s (*)%s) %s' % (self.real_t, arr, vec))
                statements.append(cast_pointer)
            # read from file
            statements.append(cgen.Statement('opesci_read_simple_binary_ptr("%s", _%s_vec, %d)' % (self.rho_file, self.rho.label, vsize)))
            statements.append(cgen.Statement('opesci_read_simple_binary_ptr("%s", _%s_vec, %d)' % (self.vp_file, self.vp.label, vsize)))
            statements.append(cgen.Statement('opesci_read_simple_binary_ptr("%s", _%s_vec, %d)' % (self.vs_file, self.vs.label, vsize)))
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
            kernel = cgen.Assign(ccode(self.beta[0][idx]), ccode(1.0/self.rho[idx]))
            statements.append(self.simple_loop(kernel))
            # beta1 (effective bouyancy in x direction)
            kernel = cgen.Assign(ccode(self.beta[1][idx]), ccode((self.beta[0][idx] + self.beta[0][idx100])/2.0))
            statements.append(self.simple_loop(kernel))
            # beta2 (effective bouyancy in y direction)
            kernel = cgen.Assign(ccode(self.beta[2][idx]), ccode((self.beta[0][idx] + self.beta[0][idx010])/2.0))
            statements.append(self.simple_loop(kernel))
            # beta3 (effective bouyancy in z direction)
            kernel = cgen.Assign(ccode(self.beta[3][idx]), ccode((self.beta[0][idx] + self.beta[0][idx001])/2.0))
            statements.append(self.simple_loop(kernel))
            # lambda
            kernel = cgen.Assign(ccode(self.lam[idx]), ccode(self.rho[idx]*(self.vp[idx]**2-2*self.vs[idx]**2)))
            statements.append(self.simple_loop(kernel))
            # mu
            kernel = cgen.Assign(ccode(self.mu[0][idx]), ccode(self.rho[idx]*(self.vs[idx]**2)))
            statements.append(self.simple_loop(kernel))
            # mu12 (effective shear modulus for shear stress sigma_xy)
            kernel = cgen.Assign(ccode(self.mu[1][idx]), ccode(1.0/(0.25*(1.0/self.mu[0][idx]+1.0/self.mu[0][idx100] + 1.0/self.mu[0][idx010]+1.0/self.mu[0][idx110]))))
            statements.append(self.simple_loop(kernel))
            # mu13 (effective shear modulus for shear stress sigma_xz)
            kernel = cgen.Assign(ccode(self.mu[2][idx]), ccode(1.0/(0.25*(1.0/self.mu[0][idx]+1.0/self.mu[0][idx100] + 1.0/self.mu[0][idx001]+1.0/self.mu[0][idx101]))))
            statements.append(self.simple_loop(kernel))
            # mu23 (effective shear modulus for shear stress sigma_yz)
            kernel = cgen.Assign(ccode(self.mu[3][idx]), ccode(1.0/(0.25*(1.0/self.mu[0][idx]+1.0/self.mu[0][idx010] + 1.0/self.mu[0][idx001]+1.0/self.mu[0][idx011]))))
            statements.append(self.simple_loop(kernel))
        return statements

    def simple_loop(self, kernel):
        """
        - helper function to generate simple nested loop over the entire domain
        (not including ghost cells) with kernel at the inner loop
        - variables defined in self.index are used as loop variables
        """
        result = kernel
        m = self.margin.value
        for d in range(self.dimension-1, -1, -1):
            result = cgen.For(cgen.InlineInitializer(cgen.Value('int', self.index[d]), m), cgen.Line('%s<%s' % (self.index[d], ccode(self.dim[d]-m))), cgen.Line('++%s' % self.index[d]), result)
        return result

    @property
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
        m = self.margin.value
        loop = [Symbol('_'+x.name) for x in self.index]  # symbols for loop

        statements = []
        for field in self.sfields+self.vfields:
            body = []
            if self.omp:
                statements.append(cgen.Pragma('omp for schedule(static,1)'))
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
                pre = [cgen.Initializer(cgen.Value(self.real_t, self.index[d].name), ccode(expr))]

                post = []
                if d == self.dimension-1:
                    # inner loop
                    # first time step
                    t0 = self.dt.value/2 if field.staggered[0] else 0
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
        return str(cgen.Module(statements))

    def simple_kernel(self, grid_field, indexes):
        """
        Generate the inner loop with all fields from stress or velocity
        :param grid_field: stress or velocity field array
        :param indexes: array with dimension, dimension var, initial margin, final margin
        - iterate through fields and replace mako template
        - return inner loop code as string
        """
        body = []
        idx = [self.time[1]] + self.index
        # This loop create the most inner loop with all fields
        for field in grid_field:
            kernel = self.transform_kernel(field)
            if self.read:
                kernel = self.resolve_media_params(kernel)
            body.append(cgen.Assign(ccode(field[idx]), ccode(kernel.xreplace({self.t+1: self.time[1], self.t: self.time[0]}))))
        body = [cgen.For(cgen.InlineInitializer(cgen.Value('int', indexes[1]), indexes[2]), cgen.Line('%s<%s' % (indexes[1], indexes[3])), cgen.Line('++%s' % indexes[1]), cgen.Block(body))]
        if not self.pluto and self.ivdep and indexes[0] == self.dimension-1:
            body.insert(0, cgen.Statement(self.compiler._ivdep))
        if not self.pluto and self.simd and indexes[0] == self.dimension-1:
            body.insert(0, cgen.Pragma('simd'))

        return body

    def fission_kernel(self, grid_field, indexes):
        """
        Generate the inner loop with all fields from stress or velocity
        :param grid_field: stress or velocity field array
        :param indexes: array with dimension, dimension var, initial margin, final margin
        - iterate through fields and for each dimension separate minus, plus and unitary strides
        on its own loop replacing it on mako template
        - return inner loop code as string
        """
        body = []
        body_tmp = []
        operator = ['=', '+=']
        idx = [self.time[1]] + self.index
        for field in grid_field:
            remainder_kernel = ()
            operator_idx = 0
            kernel = self.transform_kernel(field)
            if self.read:
                kernel = self.resolve_media_params(kernel)
            kernel_args = kernel.args
            for dim in range(self.dimension-1, -1, -1):
                dimension = self.index[dim]
                kernel_stmt_pos = kernel
                kernel_stmt_neg = kernel
                # For each dimension in each field iterate through its expressions to separate
                # positive and negative strides
                for arg in kernel_args:
                    if not (str(dimension) + " -" in str(arg)):
                        kernel_stmt_neg = kernel_stmt_neg.subs({arg: 0}, simultaneous=True)
                    if not (str(dimension) + " +" in str(arg)):
                        kernel_stmt_pos = kernel_stmt_pos.subs({arg: 0}, simultaneous=True)
                remainder_kernel += kernel_stmt_pos.args
                remainder_kernel += kernel_stmt_neg.args
                # Create the inner loop for with negative strides expressions
                if not (len(kernel_stmt_neg.args) == 0):
                    body_tmp = [cgen.Statement(ccode(field[idx]) + operator[operator_idx] + ccode(kernel_stmt_neg.xreplace({self.t+1: self.time[1], self.t: self.time[0]})))]
                    body_tmp = [cgen.For(cgen.InlineInitializer(cgen.Value('int', indexes[1]), indexes[2]), cgen.Line('%s<%s' % (indexes[1], indexes[3])), cgen.Line('++%s' % indexes[1]), cgen.Block(body_tmp))]
                    if not self.pluto and self.ivdep and indexes[0] == self.dimension-1:
                        body_tmp.insert(0, cgen.Statement(self.compiler._ivdep))
                    if not self.pluto and self.simd and indexes[0] == self.dimension-1:
                        body_tmp.insert(0, cgen.Pragma('simd'))
                    body = body + body_tmp
                    operator_idx = 1
                # Create the inner loop for with positive strides expressions
                if not (len(kernel_stmt_pos.args) == 0):
                    body_tmp = [cgen.Statement(ccode(field[idx]) + operator[operator_idx] + ccode(kernel_stmt_pos.xreplace({self.t+1: self.time[1], self.t: self.time[0]})))]
                    body_tmp = [cgen.For(cgen.InlineInitializer(cgen.Value('int', indexes[1]), indexes[2]), cgen.Line('%s<%s' % (indexes[1], indexes[3])), cgen.Line('++%s' % indexes[1]), cgen.Block(body_tmp))]
                    if not self.pluto and self.ivdep and indexes[0] == self.dimension-1:
                        body_tmp.insert(0, cgen.Statement(self.compiler._ivdep))
                    if not self.pluto and self.simd and indexes[0] == self.dimension-1:
                        body_tmp.insert(0, cgen.Pragma('simd'))
                    body = body + body_tmp
                    operator_idx = 1
            # Create the inner loop for unit strided array access
            kernel_stmt = kernel
            for arg in remainder_kernel:
                kernel_stmt = kernel_stmt.subs({arg: 0}, simultaneous=True)
            body_tmp = [cgen.Statement(ccode(field[idx]) + '+=' + ccode(kernel_stmt.xreplace({self.t+1: self.time[1], self.t: self.time[0]})))]
            body_tmp = [cgen.For(cgen.InlineInitializer(cgen.Value('int', indexes[1]), indexes[2]), cgen.Line('%s<%s' % (indexes[1], indexes[3])), cgen.Line('++%s' % indexes[1]), cgen.Block(body_tmp))]
            if not self.pluto and self.ivdep and indexes[0] == self.dimension-1:
                body_tmp.insert(0, cgen.Statement(self.compiler._ivdep))
            if not self.pluto and self.simd and indexes[0] == self.dimension-1:
                body_tmp.insert(0, cgen.Pragma('simd'))
            body = body + body_tmp

        return body

    @property
    def stress_loop(self):
        """
        generate code for stress field update loop
        - loop through stress fields to generate code of computation kernel
        - generate inner loop by inserting kernel into Mako template
        - recursive insertion to generate nested loop
        return generated code as string
        """
        return self.generate_loop(self.sfields)

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
    def velocity_loop(self):
        """
        generate code for velocity field update loop
        - loop through velocity fields to generate code of computation kernel
        - generate inner loop by inserting kernel into Mako template
        - recursive insertion to generate nested loop
        return generated code as string
        """
        return self.generate_loop(self.vfields)

    @property
    def stress_bc(self):
        return self.stress_bc_getter()

    def stress_bc_getter(self, init=False):
        """
        generate code for updating stress field boundary ghost cells
        - generate inner loop by inserting boundary code (saved in field.bc)
        - recursive insertion to generate nested loop
        - loop through all stress fields and sides
        - if init=True (initialisation), no need to generate code to overwrite Txx, Tyy, Tzz
        return generated code as string
        """
        result = []
        body = []
        if self.eval_const:
            self.create_const_dict()
        for field in self.sfields:
            # normal stress, not shear stress
            normal = field.direction[0] == field.direction[1]
            for d in range(self.dimension):
                if init and normal and (not field.direction[0] == d+1):
                    # skip re-calc Txx, Tyy at z plane etc, when initialisation
                    continue
                for side in range(2):
                    # skip if this boundary calculation is not needed
                    if field.bc[d+1][side] == []:
                        continue
                    if self.omp:
                        result += [cgen.Pragma('omp for schedule(static,1)')]
                    body = []
                    for d2 in range(self.dimension-1, -1, -1):
                        # loop through other dimensions
                        if not d2 == d:
                            i = self.index[d2]
                            if normal:
                                if field.direction[0] == d+1:
                                    # Txx in x plane
                                    i0 = 0
                                    i1 = self.dim[d2]
                                else:
                                    # Txx in y, z plane
                                    i0 = self.margin.value+1
                                    i1 = self.dim[d2]-self.margin.value-1
                            else:
                                i0 = 0
                                i1 = self.dim[d2]

                            if not body:
                                # inner loop, populate ghost cell calculation
                                bc_list = self.transform_bc(field, d+1, side)
                                if self.read:
                                    body = [cgen.Statement(ccode_eq(self.resolve_media_params(bc)).replace('[_t + 1]', '[_t1]').replace('[_t]', '[_t0]')) for bc in bc_list]
                                else:
                                    body = [cgen.Statement(ccode_eq(bc).replace('[_t + 1]', '[_t1]').replace('[_t]', '[_t0]')) for bc in bc_list]
                                body = [cgen.For(cgen.InlineInitializer(cgen.Value('int', i), i0), cgen.Line('%s<%s' % (i, i1)), cgen.Line('++%s' % i), cgen.Block(body))]
                                if self.ivdep:
                                    body.insert(0, cgen.Pragma('ivdep'))
                                if self.simd:
                                    body.insert(0, cgen.Pragma('simd'))
                            else:
                                body = [cgen.For(cgen.InlineInitializer(cgen.Value('int', i), i0), cgen.Line('%s<%s' % (i, i1)), cgen.Line('++%s' % i), cgen.Block(body))]
                    result += body
        return cgen.Module(result)

    @property
    def velocity_bc(self):
        """
        generate code for updating stress field boundary ghost cells
        - generate inner loop by inserting boundary code (saved in field.bc)
        - recursive insertion to generate nested loop
        - loop through all velocity fields and sides
        return generated code as string
        """
        result = []
        if self.eval_const:
            self.create_const_dict()
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
                        result += [cgen.Pragma('omp for schedule(static,1)')]
                    body = []
                    for d2 in range(self.dimension-1, -1, -1):
                        # loop through other dimensions
                        if not d2 == d:
                            i = self.index[d2]
                            i0 = 1
                            i1 = self.dim[d2]-1
                            if not body:
                                # inner loop, populate ghost cell calculation
                                # body = field.bc[d+1][side]
                                bc_list = self.transform_bc(field, d+1, side)
                                if self.read:
                                    body = [cgen.Statement(ccode_eq(self.resolve_media_params(bc)).replace('[_t + 1]', '[_t1]').replace('[_t]', '[_t0]')) for bc in bc_list]
                                else:
                                    body = [cgen.Statement(ccode_eq(bc).replace('[_t + 1]', '[_t1]').replace('[_t]', '[_t0]')) for bc in bc_list]
                                body = [cgen.For(cgen.InlineInitializer(cgen.Value('int', i), i0), cgen.Line('%s<%s' % (i, i1)), cgen.Line('++%s' % i), cgen.Block(body))]
                                if self.ivdep:
                                    body.insert(0, cgen.Pragma('ivdep'))
                                if self.simd:
                                    body.insert(0, cgen.Pragma('simd'))
                            else:
                                body = [cgen.For(cgen.InlineInitializer(cgen.Value('int', i), i0), cgen.Line('%s<%s' % (i, i1)), cgen.Line('++%s' % i), cgen.Block(body))]

                    result += body

        return cgen.Module(result)

    @property
    def initialise_bc(self):
        """
        - generate code for initialisation of boundary ghost cells
        - generate generic boundary cell code
        replace array indices [t] with [0]
        - return generated code as string
        """
        t1 = "'["+str(self.t)+"1]'"
        rep = "'[0]'"
        result = [cgen.replace_in_code(self.stress_bc_getter(init=True), t1, rep)]
        result += [cgen.replace_in_code(self.velocity_bc, t1, rep)]
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
    def output_step(self):
        """
        - generate code for output at each time step
        - typically output selected fields in vtk format
        - return generated code as string
        """
        if self.output_vts:
            return self.vfields[0].vtk_save_field()
        return ''

    @property
    def define_convergence(self):
        """Code fragment that defines convergence norms"""
        result = []
        for f in self.fields:
            result.append(cgen.Value(self.real_t, '%s_l2' % ccode(f.label)))
        return str(cgen.Module(result))

    @property
    def print_convergence(self):
        """Code fragment that prints convergence norms"""
        statements = [cgen.Statement('printf("%s %s\\n", conv.%s_l2)' % (ccode(f.label), '\t%.10f', ccode(f.label))) for f in self.fields]
        return str(cgen.Module(statements))

    @property
    def converge_test(self):
        """
        - generate code for convergence test
        - convergence test implemented by calculating L2 norm
        of the simulation against analytical solution
        - L2 norm of each field is calculated and output with printf()
        - return generated code as string
        """
        result = []
        if not self.converge:
            return str(cgen.Module(result))

        m = self.margin.value
        ti = self.ntsteps.value % 2  # last updated grid
        loop = [Symbol('_'+x.name) for x in self.index]  # symbols for loop

        for i in range(len(self.spacing)):
            result.append(cgen.Statement('printf("%d\\n")' % self.spacing[i].value))

        for field in self.sfields+self.vfields:
            body = []
            l2 = ccode(field.label)+'_l2'
            idx = [ti] + loop
            result.append(cgen.Initializer(cgen.Value(self.real_t, l2), 0.0))
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
                pre = [cgen.Initializer(cgen.Value(self.real_t, self.index[d].name), ccode(expr))]

                if d == self.dimension-1:
                    # inner loop
                    tn = self.dt.value*self.ntsteps.value \
                        if not field.staggered[0] \
                        else self.dt.value*self.ntsteps.value \
                        + self.dt.value/2.0
                    body = [cgen.Statement('%s += %s' % (l2, ccode((field[idx] - (field.sol.subs(self.t, tn)))**2.0)))]
                body = pre+body
                body = [cgen.For(cgen.InlineInitializer(cgen.Value('int', i), i0), cgen.Line('%s<%s' % (i, i1)), cgen.Line('++%s' % i), cgen.Block(body))]

            result += body
            volume = 1.0
            for i in range(len(self.spacing)):
                volume *= self.spacing[i].value
            l2_value = 'pow(' + l2 + '*' + ccode(volume) + ', 0.5)'
            result.append(cgen.Statement('conv->%s = %s' % (l2, l2_value)))

        return str(cgen.Module(result))

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
