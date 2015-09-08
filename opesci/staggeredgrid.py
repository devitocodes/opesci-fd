from grid import Grid
from variable import Variable
from fields import Field, VField, get_ops_expr
from codeprinter import ccode, render
from compilation import get_package_dir

from sympy import Symbol, Rational, solve, expand
from mako.lookup import TemplateLookup
import mmap
from os import path

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
    """
    template_base = 'staggered3d_tmpl.cpp'

    template_keys = ['io', 'time_stepping', 'define_constants', 'declare_fields',
                     'define_fields', 'store_fields', 'load_fields',
                     'initialise', 'initialise_bc', 'stress_loop',
                     'velocity_loop', 'stress_bc', 'velocity_bc', 'output_step',
                     'define_convergence', 'converge_test', 'print_convergence']

    _switches = ['omp', 'ivdep', 'simd', 'double', 'expand', 'eval_const',
                 'output_vts', 'converge','polly']

    def __init__(self, dimension, domain_size=None, grid_size=None,
                 time_step=None, stress_fields=None, velocity_fields=None,
                 omp=True, ivdep=True, simd=False, double=False, io=False,
                 expand=True, eval_const=True, output_vts=False, converge=False,polly=False):
        self.dimension = dimension

        template_dir = path.join(get_package_dir(), "templates")
        staggered_dir = path.join(get_package_dir(), "templates/staggered")
        self.lookup = TemplateLookup(directories=[template_dir, staggered_dir])

        # List of associated fields
        self.sfields = []
        self.vfields = []

        # Switches
        self.polly = polly
        self.omp = omp
        self.ivdep = ivdep
        self.simd = simd
        self.double = double
        self.expand = expand
        self.eval_const = eval_const
        self.real_t = 'double' if self.double else 'float'
        self.output_vts = output_vts
        self.converge = converge

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

        # Optional further grid settings
        if stress_fields:
            self.set_stress_fields(stress_fields)
        if velocity_fields:
            self.set_velocity_fields(velocity_fields)
        if domain_size:
            self.set_domain_size(domain_size)
        if grid_size:
            self.set_grid_size(grid_size)

        self.t = Symbol('t')
        self.dt = Variable('dt', 0.01, self.real_t, True)
        self.ntsteps = Variable('ntsteps', 100, 'int', True)
        self.alignment = mmap.PAGESIZE  # default alignment for malloc

        # user defined variables
        # use dictionary because of potential duplication
        self.defined_variable = {}

        self._update_spacing()

    @property
    def fields(self):
        return self.sfields + self.vfields

    @property
    def io(self):
        """Flag whether to include I/O headers"""
        return self.read or self.output_vts

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
        if self.polly:
            index = [t1]+self.index#polly
        else:
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
            #print eq
            #print field[index]
            kernel = solve(eq, field[index])[0]
           # print kernel
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
        # if self.polly:
        #     index = index[1:]       

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
        store = 0
        add = 0
        mul = 0
        arrays = []  # to store name of arrays loaded
        for field in self.vfields:
            store += 1  # increment STORE by 1 (assignment)
            expr = field.fd_align
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
        - get the arithmetic intensity of velocity kernel
        - get the number of different operations of the stress field kernel
        - types of operations are ADD (inc -), MUL (inc /), LOAD, STORE
        - #LOAD = number of unique fields in the kernel
        - return tuple (#ADD, #MUL, #LOAD, #STORE)
        - arithmetic intensity AI = (ADD+MUL)/[(LOAD+STORE)*word size]
        - weighted AI, AI_w = (ADD+MUL)/(2*Max(ADD,MUL)) * AI
        """
        store = 0
        add = 0
        mul = 0
        arrays = []  # to store name of arrays loaded
        for field in self.sfields:
            store += 1  # increment STORE by 1 (assignment)
            expr = field.fd_align
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

    # ------------------- sub-routines for output -------------------- #

    @property
    def define_constants(self):
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

    @property
    def define_fields(self):
        """Code fragment that defines field arrays"""
        return '\n'.join(['%s *%s;' % (self.real_t, ccode(f.label))
                      for f in self.fields])


    @property
    def store_fields(self):
        """Code fragment that stores field arrays to 'grid' struct"""
        return '\n'.join(['grid->%s = (%s*) %s;' %
                          (ccode(f.label), self.real_t, ccode(f.label))
                          for f in self.fields])

    @property
    def load_fields(self):
        """Code fragment that loads field arrays from 'grid' struct"""
        if self.polly:
            idxs = ''.join(['[%d]' % d.value for d in self.dim[1:]])
        else:
            idxs = ''.join(['[%d]' % d.value for d in self.dim])

        return '\n'.join(['%s (*%s)%s = (%s (*)%s) grid->%s;' %
                          (self.real_t, ccode(f.label), idxs,
                           self.real_t, idxs, ccode(f.label))
                          for f in self.fields])

    @property
    def declare_fields(self):
        """
        - generate code for delcaring fields
        - the generated code first declare fields as std::vector
        of size=vec_size, then cast to multidimensional array
        - return the generated code as string
        """
        result = ''
        arr = ''  # = [dim1][dim2][dim3]...
        if self.polly:
            for d in self.dim[1:]:
                arr += '[' + d.name + ']'
        else:
            for d in self.dim:
                arr += '[' + d.name + ']'

        vsize = 1
        for d in self.dim:
            vsize *= d.value
        vsize *= self.order[0]*2

        if self.polly:
            vsize /=2 
            for field in self.sfields + self.vfields:
                vec = '_' + ccode(field.label) + '_vec_0'
                # alloc aligned memory (on windows and linux)
                result += self.real_t + ' *' + vec + ';\n'
                result += '#ifdef _MSC_VER\n'
                result += vec + ' = (' + self.real_t + '*) _aligned_malloc(' + str(vsize) \
                    + '*sizeof(' + self.real_t + '), ' + str(self.alignment) + ');\n'
                result += '#else\n'
                result += 'posix_memalign((void **)(&' + vec + '), ' + str(self.alignment) \
                    + ', ' + str(vsize) + '*sizeof(' + self.real_t + '));\n'
                result += '#endif\n'
                # cast pointer to multidimensional array
                result += self.real_t + ' (*' + ccode(field.label) + '_0)' + arr \
                    + '= (' + self.real_t + ' (*)' + arr + ') ' + vec + ';\n'
            
            result += '\n'

            for field in self.sfields + self.vfields:
                vec = '_' + ccode(field.label) + '_vec_1'
                # alloc aligned memory (on windows and linux)
                result += self.real_t + ' *' + vec + ';\n'
                result += '#ifdef _MSC_VER\n'
                result += vec + ' = (' + self.real_t + '*) _aligned_malloc(' + str(vsize) \
                    + '*sizeof(' + self.real_t + '), ' + str(self.alignment) + ');\n'
                result += '#else\n'
                result += 'posix_memalign((void **)(&' + vec + '), ' + str(self.alignment) \
                    + ', ' + str(vsize) + '*sizeof(' + self.real_t + '));\n'
                result += '#endif\n'
                # cast pointer to multidimensional array
                result += self.real_t + ' (*' + ccode(field.label) + '_1)' + arr \
                    + '= (' + self.real_t + ' (*)' + arr + ') ' + vec + ';\n'
            
            result += '\n'

            for field in self.sfields + self.vfields:
                result += self.real_t + ' (*' + ccode(field.label) + ')[' \
                        + str(self.dim[1])+'][' + str(self.dim[2])+'];\n'
            for field in self.sfields + self.vfields:
                result += self.real_t + ' (*' + ccode(field.label) + '_old)[' \
                        + str(self.dim[1])+'][' + str(self.dim[2])+'];\n'            

        else :   
            for field in self.sfields + self.vfields:
                vec = '_' + ccode(field.label) + '_vec'
                # alloc aligned memory (on windows and linux)
                result += self.real_t + ' *' + vec + ';\n'
                result += '#ifdef _MSC_VER\n'
                result += vec + ' = (' + self.real_t + '*) _aligned_malloc(' + str(vsize) \
                    + '*sizeof(' + self.real_t + '), ' + str(self.alignment) + ');\n'
                result += '#else\n'
                result += 'posix_memalign((void **)(&' + vec + '), ' + str(self.alignment) \
                    + ', ' + str(vsize) + '*sizeof(' + self.real_t + '));\n'
                result += '#endif\n'
                # cast pointer to multidimensional array
                result += self.real_t + ' (*' + ccode(field.label) + ')' + arr \
                    + '= (' + self.real_t + ' (*)' + arr + ') ' + vec + ';\n'

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
            loop = [self.rho, self.vp, self.vs] + self.beta + [self.lam] + self.mu
            for field in loop:
                vec = '_' + ccode(field.label) + '_vec'
                # alloc aligned memory (on windows and linux)
                result += self.real_t + ' *' + vec + ';\n'
                result += '#ifdef _MSC_VER\n'
                result += vec + ' = (' + self.real_t + '*) _aligned_malloc(' + str(vsize) \
                    + '*sizeof(' + self.real_t + '), ' + str(self.alignment) + ');\n'
                result += '#else\n'
                result += 'posix_memalign((void **)(&' + vec + '), ' + str(self.alignment) \
                    + ', ' + str(vsize) + '*sizeof(' + self.real_t + '));\n'
                result += '#endif\n'
                # cast pointer to multidimensional array
                result += self.real_t + ' (*' + ccode(field.label) + ')' + arr \
                    + '= (' + self.real_t + ' (*)' + arr + ') ' + vec + ';\n'

            # read from file
            result += 'opesci_read_simple_binary_ptr("' + self.rho_file + '",_' \
                + ccode(self.rho.label) + '_vec, ' + str(vsize) + ');\n'
            result += 'opesci_read_simple_binary_ptr("' + self.vp_file + '",_' \
                + ccode(self.vp.label) + '_vec, ' + str(vsize) + ');\n'
            result += 'opesci_read_simple_binary_ptr("' + self.vs_file + '",_' \
                + ccode(self.vs.label) + '_vec, ' + str(vsize) + ');\n'
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

        tmpl = self.lookup.get_template('generic_loop.txt')
        result = ''
        m = self.margin.value
        loop = [Symbol('_'+x.name) for x in self.index]  # symbols for loop

        #refine later
        if self.polly:
            for f in self.sfields+self.vfields:
                result += str(f) + " = " + str(f)+ "_0;\n"


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
                    if(self.polly):
                        body = ccode(field[loop]) + '=' \
                            + ccode(sol) + ';\n'

                    else:
                        body = ccode(field[[0]+loop]) + '=' \
                            + ccode(sol) + ';\n'
                body = pre + body + post
                dict1 = {'i': i, 'i0': i0, 'i1': i1, 'body': body}
                body = render(tmpl, dict1)

            result += body
        return result

    @property
    def initialise_bc(self):
        """
        - generate code for initialisation of boundary ghost cells
        - generate generic boundary cell code
        replace array indices [t] with [0]
        - return generated code as string
        """
        if(self.polly):
            self.polly = False;
            result = self.stress_bc.replace('[t1]', '')
            result += self.velocity_bc.replace('[t1]', '')
            self.polly = True;
        else:
            result = self.stress_bc.replace('[t1]', '[0]')
            result += self.velocity_bc.replace('[t1]', '[0]')
        return result

    @property
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
                    print field[idx]
               #     print ', '.join("%s: %s" % item for item in vars(field).items())
                    body += ccode(field[idx]) + '=' \
                        + ccode(field.fd_align.xreplace({self.t+1:
                                                        self.time[1],
                                                        self.t:
                                                        self.time[0]})) \
                        + ';\n'
            dict1 = {'i': i, 'i0': i0, 'i1': i1, 'body': body}
            body = render(tmpl, dict1)
            if self.ivdep and d == self.dimension-1:
                    body = '%s\n' % self.compiler._ivdep + body
            if self.simd and d == self.dimension-1:
                    body = '#pragma simd\n' + body

        if self.omp:
            body = '#pragma omp for\n' + body

        return body

    @property
    def velocity_loop(self):
        """
        generate code for velocity field update loop
        - loop through velocity fields to generate code of computation kernel
        - generate inner loop by inserting kernel into Mako template
        - recursive insertion to generate nested loop
        return generated code as string
        """

        body = ''
        tmpl = self.lookup.get_template('generic_loop.txt')
        m = self.margin.value
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
                    body = '%s\n' % self.compiler._ivdep + body
            if self.simd and d == self.dimension-1:
                    body = '#pragma simd\n' + body

        if self.omp:
            body = '#pragma omp for\n' + body

        return body

    @property
    def stress_bc(self):
        """
        generate code for updating stress field boundary ghost cells
        - generate inner loop by inserting boundary code (saved in field.bc)
        - recursive insertion to generate nested loop
        - loop through all stress fields and sides
        return generated code as string
        """

        result = ''
        _ti = Symbol('_ti')
        if self.polly:
            result += "if("+ccode(_ti)+"%2==0){\n"
            for s in self.sfields :
                result += str(s) + " = " + str(s) + "_1;\n"
            result += "\n"
            result += "}else{\n"
            for s in self.sfields :
                result += str(s) + " = " + str(s) + "_0;\n"
            result += "}\n"

        tmpl = self.lookup.get_template('generic_loop.txt')
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

                                if self.polly:
                                    body = render(tmpl, dict1).replace('[t]',
                                                                   '')
                                else:    
                                    body = render(tmpl, dict1).replace('[t]',
                                                                   '[t1]')
                                if self.ivdep:
                                    body = '%s\n' % self.compiler._ivdep + body
                                if self.simd:
                                    body = '#pragma simd\n' + body
                            else:
                                dict1 = {'i': i, 'i0': i0,
                                         'i1': i1, 'body': body}
                                if self.polly:
                                    body = render(tmpl, dict1).replace('[t]',
                                                                       '')
                                else:
                                    body = render(tmpl, dict1).replace('[t]',
                                                                       '[t1]')

                    result += body

        if self.polly:
            result +=  'if('
            result += ccode(_ti)
            result += '%2==0){\n'
            for t in self.sfields:
                result += '\t\t' +ccode(t.label) + ' = ' + '%s_1;\n'%t.label
            
            for t in self.vfields:
                result += '\t\t' +ccode(t.label) + ' = ' + '%s_1;\n'%t.label
            for t in self.vfields:
                result += '\t\t' +ccode(t.label) + '_old = ' + '%s_0;\n'%t.label
            result += '}else{'
            for t in self.sfields:
                result += '\t\t' +ccode(t.label) + ' = ' + '%s_0;\n'%t.label
            for t in self.vfields:
                result += '\t\t' +ccode(t.label) + ' = ' + '%s_0;\n'%t.label
            for t in self.vfields:
                result += '\t\t' +ccode(t.label) + '_old = ' + '%s_1;\n'%t.label
            result += '}\n'

        return result

    @property
    def velocity_bc(self):
        """
        generate code for updating stress field boundary ghost cells
        - generate inner loop by inserting boundary code (saved in field.bc)
        - recursive insertion to generate nested loop
        - loop through all velocity fields and sides
        return generated code as string
        """


        result = ''
        _ti = Symbol('_ti')
        if self.polly:
            result += "if("+ccode(_ti)+"%2==0){\n"
            for s in self.vfields :
                result += str(s) + " = " + str(s) + "_1;\n"
            result += "\n"
            result += "}else{\n"
            for s in self.vfields :
                result += str(s) + " = " + str(s) + "_0;\n"
            result += "}\n"

        tmpl = self.lookup.get_template('generic_loop.txt')
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
                                if self.polly:
                                    body = render(tmpl, dict1).replace('[t]',
                                                               '')
                                else:
                                    body = render(tmpl, dict1).replace('[t]',
                                                                   '[t1]')
                                if self.ivdep:
                                    body = '%s\n' % self.compiler._ivdep + body
                                if self.simd:
                                    body = '#pragma simd\n' + body
                            else:
                                dict1 = {'i': i, 'i0': i0,
                                         'i1': i1, 'body': body}
                                if self.polly:
                                    body = render(tmpl, dict1).replace('[t]',
                                                               '')
                                else:
                                    body = render(tmpl, dict1).replace('[t]',
                                                                   '[t1]')

                    result += body

        return result

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

        result = ''
        tmpl = self.lookup.get_template('time_stepping.txt')
        _ti = Symbol('_ti')
        body = ''

        if self.polly:
            body =  'if('
            body += ccode(_ti)
            body += '%2==0){\n'
            for t in self.sfields:
                body += '\t\t' +ccode(t.label) + ' = ' + '%s_0;\n'%t.label
            for t in self.sfields:
                body += '\t\t' +ccode(t.label) + '_old = ' + '%s_1;\n'%t.label

            for t in self.vfields:
                body += '\t\t' +ccode(t.label) + ' = ' + '%s_0;\n'%t.label
            for t in self.vfields:
                body += '\t\t' +ccode(t.label) + '_old = ' + '%s_1;\n'%t.label
            body += '}else{'
            for t in self.sfields:
                body += '\t\t' +ccode(t.label) + ' = ' + '%s_1;\n'%t.label
            for t in self.sfields:
                body += '\t\t' +ccode(t.label) + '_old = ' + '%s_0;\n'%t.label

            for t in self.vfields:
                body += '\t\t' +ccode(t.label) + ' = ' + '%s_1;\n'%t.label
            for t in self.vfields:
                body += '\t\t' +ccode(t.label) + '_old = ' + '%s_0;\n'%t.label
            body += '}\n'

        else:
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

    @property
    def output_step(self):
        """
        - generate code for output at each time step
        - typically output selected fields in vtk format
        - return generated code as string
        """
        result = ''
        if self.output_vts:
            result += self.vfields[0].vtk_save_field()
        return result

    @property
    def define_convergence(self):
        """Code fragment that defines convergence norms"""
        return '\n'.join(['%s %s_l2;' % (self.real_t, ccode(f.label))
                          for f in self.fields])

    @property
    def print_convergence(self):
        """Code fragment that prints convergence norms"""
        return '\n'.join(['printf("%s %s\\n", conv.%s_l2);' %
                          (ccode(f.label), '\t%.10f', ccode(f.label))
                          for f in self.fields])

    @property
    def converge_test(self):
        """
        - generate code for convergence test
        - convergence test implemented by calculating L2 norm
        of the simulation against analytical solution
        - L2 norm of each field is calculated and output with printf()
        - return generated code as string
        """
        result = ''
        if not self.converge:
            return result
        tmpl = self.lookup.get_template('generic_loop.txt')
        m = self.margin.value
        ti = self.ntsteps.value % 2  # last updated grid
        loop = [Symbol('_'+x.name) for x in self.index]  # symbols for loop

        for i in range(len(self.spacing)):
            result += 'printf("' + str(self.spacing[i].value) + '\\n");\n'

        for field in self.sfields+self.vfields:
            body = ''
            l2 = ccode(field.label)+'_l2'

            if self.polly:
                idx = loop
            else:
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
            result += 'conv->%s = %s;\n' % (l2, l2_value)

        return result
