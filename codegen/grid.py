from sympy import Symbol, symbols, factorial, Matrix, Rational, Indexed, simplify, IndexedBase, solve, Eq
from mako.template import Template
from mako.lookup import TemplateLookup
from mako.runtime import Context
from StringIO import StringIO

from sympy.printing.ccode import CCodePrinter

hf = Rational(1,2) # 1/2


class MyCPrinter(CCodePrinter):
	def __init__(self, settings={}):
		CCodePrinter.__init__(self, settings)

	def _print_Indexed(self, expr):
		# Array base and append indices
		output = self._print(expr.base.label) + ''.join(['[' + self._print(x) + ']' for x in expr.indices])
		return output

	def _print_Rational(self, expr):
		p, q = int(expr.p), int(expr.q)
		return '%d.0F/%d.0F' % (p, q) # float precision by default

	def _print_Mod(self, expr):
		# print mod as % in C++
		args = map(ccode, expr.args)
		args = ['('+x+')' for x in args]
		result = '%'.join(args)
		return result

def ccode(expr, **settings):
	return MyCPrinter(settings).doprint(expr, None)

def ccode_eq(eq, **settings):
	# print assignment code from equation
	return MyCPrinter(settings).doprint(eq.lhs, None) + ' = ' + MyCPrinter(settings).doprint(eq.rhs, None)

def render(tmpl, dict1):
	buf = StringIO()
	ctx = Context(buf, **dict1)
	tmpl.render_context(ctx)
	return buf.getvalue()

def IndexedBases(s):
	l = s.split();
	bases = [IndexedBase(x) for x in l]
	return tuple(bases)

def tc(dx, n):
    # return coefficient of power n Taylor series term
    return (dx**n)/factorial(n)

def Taylor(dx, n):
    # return Matrix of Taylor Coefficients M, such that M * D = R
    # where D is list of derivatives at x [f, f', f'' ..]
    # R is value at neighbour grid point [.. f(x-dx), f(x), f(x+dx) ..]
    l = []
    for i in range(-n+1,n):
    	ll = [tc((i*dx),j) for j in range(2*n-1)]
    	l.append(ll)
    return Matrix(l)

def Taylor_half(dx, n):
    # return Matrix of Taylor Coefficients M, such that M * D = R
    # where D is list of derivatives at x [f, f', f'' ..]
    # R is value at neighbour grid point [.. f(x-dx/2), f(x), f(x+dx/2), f(x+3/2*dx) ..]
    l = []
    for i in range(-n*2+1,n*2,2):
    	ll = [tc((i*dx/2),j) for j in range(2*n)]
    	l.append(ll)
    return Matrix(l)

def Deriv(U, i, k, d, n):
	# get the FD approximation for nth derivative in terms of grid U
	# i is list of indices of U, e.g. [x,y,z,t] for 3D
	# k = which dimension to expand, k=0 for x, k=1 for t etc
	M = Taylor(d, n)
	s = [0]*len(i)
	s[k] = 1 # switch on kth dimension
	# generate matrix of RHS, i.e. [ ... U[x-1], U[x], U[x+1] ... ]
	if len(i)==1:
		RX = Matrix([U[i[0]+s[0]*x] for x in range(-n+1,n)])
	elif len(i)==2:
		RX = Matrix([U[i[0]+s[0]*x,i[1]+s[1]*x] for x in range(-n+1,n)])
	elif len(i)==3:
		RX = Matrix([U[i[0]+s[0]*x,i[1]+s[1]*x,i[2]+s[2]*x] for x in range(-n+1,n)])
	elif len(i)==4:
		RX = Matrix([U[i[0]+s[0]*x,i[1]+s[1]*x,i[2]+s[2]*x,i[3]+s[3]*x] for x in range(-n+1,n)])
	else:
		raise NotImplementedError(">4 dimensions, need to fix")

	return M.inv() * RX

def Deriv_half(U, i, k, d, n):
	# get the FD approximation for nth derivative in terms of grid U
	# i is list of indices of U, e.g. [x,y,z,t] for 3D
	# k = which dimension to expand, k=0 for x, k=1 for t etc
	M = Taylor_half(d, n)
	s = [0]*len(i)
	s[k] = 1 # switch on kth dimension
	# generate matrix of RHS, i.e. [ ... U[x-1], U[x], U[x+1] ... ]
	if len(i)==1:
		RX = Matrix([U[i[0]+s[0]*x*hf] for x in range(-n*2+1,n*2,2)])
	elif len(i)==2:
		RX = Matrix([U[i[0]+s[0]*x*hf,i[1]+s[1]*x*hf] for x in range(-n*2+1,n*2,2)])
	elif len(i)==3:
		RX = Matrix([U[i[0]+s[0]*x*hf,i[1]+s[1]*x*hf,i[2]+s[2]*x*hf] for x in range(-n*2+1,n*2,2)])
	elif len(i)==4:
		RX = Matrix([U[i[0]+s[0]*x*hf,i[1]+s[1]*x*hf,i[2]+s[2]*x*hf,i[3]+s[3]*x*hf] for x in range(-n*2+1,n*2,2)])
	else:
		raise NotImplementedError(">4 dimensions, need to fix")

	result = M.inv() * RX
	return result

def is_half(expr):
	d = {x:0 for x in expr.free_symbols}
	return not expr.subs(d).is_Integer

def shift_grid(expr):
	if expr.is_Symbol:
		return expr
	if expr.is_Number:
		return expr
	if isinstance(expr,Indexed):
		b = expr.base
		idx = [x-hf if is_half(x) else x for x in list(expr.indices)]
		t = Indexed(b,*idx)
		return t
	args = tuple([shift_grid(arg) for arg in expr.args])
	expr2 = expr.func(*args)
	return expr2

def shift_index(expr, k, s):
	if expr.is_Symbol:
		return expr
	if expr.is_Number:
		return expr
	if isinstance(expr,Indexed):
		b = expr.base
		idx = list(expr.indices)
		idx[k] += s
		t = Indexed(b,*idx)
		return t
	args = tuple([shift_index(arg,k,s) for arg in expr.args])
	expr2 = expr.func(*args)
	return expr2

class Field(IndexedBase):

	def set(self, dimension, staggered):
		self.dimension = dimension
		self.staggered = staggered
		self.lookup = TemplateLookup(directories=['templates/staggered/'])
		self.d = [[None]*4 for x in range(dimension+1)] # list of list to store derivative expressions
		self.bc = [[None]*2 for x in range(dimension+1)] # list of list to store boundary ghost cell code

	def set_analytic_solution(self, function):
		self.sol = function

	def calc_derivative(self, l, k, d, n):
		self.d[k][n] = Deriv_half(self, l, k, d, n)[1]

	def align(self, expr):
		# align staggered cells to whole number indices
		if expr.is_Symbol or expr.is_Number:
			return expr
		if isinstance(expr,Indexed):
			b = expr.base
			# align the index if expression field staggered different from this field
			idx = []
			for k in range(len(expr.indices)):
				if self.staggered[k] == b.staggered[k]:
					idx += [expr.indices[k]]
				elif self.staggered[k]:
					idx += [expr.indices[k]+hf]
				else:
					idx += [expr.indices[k]-hf]
			tmp = Indexed(b,*idx)
			return tmp
		# recursive call for all arguments of expr
		args = tuple([self.align(arg) for arg in expr.args])
		result = expr.func(*args)
		return result

	def set_fd(self,fd):
		self.fd = fd
		t = self.align(fd)
		self.fd_align = t

	def vtk_save_field(self):
		tmpl = self.lookup.get_template('save_field.txt')
		result = ''
		dict1 = {'filename':ccode(self.label)+'_','field':ccode(self.label)}
		result = render(tmpl, dict1)

		return result

	def set_dt(self,dt):
		self.dt = dt

class VField(Field):

	def set(self, dimension, direction):
		self.direction = direction
		staggered = [False] * (dimension+1)
		staggered[0] = True
		staggered[direction] = True
		Field.set(self,dimension,staggered)

	def associate_stress_fields(self,sfields):
		self.sfields = sfields

	def set_free_surface(self, indices, d, b, side):
		# boundary at dimension[d] = b
		field = self.sfields[d] # use this stress field to solve for ghost cell expression
		idx = list(indices)
		if not field.staggered[d]:
			eq = Eq(field.dt)
			shift = hf
			t = b - hf
		else:
			eq = Eq(field.dt.subs(indices[d],indices[d]-hf),field.dt.subs(indices[d],indices[d]+hf))
			shift = 1
			t = b

		idx[d] -= ((-1)**side)*shift

		lhs = self[idx]
		rhs = solve(eq,lhs)[0]
		lhs = lhs.subs(indices[d],t)
		rhs = self.align(rhs.subs(indices[d],t))

		self.bc[d][side] = ccode(lhs) + ' = ' + ccode(rhs) + ';\n'


class SField(Field):

	def set(self, dimension, direction):
		self.direction = direction
		staggered = [False] * (dimension+1)
		if direction[0]==direction[1]:
			# compression stress, not staggered
			Field.set(self,dimension,staggered)
		else:
			# sheer stress, staggered
			for i in range(len(direction)):
				staggered[direction[i]] = True
			Field.set(self,dimension,staggered)

	def set_free_surface(self, indices, d, b, side):
		# boundary at dimension[d] = b
		if d not in self.direction:
			# x boundary for Tyy etc are not needed
			self.bc[d][side] = '// nothing\n'
			return

		idx = list(indices) # ghost cell
		idx2 = list(indices) # cell inside domain
		
		if not self.staggered[d]:
			idx[d] = b
			idx2[d] = b
			eq1 = Eq(self[idx])
		else:
			idx[d] = b - (1-side)
			idx2[d] = idx[d] + (-1)**side
			eq1 = Eq(self[idx], -self[idx2])

		idx[d] -= (-1)**side
		idx2[d] += (-1)**side
		eq2 = Eq(self[idx], -self[idx2])
		self.bc[d][side] = ccode_eq(eq1) +';\n' + ccode_eq(eq2) + ';\n'

class Variable(Symbol):
	""" wrapper for Symbol to store extra information """
	def __new__(cls, name, *args):
		return Symbol.__new__(cls,name)

	def __init__(self, name, value=0, type='int', constant=False):
		self.type = type
		self.constant = constant
		self.value = value


class StaggeredGrid:
	""" description of staggered grid for finite difference method """
	
	def __init__(self, dimension):
		self.dimension = dimension
		self.lookup = TemplateLookup(directories=['templates/staggered/'])
		self.omp = True # switch for inserting #pragma omp for
		self.size = [1.0] * dimension # default domain size
		self.spacing = [Variable('dx'+str(k+1), 0.1, 'float', True)  for k in range(dimension)] # spacing symbols, dx1, dx2, ...
		self.index = [Symbol('x'+str(k+1))  for k in range(dimension)] # indices symbols, x1, x2 ...

		self.order = (1,2,2,2)

		self.t = Symbol('t')
		self.dt = Variable('dt', 0.01, 'float', True)
		self.tp = Variable('tp', self.order[0]*2, 'int', True) # periodicity for time stepping
		self.margin = Variable('margin', 2, 'int', True)
		self.ntsteps= Variable('ntsteps', 100, 'int', True)

		self.defined_variable = {} # user defined variables, use dictionary because of potential duplication when using list

		# add time variables for time stepping
		self.time = []
		for k in range(self.order[0]*2):
			name = 't' + str(k)
			v = Variable(name, 0, 'int', False)
			self.time.append(v)

		self._update_domain_size()

		#self.float_symbols = {self.h[0]:0.1,self.h[1]:0.1,self.h[2]:0.1,self.dt:1.0} # dictionary to hold symbols and their values
		#self.int_symbols = {self.ntsteps:1,self.dim[0]:0,self.dim[1]:0,self.dim[2]:0}

	def _update_domain_size(self):
		# set dimension symbols, dim1, dim2, ...
		self.dim = [Variable('dim'+str(k+1), int(self.size[k]/self.spacing[k].value)+1+self.margin.value*2, 'int', True)  for k in range(self.dimension)]
		expr = 2*self.order[0]
		for d in self.dim:
			expr *= d
		self.vec_size = Variable('vec_size', expr, 'int', True)

	def set_domain_size(self, size):
		self.size = size
		self._update_domain_size()

	def set_spacing(self, spacing):
		self.spacing = [Variable('dx'+str(k+1), spacing[k], 'float', True)  for k in range(self.dimension)] # spacing symbols, dx1, dx2, ...
		self._update_domain_size()

	def set_index(self, indices):
		self.index = indices

	def set_variable(self, var, value=0, type='int', constant=False):
		if isinstance(var, Symbol):
			var = var.name
		self.defined_variable[var] = Variable(var,value,type,constant)

	def get_time_step_limit(self):
		# Vp = sqrt((lambda+2*mu)/rho)
		l = self.defined_variable['lambda'].value
		m = self.defined_variable['mu'].value
		r = self.defined_variable['rho'].value
		Vp = ((l + 2*m)/r)**0.5
		h = min([sp.value for sp in self.spacing])
		return 0.5*h/Vp

	def set_time_step(self, dt, tmax):
		self.dt.value = dt
		# self.float_symbols[Symbol('tmax')] = tmax
		self.ntsteps.value = int(tmax/dt)

	def set_stress_fields(self, sfields):
		num = self.dimension+self.dimension*(self.dimension-1)/2
		if not len(sfields)==num:
			raise Exception('wrong number of stress fields: '+str(num)+' fields required.')
		self.sfields = sfields

	def set_velocity_fields(self, vfields):
		num = self.dimension
		if not len(vfields)==num:
			raise Exception('wrong number of velocity fields: '+str(num)+' fields required.')
		self.vfields = vfields

	def calc_derivatives(self):
		l = [self.dt] + self.spacing
		for field in self.sfields+self.vfields:
			for k in range(self.dimension+1):
				# loop through dimensions
				h = l[k]
				for o in range(1,self.order[k]+1):
					# loop through order of derivatives
					field.calc_derivative([self.t]+self.index,k,h,o)

	def solve_fd(self,eqs):
		t = self.t
		index = [t+hf] + self.index
		self.eqs = eqs
		for field, eq in zip(self.vfields+self.sfields, eqs):
			field.set_fd(solve(eq,field[index])[0].subs({t:t+hf}))

	def associate_fields(self):
		# match up stress fields with velocity fields
		for v in self.vfields:
			fields = [None]
			for d in range(1,self.dimension+1):
				lookfor = tuple(sorted([v.direction,d]))
				fields.append([f for f in self.sfields if f.direction==lookfor][0])
			if not len(fields)==self.dimension+1:
				print len(fields)
				raise Exception('error in field directions')
			v.associate_stress_fields(fields)

	def set_free_surface_boundary(self, dimension, side):
		self.associate_fields()
		index = [self.t] + self.index
		for field in self.sfields+self.vfields:
			if side==0:
				field.set_free_surface(index, dimension, self.margin.value, side)
			else:
				field.set_free_surface(index, dimension, self.dim[dimension-1]-self.margin.value-1, side)

	############## sub-routines for output ##############

	def define_variables(self):
		result = ''
		variables = self.dim + self.spacing + self.time + [self.tp, self.dt, self.margin, self.ntsteps, self.vec_size] + self.defined_variable.values()
		for v in variables:
			line = ''
			if v.constant:
				line += 'const '
			line += v.type + ' ' + v.name + ' = ' + str(v.value) + ';\n'
			result += line
		return result

	def declare_fields(self):
		result = ''
		arr = '' # =[dim1][dim2][dim3]...
		for d in self.dim:
			arr += '[' + d.name + ']'
		for field in self.sfields + self.vfields:
			vec = '_' + ccode(field.label) + '_vec'
			result += 'std::vector<float> ' + vec + '(' + self.vec_size.name + ');\n'
			result += 'float (*' + ccode(field.label) + ')' + arr + '= (float (*)' + arr + ') ' + vec + '.data();\n'
		return result

	def initialise(self):

		tmpl = self.lookup.get_template('generic_loop.txt')
		result = ''
		m = self.margin.value
		loop = [Symbol('_'+x.name) for x in self.index] # symbols for loop

		for field in self.sfields+self.vfields:
			body = ''
			if self.omp:
				result += '#pragma omp for\n'
			# populate xvalue, yvalue zvalue code
			for d in range(self.dimension-1,-1,-1):
				i = loop[d]
				i0 = m
				if field.staggered[d+1]:
					i1 = ccode(self.dim[d]-m-1)
					expr = self.spacing[d]*(loop[d] - self.margin.value + 0.5)
				else:
					i1 = ccode(self.dim[d]-m)
					expr = self.spacing[d]*(loop[d] - self.margin.value)
				pre = 'float ' + self.index[d].name + '= ' + ccode(expr) + ';\n'
				post = ''
				if d==self.dimension-1:
					# inner loop
					t0 = self.dt.value/2 if field.staggered[0] else 0 # first time step
					body = ccode(field[[0]+loop]) + '=' + ccode(field.sol.subs(self.t, t0)) +';\n'
				body = pre + body
				dict1 = {'i':i,'i0':i0,'i1':i1,'body':body}
				body = render(tmpl, dict1)

			result += body
		return result

	def initialise_boundary(self):
		result = self.stress_bc().replace('[t1]','[0]')
		result += self.velocity_bc().replace('[t1]','[0]')
		return result

	def stress_loop(self):
		tmpl = self.lookup.get_template('generic_loop.txt')
		m = self.margin.value
		body = ''
		for d in range(self.dimension-1,-1,-1):
			i = self.index[d]
			i0 = m
			i1 = ccode(self.dim[d]-m)
			if d==self.dimension-1:
				# inner loop
				idx = [self.time[1]] + self.index
				for field in self.sfields:
					body += ccode(field[idx]) + '=' + ccode(field.fd_align.xreplace({self.t+1:self.time[1], self.t:self.time[0]})) + ';\n'
			dict1 = {'i':i,'i0':i0,'i1':i1,'body':body}
			body = render(tmpl, dict1)

		if self.omp:
			body = '#pragma omp for\n' + body

		return body

	def velocity_loop(self):
		tmpl = self.lookup.get_template('generic_loop.txt')
		m = self.margin.value
		body = ''
		for d in range(self.dimension-1,-1,-1):
			i = self.index[d]
			i0 = m
			i1 = ccode(self.dim[d]-m)
			if d==self.dimension-1:
				# inner loop
				idx = [self.time[1]] + self.index
				for field in self.vfields:
					body += ccode(field[idx]) + '=' + ccode(field.fd_align.xreplace({self.t+1:self.time[1], self.t:self.time[0]})) + ';\n'
			dict1 = {'i':i,'i0':i0,'i1':i1,'body':body}
			body = render(tmpl, dict1)

		if self.omp:
			body = '#pragma omp for\n' + body

		return body

	def stress_bc(self):
		tmpl = self.lookup.get_template('generic_loop.txt')
		result = ''
		for field in self.sfields:
			for d in range(self.dimension):
				for side in range(2):
					if self.omp:
						result += '#pragma omp for\n'
					body = ''
					for d2 in range(self.dimension-1,-1,-1):
						# loop through other dimensions
						if not d2==d:
							i = self.index[d2]
							i0 = 0
							i1 = self.dim[d2]
							if body=='':
								# inner loop, populate ghost cell calculation
								body = field.bc[d+1][side]
							dict1 = {'i':i,'i0':0,'i1':i1,'body':body}
							body = render(tmpl, dict1).replace('[t]','[t1]')

					result += body

		return result

	def velocity_bc(self):
		tmpl = self.lookup.get_template('generic_loop.txt')
		result = ''
		for d in range(self.dimension):
			# update the staggered field first because other fields depends on it
			sequence = [f for f in self.vfields if f.staggered[d+1]] + [f for f in self.vfields if not f.staggered[d+1]]
			for field in sequence:
				for side in range(2):
					if self.omp:
						result += '#pragma omp for\n'
					body = ''
					for d2 in range(self.dimension-1,-1,-1):
						# loop through other dimensions
						if not d2==d:
							i = self.index[d2]
							i0 = 1
							i1 = self.dim[d2]-1
							if body=='':
								# inner loop, populate ghost cell calculation
								body = field.bc[d+1][side]
							dict1 = {'i':i,'i0':i0,'i1':i1,'body':body}
							body = render(tmpl, dict1).replace('[t]','[t1]')

					result += body

		return result

	def time_stepping(self):
		# generate time index for time stepping
		result = ''
		tmpl = self.lookup.get_template('time_stepping.txt')
		_ti = Symbol('_ti')
		body = ''

		for i in range(len(self.time)):
			lhs = self.time[i].name
			if i==0:
				rhs = ccode(_ti % self.tp)
			else:
				rhs = ccode((self.time[i-1]+1) % self.tp)
			body += lhs + ' = ' + rhs + ';\n'

		dict1 = {'body':body}
		result = render(tmpl, dict1)
		return result

	def output_step(self):
		result = ''
		#result += self.vfields[0].vtk_save_field()
		return result

	def converge_test(self):
		tmpl = self.lookup.get_template('generic_loop.txt')
		result = ''
		m = self.margin.value
		ti = self.ntsteps.value % 2 # last updated grid
		loop = [Symbol('_'+x.name) for x in self.index] # symbols for loop

		for field in self.sfields+self.vfields:
			body = ''
			l2 = ccode(field.label)+'_l2'
			idx = [ti] + loop
			result += 'float ' + l2 + ' = 0.0;\n'
			# populate xvalue, yvalue zvalue code
			for d in range(self.dimension-1,-1,-1):
				i = loop[d]
				i0 = m
				if field.staggered[d+1]:
					i1 = ccode(self.dim[d]-m-1)
					expr = self.spacing[d]*(loop[d] - self.margin.value + 0.5)
				else:
					i1 = ccode(self.dim[d]-m)
					expr = self.spacing[d]*(loop[d] - self.margin.value)
				pre = 'float ' + self.index[d].name + '= ' + ccode(expr) + ';\n'
				if d==self.dimension-1:
					# inner loop
					tn = self.dt.value*self.ntsteps.value if not field.staggered[0] else self.dt.value*self.ntsteps.value + self.dt.value/2.0
					body = l2 + '+=' + ccode((field[idx]-(field.sol.subs(self.t,tn)))**2.0) + ';\n'
				body = pre + body
				dict1 = {'i':i,'i0':i0,'i1':i1,'body':body}
				body = render(tmpl, dict1)

			result += body
			result += 'printf("' + l2 + ' = %.10f\\n", ' + l2 + ');\n'

		return result
