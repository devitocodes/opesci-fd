from sympy import Symbol, symbols, factorial, Matrix, Rational, Indexed, simplify, IndexedBase, solve, Eq
from mako.template import Template
from mako.lookup import TemplateLookup
from mako.runtime import Context
from StringIO import StringIO

from sympy.printing.ccode import CCodePrinter

hf = Rational(1,2) # 1/2

def IndexedBases(s):
	l = s.split();
	bases = [IndexedBase(x) for x in l]
	return tuple(bases)

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

def ccode(expr, **settings):
	return MyCPrinter(settings).doprint(expr, None)

def render(tmpl, dict1):
	buf = StringIO()
	ctx = Context(buf, **dict1)
	tmpl.render_context(ctx)
	return buf.getvalue()

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

	def __new__(cls, name, *args, **kargs):
		return IndexedBase.__new__(cls, name)

	def __init__(self, name, dimension, staggered, *args, **kargs):
		self.lookup = TemplateLookup(directories=['templates/staggered/'])
		self.dimension = dimension
		self.staggered = staggered
		self.d = [[None]*4 for x in range(dimension+1)] # list of list to store derivative expressions
		self.bc = [[None]*2 for x in range(dimension+1)] # list of list to store boundary condition

	def set_analytic_solution(self, function):
		self.sol = function

	def calc_derivative(self, l, k, d, n):
		self.d[k][n] = Deriv_half(self, l, k, d, n)[1]

	def recenter(self, expr):
		result = expr
		for k in range(len(self.offset)):
			if self.offset[k]:
				result = shift_index(result,k,hf)
		return result

	def set_fd(self,fd):
		self.fd = fd
		t = self.recenter(fd)
		self.shift_fd = shift_grid(t)

	def save_field(self):
		tmpl = self.lookup.get_template('save_field.txt')
		result = ''
		dict1 = {'filename':ccode(self.name.label)+'_','field':ccode(self.name.label)}
		result = render(tmpl, dict1)

		return result

	def associate_stress_fields(self,sfields):
		self.sfields = sfields

	def set_dt(self,dt):
		self.dt = dt

class VField(Field):
	def __new__(cls, name, *args, **kargs):

		return Field.__new__(cls, name)

	def __init__(self, name, dimension, direction):
		self.direction = direction
		staggered = [False] * (dimension+1)
		staggered[0] = True
		staggered[direction] = True
		Field.__init__(self, name, dimension, staggered)

	def set_free_surface(self, indices, d, b, side):
		# boundary at dimension[d] = b
		field = self.sfields[d-1]
		idx = list(indices)
		if not field.offset[d]:
			eq = Eq(field.dt.subs(indices[d],b))
			shift = hf
		else:
			eq = Eq(field.dt.subs(indices[d],b-hf),field.dt.subs(indices[d],b+hf))
			shift = 1

		if side==0:
			idx[d] = b-shift
		else:
			idx[d] = b+shift

		lhs = self.name[idx]
		rhs = solve(eq,lhs)[0]
		if shift == 1:
			rhs = self.recenter(rhs) # shift if not shifted already
		self.bc[d][side] = ccode(shift_grid(lhs)) + '=' + ccode(shift_grid(rhs)) + ';\n\t\t\t\t'


class SField(Field):
	def __new__(cls, name, *args, **kargs):

		return Field.__new__(cls, name)

	def __init__(self, name, dimension, direction):
		self.direction = direction
		staggered = [False] * (dimension+1)
		for i in range(len(direction)):
			staggered[direction[i]] = True

		Field.__init__(self, name, dimension, staggered)

	def set_free_surface(self, indices, d, b, side):
		# boundary at dimension[d] = b
		result = ''
		if d == 1:
			ch = 'x'
		elif d == 2:
			ch = 'y'
		else:
			ch = 'z'

		if ch not in ccode(self.name.label):
			self.bc[d][side] = '// nothing'
			return

		idx = list(indices)
		idx2 = list(indices)
		if not self.offset[d]:
			idx[d] = b
			result += ccode(self.name[idx]) +' = 0;\n\t\t\t\t'
			if side==0:
				idx[d] = b-1
				idx2[d] = b+1
			else:
				idx[d] = b+1
				idx2[d] = b-1
			result += ccode(self.name[idx]) +' = ' + ccode(-self.name[idx2])+';\n\t\t\t\t'
		else:
			if side==0:
				idx[d] = b-1
				idx2[d] = b
				result += ccode(self.name[idx]) +' = ' + ccode(-self.name[idx2])+';\n\t\t\t\t'
				idx[d] = b-2
				idx2[d] = b+1
				result += ccode(self.name[idx]) +' = ' + ccode(-self.name[idx2])+';\n\t\t\t\t'
			else:
				idx[d] = b
				idx2[d] = b-1
				result += ccode(self.name[idx]) +' = ' + ccode(-self.name[idx2])+';\n\t\t\t\t'
				idx[d] = b+1
				idx2[d] = b-2
				result += ccode(self.name[idx]) +' = ' + ccode(-self.name[idx2])+';\n\t\t\t\t'

		self.bc[d][side] = result


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
		self.size = [1.0] * dimension # default domain size
		self.spacing = [Variable('dx'+str(k+1), 0.1, 'float', True)  for k in range(dimension)] # spacing symbols, dx1, dx2, ...
		self.index = [Symbol('x'+str(k+1))  for k in range(dimension)] # indices symbols, x1, x2 ...

		self.t = Symbol('t')
		self.dt = Variable('dt', 0.01, 'float', True)
		self.margin = Variable('margin', 2, 'int', True)
		self.ntsteps= Variable('ntsteps', 100, 'int', True)

		self.order = [1,2,2,2]

		self.defined_variable = {} # user defined variables

		self._update_domain_size()

		#self.float_symbols = {self.h[0]:0.1,self.h[1]:0.1,self.h[2]:0.1,self.dt:1.0} # dictionary to hold symbols and their values
		#self.int_symbols = {self.ntsteps:1,self.dim[0]:0,self.dim[1]:0,self.dim[2]:0}

	def _update_domain_size(self):
		# set dimension symbols, dim1, dim2, ...
		self.dim = [Variable('dim'+str(k+1), int(self.size[k]/self.spacing[k].value)+1+self.margin.value*2, 'int', True)  for k in range(self.dimension)]

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
		self.sfields = sfields

	def set_velocity_fields(self, vfields):
		self.vfields = vfields

	def calc_derivatives(self):
		l = [self.dt] + self.spacing
		for field in self.sfields+self.vfields:
			for k in range(self.dimension+1):
				h = l[k]
				for o in range(1,self.order[k]+1):
					field.calc_derivative([self.t]+self.index,k,h,o)

	def solve_fd(self,eqs):
		t, x, y, z = self.index
		self.eqs = eqs
		for field, eq in zip(self.vfields+self.sfields, eqs):
			field.set_fd(solve(eq,field.name[t+hf,x,y,z])[0].subs({t:t+hf}))

	def set_free_surface_boundary(self, d, side):
		for field in self.sfields:
			if side==0:
				field.set_free_surface_stress(self.index, d, self.margin, side)
			else:
				field.set_free_surface_stress(self.index, d, self.dim[d-1]-self.margin-1, side)

		for field in self.vfields:
			if side==0:
				field.set_free_surface_velocity(self.index, d, self.margin, side)
			else:
				field.set_free_surface_velocity(self.index, d, self.dim[d-1]-self.margin-1, side)

	def define_const(self):
		result = ''
		for v in self.float_symbols:
			result += 'const float ' + v.name + ' = ' + str(self.float_symbols[v]) + ';\n'
		for v in self.int_symbols:
			result += 'const int ' + v.name + ' = ' + str(self.int_symbols[v]) + ';\n'
		return result

	def declare_fields(self):
		result = ''
		for field in self.sfields + self.vfields:
			vec = '_' + ccode(field.name.label) + '_vec'
			result += 'std::vector<float> ' + vec + '(2*dimx*dimy*dimz);\n'
			result += 'float (*' + ccode(field.name.label) + ')[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) ' + vec + '.data();\n'
		return result

	def initialise(self):
		tmpl = self.lookup.get_template('init_loop_3d.txt')
		i, j, k = symbols('i j k')
		ijk = [i,j,k]
		m = self.margin
		t = self.index[0]
		xyzvalue = [None]*3
		ijk0 = [None]*3
		ijk1 = [None]*3

		result = ''

		for field in self.sfields+self.vfields:
			# populate xvalue, yvalue zvalue code
			for d in range(3):
				if not field.offset[d+1]:
					xyzvalue[d] = ccode((ijk[d]-m)*self.h[d])
					ijk0[d] = ccode(m)
					ijk1[d] = ccode(self.dim[d]-m)
				else:
					xyzvalue[d] = ccode((ijk[d]-m+0.5)*self.h[d])
					ijk0[d] = ccode(m) # same as first case
					ijk1[d] = ccode(self.dim[d]-m-1)

			t0 = self.float_symbols[self.dt]/2.0 if field.offset[0] else 0.0
			body = ccode(field.name[0,i,j,k]) + '=' + ccode(field.func.subs(t,t0)) + ';'
			dict1 = {'i':'i','j':'j','k':'k','i0':ijk0[0],'i1':ijk1[0],'j0':ijk0[1],'j1':ijk1[1],'k0':ijk0[2],'k1':ijk1[2],'xvalue':xyzvalue[0],'yvalue':xyzvalue[1],'zvalue':xyzvalue[2],'body':body}
			result += render(tmpl, dict1)

		return result

	def initialise_boundary(self):
		result = self.stress_bc().replace('[t1]','[0]')
		result += self.velocity_bc().replace('[t1]','[0]')
		return result

	def stress_loop(self):
		tmpl = self.lookup.get_template('update_loop_3d.txt')
		t, x, y, z = self.index
		t1 = Symbol('t1')
		m = self.margin
		i1 = ccode(self.dim[0]-m)
		j1 = ccode(self.dim[1]-m)
		k1 = ccode(self.dim[2]-m)
		body = ''
		for field in self.sfields:
			body += ccode(field.name[t1,x,y,z]) + '=' + ccode(field.shift_fd) + ';\n\t\t\t'
		dict1 = {'i':'x','j':'y','k':'z','i0':m,'i1':i1,'j0':m,'j1':j1,'k0':m,'k1':k1,'body':body}
		result = render(tmpl, dict1)
		return result

	def velocity_loop(self):
		tmpl = self.lookup.get_template('update_loop_3d.txt')
		t, x, y, z = self.index
		t1 = Symbol('t1')
		m = self.margin
		i1 = ccode(self.dim[0]-m)
		j1 = ccode(self.dim[1]-m)
		k1 = ccode(self.dim[2]-m)
		body = ''
		for field in self.vfields:
			body += ccode(field.name[t1,x,y,z]) + '=' + ccode(field.shift_fd.replace(t+1,t1)) + ';\n\t\t\t'
		dict1 = {'i':'x','j':'y','k':'z','i0':m,'i1':i1,'j0':m,'j1':j1,'k0':m,'k1':k1,'body':body}
		result = render(tmpl, dict1)
		return result

	def _get_ij(self, d):
		if d==1:
			i = 'y'
			j = 'z'
			i1 = self.dim[1]
			j1 = self.dim[2]
		elif d==2:
			i = 'x'
			j = 'z'
			i1 = self.dim[0]
			j1 = self.dim[2]
		else:
			i = 'x'
			j = 'y'
			i1 = self.dim[0]
			j1 = self.dim[1]						
		return i,j,i1,j1

	def stress_bc(self):
		tmpl = self.lookup.get_template('ghost_cell_3d.txt')
		result = ''
		for field in self.sfields:
			for d in range(1,4):
				for side in range(2):
					i,j,i1,j1 = self._get_ij(d)
					dict1 = {'i':i,'j':j,'i0':0,'i1':i1,'j0':0,'j1':j1,'body':field.bc[d][side]}
					result += render(tmpl, dict1).replace('[t]','[t1]')


		return result

	def velocity_bc(self):
		tmpl = self.lookup.get_template('ghost_cell_3d.txt')
		result = ''
		for d in range(1,4):
			# update the staggered field first because other fields depends on it
			sequence = [f for f in self.vfields if f.offset[d]] + [f for f in self.vfields if not f.offset[d]]
			for field in sequence:
				for side in range(2):
					i,j,i1,j1 = self._get_ij(d)
					dict1 = {'i':i,'j':j,'i0':1,'i1':i1-1,'j0':1,'j1':j1-1,'body':field.bc[d][side]}
					result += render(tmpl, dict1).replace('[t]','[t1]')

		return result

	def output_step(self):
		result = ''
		result += self.vfields[0].save_field()
		return result

	def converge_test(self):
		tmpl = self.lookup.get_template('init_loop_3d.txt')
		i, j, k = symbols('i j k')
		ijk = [i,j,k]
		m = self.margin
		t = self.index[0]
		xyzvalue = [None]*3
		ijk0 = [None]*3
		ijk1 = [None]*3
		dt = self.float_symbols[self.dt]
		ntsteps = self.int_symbols[self.ntsteps]
		ti = 0 if ntsteps%2 == 0 else 1

		result = ''

		for field in self.sfields+self.vfields:
			# populate xvalue, yvalue zvalue code
			l2 = ccode(field.name.label)+'_l2'
			result += 'float ' + l2 + ' = 0.0;\n\t\t'
			for d in range(3):
				if not field.offset[d+1]:
					xyzvalue[d] = ccode((ijk[d]-m)*self.h[d])
					ijk0[d] = ccode(m)
					ijk1[d] = ccode(self.dim[d]-m)
				else:
					xyzvalue[d] = ccode((ijk[d]-m+0.5)*self.h[d])
					ijk0[d] = ccode(m) # same as first case
					ijk1[d] = ccode(self.dim[d]-m-1)

			tn = dt*ntsteps if not field.offset[0] else dt*ntsteps + dt/2.0
			body = l2 + '+=' + ccode((field.name[ti,i,j,k]-(field.func.subs(t,tn)))**2.0) + ';'
			dict1 = {'i':'i','j':'j','k':'k','i0':ijk0[0],'i1':ijk1[0],'j0':ijk0[1],'j1':ijk1[1],'k0':ijk0[2],'k1':ijk1[2],'xvalue':xyzvalue[0],'yvalue':xyzvalue[1],'zvalue':xyzvalue[2],'body':body}
			result += render(tmpl, dict1)
			result += 'printf("' + l2 + ' = %.10f\\n", ' + l2 + ');\n\t\t'

		return result
