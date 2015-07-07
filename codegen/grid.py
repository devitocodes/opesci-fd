from sympy import Symbol, symbols, factorial, Matrix, Rational, Indexed, simplify, IndexedBase
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

def print_myccode(expr, **settings):
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

class StaggeredGrid2D:
	"""description of staggered grid for finite difference method"""

	lookup = TemplateLookup(directories=['templates/staggered/'])
	functions = {} # empty dictionary for mapping fields with functions
	accuracy_time = 1
	accuracy_space = 2
	Dt = {} # dictionary for mapping fields to their time derivatives
	Dx = {} # dictionary for mapping fields to their first derivatives
	Dy = {} # dictionary for mapping fields to their first derivatives
	Dx_2 = {} # first order spatial derivatives
	Dy_2 = {} # first order spatial derivatives
	fd = {} # dictionary for mapping fields to fd expression for t+1
	fd_shifted = {} # dictionary for mapping code
	bc = [[{},{}],[{},{}]] # list of list of dictionary
	index = list(symbols('t x y'))

	def __init__(self):
		return

	def assign_grid_size(self, h, dt):
		self.h = h
		self.dt = dt

	def assign_dimensions(self, dim):
		self.dim = dim

	def assign_stress_fields(self, Txx, Tyy, Txy):
		self.Txx = Txx
		self.Tyy = Tyy
		self.Txy = Txy

	def assign_velocity_fields(self, U,V):
		self.U = U
		self.V = V

	def assign_function(self, field, function):
		self.functions[field] = function

	def calc_derivatives(self):
		for field in [self.U, self.V, self.Txx, self.Tyy, self.Txy]:
			self.Dt[field] = Deriv_half(field, self.index, 0, self.dt, self.accuracy_time)[1]
			self.Dx[field] = Deriv_half(field, self.index, 1, self.h[0], self.accuracy_space)[1]
			self.Dy[field] = Deriv_half(field, self.index, 2, self.h[1], self.accuracy_space)[1]
			self.Dx_2[field] = Deriv_half(field, self.index, 1, self.h[0], 1)[1]
			self.Dy_2[field] = Deriv_half(field, self.index, 2, self.h[1], 1)[1]

	def assign_fd(self, field, fd):
		self.fd[field] = fd
		self.fd_shifted[field] = shift_grid(fd)

	def assign_bc(self, field, dim, side, expr):
		self.bc[dim][side][field] = expr

	def eval_function(self, field, t):
		return self.functions[field].subs(self.t,t)

	def _update(self):
		self.margin = 2
		self.l_x = self.margin
		self.l_y = self.margin
		self.h_x_long = print_myccode(self.dim[0]-self.margin)
		self.h_x_short = print_myccode(self.dim[0]-self.margin-1)
		self.h_y_long = print_myccode(self.dim[1]-self.margin)
		self.h_y_short = print_myccode(self.dim[1]-self.margin-1)


	def initialise(self):
		self._update()
		tmpl = self.lookup.get_template('generic_loop_2d.txt')
		i, j = symbols('i j')
		m = self.margin
		t = self.index[0]

		x_value_whole = print_myccode((i-m)*self.h[0])
		y_value_whole = print_myccode((j-m)*self.h[1])
		x_value_half = print_myccode((i-m+0.5)*self.h[0])
		y_value_half = print_myccode((j-m+0.5)*self.h[1])

		# Txx
		body = print_myccode(self.Txx[0,i,j]) + '=' + print_myccode(self.functions[self.Txx].subs(t,0)) + ';'
		dict1 = {'i':'i','j':'j','l_i':self.l_x,'h_i':self.h_x_long,'l_j':self.l_y,'h_j':self.h_y_long,'x_value':x_value_whole,'y_value':y_value_whole,'body':body}
		result = render(tmpl, dict1)
		# Tyy
		body = print_myccode(self.Tyy[0,i,j]) + '=' + print_myccode(self.functions[self.Tyy].subs(t,0)) + ';'
		dict1.update({'body':body})
		result += render(tmpl, dict1)
		# Txy
		body = print_myccode(self.Txy[0,i,j]) + '=' + print_myccode(self.functions[self.Txy].subs(t,0)) + ';'
		dict1.update({'h_i':self.h_x_short,'h_j':self.h_y_short,'x_value':x_value_half,'y_value':y_value_half,'body':body})
		result += render(tmpl, dict1)
		# U
		body = print_myccode(self.U[0,i,j]) + '=' + print_myccode(self.functions[self.U].subs(t,self.dt/2)) + ';'
		dict1.update({'h_i':self.h_x_short,'h_j':self.h_y_long,'x_value':x_value_half,'y_value':y_value_whole,'body':body})
		result += render(tmpl, dict1)
		# V
		body = print_myccode(self.V[0,i,j]) + '=' + print_myccode(self.functions[self.V].subs(t,self.dt/2)) + ';'
		dict1.update({'h_i':self.h_x_long,'h_j':self.h_y_short,'x_value':x_value_whole,'y_value':y_value_half,'body':body})
		result += render(tmpl, dict1)
		result += self.initialise_boundary()

		return result

	def initialise_boundary(self):
		tmpl = self.lookup.get_template('ghost_stress_x.txt')
		dict1 = {'dimx':self.dim[0],'dimy':self.dim[1],'t':0}
		result = render(tmpl, dict1)

		tmpl = self.lookup.get_template('ghost_stress_y.txt')
		result += render(tmpl, dict1)
		result += self._velocity_bc(0)
		return result

	def stress_loop(self):
		self._update()
		tmpl = self.lookup.get_template('update_loop_2d.txt')
		t, x, y = self.index
		t1 = Symbol('t1')

		body = print_myccode(self.Txx[t1,x,y]) + '=' + print_myccode(self.fd_shifted[self.Txx]) + ';\n\t\t\t'
		body += print_myccode(self.Tyy[t1,x,y]) + '=' + print_myccode(self.fd_shifted[self.Tyy]) + ';\n\t\t\t'
		body += print_myccode(self.Txy[t1,x,y]) + '=' + print_myccode(self.fd_shifted[self.Txy]) + ';'
		dict1 = {'i':'x','j':'y','l_i':self.l_x,'h_i':self.h_x_long,'l_j':self.l_y,'h_j':self.h_y_long,'body':body}
		result = render(tmpl, dict1)
		return result

	def velocity_loop(self):
		self._update()
		tmpl = self.lookup.get_template('update_loop_2d.txt')
		t, x, y = self.index
		t1 = Symbol('t1')

		body = print_myccode(self.U[t1,x,y]) + '=' + print_myccode(simplify(self.fd_shifted[self.U]-self.U[t,x,y]).subs(t,t1)+self.U[t,x,y]) + ';\n\t\t\t'
		body += print_myccode(self.V[t1,x,y]) + '=' + print_myccode(simplify(self.fd_shifted[self.V]-self.V[t,x,y]).subs(t,t1)+self.V[t,x,y]) + ';'
		dict1 = {'i':'x','j':'y','l_i':self.l_x,'h_i':self.h_x_long,'l_j':self.l_y,'h_j':self.h_y_long,'body':body}
		result = render(tmpl, dict1)
		return result

	def stress_bc(self):
		tmpl = self.lookup.get_template('ghost_stress_x.txt')
		dict1 = {'dimx':self.dim[0],'dimy':self.dim[1],'t':'t1'}
		result = render(tmpl, dict1)

		tmpl = self.lookup.get_template('ghost_stress_y.txt')
		result += render(tmpl, dict1)
		return result

	def velocity_bc(self):
		t1 = Symbol('t1')
		return self._velocity_bc(t1)

	def _velocity_bc(self, t1):
		tmpl = self.lookup.get_template('ghost_velocity_x.txt')
		t, x, y = self.index

		U0 = print_myccode(self.U[t1,self.margin-1,y]) + '=' + print_myccode(shift_grid(self.bc[0][0][self.U].subs({t:t1,x:1+hf}))) +';\n\t\t\t'
		U1 = print_myccode(self.U[t1,self.dim[0]-self.margin-1,y]) + '=' + print_myccode(shift_grid(self.bc[0][1][self.U].subs({t:t1,x:self.dim[0]-2-hf}))) +';'
		dict1 = {'dimx':self.dim[0],'dimy':self.dim[1],'body':U0+U1}
		result = render(tmpl, dict1)
		V0 = print_myccode(self.V[t1,self.margin-1,y]) + '=' + print_myccode(shift_grid(self.bc[0][0][self.V].subs({t:t1,x:1}))) +';\n\t\t\t'
		V1 = print_myccode(self.V[t1,self.dim[0]-self.margin,y]) + '=' + print_myccode(shift_grid(self.bc[0][1][self.V].subs({t:t1,x:self.dim[0]-2}))) +';'
		dict1 = {'dimx':self.dim[0],'dimy':self.dim[1],'body':V0+V1}
		result += render(tmpl, dict1)

		tmpl = self.lookup.get_template('ghost_velocity_y.txt')
		V0 = print_myccode(self.V[t1,x,self.margin-1]) + '=' + print_myccode(shift_grid(self.bc[1][0][self.V].subs({t:t1,y:1+hf}))) +';\n\t\t\t'
		V1 = print_myccode(self.V[t1,x,self.dim[1]-self.margin-1]) + '=' + print_myccode(shift_grid(self.bc[1][1][self.V].subs({t:t1,y:self.dim[1]-2-hf}))) +';'
		dict1 = {'dimx':self.dim[0],'dimy':self.dim[1],'body':V0+V1}
		result += render(tmpl, dict1)
		U0 = print_myccode(self.U[t1,x,self.margin-1]) + '=' + print_myccode(shift_grid(self.bc[1][0][self.U].subs({t:t1,y:1}))) +';\n\t\t\t'
		U1 = print_myccode(self.U[t1,x,self.dim[1]-self.margin]) + '=' + print_myccode(shift_grid(self.bc[1][1][self.U].subs({t:t1,y:self.dim[1]-2}))) +';'
		dict1 = {'dimx':self.dim[0],'dimy':self.dim[1],'body':U0+U1}
		result += render(tmpl, dict1)

		return result

	def converge_test(self):
		tmpl = self.lookup.get_template('generic_loop_2d_2.txt')
		i, j = symbols('i j')
		m = self.margin
		t = self.index[0]
		t1, tf1, tf2 = symbols('t1, tf1, tf2')

		x_value_whole = print_myccode((i-m)*self.h[0])
		y_value_whole = print_myccode((j-m)*self.h[1])
		x_value_half = print_myccode((i-m+0.5)*self.h[0])
		y_value_half = print_myccode((j-m+0.5)*self.h[1])

		# Txx
		body = 'Txx_diff += pow(' + print_myccode(self.Txx[t1,i,j]) + '-(' + print_myccode(self.functions[self.Txx].subs(t,tf1)) + '),2);'
		dict1 = {'i':'i','j':'j','l_i':self.l_x,'h_i':self.h_x_long,'l_j':self.l_y,'h_j':self.h_y_long,'x_value':x_value_whole,'y_value':y_value_whole,'body':body}
		result = render(tmpl, dict1)
		# Tyy
		body = 'Tyy_diff += pow(' + print_myccode(self.Tyy[t1,i,j]) + '-(' + print_myccode(self.functions[self.Tyy].subs(t,tf1)) + '),2);'
		dict1.update({'body':body})
		result += render(tmpl, dict1)
		# Txy
		body = 'Txy_diff += pow(' + print_myccode(self.Txy[t1,i,j]) + '-(' + print_myccode(self.functions[self.Txy].subs(t,tf1)) + '),2);'
		dict1.update({'h_i':self.h_x_short,'h_j':self.h_y_short,'x_value':x_value_half,'y_value':y_value_half,'body':body})
		result += render(tmpl, dict1)
		# U
		body = 'U_diff += pow(' + print_myccode(self.U[t1,i,j]) + '-(' + print_myccode(self.functions[self.U].subs(t,tf2)) + '),2);'
		dict1.update({'h_i':self.h_x_short,'h_j':self.h_y_long,'x_value':x_value_half,'y_value':y_value_whole,'body':body})
		result += render(tmpl, dict1)
		# V
		body = 'V_diff += pow(' + print_myccode(self.V[t1,i,j]) + '-(' + print_myccode(self.functions[self.V].subs(t,tf2)) + '),2);'
		dict1.update({'h_i':self.h_x_long,'h_j':self.h_y_short,'x_value':x_value_whole,'y_value':y_value_half,'body':body})
		result += render(tmpl, dict1)

		return result


class Field:
	d = [[None]*4]*4 # list of list to store derivative expressions
	def __init__(self, name, offset):
		self.name = IndexedBase(name)
		self.offset = offset

	def set_analytic_func(self, function):
		self.func = function

	def calc_derivative(self, l, k, d, n):
		self.d[k][n] = Deriv_half(self.name, l, k, d, n)[1]

	def calc_fd(self, eq):
		self.fd = 

class StaggeredGrid3D:
	"""description of staggered grid for finite difference method"""

	lookup = TemplateLookup(directories=['templates/staggered/'])
	solution = {} # empty dictionary for mapping fields with functions
	order = [1,2,2,2]

	Dt = {} # dictionary for mapping fields to their time derivatives
	Dx = {} # dictionary for mapping fields to their first derivatives
	Dy = {} # dictionary for mapping fields to their first derivatives
	Dx_2 = {} # first order spatial derivatives
	Dy_2 = {} # first order spatial derivatives
	fd = {} # dictionary for mapping fields to fd expression for t+1
	fd_shifted = {} # dictionary for mapping code
	bc = [[{},{}],[{},{}]] # list of list of dictionary
	index = list(symbols('t x y z'))

	def __init__(self, h, dt, dim):
		self.h = symbols(h)
		self.dt = Symbol(dt)
		self.dim = symbols(dim)

	def set_stress_fields(self, sfields):
		self.sfields = sfields

	def set_velocity_fields(self, vfields):
		self.vfields = vfields

	def calc_derivatives(self):
		l = [self.dt] + list(self.h)
		for field in self.sfields+self.vfields:
			for k in range(3):
				field.calc_derivative(self.index,k,l[k],1)
				field.calc_derivative(self.index,k,l[k],2)

	def assign_fd(self, field, fd):
		self.fd[field] = fd
		self.fd_shifted[field] = shift_grid(fd)

	def assign_bc(self, field, dim, side, expr):
		self.bc[dim][side][field] = expr

	def eval_function(self, field, t):
		return self.functions[field].subs(self.t,t)

	def _update(self):
		self.margin = 2
		self.l_x = self.margin
		self.l_y = self.margin
		self.h_x_long = print_myccode(self.dim[0]-self.margin)
		self.h_x_short = print_myccode(self.dim[0]-self.margin-1)
		self.h_y_long = print_myccode(self.dim[1]-self.margin)
		self.h_y_short = print_myccode(self.dim[1]-self.margin-1)


	def initialise(self):
		self._update()
		tmpl = self.lookup.get_template('generic_loop_2d.txt')
		i, j = symbols('i j')
		m = self.margin
		t = self.index[0]

		x_value_whole = print_myccode((i-m)*self.h[0])
		y_value_whole = print_myccode((j-m)*self.h[1])
		x_value_half = print_myccode((i-m+0.5)*self.h[0])
		y_value_half = print_myccode((j-m+0.5)*self.h[1])

		# Txx
		body = print_myccode(self.Txx[0,i,j]) + '=' + print_myccode(self.functions[self.Txx].subs(t,0)) + ';'
		dict1 = {'i':'i','j':'j','l_i':self.l_x,'h_i':self.h_x_long,'l_j':self.l_y,'h_j':self.h_y_long,'x_value':x_value_whole,'y_value':y_value_whole,'body':body}
		result = render(tmpl, dict1)
		# Tyy
		body = print_myccode(self.Tyy[0,i,j]) + '=' + print_myccode(self.functions[self.Tyy].subs(t,0)) + ';'
		dict1.update({'body':body})
		result += render(tmpl, dict1)
		# Txy
		body = print_myccode(self.Txy[0,i,j]) + '=' + print_myccode(self.functions[self.Txy].subs(t,0)) + ';'
		dict1.update({'h_i':self.h_x_short,'h_j':self.h_y_short,'x_value':x_value_half,'y_value':y_value_half,'body':body})
		result += render(tmpl, dict1)
		# U
		body = print_myccode(self.U[0,i,j]) + '=' + print_myccode(self.functions[self.U].subs(t,self.dt/2)) + ';'
		dict1.update({'h_i':self.h_x_short,'h_j':self.h_y_long,'x_value':x_value_half,'y_value':y_value_whole,'body':body})
		result += render(tmpl, dict1)
		# V
		body = print_myccode(self.V[0,i,j]) + '=' + print_myccode(self.functions[self.V].subs(t,self.dt/2)) + ';'
		dict1.update({'h_i':self.h_x_long,'h_j':self.h_y_short,'x_value':x_value_whole,'y_value':y_value_half,'body':body})
		result += render(tmpl, dict1)
		result += self.initialise_boundary()

		return result

	def initialise_boundary(self):
		tmpl = self.lookup.get_template('ghost_stress_x.txt')
		dict1 = {'dimx':self.dim[0],'dimy':self.dim[1],'t':0}
		result = render(tmpl, dict1)

		tmpl = self.lookup.get_template('ghost_stress_y.txt')
		result += render(tmpl, dict1)
		result += self._velocity_bc(0)
		return result

	def stress_loop(self):
		self._update()
		tmpl = self.lookup.get_template('update_loop_2d.txt')
		t, x, y = self.index
		t1 = Symbol('t1')

		body = print_myccode(self.Txx[t1,x,y]) + '=' + print_myccode(self.fd_shifted[self.Txx]) + ';\n\t\t\t'
		body += print_myccode(self.Tyy[t1,x,y]) + '=' + print_myccode(self.fd_shifted[self.Tyy]) + ';\n\t\t\t'
		body += print_myccode(self.Txy[t1,x,y]) + '=' + print_myccode(self.fd_shifted[self.Txy]) + ';'
		dict1 = {'i':'x','j':'y','l_i':self.l_x,'h_i':self.h_x_long,'l_j':self.l_y,'h_j':self.h_y_long,'body':body}
		result = render(tmpl, dict1)
		return result

	def velocity_loop(self):
		self._update()
		tmpl = self.lookup.get_template('update_loop_2d.txt')
		t, x, y = self.index
		t1 = Symbol('t1')

		body = print_myccode(self.U[t1,x,y]) + '=' + print_myccode(simplify(self.fd_shifted[self.U]-self.U[t,x,y]).subs(t,t1)+self.U[t,x,y]) + ';\n\t\t\t'
		body += print_myccode(self.V[t1,x,y]) + '=' + print_myccode(simplify(self.fd_shifted[self.V]-self.V[t,x,y]).subs(t,t1)+self.V[t,x,y]) + ';'
		dict1 = {'i':'x','j':'y','l_i':self.l_x,'h_i':self.h_x_long,'l_j':self.l_y,'h_j':self.h_y_long,'body':body}
		result = render(tmpl, dict1)
		return result

	def stress_bc(self):
		tmpl = self.lookup.get_template('ghost_stress_x.txt')
		dict1 = {'dimx':self.dim[0],'dimy':self.dim[1],'t':'t1'}
		result = render(tmpl, dict1)

		tmpl = self.lookup.get_template('ghost_stress_y.txt')
		result += render(tmpl, dict1)
		return result

	def velocity_bc(self):
		t1 = Symbol('t1')
		return self._velocity_bc(t1)

	def _velocity_bc(self, t1):
		tmpl = self.lookup.get_template('ghost_velocity_x.txt')
		t, x, y = self.index

		U0 = print_myccode(self.U[t1,self.margin-1,y]) + '=' + print_myccode(shift_grid(self.bc[0][0][self.U].subs({t:t1,x:1+hf}))) +';\n\t\t\t'
		U1 = print_myccode(self.U[t1,self.dim[0]-self.margin-1,y]) + '=' + print_myccode(shift_grid(self.bc[0][1][self.U].subs({t:t1,x:self.dim[0]-2-hf}))) +';'
		dict1 = {'dimx':self.dim[0],'dimy':self.dim[1],'body':U0+U1}
		result = render(tmpl, dict1)
		V0 = print_myccode(self.V[t1,self.margin-1,y]) + '=' + print_myccode(shift_grid(self.bc[0][0][self.V].subs({t:t1,x:1}))) +';\n\t\t\t'
		V1 = print_myccode(self.V[t1,self.dim[0]-self.margin,y]) + '=' + print_myccode(shift_grid(self.bc[0][1][self.V].subs({t:t1,x:self.dim[0]-2}))) +';'
		dict1 = {'dimx':self.dim[0],'dimy':self.dim[1],'body':V0+V1}
		result += render(tmpl, dict1)

		tmpl = self.lookup.get_template('ghost_velocity_y.txt')
		V0 = print_myccode(self.V[t1,x,self.margin-1]) + '=' + print_myccode(shift_grid(self.bc[1][0][self.V].subs({t:t1,y:1+hf}))) +';\n\t\t\t'
		V1 = print_myccode(self.V[t1,x,self.dim[1]-self.margin-1]) + '=' + print_myccode(shift_grid(self.bc[1][1][self.V].subs({t:t1,y:self.dim[1]-2-hf}))) +';'
		dict1 = {'dimx':self.dim[0],'dimy':self.dim[1],'body':V0+V1}
		result += render(tmpl, dict1)
		U0 = print_myccode(self.U[t1,x,self.margin-1]) + '=' + print_myccode(shift_grid(self.bc[1][0][self.U].subs({t:t1,y:1}))) +';\n\t\t\t'
		U1 = print_myccode(self.U[t1,x,self.dim[1]-self.margin]) + '=' + print_myccode(shift_grid(self.bc[1][1][self.U].subs({t:t1,y:self.dim[1]-2}))) +';'
		dict1 = {'dimx':self.dim[0],'dimy':self.dim[1],'body':U0+U1}
		result += render(tmpl, dict1)

		return result

	def converge_test(self):
		tmpl = self.lookup.get_template('generic_loop_2d_2.txt')
		i, j = symbols('i j')
		m = self.margin
		t = self.index[0]
		t1, tf1, tf2 = symbols('t1, tf1, tf2')

		x_value_whole = print_myccode((i-m)*self.h[0])
		y_value_whole = print_myccode((j-m)*self.h[1])
		x_value_half = print_myccode((i-m+0.5)*self.h[0])
		y_value_half = print_myccode((j-m+0.5)*self.h[1])

		# Txx
		body = 'Txx_diff += pow(' + print_myccode(self.Txx[t1,i,j]) + '-(' + print_myccode(self.functions[self.Txx].subs(t,tf1)) + '),2);'
		dict1 = {'i':'i','j':'j','l_i':self.l_x,'h_i':self.h_x_long,'l_j':self.l_y,'h_j':self.h_y_long,'x_value':x_value_whole,'y_value':y_value_whole,'body':body}
		result = render(tmpl, dict1)
		# Tyy
		body = 'Tyy_diff += pow(' + print_myccode(self.Tyy[t1,i,j]) + '-(' + print_myccode(self.functions[self.Tyy].subs(t,tf1)) + '),2);'
		dict1.update({'body':body})
		result += render(tmpl, dict1)
		# Txy
		body = 'Txy_diff += pow(' + print_myccode(self.Txy[t1,i,j]) + '-(' + print_myccode(self.functions[self.Txy].subs(t,tf1)) + '),2);'
		dict1.update({'h_i':self.h_x_short,'h_j':self.h_y_short,'x_value':x_value_half,'y_value':y_value_half,'body':body})
		result += render(tmpl, dict1)
		# U
		body = 'U_diff += pow(' + print_myccode(self.U[t1,i,j]) + '-(' + print_myccode(self.functions[self.U].subs(t,tf2)) + '),2);'
		dict1.update({'h_i':self.h_x_short,'h_j':self.h_y_long,'x_value':x_value_half,'y_value':y_value_whole,'body':body})
		result += render(tmpl, dict1)
		# V
		body = 'V_diff += pow(' + print_myccode(self.V[t1,i,j]) + '-(' + print_myccode(self.functions[self.V].subs(t,tf2)) + '),2);'
		dict1.update({'h_i':self.h_x_long,'h_j':self.h_y_short,'x_value':x_value_whole,'y_value':y_value_half,'body':body})
		result += render(tmpl, dict1)

		return result