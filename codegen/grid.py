from sympy import Symbol, symbols, factorial, Matrix, Rational
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



class StaggeredGrid2D:
	"""description of staggered grid for finite difference method"""

	lookup = TemplateLookup(directories=['templates/staggered/'])
	functions = {} # empty dictionary for mapping fields with functions
	accuracy_time = 1
	accuracy_space = 2
	Dt = {} # dictionary for mapping fields to their time derivatives
	Dx = {} # dictionary for mapping fields to their first derivatives
	Dy = {} # dictionary for mapping fields to their first derivatives
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


	def eval_function(self, field, t):
		return self.functions[field].subs(self.t,t)

	def initialise(self):
		tmpl = self.lookup.get_template('generic_loop_2d.txt')
		i, j = symbols('i j')
		margin = 2
		l_i = margin
		l_j = margin
		t = self.index[0]

		h_i_long = print_myccode(self.dim[0]-margin)
		h_i_short = print_myccode(self.dim[0]-margin-1)
		h_j_long = print_myccode(self.dim[1]-margin)
		h_j_short = print_myccode(self.dim[1]-margin-1)

		x_value_whole = print_myccode((i-margin)*self.h[0])
		y_value_whole = print_myccode((j-margin)*self.h[1])
		x_value_half = print_myccode((i-margin+0.5)*self.h[0])
		y_value_half = print_myccode((j-margin+0.5)*self.h[1])

		# Txx
		body = print_myccode(self.Txx[0,i,j]) + '=' + print_myccode(self.functions[self.Txx].subs(t,0))
		dict1 = {'i':'i','j':'j','l_i':l_i,'h_i':h_i_long,'l_j':l_j,'h_j':h_j_long,'x_value':x_value_whole,'y_value':y_value_whole,'body':body}
		result = render(tmpl, dict1)
		# Tyy
		body = print_myccode(self.Tyy[0,i,j]) + '=' + print_myccode(self.functions[self.Tyy].subs(t,0))
		dict1.update({'body':body})
		result += render(tmpl, dict1)
		# Txy
		body = print_myccode(self.Txy[0,i,j]) + '=' + print_myccode(self.functions[self.Txy].subs(t,0))
		dict1.update({'h_i':h_i_short,'h_j':h_j_short,'x_value':x_value_half,'y_value':y_value_half,'body':body})
		result += render(tmpl, dict1)
		# U
		body = print_myccode(self.U[0,i,j]) + '=' + print_myccode(self.functions[self.U].subs(t,self.dt/2))
		dict1.update({'h_i':h_i_short,'h_j':h_j_long,'x_value':x_value_half,'y_value':y_value_whole,'body':body})
		result += render(tmpl, dict1)
		# V
		body = print_myccode(self.V[0,i,j]) + '=' + print_myccode(self.functions[self.V].subs(t,self.dt/2))
		dict1.update({'h_i':h_i_long,'h_j':h_j_short,'x_value':x_value_whole,'y_value':y_value_half,'body':body})
		result += render(tmpl, dict1)

		return result

	def initialise_boundary(self):
		tmpl = self.lookup.get_template('ghost_stress_x.txt')
		dict1 = {'dimx':self.dim[0],'dimy':self.dim[1],'t':0}
		result = render(tmpl, dict1)

		tmpl = self.lookup.get_template('ghost_stress_y.txt')
		result += render(tmpl, dict1)
		return result
