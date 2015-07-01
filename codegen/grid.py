from sympy import Symbol
from mako.template import Template
from mako.lookup import TemplateLookup
from mako.runtime import Context
from StringIO import StringIO


def render(tmpl, dict1):
	buf = StringIO()
	ctx = Context(buf, **dict1)
	tmpl.render_context(ctx)
	return buf.getvalue()

class StaggeredGrid2D:
	"""description of staggered grid for finite difference method"""

	lookup = TemplateLookup(directories=['templates/staggered/'])
	functions={} # empty dictionary for mapping fields with functions


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

	def initialise(self):
		tmpl = self.lookup.get_template('generic_loop_2d.txt')
		dict1 = {'i':'x','j':'y','l_i':0,'h_i':10,'l_j':0,'h_j':10,'body':'test'}
		return render(tmpl, dict1)

	def eval_function(self, field, t):
		return self.functions[field].subs(Symbol('t'),t)
