from sympy import factorial, Symbol, symbols, Matrix, IndexedBase, Idx, Eq, solve
from sympy.printing.ccode import CCodePrinter
from matplotlib import animation
import numpy as np
import matplotlib.pyplot as plt

class MyCPrinter(CCodePrinter):
	def __init__(self, settings={}):
		CCodePrinter.__init__(self, settings)

	def _print_Indexed(self, expr):
		# Array base and append indices
		output = self._print(expr.base.label) + ''.join(['[' + self._print(x) + ']' for x in expr.indices])
		return output

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

def print_myccode(expr, assign_to=None, **settings):

    return MyCPrinter(settings).doprint(expr, assign_to)

def main():
	dx, dt, x, y, z, t, c = symbols('dx dt x y z t c')
	U = IndexedBase('U')
	n = 2
	Uxx = Deriv(U,[x,y,z,t],0,dx,n)[2]
	Utt = Deriv(U,[x,y,z,t],3,dt,n)[2]
	eq = Eq(Utt, (c**2)*Uxx)
	code = print_myccode(U[x,y,z,t+1]) + "=" + print_myccode(solve(eq, U[x,y,z,t+1])[0])
	print(code)

if __name__ == "__main__":
    main()