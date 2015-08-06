from grid_test import *
import sys

def default():
	domain_size = (1.0,1.0,1.0)
	spacing = (0.01,0.01,0.01)
	dt = 0.002
	tmax = 2.0
	run_test(domain_size, spacing, dt, tmax, o_step=False, o_converge=True, omp=True, simd=False, ivdep=True, vtk=False, filename='src/tests/test3d.cpp')

def cx1():
	domain_size = (1.0,1.0,1.0)
	spacing = (0.005,0.005,0.005)
	dt = 0.001
	tmax = 5.0
	run_test(domain_size, spacing, dt, tmax, o_step=False, o_converge=False, omp=True, simd=False, ivdep=True, vtk=False, filename='src/tests/test3d_ivdep.cpp')
	run_test(domain_size, spacing, dt, tmax, o_step=False, o_converge=False, omp=True, simd=True, ivdep=False, vtk=False, filename='src/tests/test3d_simd.cpp')

def main():
	if len(sys.argv)>2:
		print 'too many arguments!'
		return
	if len(sys.argv)==1:
		default()
		return
	if sys.argv[1]=='cx1':
		cx1()

if __name__ == "__main__":
	main()