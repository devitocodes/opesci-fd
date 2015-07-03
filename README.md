# propagator

good place to start is probably the notebook codegen/refactor.ipynb, running it will generate test.cpp, to be compiled with C++ with OpenMP.

It outputs L2 differences between numerical and analytical solutions for the 5 fields.

the symbolic library is grid.py
fdlib.py is the old python library, I'm gradually migrating towards grid.py
