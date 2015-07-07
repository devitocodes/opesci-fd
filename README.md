# propagator

FD simulation for 3D elastic wave equation

running the python script codegen/grid3d.py (this script is the same as the notebook grid3d.ipynb), this will generate test.cpp, which can be compiled with g++ -fopenmp

The outputs are L2 norms between numerical and analytical solutions for the stress and velocity fields

symbolic subroutines are defined in grid.py

