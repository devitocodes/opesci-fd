# propagator

FD simulation for 3D elastic wave equation.

running the python script codegen/grid3d.py (this script is the same as the notebook grid3d.ipynb), this will generate test3d.cpp. This can be compiled manually with g++ -fopenmp for instance, or build the make file according to the .travis.yml script.

Testing with eigen waves on unit cube is implemented in grid_test.py. Run run_test() with relevant parameters to generate the source code. Some test cases with common settings are defined in grid3d.py, with explanation.

The outputs are grid sizes (dimx, dimy, dimz), followed by L2 norms between numerical and analytical solutions for the stress and velocity fields. To switch off this output, change output_convergence from True to False in main().

symbolic subroutines are defined in grid.py
