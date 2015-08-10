# propagator

FD simulation for 3D elastic wave equation.

running the python script codegen/grid3d.py (this script is the same as the notebook grid3d.ipynb), this will generate test3d.cpp. This can be compiled manually with g++ -fopenmp for instance, or build the make file according to the .travis.yml scrip.

To enable vtk output, go to grid3d.py and change vtk from False to True in main(), and change output_step to True. Compile with cmake in this case.

Other grid parameters can be set in the grid3d.py script as well.

The outputs are grid spacing (dx, dy, dz), followed by L2 norms between numerical and analytical solutions for the stress and velocity fields. To switch off this output, change output_convergence from True to False in main().

symbolic subroutines are defined in grid.py
