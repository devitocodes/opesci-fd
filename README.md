# Opesci-FD

Opesci-FD is a software package to automatically generate Finite
Difference models from a high-level description of the model
equations. It allows the rapid development, analysis and optimisation
of propagator codes for use in seismic imaging.

### Installation

Opesci-FD can be installed with:
```
pip install --local git+https://github.com/opesci/opesci-fd.git
```
This will install the latest version of the `opesci` python package
locally. To get the latests updates for your local copy simply add
`--upgrade` to the above command. For a developer checkout of the code
run:
```
git clone https://github.com/opesci/opesci-fd.git
cd opesci-fd; pip install -r requirements.txt
export PYTHONPATH=`pwd`:$PYTHONPATH
python setup.py build_clib --build-clib=`pwd`
```

### Getting started

An example of the high-level model description is provided under
`tests/eigenwave3d.py`. This script generates a stencil code that
models the propagation of an eigenwave on a unit cube by solving the
3D elastic wave equation on a staggered grid. To generate the model
code run:
```
python tests/eigenwave3d.py
```

This will generate the mode source code in `tests/src/eigenwave3d.cpp`,
as well as the compiled executable `tests/src/eigenwave3d`. The
generated binary can now be executed with:
```
cd tests/src
./eigenwave3d
```

The generated source file can also be compiled by hand using the
provided CMake file:
```
mkdir tests/build
cd tests/build
cmake ../src
make
bin/eigenwave3d
```

The outputs are grid spacings (dx, dy, dz), followed by L2 norms
between numerical and analytical solutions for the stress and velocity
fields. To switch off this output, change output_convergence from True
to False in the Python test definition. Further parameter switches for
controlling model input, the data type used (single of double
precision or explicit vectorisation are also provided.
