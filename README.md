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

This will generate the mode source code in `tests/src/eigenwave3d.cpp`.
The source code can be compiled and executed either manually or 
automatically.

##### Automatic compilation and execution

Opesci-FD provides automatic compilation and execution that allows 
developers to test their code directly from the Python environment. 
To compile the generted souce code:
```
python tests/eigenwave3d.py --compiler <cc>
```
where `<cc>` is either `g++`,`clang` or `icpc`, indicating which compiler to
use.Make sure your [clang](http://clang-omp.github.io/) compiler support openmp 
for multithreads program.
To compile and execute the above test case in parallel run:
```
python tests/eigenwave3d.py --compiler <cc> --execute --nthreads <nt>
```
where `<nt>` is the number of threads to use during execution. 
For additional options please see:
```
python tests/eigenwave3d.py --help
```

##### Manual compilation

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
