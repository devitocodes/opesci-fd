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

##### Profiling

If the PAPI library is found on your system during the initial build,
Opesci-FD can also provide profiling information, such as the achieved
number of MFlops/s during automated runs. To enable this feature
simply add the `--profiling` flag to the example command above. The
user can also supply a list of PAPI event names, for example
`PAPI_TOT_CYC` or `PAPI_FP_OPS`, via the `--papi-events` flag:
```
python tests/eigenwave3d.py -c g++ -x --n 4 --profiling --papi-events PAPI_TOT_CYC PAPI_FP_OPS
```

Please note that the availability of PAPI events is highly dependent
on the hardware you are running on and the local PAPI install.

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


##### Auto-tuning for pluto

Here is a script to test the best tile size to get best optimisation effect:
```
python tests/pluto_tile_test.py -b -s -l
```
The results will be in "results" folder.
The list "sizes" contains all the tile sizes will be tested. 
This script depends on [pybench.py](https://github.com/firedrakeproject/pybench)

