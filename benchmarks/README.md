# Propagator benchmarks

Set of benchmarking tools for evaluating the performance of generated
propagator codes. It provides two scripts, one for executing
benchmarks on a given target system and one for plotting the recorded
results. The scripts utilise the
[pybench](https://github.com/firedrakeproject/pybench) package
developed by the [Firedrake](http://www.firedrakeproject.org)
group. To install this dependency via pip:
```
pip install git+https://github.com/firedrakeproject/pybench
```
Alternatively, clone the package directly and add to `PYTHONPATH`:
```
git clone https://github.com/firedrakeproject/pybench
export PYTHONPATH=$PYTHONPATH:$PWD/pybench
```

## Benchmarking
To run benchmarks with various parameters and record results run:
```
python prop_bench.py -b -l -s -- <param1>=<val1> <param2>=<val2>
```
The recorded timings will be stored in a `results` directory. The
result files are indexed by the different parameter values, which
allows easy sweeping of a parameter space, for example:
```
for CC in icpc g++ clang; do
    for OPT in 2 3; do
        python prop_bench.py -b -l -s -- compiler=$CC opt_level $OPT
    done
done
```

## Plotting
Note that plotting the results might be not carried out on the machine
the benchmarks were run on, in which case the `results` directory
simply needs to be copied over:
```
mkdir machinex
scp -r me@machinex:path/to/propagator/benchmarks/results machinex
```
A bar chart comparison of the different parameter values (compiler
flags in the above example) can now be plotted with:
```
python prop_plot.py -i machinex/results -o machinex/plots --compiler icpc g++ clang --opt_level 2 3
```
The plotting script can also plot a parallel scaling curve when given multiple values for the `-p/--parallel` argument, for example:
```
python prop_plot.py -i machinex/results -o machinex/plots -p 1 2 4 8 --compiler icpc --opt_level 2 3
```
Note that the number of threads used during benchmarking is set via
the the `nthreads` parameter in `prop_bench.py`. For a list of other
parameters, please refer to the code.
