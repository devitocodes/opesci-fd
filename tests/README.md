### Command

```
python eigenwave_bench.py -b -l -s -- compiler=g++ opt_level=3  nthreads=4

python eigenwave_plot.py roofline --compiler g++ --opt_level 3 --parallel 4  --ai 0.85 --precision single

```

### Compare two or more codes

```

python eigenwave_plot.py roofline --basename <code1> <code2> <code3>  --ai <ai_code1> <ai_code2> <ai_code3> --compiler g++ --opt_level 3 --parallel 4  --precision single

```
