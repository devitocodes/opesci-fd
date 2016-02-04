import cgen
import os

def common_include():
    libraries = ['cassert', 'cstdlib', 'cmath', 'iostream', 'fstream', 'vector', 'cstdio', 'string']
    statements = [cgen.Include(s) for s in libraries]
    return statements

def io_include():
    libraries = ['opesciIO.h', 'opesciHandy.h']
    statements = [cgen.Include(s, False) for s in libraries]
    return statements

def profiling_include():
    libraries = ['opesciProfiling.h']
    statements = [cgen.Include(s, False) for s in libraries]
    return statements

def pluto_include():
    libraries  = ['omp.h', 'math.h']
    statements = [cgen.Include(s) for s in libraries]
    statements += [cgen.Define("ceild(n,d)", "ceil(((double)(n))/((double)(d)))")]
    statements += [cgen.Define("floord(n,d)", "floor(((double)(n))/((double)(d)))")]
    statements += [cgen.Define("max(x,y)", "((x) > (y)? (x) : (y))")]
    statements += [cgen.Define("min(x,y)", "((x) < (y)? (x) : (y))")]
    return statements

def copyright():
    copyright_directory = os.path.dirname(os.path.realpath(__file__))
    copyright_filename = 'copyright.txt'
    content = []
    with open(copyright_directory+"/"+copyright_filename) as f:
        content = f.readlines()
    statements = [cgen.Line(s) for s in content]
    return statements
