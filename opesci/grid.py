from compilation import GNUCompiler, IntelCompiler
from codeprinter import ccode

from StringIO import StringIO
from mako.runtime import Context
from ctypes import cdll, Structure, POINTER, c_float, pointer


class Grid(object):
    """Base class for grid objects that provides the code generation,
    compilation and execution infrastructure"""

    # mako.lookup.TemplateLookup object
    lookup = None

    # Name of the base template
    template_base = None

    # List of keys used the template that corresponds to an equivalent
    # property field in the derived grid class
    template_keys = []

    # Default compiler is GNU
    _compiler = GNUCompiler()

    # Placeholders for generated code and associated files
    src_code = None
    src_file = None
    src_lib = None

    # Dynamic execution objects
    _library = None
    _arg_grid = None
    _arg_conv = None

    def __del__(self):
        # Correctly close compiled kernel library
        if self._library is not None:
            cdll.LoadLibrary('libdl.so').dlclose(self._library._handle)

    def _load_library(self, src_lib):
        """Load a compiled dynamic binary using ctypes.cdll"""
        libname = src_lib or self.src_lib
        try:
            self._library = cdll.LoadLibrary(libname)
        except OSError as e:
            print "Library load error: ", e
            raise Exception("Failed to load %s" % libname)

    @property
    def compiler(self):
        return self._compiler

    @compiler.setter
    def compiler(self, compiler):
        if compiler in ['g++', 'gnu']:
            self.compiler = GNUCompiler()
        elif compiler in ['icpc', 'intel']:
            self.compiler = IntelCompiler()
        else:
            raise ValueError("Unknown compiler.")

    def generate(self, filename, compiler=None):
        if compiler:
            compiler = compiler

        # Generate a dictionary that maps template keys to code fragments
        template = self.lookup.get_template(self.template_base)
        template_keys = dict([(k, getattr(self, k)) for k in self.template_keys])

        # Render code from provided template
        buf = StringIO()
        ctx = Context(buf, **template_keys)
        template.render_context(ctx)
        self.src_code = buf.getvalue()

        # Generate compilable source code
        self.src_file = filename
        with file(self.src_file, 'w') as f:
            f.write(self.src_code)

        print "Generated:", self.src_file

    def compile(self, filename, compiler=None, shared=True):
        if compiler:
            compiler = compiler

        # Generate code if this hasn't been done yet
        if self.src_file is None:
            self.generate(filename)

        # Compile source file
        out = self.compiler.compile(self.src_file, shared=shared)
        if shared:
            self.src_lib = out

    def execute(self, filename, compiler='g++'):
        # Compile code if this hasn't been done yet
        if self.src_lib is None:
            self.compile(filename, compiler=compiler, shared=True)

        # Load compiled binary
        self._load_library(src_lib=self.src_lib)

        # Define OpesciGrid struct
        class OpesciGrid(Structure):
            _fields_ = [(ccode(f.label), POINTER(c_float)) for f in self.fields]

        # Generate the grid argument
        self._arg_grid = OpesciGrid()
        self._arg_grid.values = [POINTER(c_float)() for f in self.fields]

        # Load opesci_run, define it's arguments and run
        print "Executing core computation..."
        opesci_execute = self._library.opesci_execute
        opesci_execute.argtypes = [POINTER(OpesciGrid)]
        opesci_execute(pointer(self._arg_grid))

    def convergence(self):
        """Compute L2 norms for convergence testing"""
        if self._library is None:
            self._load_library()
        if self._arg_grid is None:
            raise RuntimeError("""Convergence could not find grid argument!
You need to you run grid.execute() first!""")

        # Define OpesciConvergence struct
        class OpesciConvergence(Structure):
            _fields_ = [('%s_l2' % ccode(f.label), c_float)
                        for f in self.fields]

        # Generate the convergence argument
        arg_conv = OpesciConvergence()
        arg_conv.values = [c_float(0.) for f in self.fields]

        # Load opesci_convergence, define it's arguments and run
        print "Convergence:"
        opesci_convergence = self._library.opesci_convergence
        opesci_convergence.argtypes = [POINTER(self._arg_grid.__class__),
                                       POINTER(OpesciConvergence)]
        opesci_convergence(pointer(self._arg_grid), pointer(arg_conv))
        for field, _ in arg_conv._fields_:
            print "%s: %.10f" % (field, getattr(arg_conv, field))
