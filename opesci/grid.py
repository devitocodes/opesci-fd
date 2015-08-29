from compilation import GNUCompiler
from codeprinter import ccode

from StringIO import StringIO
from mako.runtime import Context
from ctypes import cdll, Structure, POINTER, c_float, pointer

class Grid:
    """Base class for grid objects that provides the code generation,
    compilation and execution infrastructure"""

    # mako.lookup.TemplateLookup object
    lookup = None

    # Name of the base template
    template_base = None

    # List of keys used the template that corresponds to an equivalent
    # property field in the derived grid class
    template_keys = []

    # Placeholders for generated code and associated files
    src_code = None
    src_file = None
    src_lib = None

    def generate(self, filename):
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

    def compile(self, filename, compiler='g++', shared=True):
        # Generate code if this hasn't been done yet
        if self.src_file is None:
            self.generate(filename)

        # Compile source file with appropriate compiler
        if compiler in ['g++', 'gnu']:
            self._compiler = GNUCompiler()
            out = self._compiler.compile(self.src_file, shared=shared)
        if shared:
            self.src_lib = out

    def execute(self, filename, compiler='g++'):
        # Compile code if this hasn't been done yet
        if self.src_lib is None:
            self.compile(filename, compiler=compiler, shared=True)

        # Load compiled binary
        try:
            library = cdll.LoadLibrary(self.src_lib)
        except OSError as e:
            print "Library load error: ", e
            raise Exception("Failed to load %s" % self.src_lib)

        # Define our custom OpesciGrid struct
        class OpesciGrid(Structure):
            _fields_ = [(ccode(f.label), POINTER(c_float)) for f in self.fields]

        # Execute opesci_run function with grid struct
        grid_arg = OpesciGrid()
        grid_arg.values = [POINTER(c_float)() for f in self.fields]

        # Load opesci_run, define it's arguments and run
        print "Executing core computation..."
        opesci_execute = library.opesci_execute
        opesci_execute.argtypes = [POINTER(OpesciGrid)]
        opesci_execute(pointer(grid_arg))

        # Load opesci_convergence, define it's arguments and run
        print "Computing convergence:"
        opesci_convergence = library.opesci_convergence
        opesci_convergence.argtypes = [POINTER(OpesciGrid)]
        opesci_convergence(pointer(grid_arg))

        # Close compiled kernel library
        cdll.LoadLibrary('libdl.so').dlclose(library._handle)
