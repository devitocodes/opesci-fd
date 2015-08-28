from compilation import GNUCompiler

from StringIO import StringIO
from mako.runtime import Context

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

    def compile(self, filename, compiler='g++'):
        # Generate code if this hasn't been done yet
        if self.src_file is None:
            self.generate(filename)

        # Compile cource file with appropriate compiler
        if compiler in ['g++', 'gnu']:
            self._compiler = GNUCompiler()
            self._compiler.compile(self.src_file)
