from sympy.printing.ccode import CCodePrinter
from mpmath.libmp import to_str, prec_to_dps
from mako.runtime import Context
from StringIO import StringIO


class CodePrinter(CCodePrinter):
    """
    CodePrinter extends sympy.printing.ccode.CCodePrinter
    """
    def __init__(self, settings={}):
        CCodePrinter.__init__(self, settings)

    def _print_Indexed(self, expr):
        """
        override method in CCodePrinter
        Print field as C style multidimensional array
        e.g. U[t,x,y,z] -> U[t][x][y][z]
        """
        output = self._print(expr.base.label) \
            + ''.join(['[' + self._print(x) + ']' for x in expr.indices])
        return output

    def _print_Rational(self, expr):
        """
        override method in CCodePrinter
        print fractional number as float/float
        (default was long double/long double)
        """
        p, q = int(expr.p), int(expr.q)
        return '%d.0F/%d.0F' % (p, q)  # float precision by default

    def _print_Mod(self, expr):
        """
        override method in CCodePrinter
        print mod using % operator in C++
        """
        args = map(ccode, expr.args)
        args = ['('+x+')' for x in args]
        result = '%'.join(args)
        return result

    def _print_Float(self, expr):
        """
        override method in StrPrinter
        always printing floating point numbers in scientific notation
        """
        prec = expr._prec
        if prec < 5:
            dps = 0
        else:
            dps = prec_to_dps(expr._prec)
        if self._settings["full_prec"] is True:
            strip = False
        elif self._settings["full_prec"] is False:
            strip = True
        elif self._settings["full_prec"] == "auto":
            strip = self._print_level > 1
        rv = to_str(expr._mpf_, dps, strip_zeros=strip, max_fixed=-2, min_fixed=2)
        if rv.startswith('-.0'):
            rv = '-0.' + rv[3:]
        elif rv.startswith('.0'):
            rv = '0.' + rv[2:]
        return rv


def ccode(expr, **settings):
    """
    generate C++ code from expression expr
    calling CodePrinter class
    """
    return CodePrinter(settings).doprint(expr, None)


def ccode_eq(eq, **settings):
    """
    genereate C++ assignment from equation eq
    assigning RHS to LHS
    """
    return CodePrinter(settings).doprint(eq.lhs, None) \
        + ' = ' + CodePrinter(settings).doprint(eq.rhs, None)


def render(tmpl, dict1):
    """
    render Mako template using dictionary dict1
    returns the result as string
    """
    buf = StringIO()
    ctx = Context(buf, **dict1)
    tmpl.render_context(ctx)
    return buf.getvalue()
