from cgen import *


class IfDef(Module):
    """
    Class to represent IfDef-Else-EndIf construct for the C preprocessor.
    While Cgen has classes for #define, #include, #pragma etc., it has nothing for IfDef.
    :param condition: the condition in IfDef
    :param iflines: the block of code inside the if [an array of type Generable]
    :param elselines: the block of code inside the else [an array of type Generable]
    """
    def __init__(self, condition, iflines, elselines):
        ifdef_line = Line('#ifdef %s' % condition)
        else_line = Line('#else')
        endif_line = Line('#endif')
        lines = [ifdef_line]+iflines+[else_line]+elselines+[endif_line]
        super(IfDef, self).__init__(lines)


class InlineInitializer(Initializer):
    """
    Class to represent Initializers (int i=0) without the semicolon in the end(e.g. in a for statement)
    Usage: same as cgen.Initializer
    Result: same as cgen.Initializer except for the lack of a semi-colon at the end
    """
    def generate(self):
        result = super(InlineInitializer, self).generate()
        for v in result:
            if v.endswith(';'):
                yield v[:-1]
            else:
                yield v


def replace_in_code(code, str_from, str_to):
    """
    Helper function to find and replace strings inside code blocks
    :param code: The code block in which to carry out the find/replace. This can be array of Generables or a Block or a Loop
    :param str_from: The string to find (and to be replaced)
    :param str_to: The string to replace all instances of str_from with
    """
    if isinstance(code, Block):
        replace_in_code(code.contents, str_from, str_to)
        return code
    if isinstance(code, Loop):
        replace_in_code(code.body, str_from, str_to)
        return code
    for code_element in code:
        if isinstance(code_element, Statement) or isinstance(code_element, Line):
            code_element.text.replace(str_from, str_to)
        if isinstance(code_element, Block):
            replace_in_code(code_element.contents, str_from, str_to)
        if isinstance(code_element, Loop):
            replace_in_code(code_element.body, str_from, str_to)
    return code
