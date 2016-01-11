from cgen import *



class IfDef(Module):
    def __init__(self, condition, iflines, elselines):
        
        ifdef_line = Line('#ifdef %s'%condition)
        else_line = Line('#else')
        endif_line = Line('#endif')
        lines = [ifdef_line]+iflines+[else_line]+elselines+[endif_line]
        
        super(IfDef, self).__init__(lines)

class InlineInitializer(Initializer):
    def generate(self):
        result = super(InlineInitializer, self).generate()
        for v in result:
            if v.endswith(';'):
                yield v[:-1]
            else:
                yield v

def replace_in_code(code, str_from, str_to):
    if isinstance(code, Block):
        replace_in_code(code.contents, str_from, str_to)
        return code
    if isinstance(code,Loop):
        replace_in_code(code.body, str_from, str_to)
        return code
    for code_element in code:
        if isinstance(code_element, Statement) or isinstance(code_element,Line):
            code_element.text.replace(str_from, str_to)
        if isinstance(code_element, Block):
            replace_in_code(code_element.contents, str_from, str_to)
        if isinstance(code_element, Loop):
            replace_in_code(code_element.body, str_from, str_to)
    return code
        