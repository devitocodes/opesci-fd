from cgen import *


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
            code_element.text = code_element.text.replace(str_from, str_to)
        if isinstance(code_element, Block):
            replace_in_code(code_element.contents, str_from, str_to)
        if isinstance(code_element, Loop):
            replace_in_code(code_element.body, str_from, str_to)
    return code
