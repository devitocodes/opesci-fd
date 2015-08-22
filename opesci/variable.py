from sympy import Symbol

__all__ = ['Variable']


class Variable(Symbol):
    """
    wrapper class of Symbol to store extra information
    """
    def __new__(cls, name, *args):
        return Symbol.__new__(cls, name)

    def __init__(self, name, value=0, type='int', constant=False):
        """
        :param type: type string of the variable to be used in generated code
        :param value: associated value of the variable
        :param constant: if the variable is a constant
        """
        self.type = type
        self.constant = constant
        self.value = value
