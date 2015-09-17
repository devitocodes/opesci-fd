from sympy import Symbol

__all__ = ['DDerivative']


class DDerivative(Symbol):
    """
    wrapper class to represent derivatives
    Sympy already have a "Derivative" class, thus double D
    """
    def __new__(cls, name, *args):
        return Symbol.__new__(cls, name)

    def __init__(self, name, var, order, max_accuracy):
        """
        :param name: name string, should use "D_field_var_order"
        :param var: the dependent variable (symbol)
        :param order: the order of the derivative, i.e. order=1 for 1st derivatives
        :param max_accuracy: the maximum order of accuracy stored
        """
        self.name = name
        self.var = var
        self.order = order
        self.max_accuracy = max_accuracy
        # store fd approxmiate of different accuracy
        self.fd = [None] * (max_accuracy + 1)
