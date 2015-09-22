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
        self.fd_1side = []

    def set_fd_1side(self, fds):
        """
        set the non-symmetric FD approximation of this Derivative object
        :param list fds: list of non-symmetric FD approximation, where fds[x] indicates FD approximation using left=x, right=order-x
        """
        if len(self.fd_1side) == len(fds):
            for i in range(len(self.fd_1side)):
                if not fds[i] is None:
                    self.fd_1side[i] = fds[i]
        else:
            self.fd_1side = fds
