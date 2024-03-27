import sympy as sym
from sympy.utilities.lambdify import lambdify
import numpy as np
import matplotlib.pyplot as plt

err_abs = lambda f, g: np.abs(f - g)

class FunctionApprox:

    def __init__(self, sym_expr, interval, real_alfas, x0, prev_x):
        self.sym_expr = sym_expr
        self.interval = interval
        self.real_alfas = real_alfas
        self.x_sym = sym.Symbol('x')
        self.compute()
        self.def_get_extr()
        self.derivative()
        self.x0 = x0
        self.prev_x = prev_x
        
    def lambda_map(self, expr):
        return lambdify(self.x_sym, expr, np)

    def derivative(self, multeplicity=1):
        self.derivate = lambdify(self.x_sym, sym.diff(self.sym_expr, self.x_sym, multeplicity), np)

    def compute(self):
        self.func = self.lambda_map(self.sym_expr)

    def def_get_extr(self):
        self.a, self.b = self.interval

    def solve_radix(self, method):
        alfa, self.i, self.vec_xk = method

    def plot_error(self, method):
        self.solve_radix(method)
        plt.semilogy(np.arange(self.i), err_abs(np.array(self.vec_xk), self.real_alfas), "o-")
    
