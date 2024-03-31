import sympy as sym
from sympy.utilities.lambdify import lambdify
import numpy as np
import matplotlib.pyplot as plt
import IPython.display as notebook

err_abs = lambda f, g: np.abs(f - g)

class FunctionApprox:

    def __init__(self, sym_expr, interval=None, real_alfas=None, x0=None, prev_x=None, multiplicity=1, symbols=[sym.symbols('x')], notation='f'):
        self.sym_expr = sym_expr(*symbols)
        self.symbols = symbols
        self.interval = interval
        self.real_alfas = real_alfas
        self.compute()
        if interval is not None:
            self.def_get_extr() 
        self.derivative(multiplicity)
        self.x0 = x0
        self.prev_x = prev_x
        self.notation = notation
        
    def lambda_map(self, expr):
        return lambdify(self.symbols, expr, np)

    def derivative(self, multiplicity):
        self.derivate = lambdify(self.symbols, sym.diff(self.sym_expr, *self.symbols, multiplicity), np)

    def partial_derivative(self, x_n: str, multiplicity=1, display=False):
        p_diff = sym.diff(self.sym_expr, sym.symbols(x_n), multiplicity)
        
        if display:
            notebook.display(notebook.Latex(r"$\frac{\partial "+self.notation+f"{self.symbols}"+r"}{\partial  "+x_n+r"}= "+f"{sym.latex(p_diff)}"+r"$"))
        self.derivate = lambdify(sym.symbols(x_n), p_diff, np)
        return sym.latex(p_diff)

    def compute(self):
        self.func = self.lambda_map(self.sym_expr)

    def def_get_extr(self):
        self.a, self.b = self.interval

    def solve_radix(self, method):
        zero, self.i, self.vec_xk = method

    def display_func(self):
        notebook.display(notebook.Latex(f"${self.notation}{self.symbols}={sym.latex(self.sym_expr)}$"))

    def plot_error(self, method, plotter=plt):
        self.solve_radix(method)
        plotter.semilogy(np.arange(self.i), err_abs(np.array(self.vec_xk), self.real_alfas), "o-")

