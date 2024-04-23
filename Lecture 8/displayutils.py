import sympy as sym
from sympy.utilities.lambdify import lambdify
import numpy as np
import matplotlib.pyplot as plt
import IPython.display as nb


def display_vector(vector, notation='v', round_n=5):
    str_build: str = r"$ \mathbf{"+ notation + r"} = \begin{bmatrix}"
    rows = len(vector)
    for i in range(rows):
        str_build += str(round(vector[i], round_n))
        str_build += "" if i == (rows - 1)  else r" \\ "
    
    str_build += r"\end{bmatrix}$"
    nb.display(nb.Latex(str_build))

def plot_alfa(x, f, title="", radix=None, log=False, plotter=plt):
    plotter.plot(x, f(x), x, x * 0, 'r')
    plotter.scatter(radix if radix is not None else 0, 0)
    plotter.grid(True)

def display_function(func_expr, varsym=sym.symbols("x")):
    nb.display(nb.Latex(r"$f"+f"{varsym}"+r"="+f"{sym.latex(func_expr)}"+r"$"))


def display_matrix(matrix, notation="A(X)", round_n=5):
    str_build: str = r"$ "+ notation + r"= \begin{bmatrix}"
    rows, cols = matrix.shape
    for i in range(rows):
        for j in range(cols):
            str_build += str(round(matrix[i, j], round_n))
            str_build += "" if j == (cols - 1)  else r" & "
        str_build += r" \\ "
    
    str_build += r"\end{bmatrix}$"
    nb.display(nb.Latex(str_build))

def plot_contour(x, y, f, i, label):
    col=('grey', 'blue')
    X, Y = np.meshgrid(x, y)
    ct = plt.contour(X, Y, f(X, Y), levels=0, colors=col[i])
    plt.clabel(ct, inline=True, colors=col[i], fmt=label)
    plt.grid(axis='x')

def plot_3D_projection(x, y, f):
    X, Y = np.meshgrid(x, y)
    ax = plt.figure().add_subplot(111, projection='3d')
    Z = f(X, Y)
    ax.plot_surface(X, Y, Z, color='blue', alpha = 0.7,
                linewidth = 0.3, edgecolor = 'black')
    ax.grid(color='black')
    return ax