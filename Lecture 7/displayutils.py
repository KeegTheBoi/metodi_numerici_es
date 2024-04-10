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


def display_matrix(matrix, notation="A(X)"):
    str_build: str = r"$ "+ notation + r"= \begin{bmatrix}"
    rows, cols = matrix.shape
    for i in range(rows):
        for j in range(cols):
            str_build += str(matrix[i, j])
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