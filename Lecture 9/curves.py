import numpy as np
import IPython.display as nb
import matplotlib.pyplot as plt


def create_mesh(x, y, f, mat):
    A, b = mat
    X, Y = np.meshgrid(x, y)
    rows, cols = X.shape
    Z = np.array([[f((X[i, j], Y[i, j]), A, b) for j in range(cols)] for i in range(rows)])
    return X, Y, Z

def plot_contour(X, Y, Z, f, mat, sol):
    A, b = mat
    ct = plt.contour(X, Y, Z, levels=f(sol, A, b).flatten())
    plt.plot(*sol, 'ro')
    plt.legend(['x at each curve'])

def chached_countor(sol, mesh, Z, f, mat):
    X, Y = mesh
    plot_contour(X, Y, Z, f, mat, sol)

def plot_3D_projection(X, Y, Z, f, mat):
    ax = plt.figure().add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, color='blue', alpha = 0.7,
                linewidth = 0.3, edgecolor = 'black')
    ax.grid(color='black')
    return ax