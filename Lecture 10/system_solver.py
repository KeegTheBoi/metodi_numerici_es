from SolveTriangular import *
import numpy.linalg as lin
import scipy

def LU_solve(A, b):
    PT, L, U = scipy.linalg.lu(A)
    P = PT.copy().T
    return Usolve(U, Lsolve(L, P @ b))

def cholensky_solve(A, b):
    L = lin.cholesky(A)
    return Usolve(L.T.copy(), Lsolve(L, b))

def QR_solve(A, b):
    Q, R = scipy.linalg.qr(A)
    return Usolve(R, Q.T @ b)
    