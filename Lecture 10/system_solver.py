from SolveTriangular import *
import numpy.linalg as lin
import scipy

def LU_solve(A, b):
    PT, L, U = scipy.linalg.lu(A)
    return Usolve(U, Lsolve(L, PT.T @ b))

def cholensky_solve(A, b):
    L = lin.cholesky(A)
    return Usolve(L, Lsolve(L, b))

def QR_solve(A, b):
    Q, R = scipy.linalg.qr(A)
    return Usolve(R, Q.T @ b)
    
