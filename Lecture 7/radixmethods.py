import math
import numpy as np
from itertools import takewhile
from funcutils import FunctionApprox

sign = lambda x: math.copysign(1, x)
reg_mid = lambda a, b, f: a + (b - a) / 2
false_mid = lambda a, b, f: a - f(a) * (b-a) / (f(b) - f(a))
secure_bound = lambda a, b, tolx, i, maxi, fxk, tolf: abs(b - a) > tolx
over_flow_check = lambda a, b, tolx, i, maxi, fxk, tolf:  \
    i < maxi and secure_bound(a, b, tolx, i, maxi, fxk, tolf) and abs(fxk) > tolf


def iterate(func, start):
    while True:
        yield start
        start = func(*start)


def non_linear(f, a, b, tolx, maxit, hypo_mid, arrest_cond, tolf=None):
    
    fa=f(a)
    fb=f(b)
    if sign(fa)*sign(fb)>=0:
        print("Non è possibile applicare il metodo di bisezione \n")
        return None, None,None

    i = 0
    v_xk = []

    fxk = 10 #any positive number greater than tolerance
    
    while arrest_cond(a, b, tolx, i, maxit, fxk, tolf):
        xk = hypo_mid(a, b, f)

        v_xk.append(xk)
        i += 1
        fxk=f(xk)
        if fxk==0:
            return xk, i, v_xk
        if sign(fa)*sign(fxk)>0:  #continua su [xk,b]
            a = xk
            fa=fxk
        elif sign(fxk)*sign(fb)>0:   #continua su [a,xk]
            b = xk
            fb=fxk
            
    return xk, i, v_xk



def linear(f, x0, tolx, tolf, nmax, delta, dx):
    
    def stopping_criteria(k, x):
        return k < nmax and abs(f(x)) >= tolf and abs(delta(x, f, dx)) >= tolx * abs(x) and abs(dx(x)) > tolf
        
    iter = iterate(lambda xk: xk - delta(xk, f, dx), x0)
    arr_xk = list(dict(takewhile(lambda param: (stopping_criteria(*param)), enumerate(iter))).values()) + [next(iter)]

    return arr_xk[-1], len(arr_xk), arr_xk

#newton_raphson

#general
def newton_raphson(f, J, x0, tolx, tolf, nmax, s):       
    
    def stop_criteria(i, X):
        X = X[0]
        return i < nmax and np.linalg.norm(s(f, J, X)) / np.linalg.norm(X) > tolx and np.linalg.norm(f(X)) >= tolf and np.linalg.det(J(X)) != 0

    iter = iterate(lambda xk, err: (xk + s(f, J, xk), np.linalg.norm(s(f, J, xk)) / np.linalg.norm(xk)), (x0, 1)) 
    arr = list(dict(takewhile(lambda param: stop_criteria(*param), \
                              enumerate(iter))).values()) + [next(iter)] 
    return arr[-1][0], len(arr), [e for _, e in arr]
    
#variants
def newton_raphson_teoretical(f, J, x0, tolx, tolf, nmax):
    return newton_raphson(f, J, x0, tolx, tolf, nmax, lambda f, J, X: np.linalg.inv(J(X)) @ (-f(X)))

def newton_raphson_optimized(f, J, x0, tolx, tolf, nmax):
    return newton_raphson(f, J, x0, tolx, tolf, nmax, lambda f, J, X: -np.linalg.solve(J(X), f(X)))

def newton_raphson_chord(f, J, x0, tolx, tolf, nmax):
    return newton_raphson(f, J, x0, tolx, tolf, nmax, lambda f, J, X: -np.linalg.solve(J(x0), f(X)))
    


def bisection(f, a, b, maxit, tolx, tolf=None):
    return non_linear(f, a, b, tolx, maxit, reg_mid, secure_bound)

def regula_falsi(f, a, b, maxit, tolx, given_tolf):
    return non_linear(f, a, b, tolx, maxit, false_mid, over_flow_check, tolf=given_tolf)

def newton(f, x0, tolx, tolf, nmax, fa: FunctionApprox, multiplicity=1):
    dx = fa.derivate
    return linear(f, x0, tolx, tolf, nmax, lambda x, f, dx: multiplicity * f(x) / dx(x), dx)

def newton_mod_2(f, x0, tolx, tolf, nmax, fa: FunctionApprox):
    return newton(f, x0, tolx, tolf, nmax, fa, multiplicity=2)

def chord(f, x0, tolx, tolf, nmax, fa: FunctionApprox):
    dx = fa.derivate
    if fa.b - fa.a == 0:
        raise ValueError
    m = (fa.func(fa.b) - fa.func(fa.a)) / (fa.b - fa.a)
    return linear(f, x0, tolx, tolf, nmax, lambda x, f, dx: f(x) / m, dx)

def secant(f, x0, tolx, tolf, nmax, fa: FunctionApprox):
    dx = fa.derivate
    m = lambda x: (fa.func(x) - fa.func(fa.prev_x)) / (x - fa.prev_x)
    return linear(f, x0, tolx, tolf, nmax, lambda x, f, dx: f(x) / m(x), dx)