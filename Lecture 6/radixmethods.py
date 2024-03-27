import math
import numpy as np
from more_itertools import iterate
from itertools import takewhile
from funcutils import FunctionApprox

sign = lambda x: math.copysign(1, x)
reg_mid = lambda a, b, f: a + (b - a) / 2
false_mid = lambda a, b, f: a - f(a) * (b-a) / (f(b) - f(a))
secure_bound = lambda a, b, tolx, i, maxi, fxk, tolf: abs(b - a) > tolx
over_flow_check = lambda a, b, tolx, i, maxi, fxk, tolf:  \
    i < maxi and secure_bound(a, b, tolx, i, maxi, fxk, tolf) and abs(fxk) > tolf


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
    
    stream = takewhile(lambda xk: stopping_criteria(*xk), enumerate(iterate(lambda xk: xk - delta(xk, f, dx), x0)))
    arr_xk = [v for _, v in stream]
    if arr_xk == []:
        print("nothing has been found")
        
    return None if arr_xk == [] else arr_xk[-1], len(arr_xk), arr_xk


def bisection(f, a, b, maxit, tolx, tolf=None):
    return non_linear(f, a, b, tolx, maxit, reg_mid, secure_bound)

def regula_falsi(f, a, b, maxit, tolx, given_tolf):
    return non_linear(f, a, b, tolx, maxit, false_mid, over_flow_check, tolf=given_tolf)

def newton(f, x0, tolx, tolf, nmax, fa: FunctionApprox):
    dx = fa.derivate
    return linear(f, x0, tolx, tolf, nmax, lambda x, f, dx: f(x) / dx(x), dx)

def chord(f, x0, tolx, tolf, nmax, fa: FunctionApprox):
    dx = fa.derivate
    m = (fa.func(fa.b) - fa.func(fa.a)) / (fa.b - fa.a)
    return linear(f, x0, tolx, tolf, nmax, lambda x, f, dx: f(x) / m, dx)

def secant(f, x0, tolx, tolf, nmax, fa: FunctionApprox):
    dx = fa.derivate
    m = lambda x: (fa.func(x) - fa.func(fa.prev_x)) / (x - fa.prev_x)
    return linear(f, x0, tolx, tolf, nmax, lambda x, f, dx: f(x) / m(x), dx)