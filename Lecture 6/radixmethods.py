import math
import numpy as np

sign = lambda x: math.copysign(1, x)
reg_mid = lambda a, b, f: a + (b - a) / 2
false_mid = lambda a, b, f: a - f(a) * (b-a) / (f(b) - f(a))
secure_bound = lambda a, b, tolx, i, maxi, fxk, tolf: abs(b - a) > tolx
over_flow_check = lambda a, b, tolx, i, maxi, fxk, tolf:  i < maxi and secure_bound(a, b, tolx, i, maxi, fxk, tolf) and abs(fxk) > tolf


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


def bisection(f, a, b, maxit, tolx, tolf=None):
    return non_linear(f, a, b, tolx, maxit, reg_mid, secure_bound)

def regula_falsi(f, a, b, maxit, tolx, given_tolf):
    return non_linear(f, a, b, tolx, maxit, false_mid, over_flow_check, tolf=given_tolf)
