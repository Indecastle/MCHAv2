import sympy as sm
from sympy import Derivative, diff

#n = None

def maper(expr):
    if (expr.func == sm.Mul or expr.func == sm.Symbol):
        return [expr]
    return list(expr.args)


def maper2(expr):
    if(isinstance(expr, int) or isinstance(expr, float)):
        return float(expr)
    if (expr.func == sm.Symbol):
        return 1.0
    elif (expr == 0):
        return 0.0
    return float(expr.args[0])


def maper3(express, n):
    for e in express:
        e.sort(key=comparator)
    for e in express:
        if (len(e) != n+1):
            for i, v in enumerate(e):
                if comparator(v) != i:
                    e.insert(i, 0)
            if (len(e) != n+1):
                e.extend([0]*(n+1-len(e)))


def comparator(expr):
    if (expr.args == () ):
        text = str(expr)
    else:
        text = str(expr.args[1])
    return int(text[1:])


def TDMA(a,b,c,f):
    a, b, c, f = map(lambda k_list: list(map(float, k_list)), (a, b, c, f))
    a[0] = 0
    n = len(f)
    c[n-1] = 0
    alpha = [0]*n
    beta = [0]*n
    x = [0]*n
    u = [0]*n
    
    x[0] = b[0]
    alpha[0] = -c[0]/x[0]
    beta[0] = f[0]/x[0]
    
    for i in range(1,n-1):
        x[i] = b[i] + a[i]*alpha[i-1]
        alpha[i] = -c[i]/x[i]
        beta[i] = (f[i] - a[i]*beta[i-1])/x[i]
    
    x[n-1] = b[n-1] + a[n-1]*alpha[n-2]
    beta[n-1] = (f[n-1] - a[n-1]*beta[n-2])/x[n-1]
    
    u[n-1] = beta[n-1]
    for i in range(n-2,-1,-1):
        u[i] = alpha[i]*u[i+1] + beta[i]
    
    return u
