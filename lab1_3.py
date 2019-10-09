import math
import numpy as np
import sympy as sm
from sympy import Derivative
from itertools import groupby
import matplotlib.pyplot as plt
from tabulate import tabulate
from pprint import pprint

from IPython.display import display, Latex, clear_output
from sympy import latex

sm.init_printing(use_latex='mathjax')
x, y, c = sm.symbols('x, y, c', real=True)
C1, C2 = sm.symbols("C1, C2", real=True)
c = sm.symbols("c")
u = sm.Function('u')
ux = u(x)
ux1 = Derivative(u(x), x)
ux2 = Derivative(u(x), (x, 2))

h, uup1, uum1, uu0 = sm.symbols('h, u_i+1, u_i-1, u_i')
uux1 = (uup1 - uum1)/(2*h)
uux2 = (uup1 - 2*uu0 + uum1)/(h**2)


zad1 = {
    'p': sm.sqrt(1 + x**2)*0.4,
    'q': 4*(1 + x**2),
    'f': 20*sm.exp(-x),
    'a': 0, 'UA': 0,
    'b': 2.5, 'UB': 0,
    'E': 0.05
}

#apro = mainEq.subs({ux: uu0, ux2: uux2, ux1: uux1 })
apro = sm.Eq(uux2 + zad1['p']*uux1 + zad1['q']*uu0,  zad1['f'])


n = 10
n -= 1
count=10

hh = (zad1['b'] - zad1['a']) / n
us = sm.symbols(' '.join([f"u{i}" for i in range(n+1)]))
list_x = [zad1['a'] + hh*i for i in range(n+1)]
display(hh)
display(us)
display(list_x)


tyts = []
tyts.append(sm.Eq(us[0] + (-us[2]+4*us[1]-3*us[0])/(2*hh), zad1['UA']))
tyts.extend([apro.subs({x: list_x[i], h: hh, uum1: us[i-1], uu0: us[i], uup1: us[i+1]})  for i in range(1,n)])
tyts.append(sm.Eq(us[n] - 0.5*(3*us[n]-4*us[n-1]+us[n-2])/(2*hh), zad1['UB']))
display(*tyts)


def comparator(expr):
    text = str(expr.args[1])
    return int(text[1:])
express = list(map(lambda t: list(t.lhs.args), tyts))

for e in express:
    e.sort(key=comparator)

expr_values = list(map(lambda t: list(map(lambda t2: float(t2.args[0]), t)), express))
print(tabulate(expr_values, floatfmt=".3f"))



ptable = np.eye(count)
ptable[0,0:3] = expr_values[0]

for i, e in enumerate(expr_values[:-2], 1):
    ptable[i,i-1:i+2] = expr_values[i]
ptable[count-1,0:3] = expr_values[count-1]
print(ptable)
print()
print()
print()
print()
print(latex(ux2))
