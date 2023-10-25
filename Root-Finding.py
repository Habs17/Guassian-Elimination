import math
from math import e
from numpy import log as ln
def secant_method(fn, fn2, n): #this is in the form of x2 = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0))
    if n == 5:
        return
    func1 = ln(fn-.95)+.85*math.cos(fn-.95)
    func2 = ln(fn2-.95)+.85*math.cos(fn2-.95)
    ans = fn - (func1)*((fn-fn2)/((func1)-(func2)))
    sol = func1
    print(ans, end = ' ')
    print(sol)
    return secant_method(ans,fn,n+1)

def newtons_method(fn,n):
    if n == 5:
        return
    func = ln(fn-.95)+.85*math.cos(fn-.95)
    deriv = (1/(fn-(19/20)))-((17*math.sin(fn-(19/20)))/20)
    ans = fn - (func)/(deriv)
    sol = func
    print(ans,end = ' ')
    print(sol)
    return newtons_method(ans,n+1)

def bisection_method(x0,x1,n):
    fn = ln(x0-0.95)+.85*math.cos(x0-0.95)
    fn1 = ln(x1-0.95)+.85*math.cos(x1-0.95)
    if abs(fn-fn1) < .00000001:
        print(n)
        return
    c0 = (x0+x1)/2
    fnc0 = ln(c0-0.95)+.85*math.cos(c0-0.95)
    if fnc0 > 0:
        print(x0, ' ', c0)
        bisection_method(x0,c0,n+1)
    else:
        print(c0, ' ', x1)
        bisection_method(x0,c0,n+1)

secant_method(2, 1.25, 0)
