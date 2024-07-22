from sympy import *
def newton_raphson(estimate: float, function: str, symbol: str) -> None:
    """Implementation of the Newton-Raphson method"""
    x = symbols(symbol)
    return ((sympify(function).subs(x,estimate))/(sympify(diff(function),x).subs(x,estimate)))
