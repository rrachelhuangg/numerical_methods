"""Newton_Raphson method implementations for single and multivariable real functions"""
from sympy import *
import numpy as np

def newton_raphson(functions, initial_approximations, tolerance: float, max_iterations: int, symbols) -> list[float]:
    """Implementation of the Newton-Raphson method"""
    root_approximations = []
    for approximation in initial_approximations:
        estimate = float(approximation)
        iter_estimate, x, function = estimate, str(np.asarray(symbols)[0]), functions[0]
        if not max_iterations:
            if tolerance:
                while True:
                    iter_estimate = estimate - ((sympify(function).subs(x,estimate))/(sympify(diff(function),x).subs(x,estimate)))
                    if abs(iter_estimate-estimate) < tolerance:
                        break
                    estimate = iter_estimate
            elif not tolerance:
                max_iterations = 10
                for i in range(max_iterations):
                    iter_estimate = estimate - ((sympify(function).subs(x,estimate))/(sympify(diff(function),x).subs(x,estimate)))
                    estimate = iter_estimate
        if max_iterations:
            for i in range(max_iterations):
                iter_estimate = estimate - ((sympify(function).subs(x,estimate))/(sympify(diff(function),x).subs(x,estimate)))
                if tolerance is not None and abs(iter_estimate-estimate) < tolerance:
                    break
                estimate = iter_estimate
        root_approximations+=[float(f"{iter_estimate:.6f}")]
    return root_approximations

def newton_raphson_multi(functions, initial_approximations, tolerance: float, max_iterations: int, symbols):
    """Implementation of the Newton-Raphson method for systems of 
    equations with multiple variables."""
    root_approximations = []
    for i in range(0, len(initial_approximations), len(symbols)):
        estimate, symbols = np.asarray(initial_approximations[i:i+len(symbols)]), np.asarray(symbols)
        iter_estimate = estimate
        F = [function for function in functions]
        F = Matrix([F])
        J = F.jacobian(Matrix(symbols))
        if not max_iterations:
            if tolerance:
                stop=False
                while stop==False:
                    point = {sym:estimate[i] for (i, sym) in enumerate(symbols)}
                    iter_estimate, estimate, stop = Matrix(list(estimate))-((J.subs(point).inv())*(F.subs(point)).T), Matrix(list(estimate)), True
                    checks = [abs(estimate[i]-iter_estimate[i]) for i in range(len(estimate))]
                    for c in checks:
                        if c > tolerance:
                            stop=False
                    estimate = iter_estimate
            elif not tolerance:
                max_iterations = 10
                for i in range(max_iterations):
                    point = {sym:estimate[i] for (i, sym) in enumerate(symbols)}
                    iter_estimate = Matrix(list(estimate))-((J.subs(point).inv())*(F.subs(point)).T)
                    estimate = iter_estimate
        if max_iterations:
            for i in range(max_iterations):
                point = {sym:estimate[i] for (i, sym) in enumerate(symbols)}
                iter_estimate = Matrix(list(estimate))-((J.subs(point).inv())*(F.subs(point)).T)
                if tolerance is not None:
                    estimate, stop = Matrix(list(estimate)), True
                    checks = [abs(estimate[i]-iter_estimate[i]) for i in range(len(estimate))]
                    for c in checks:
                        if c > tolerance:
                            stop=False
                    if stop==True:
                        break
                estimate = iter_estimate
        root_approximations+=[tuple([f"{val:.6f}" for val in estimate])]
    return root_approximations