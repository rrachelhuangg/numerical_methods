from sympy import *
import numpy as np

def newton_raphson(functions, initial_guesses, tolerance: float, max_iterations: int) -> list[float]:
    """Implementation of the Newton-Raphson method"""
    estimate, iter_estimate, x, function = float(initial_guesses[0]), float(initial_guesses[0]), "x", functions[0]
    if not max_iterations:
        if tolerance:
            while abs(iter_estimate-estimate) >= tolerance:
                iter_estimate = estimate - ((sympify(function).subs(x,estimate))/(sympify(diff(function),x).subs(x,estimate)))
                estimate = iter_estimate
        elif not tolerance:
            max_iterations = 10
            for i in range(max_iterations):
                iter_estimate = estimate - ((sympify(function).subs(x,estimate))/(sympify(diff(function),x).subs(x,estimate)))
                estimate = iter_estimate
    if max_iterations:
        for i in range(max_iterations):
            iter_estimate = estimate - ((sympify(function).subs(x,estimate))/(sympify(diff(function),x).subs(x,estimate)))
            estimate = iter_estimate
            if tolerance is not None and abs(iter_estimate-estimate) < tolerance:
                break
    print("iter_estimate: ", [iter_estimate])
    return [iter_estimate]

def newton_raphson_multi(functions, initial_guesses, symbols, tolerance, max_iterations):
    """Implementation of the Newton-Raphson method for systems of 
    equations with multiple variables."""
    #use numpy instead?
    estimate, iter_estimate, symbols = np.asarray(initial_guesses), np.asarray(initial_guesses), np.asarray(symbols)
    functions_matrix = []
    for function in functions:
        functions_matrix += [function]
    functions_matrix = Matrix([functions_matrix])
    symbols_matrix = Matrix(symbols)
    if max_iterations:
        for i in range(max_iterations):
            jacobian_matrix = functions_matrix.jacobian(symbols_matrix)
            point = {}
            for i, sym in enumerate(symbols):
                point[sym] = estimate[i]
            jacobian_evaluated = jacobian_matrix.subs(point)
            inverse = jacobian_evaluated.inv()
            functions_matrix_evaluated = functions_matrix.subs(point)
            mult = inverse*functions_matrix_evaluated.T
            est = Matrix(list(estimate))
            iter_estimate = est-mult
            estimate = iter_estimate
    print("iter_estimate: ", list(estimate))
    return [list(estimate)]







            




