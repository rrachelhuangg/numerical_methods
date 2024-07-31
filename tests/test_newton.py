"""Test functions for the Newton-Raphson method"""
import pytest
from methods import newton_raphson, newton_raphson_multi

@pytest.mark.parametrize(
    'functions, initial_approximations, tolerance, max_iterations, symbols, expected_approximations',
    [
        (['x**3-4*x**2+1'], [0.5], None, 3, ('x',), 'Root approximations:  [0.537402]'),
        (['x**3-x'], [0.5, -0.5, 0.1], None, 5, ('x',), 'Root approximations:  [-1.0, 1.0, 0.0]'),
        (['x**4-5*x**2+4'], [-1.5, 1.5, 0.5, -0.5], None, 11, ('x',), 'Root approximations:  [-2.000189, 2.000189, 1.0, -1.0]'),
        (['x**4-5*x**2+4'], [-1.5, 1.5, 0.5, -0.5], 0.1, None, ('x',), 'Root approximations:  [-2.011062, 2.011062, 0.999995, -0.999995]'),
        (['x**4-5*x**2+4'], [-1.5, 1.5, 0.5, -0.5], 0.1, 3, ('x',), 'Root approximations:  [-7.286189, 7.286189, 0.999995, -0.999995]'),
        (['sqrt(x)-3'], [1], None, 10, ('x',), 'Root approximations:  [9.0]'),
        (['(1)/(x**2+1)'], [1], None, 10, ('x',), 'Root approximations:  [87.939691]')
    ]
)
def test_single_var_newton(functions: list[str], initial_approximations: list[float], tolerance: float, max_iterations: int, symbols, expected_approximations):
    """Test the Newton_Raphson method on single variable real functions with a single initial approximation"""
    assert newton_raphson(functions, initial_approximations, tolerance, max_iterations, symbols) == expected_approximations

@pytest.mark.parametrize(
    'functions, initial_approximations, tolerance, max_iterations, symbols, expected_approximations',
    [
        (['x+y**2','x-y**2'], [1, 1], None, 3, ('x','y',), "Root approximations:  [('0.000000', '0.125000')]"),
        (['cos(x)-y','sin(y)+x-1'], [0, 1], None, 3, ('x','y',), "Root approximations:  [('0.166039', '0.986247')]"),
        (['x+y+z','2*x+y-z', 'x+y+z**2'], [1, 1, 1], None, 3, ('x','y','z',), "Root approximations:  [('2.000000', '-3.000000', '1.000000')]"),
        (['x**2+y-z-3','x-y**2+z-1', 'x+y+z**2-3'], [1, 1, 1, 2, 2, 2], None, 3, ('x','y','z',), "Root approximations:  [('1.582129', '1.079010', '0.582132'), ('1.581883', '1.080814', '0.583169')]"),
        (['x**2+y-z-3','x-y**2+z-1', 'x+y+z**2-3'], [1, 1, 1, 2, 2, 2], 0.01, 3, ('x','y','z',), "Root approximations:  [('1.585572', '1.080462', '0.577934'), ('1.581883', '1.080814', '0.583169')]"),
        (['x**2+y','y**2+x'], [1, 1, -1, -1], None, 10, ('x','y',), "Root approximations:  [('0.000000', '0.000000'), ('-1.000000', '-1.000000')]"),
        (['(1)/(x**2+y**2+4)', '(1)/(x**2-y**2+4)'], [1, -1], None, 10, ('x', 'y',), "Root approximations:  [('159.530531', '-57.665039')]"),
        (['x**2+y-z-3','x-y**2+z-1', 'x+y+z**2-3'], [1, 1, 1, 2, 2, 2], 0.01, None, ('x','y','z',), "Root approximations:  [('1.582129', '1.079010', '0.582132'), ('1.582125', '1.079006', '0.582126')]"),

    ]
)
def test_multi_var_newton(functions: list[str], initial_approximations: list[float], tolerance: float, max_iterations: int, symbols, expected_approximations):
    """Test the Newton_Raphson method on single variable real functions with multiple initial approximations"""
    assert newton_raphson_multi(functions, initial_approximations, tolerance, max_iterations, symbols) == expected_approximations

@pytest.mark.parametrize(
    'functions, initial_approximations, tolerance, max_iterations, symbols, expected_approximations',
    [
        (['x**2+y**2+z**2','x**2+y**2+z**2','x**2+y**2-z**2'], [0, 0, 0], None, 5, ('x','y','z'), 'Jacobian matrix not invertible. Please try other initial approximations.'),
        (['(1)/(x**2+y**2)', '(1)/(x**2-y**2)'], [1, -1], None, 10, ('x', 'y',), 'Jacobian matrix not invertible. Please try other initial approximations.')
    ]
)
def test_non_inv_Jac_newton(functions: list[str], initial_approximations: list[float], tolerance: float, max_iterations: int, symbols, expected_approximations):
    """Test the Newton_Raphson method on single variable real functions with multiple initial approximations"""
    assert newton_raphson_multi(functions, initial_approximations, tolerance, max_iterations, symbols) == expected_approximations
