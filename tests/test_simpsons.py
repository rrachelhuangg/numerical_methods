"""Test functions for Simpson's rule."""
import pytest
from methods import simple_simpsons, composite_simpsons, simpsons_error

@pytest.mark.parametrize(
    'function, interval_start, interval_end, n_subintervals, calculate_error, toggle, expected_approximation',
    [
        ('x**3', 2, 10, 4, False, False, 'Approximation:  2496.0'),
        ('3**x', -2, 2, 4, False, False, 'Approximation:  8.148148'),
        ('(1)/(x**3+1)', 1, 7, 6, False, False, 'Approximation:  0.371542'),
        ('cos(1+sqrt(x))',0, 4, 8, False, False, 'Approximation:  -2.471601'),
        ('(1)/(1+x**2)', 0, 1, 4, False, False, 'Approximation:  0.785392'),
        ('1/x', 1, 2, 2, False, False, 'Approximation:  0.694444'),
        ('x**2', 0, 4, 4, False, False, 'Approximation:  21.333333')
    ]
)
def test_simple_simpsons(function, interval_start, interval_end, n_subintervals, calculate_error, toggle, expected_approximation):
    """Test Simpson's rule over a single interval."""
    assert simple_simpsons(function, interval_start, interval_end, n_subintervals, calculate_error, toggle) == expected_approximation

@pytest.mark.parametrize(
    'function, interval_start, interval_end, n_subintervals, calculate_error, toggle, expected_approximation',
    [
        ('x**2', 0, 4, 4, False, True, 'Approximation:  21.333332'),
        ('x**3', 2, 10, 4, False, True, 'Approximation:  2496.0'),
        ('sin(x)-cos(x)', 0, 100, 10, False, True, 'Approximation:  0.648101'),
        ('(1)/(1+x**2)', 0, 10, 10, False, True, 'Approximation:  1.459168'),
        ('sin(x)', 0, 100, 10, False, True, 'Approximation:  0.138548'),
        ('4**x', 0, 100, 10, False, True, 'Approximation:  1.1784212324565931e+60')
    ]
)
def test_composite_simpsons(function, interval_start, interval_end, n_subintervals, calculate_error, toggle, expected_approximation):
    """Test Simpson's rule over a single interval."""
    assert composite_simpsons(function, interval_start, interval_end, n_subintervals, calculate_error, toggle) == expected_approximation

@pytest.mark.parametrize(
    'function, interval_start, interval_end, n_subintervals, expected_error_bound',
    [
        ('x**4-4*x**3+6*x**2-4*x+1', 0, 2, 4, 24.0),
        ('e**x', 0, 1, 4, 5.9e-05),
        ('sin(x)', 0, 'pi', 4, 0.006641),
        ('e**x', 0, 1, 6, 1.2e-05)
    ]
)
def test_simpsons_error_calc(function, interval_start, interval_end, n_subintervals, expected_error_bound):
    """Test error bound calculation for Simpson's rule."""
    assert simpsons_error(function, interval_start, interval_end, n_subintervals) == expected_error_bound 
