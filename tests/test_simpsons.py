"""Test functions for Simpson's rule."""
import pytest
from methods import simple_simpsons, composite_simpsons

@pytest.mark.parametrize(
    'function, interval_start, interval_end, n_subintervals, expected_approximation',
    [
        ('x**3', 2, 10, 4, 2496.0),
        ('3**x', -2, 2, 4, 8.148148),
    ]
)
def test_simple_simpsons(function, interval_start, interval_end, n_subintervals, expected_approximation):
    """Test Simpson's rule over a single interval."""
    assert simple_simpsons(function, interval_start, interval_end, n_subintervals) == expected_approximation
