import click
from methods import newton_raphson

@click.group
def main():
    """Approximate the roots of a function using the Newton-Raphson method"""

@main.command()
@click.argument('function', type=str)
@click.argument('initial-guess', type=float)
@click.option(
    '--tolerance',
    default=None,
    type=float,
    help='If the absolute difference between the current approximation and the previous approximation' 
    'is less than the tolerance, iteration will be stopped.'
)
@click.option(
    '--max-iterations',
    default=None,
    type=int,
    help='The number of Newton-Raphson method iterations desired to approximate function roots.'
)

def iterate_newton(
    function: str,
    initial_guess: float,
    tolerance: float | None = None,
    max_iterations: int | None = None
):
    """Apply the Newton-Raphson method iteratively until desired accuracy is obtained."""
    estimate, iter_estimate, x = initial_guess, initial_guess, "x"
    if not max_iterations:
        if tolerance:
            while abs(iter_estimate-estimate) >= tolerance:
                iter_estimate = estimate - newton_raphson(estimate, function, x)
                estimate = iter_estimate
        elif not tolerance:
            max_iterations = 10
            for i in range(max_iterations):
                iter_estimate = estimate - newton_raphson(estimate, function, x)
                estimate = iter_estimate
    if max_iterations:
        for i in range(max_iterations):
            iter_estimate = estimate - newton_raphson(estimate, function, x)
            estimate = iter_estimate
            if tolerance is not None and abs(iter_estimate-estimate) < tolerance:
                break
    return iter_estimate

if __name__ == "__main__":
    main()
