"""CLI to support the Newton-raphson method for single and multivariable real functions"""
import click
from methods import newton_raphson, newton_raphson_multi

@click.group
def main():
    """Approximate the roots of a function using the Newton-Raphson method"""

@main.command()
@click.option(
    '--functions',
    multiple=True,
    help='Input functions.'
)
@click.option(
    '--initial-approximations',
    multiple=True,
    type=float,
    help='Initial root approximations for the input functions.'
)
@click.option(
    '--tolerance',
    default=None,
    type=float,
    help='Optional: If the absolute difference between the current approximation and the previous approximation' 
    'is less than the tolerance, iteration will be stopped.'
)
@click.option(
    '--max-iterations',
    default=None,
    type=int,
    help='Optional: The number of Newton-Raphson method iterations desired to approximate function roots. Defaults to 10 if neither'
    'a --tolerance argument or --max-iterations argument is provided.'
)
@click.option(
    '--symbols',
    multiple=True,
    help='The symbols used in the system of equations provided.'
)

def iterate_newton(
    functions,
    initial_approximations,
    tolerance,
    max_iterations,
    symbols
):
    """Apply the Newton-Raphson method iteratively until desired accuracy is obtained."""
    if len(functions) == 1:
        print(newton_raphson(functions, initial_approximations, tolerance, max_iterations, symbols))
    elif len(functions) > 1:
        print(newton_raphson_multi(functions, initial_approximations, tolerance, max_iterations, symbols))

if __name__ == '__main__':
    main()
