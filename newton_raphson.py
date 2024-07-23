import click
from methods import newton_raphson, newton_raphson_multi

@click.group
def main():
    """Approximate the roots of a function using the Newton-Raphson method"""

@main.command()
@click.option(
    '--functions',
    multiple=True,
    help='temp text'
)
@click.option(
    '--initial-guesses',
    multiple=True,
    help='temp text'
)
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
@click.option(
    '--symbols',
    default = None,
    multiple=True,
    help='temp text'
)

def iterate_newton(
    functions,
    initial_guesses,
    symbols,
    tolerance: float | None = None,
    max_iterations: int | None = None
) -> list[float]:
    """Apply the Newton-Raphson method iteratively until desired accuracy is obtained."""
    if len(functions) == 1:
        return newton_raphson(functions, initial_guesses, tolerance, max_iterations)
        #python newton_raphson.py iterate-newton --functions x**3-4*x**2+1 --initial-guesses 0.5 --max-iterations 3
    elif len(functions) > 1:
        #python newton_raphson.py iterate-newton --functions x+y**2 --functions x-y**2 --symbols x --symbols y --initial-guesses 1 --initial-guesses 1 --max-iterations 3
        # CURRENT ESTIMATE:  Matrix([[0], [1/2]])
        # CURRENT ESTIMATE:  Matrix([[0], [1/4]])
        # CURRENT ESTIMATE:  Matrix([[0], [1/8]])
        return newton_raphson_multi(functions, initial_guesses, symbols, tolerance, max_iterations)

if __name__ == "__main__":
    main()