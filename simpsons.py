"""CLI to support Simpson's rule for approximating the value of definite integrals"""
import click
from methods import simple_simpsons

@click.group
def main():
    """Approximate the value of definite integrals using Simpson's rule."""

@main.command()
@click.option(
    '--function',
    type = str,
    help='Input function to integrate.'
)
@click.option(
    '--interval-start',
    type=float,
    help='Start of the desired interval that we want to be approximated.'
)
@click.option(
    '--interval-end',
    type=float,
    help='End of the desired interval that we want to be approximated.'
)
@click.option(
    '--n-subintervals',
    type=int,
    help='The number of subintervals that the integration interval is divided into for approximation.'
    'Must be an even integer.'
)
@click.option(
    '--simpsons-type',
    default = True,
    type=bool,
    help='Optional: Will apply the simple variation of Simpson\'s rule where a single interval is used if not provided or True is provided.'
    ' Will apply the composite variation of Simpson\'s rule where multiple intervals are used if False is provided.'
)

def integrate_simpsons(
    function: str,
    interval_start: float,
    interval_end: float,
    n_subintervals: int,
    simpsons_type: bool
) -> float:
    """Apply the appropriate version of Simpson's rule for approximating the integral of the input function"""
    if simpsons_type:
        print(simple_simpsons(function, interval_start, interval_end, n_subintervals))
    else:
        print(composite_simpsons(function, interval_start, interval_end, n_subintervals))

if __name__ == '__main__':
    main()
