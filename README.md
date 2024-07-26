## Implementation of the Newton-Raphson method for approximating the roots of real-valued functions

The Newton-Raphson method is an iterative method that takes an initial estimate of the input function(s) root(s) and then uses the tangent of the function(s) to those points
to improve on the estimate of the root(s). [insert the formulas formatted nicely for single variable and multivariable functions]

### Description

- Command line interface built with Click
- Running `python newton_raphson.py iterate-newton --help` will display information about Click input formatting through the command line

### Testing

- `git clone https://github.com/rrachelhuangg/numerical_methods.git`
- `pip install sympy`
- Can either test directly using the command line or by running `python -m pytest tests` to run the provided test suite
- If testing with the command line:
  - Click statements are of the form `python newton_raphson.py iterate-newton --functions x**2+4*x --initial-guesses 1 --max-iterations 10 --tolerance 0.01 --symbols x` for the example
  function x<sup>2</sup>+4x with an initial guess of x = 1 and iteration bounds of either 10 iterations or until the error of the root approximation is less than 0.01.
  - Double asterisk ** for exponents and single asterisk * for multiplication
  - The order of click arguments does not matter
  - The tolerance and max-iterations arguments are optional, max-iterations will default to 10 if neither tolerance and max-iterations are provided as input
  - example root approximation command for a multivariable system of equations with multiple initial guesses for the roots that solve the system: `python newton_raphson.py iterate-newton --functions
x**2+y+z --functions x+y+z --functions x+y+z**2 --initial-guesses 1 --initial-guesses 1 --initial-guesses 1 --initial-guesses 0 --initial-guesses 1 --initial-guesses 0 --max-iterations 5 --symbols x --symbols y --symbols z`
    - In this three variable system of equations, two initial guesses are entered: (1, 1, 1) and (0, 1, 0). The program will automatically group the entered values together
    into points with the appropriate number of coordinates regardless of what order the `--initial-guesses` Click arguments are entered in.
    - This example testing statement should return the two root approximations (1, -2, 1) and (0, 0, 0) in the form [('1.000000', '-2.000000', '1.000000'), ('0.000000', '0.000000', '0.000000')] 
