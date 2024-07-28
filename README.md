## Implementation of the Newton-Raphson method for approximating the roots of real-valued functions

The Newton-Raphson method is an iterative method that takes an initial approximation of its input function's roots and then uses the tangent of the function to those points to improve on the approximation of the roots.

For single variable functions: x<sub>n+1</sub> = x<sub>n</sub> - f(x<sub>n</sub>)/f'(x<sub>n</sub>) where x<sub>n+1</sub> is the next iteration of the approximation for a given root x.

For multivariable functions: x<sub>n+1</sub> = x<sub>n</sub> - [J(x<sub>n</sub>)]<sup>-1</sup>*F(x<sub>n</sub>) where F is the vector-valued function, x is the vector of variables, J is the Jacobian matrix of F, and n is the iteration number. 

### Description

- Command line interface built with Click
- Running `python newton_raphson.py iterate-newton --help` will display information about Click input formatting through the command line

### Testing

- `git clone https://github.com/rrachelhuangg/numerical_methods.git`
- `cd numerical_methods`
- `python -m pytest tests` runs the provided test suite
- Example Click commands for testing through the command line (the order of Click arguments does not matter):
  - `python newton_raphson.py iterate-newton --functions 'x**3-4*x**2+1' --initial-approximations 0.5 --max-iterations 3 --symbols x`
  - `python3 newton_raphson.py iterate-newton --functions 'x**2+y+z' --functions 'x+y+z' --functions 'x+y+z**2' --initial-approximations 1 --initial-approximations 1 --initial-approximations 1 --initial-approximations 0 --initial-approximations 1 --initial-approximations 0 --max-iterations 5 --symbols x --symbols y --symbols z`
    - In this three variable system of equations, two initial approximations are entered: (1, 1, 1) and (0, 1, 0). The program will automatically group the entered values together into points with the appropriate number of coordinates regardless of what order the `--initial-approximations` Click arguments are entered in.
