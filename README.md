## Implementing the Newton-Raphson method and Simpson's Rule

### Description

- Command line interface built with Click
- `git clone https://github.com/rrachelhuangg/numerical_methods.git`
- `cd numerical_methods`
- `python -m pytest tests` runs the provided test suite

### Implementation of the Newton-Raphson method for approximating the roots of real-valued functions

Running `python newton_raphson.py iterate-newton --help` will display information about Click input formatting through the command line

The Newton-Raphson method is an iterative method that takes an initial approximation of its input function's roots and then uses the tangent of the function to those points to improve on the approximation of the roots.

For single variable functions: $x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$ where $x_{n+1}$ is the next iteration of the approximation for a given root x.

For multivariable functions: $x_{n+1} = x_n - [J(x_n)]^{-1}*F(x_n)$ where F is the vector-valued function, x is the vector of variables, J is the Jacobian matrix of F, and n is the iteration number. 

### Testing the Newton-Raphson method implementation
- Example Click commands for testing through the command line (the order of Click arguments does not matter):
  - `python newton_raphson.py iterate-newton --functions 'x**3-4*x**2+1' --initial-approximations 0.5 --max-iterations 3 --symbols x`
  - `python newton_raphson.py iterate-newton --functions 'x**2+y+z' --functions 'x+y+z' --functions 'x+y+z**2' --initial-approximations 1 --initial-approximations 1 --initial-approximations 1 --initial-approximations 0 --initial-approximations 1 --initial-approximations 0 --max-iterations 5 --symbols x --symbols y --symbols z`
    - In this three variable system of equations, two initial approximations are entered: (1, 1, 1) and (0, 1, 0). The program will automatically group the entered values together into points with the appropriate number of coordinates regardless of what order the `--initial-approximations` Click arguments are entered in.

### Implementation of Simpson's rule for approximating the value of definite integrals

Running `python simpsons.py integrate-simpsons --help` will display information about Click input formatting through the command line

Simpson's rule is used to approximate the value of definite integrals by dividing a curve into smaller curves and approximating the area under each smaller curve,
usually with a parabola. These smaller areas are then added together to approximate the area under the curve as a whole.

Formula: S<sub>n</sub> = $(\frac{∆x}{3})*[f(x_0)+4f(x_1)+2f(x_2)+...+2f(x_{n-2})+4f(x_{n-1})+f(x_n)] ≈ \int_{a}^{b}f(x)dx$, where $∆x = \frac{b-a}{n}$.

The implementation of composite Simpson's rule applies Simpson's rule to each subinterval, which can be useful when integrating functions that have high degrees of variance in local behavior.

Note: The calculation of the error bound for Simpson's rule may take a while if the fifth derivative of the provided function
is difficult to solve.

### Testing the Simpsons Rule implementation
- Example Click commands for testing through the command line:
  - `python simpsons.py integrate-simpsons --function '4**x' --interval-start 0 --interval-end 100 --n-subintervals 10`
  - `python simpsons.py integrate-simpsons --function 'sin(x)' --interval-start 0 --interval-end 100 --n-subintervals 10 --simpsons-type False`
  - `python simpsons.py integrate-simpsons --function 'sin(x)' --interval-start 0 --interval-end 3.14159 --n-subintervals 4 --simpsons-type False --calculate-error True`
