"""Test functions for the Newton-Raphson method"""
#errors to catch neatly: 1) determinant = 0, can't invert
#singlevariable newton
#python newton_raphson.py iterate-newton --functions x**3-4*x**2+1 --initial-guesses 0.5 --max-iterations 3
#multivariable newton
#python newton_raphson.py iterate-newton --functions x+y**2 --functions x-y**2 --symbols x --symbols y --initial-guesses 1 --initial-guesses 1 --max-iterations 3
# CURRENT ESTIMATE:  Matrix([[0], [1/8]])
#python newton_raphson.py iterate-newton --functions x+y+z --functions 2*x+y-z --functions x+y+z**2 --initial-guesses 1 --initial-guesses 1 --initial-guesses 1 --max-iterations 3 --symbols x --symbols y --symbols z
#iter_estimate:  [2, -3, 1]