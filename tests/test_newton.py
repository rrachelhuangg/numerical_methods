"""Test functions for the Newton-Raphson method"""
#check these answers by hand first!
#errors to catch neatly: 1) determinant = 0, can't invert
#singlevariable newton
#python newton_raphson.py iterate-newton --functions x**3-4*x**2+1 --initial-guesses 0.5 --max-iterations 3
#single variable newton multiple roots
#python newton_raphson.py iterate-newton --functions x**3-x --initial-guesses 0.5 --initial-guesses -0.5 --initial-guesses 0.1 --max-iterations 5
#[-1.00000000000000, 1.00000000000000, 0]
#python newton_raphson.py iterate-newton --functions x**4-5*x**2+4 --initial-guesses -1.5 --initial-guesses 1.5 --initial-guesses 0.5 --initial-guesses -0.5 --max-iterations 11
#[-2.000189, 2.000189, 1.0, -1.0]
#multivariable newton
#python newton_raphson.py iterate-newton --functions x+y**2 --functions x-y**2 --symbols x --symbols y --initial-guesses 1 --initial-guesses 1 --max-iterations 3
# CURRENT ESTIMATE:  Matrix([[0], [1/8]])
#python newton_raphson.py iterate-newton --functions cos(x)-y --functions sin(y)+x-1 --initial-guesses 0 --initial-guesses 1 --symbols x --symbols y --max-iterations 3
#[('0.166039', '0.986247')]
#python newton_raphson.py iterate-newton --functions x+y+z --functions 2*x+y-z --functions x+y+z**2 --initial-guesses 1 --initial-guesses 1 --initial-guesses 1 --max-iterations 3 --symbols x --symbols y --symbols z
#iter_estimate:  [2, -3, 1]
#multivariable newton multiple roots
#python newton_raphson.py iterate-newton --functions cos(x)-y --functions sin(y)+x-1 --initial-guesses 0 --initial-guesses 1 --symbols x --symbols y --max-iterations 3 --initial-guesses 1 --initial-guesses 1
#(the above took a really long time, but did finish eventually)
#[('0.166039', '0.986247'), ('0.162625', '0.992249')]
#python newton_raphson.py iterate-newton --functions x**2+y-z-3 --functions x-y**2+z-1 --functions x+y+z**2-3 --initial-guesses 1 --initial-guesses 1 --initial-guesses 1 --initial-guesses 2 --initial-guesses 2 --initial-guesses 2 --symbols x --symbols y --symbols z --max-iterations 3
#[('1.582129', '1.079010', '0.582132'), ('1.581883', '1.080814', '0.583169')]