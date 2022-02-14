#!/usr/bin/env python3
#
#  ode_solve.py
#
#  Usage:

import numpy as np
import matplotlib.pyplot as plt


def main(method,output):
    print('hi')

## Functions
# dx/dt = f(x,t)
def f(x,t):
    return x

# Euler step
# f - the function
# tn - the time to estimate at
# h - the timestep
# x0 - the initial condition
def euler_step(f, tn, h, x0):
    xn1 = x0 + h*f(x0, tn)
    return xn1


# solve to
def solve_to(f,x0,t1,t2,deltat_max):
    timesteps = int((t2-t1)/deltat_max)
    x = x0
    for step in range(timesteps):
        x = euler_step(f, step, deltat_max, x)
    return x


# solve the ode
def solve_ode(f,t,x0,deltat_max):

    x_series = [x0]
    x = x0
    for i in range(len(t)-1):
        
        x = solve_to(f,x,t[i],t[i+1],deltat_max)
        x_series.append(x)

    return x_series

# if __name__ == "__main__":
#     import sys
#     main(*sys.argv[1:])

# print(solve_to(f, 1, 0, 1, 0.0001))
# print(solve_to(f, 2.718, 1, 2, 0.0001))
# print(solve_to(f, 1, 2, 3, 0.0001))
print(solve_ode(f, [0,.1,.2,.3], 1, 0.0001))