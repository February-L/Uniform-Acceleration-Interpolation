# Uniform-Acceleration Interpolation
Interpolate curves between given points on a plane by a series of uniformly accelerated motions.

Requirements: Python version 3.5+, NumPy (reference version: 1.24)  
Tested with Python IDLE 3.10.8

# Introduction
The 'uniform-acceleration interpolation' is a continuously differentiable, short-sighted, and causality-preserving interpolation algorithm. It imagines that the points on the plane are described by a single parameter, just as the time during a motion. It then interpolates values between two points (x0, y0) and (x1, y1) on a plane by simulating a uniformly accelerated motion through (x0, y0) -> (x1, y1). Note that the acceleration is only uniform between any two points, but usually discontinuous at the points.
Let the speed at (x0, y0) be fixed to 1, the motion curve is then modified by two parameters: the initial velocity's direction, and the maximum speed. The initial velocity's direction is set as the velocity direction at (x0, y0) in the preceding motion, and the maximum speed is specified by the user to indirectly control the bending of the curve.

Determination of the first initial velocity: The initial velocity at the first point in all points can be specified manually, or be determined by analysing the first three points using certain algorithm (see the comment below the definition of variable \_choice_str_tuple in the source code).

About the maximum speed: Given the initial velocity, the algorithm determines the motion by choosing the shortest-time path between two points, which is determined by the maximum speed allowed during the motion (e.g., for infinity max speed, the time cost is 0 and the path is (almost) a straight line). By definition, the maximum speed must be larger than the initial speed, and when so, it equals the speed at the end point (x1, y1); note that the latter can be smaller than the initial speed if we imagine a retarded motion, but it's meaningless for our short-sighted algorithm.

About the short sight and the preserving of causality of the algorithm: The algorithm is called 'short-sighted' because it only cares about the point it currently stands on and the next point, which means later points will not influence the results of former interpolations (causality). This short sight also means that if we're seeking the shortest-time path, the speed at (x1, y1) must exceed the speed at (x0, y0) as long as the two points don't overlap.

Note that the default value of the default start velocity algorithm (see choose_start_v()) might not be good in some cases, so try choosing another start velocity algorithm or calculate a satisfying start velocity yourself.

# Usage
The source code has implemented 3 public functions: uni_accel_interp(), calc_start_v(), choose_start_v().

While the usage is written in the function comments, I'd still like to give two short examples on how to use these functions:  
(Assuming you have imported the functions like 'from uni_accel_interp_feb import *')

1.

    choose_start_v('sketch') # The default is 'parabola', and there's another choice called 'prefer_x'  
    (xt, yt) = uni_accel_interp(x, y, 99, 2.5) # Here between each (x, y), 99 points will be interpolated, using a uniform max_speed value 2.5

(xt, yt) will be the interpolated arrays. Note that the interpolation function regards NaNs as interruptions and copies the velocity before the interruption to the first point after the interruption.

2.

    (v0_x, v0_y) = calc_start_v(x[:3], y[:3], 2., False) # Use the default algorithm 'parabola' to get a start velocity (this one is unnormalised) at max speed 2.0  
    N_array = 9. * ((x[1:] - x[:-1])**2. + (y[1:] - y[:-1])**2.)**.5 # To be passed as an array of the numbers of interpolated points within each (x, y) interval
    N_array[0] = 0 # 0 corresponds to no interpolated points
    V_array = 0. * x + 3. # To be passed as an array of max speed values within each (x, y) interval
    V_array[0] = 2.5
    (xt, yt) = uni_accel_interp(x, y, N_array, V_array, (v0_x, v0_y))

The input start velocity will be normalised, so you can ignore its normalisation when inputting as long as its norm doesn't overflow. The arrayed forms of the arguments interp_num and max_speed can also 1-1 correspond to the points rather than the intervals, and when so, the interpolation function will ignore their first elements.
