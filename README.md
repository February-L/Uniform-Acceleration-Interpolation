# Uniform-Acceleration Interpolation
Interpolate 1D curve between given points by a series of uniformly accelerated motions.

Requirement: Python version 3.5+, NumPy (reference version: 1.24)  
Tested with Python IDLE 3.10.8

With 3 public functions implemented: uni_accel_interp(), calc_start_v(), choose_start_v()

# Introduction
The 'uniform-acceleration interpolation' is a continuously differentiable, short-sighted, and causality-preserving interpolation algorithm. It imagines that the points on the plane are described by a single parameter, just as the time during a motion. It then interpolates values between two points (x0, y0) and (x1, y1) on a plane by simulating a uniformly accelerated motion through (x0, y0) -> (x1, y1). Note that the acceleration is only uniform between any two points, but usually discontinuous at the points.
Let the speed at (x0, y0) be fixed to 1, the motion curve is then modified by two parameters: the initial velocity's direction, and the maximum speed. The initial velocity's direction is set as the velocity direction at (x0, y0) in the preceding motion, and the maximum speed is specified by the user to indirectly control the bending of the curve.

Determination of the first initial velocity: The initial velocity at the first point in all points can be specified manually, or be determined by analysing the first three points using certain algorithm (see the comment below the definition of variable \_choice_str_tuple in the source code).

About the maximum speed: Given the initial velocity, the algorithm determines the motion by choosing the shortest-time path between two points, which is determined by the maximum speed allowed during the motion (e.g. for infinity max speed, the time cost is 0 and the path is (almost) a straight line). By definition, the maximum speed must be larger than the initial speed, and when so, it equals the speed at the end point (x1, y1); note that the latter can be smaller than the initial speed if we imagine a retarded motion, but it's meaningless for our short-sighted algorithm.

About the short sight and the perserving of causality of the algorithm: The algorithm is called 'short-sighted' because it only cares about the point it currently stands on and the next point, which means later points will not influence the results of former interpolations (causality). This short sight also means that if we're seeking the shortest-time path, the speed at (x1, y1) must exceed the speed at (x0, y0) as long as the two points don't overlap.

Note that the default value of the default start velocity algorithm (see choose_start_v()) might not be good in some cases, so try choosing another start velocity algorithm or calculate a satisfying start velocity yourself.
   
