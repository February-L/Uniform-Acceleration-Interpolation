from typing import Union
from sys import stderr
from numpy import inf as the_inf, nan as the_nan, ndarray as np_ndarray, \
    ravel as np_ravel, empty as np_empty, arange as np_arange, \
    nditer as np_nditer, \
    sum as np_sum, sqrt as np_sqrt, \
    any as np_any, all as np_all, logical_not as np_logical_not
from math import sqrt, isfinite
#import numpy as np ##########
#import matplotlib.pyplot as plt ##########

'''
  * Python function uni_accel_interp() by February
  * First created on 10 January 2023.
  * Version 20210116c
  * Removing the typing module import and the type hints should make
    this function work for Python versions below 3.5.
'''

_1st_v0_choice = 0
    # Choice of start velocity when the corresponding argument is None \
    # and the number of (x, y) points is not smaller than 3
_choice_str_tuple = ('parabola', 'prefer_x', 'sketch')
''' Choices:
    0 ('parabola') = to imagine a uniformly accelerated motion through (x0, y0)
        -> (x1, y1) -> (x2, y2), set |v1| to max_speed and then try to find a
        solution, and if no solution is find, the edge of the constrained domain
        of definition of the start velocity's direction is chosen
      * Note that for some point configurations, this gives 90-degree start
        direction for (perhaps almost) all max_speed values allowed.
    1 ('prefer_x') = to give a roughly acceptable start direction depending on
        the first three points, ensuring that the x component of the start
        velocity has the same sign as (x1 - x0) (useful for functions y = y(x))
    2 ('sketch') = to give a roughly acceptable start direction depending on
        the first three points; perhaps close to hand drawing
  * Currently all choices' start velocity outcomes are usually susceptible to
    the value of max_speed. The default is 0, 'parabola'.
'''


def uni_accel_interp(
    x: np_ndarray, y: np_ndarray, interp_num: Union[int, np_ndarray] = 9,
    max_speed: Union[float, np_ndarray] = 3.,
    start_velocity: Union[tuple, list, np_ndarray, None] = None
) -> (np_ndarray, np_ndarray):
    # Interpolate f(x, y) = 0 by a series of uniformly accelerated motions.
    "interp_num specifies how many points to be interpolated between " \
    "two (x, y) points.\nmax_speed (unit 1) limits the least-time-cost " \
    "acceleration between two (x, y) points.\n* both may also 1-1 " \
    "correspond to the intervals or the points (the start is excess)\n" \
    "start_velocity with 2 float elements specifies the (to-be-normalised) " \
    "start velocity.\nReturns interpolated arrays (xt, yt) based on (x, y). " \
    "Unravels all arrays.  @February"

    x = np_ravel(x)
    y = np_ravel(y)
    point_count = x.size
    interval_count = point_count - 1

    if (point_count != y.size):
        print("Error in uni_accel_interp(): Sizes of x and y are different.",
              file = stderr)
        return None
    if (point_count < 2):
        print("Error in uni_accel_interp(): Size of x or y is "
              "smaller than 2.", file = stderr)
        return None

    if (type(interp_num) == np_ndarray):
        interp_num = np_ravel(interp_num.astype('int64', copy = False))
        if (interp_num.size == point_count):
            interp_num = interp_num[1:]
        elif (interp_num.size != interval_count):
            print("Error in uni_accel_interp(): Size of arrayed interp_num is "
                  "inconsistent with sizes of x and y.", file = stderr)
            return None
        if (np_any(interp_num < 0)):
            print("Warning in uni_accel_interp(): Certain effective interp_num "
                  "element is smaller than 0; it will be treated as 0.",
                  file = stderr)
            interp_num *= (interp_num > 0)
        xt_size = np_sum(interp_num) + point_count
    else:
        interp_num = int(interp_num)
        if (interp_num < 0):
            print("Warning in uni_accel_interp(): interp_num is smaller "
                  "than 0; it will be treated as 0.", file = stderr)
            interp_num = 0
        xt_size = (point_count - 1) * (interp_num + 1) + 1

    if (type(max_speed) == np_ndarray):
        max_speed = np_ravel(max_speed.astype('float64', copy = False))
        if (max_speed.size == point_count):
            max_speed = max_speed[1:]
        elif (max_speed.size != interval_count):
            print("Error in uni_accel_interp(): Size of arrayed max_speed is "
                  "inconsistent with sizes of x and y.", file = stderr)
            return None
        vm_v0_sq_diff = max_speed**2. - 1.
            # (max speed)^2 - (initial speed)^2, the latter is normalised to 1
        if (not np_all(vm_v0_sq_diff > 0.)): # check NaN by the way
            print("Error in uni_accel_interp(): Certain effective max_speed "
                  "element (squared) isn't larger than 1.", file = stderr)
            return None
        vm0_v0_ratio = abs(max_speed[0]) # divided by unit v0
    else:
        max_speed = float(max_speed)
        vm_v0_sq_diff = max_speed * max_speed - 1.
        if (not (vm_v0_sq_diff > 0.)): # check NaN by the way
            print("Error in uni_accel_interp(): max_speed squared is "
                  "too small: %G; " % (vm_v0_sq_diff + 1.) +
                  "it should be larger than 1.", file = stderr)
            return None
        vm0_v0_ratio = abs(max_speed) # divided by unit v0

    xt = np_empty(xt_size)
    yt = np_empty(xt_size)
    ''' Determination of start velocity --> '''
    if (start_velocity is None):
        if (point_count > 2):
            dx1 = x[1] - x[0]; dy1 = y[1] - y[0]
            dx2 = x[2] - x[0]; dy2 = y[2] - y[0]
            (v0_x, v0_y) = _v0_from_3p(dx1, dy1, dx2, dy2, vm0_v0_ratio)
        else:
            v0_x = x[1] - x[0]; v0_y = y[1] - y[0]
    else:
        v0_x = float(start_velocity[0])
        v0_y = float(start_velocity[1])
    v0_norm = sqrt(v0_x * v0_x + v0_y * v0_y)
    if (v0_norm == the_inf):
        print("Warning in uni_accel_interp(): Norm of the specified start "
              "velocity (%G, %G) is infinity; the start velocity will be "
              "set to zero." % (v0_x, v0_y), file = stderr)
        v0_x = 0.; v0_y = 0.
            # if let them be divided by v0_norm, (infinity/infinity) might occur
    elif (v0_norm > 0.):
        v0_x /= v0_norm; v0_y /= v0_norm
    elif (v0_x != 0. or v0_y != 0.): # also dealing with NaN
        print("Warning in uni_accel_interp(): Norm of the start velocity "
              "(%G, %G) is too small or not a number; the start velocity "
              "will be set to zero." % (v0_x, v0_y), file = stderr)
        v0_x = 0.; v0_y = 0.
    # 'else' is when (v0_x == 0. and v0_y == 0.)

    ''' Interpolation --> '''
    with np_nditer(
        (x[:interval_count], y[:interval_count], x[1:], y[1:],
         interp_num, vm_v0_sq_diff)
    ) as iter_np:

        prev_N = -1
        i_t = 0
        for (x0, y0, x1, y1, N, vm2_sub_v02) in iter_np:

            xt[i_t] = x0; yt[i_t] = y0
            i_t += 1
            dx = x1 - x0; dy = y1 - y0
            dot_x = dx * dx + dy * dy
            if (not isfinite(dot_x)):
                xt[i_t:(i_t + N)] = x0
                yt[i_t:(i_t + N)] = y0
                i_t += N
                continue

            t_base =  (dx * v0_x + dy * v0_y) / vm2_sub_v02
                # term pre-calculated to enchance performance & reduce overflow
            dt_min = 2. * (sqrt(
                t_base * t_base + dot_x / vm2_sub_v02
            ) - t_base)
            if (dt_min > 0. and dt_min != the_inf):
                vd_x = dx / dt_min - v0_x; half_a_x = vd_x / dt_min
                vd_y = dy / dt_min - v0_y; half_a_y = vd_y / dt_min
                    # the acceleration a = 2 * vd / dt
                    # it should minimise the time cost from (x0, y0) to (x1, y1)
            else: # dt_min can't be NaN since isfinite(dot_x) is ensured
                half_a_y = the_inf # simply for the 'if' statement below

            if (N != prev_N):
                frac_arr = np_arange(1, N + 1, dtype = 'float64') / (N + 1)
                    # probably reusable fraction array [1, 2, ..., N] / (N + 1)
                prev_N = N
            if (isfinite(half_a_y) and isfinite(half_a_x)):
                t_arr = frac_arr * dt_min
                xt[i_t:(i_t + N)] = t_arr * (t_arr * half_a_x + v0_x) + x0
                yt[i_t:(i_t + N)] = t_arr * (t_arr * half_a_y + v0_y) + y0
                v0_x += 2. * vd_x
                v0_y += 2. * vd_y
            else:
                t_arr = frac_arr**2.
                    # not the actual time array; just reusing the variable here
                xt[i_t:(i_t + N)] = t_arr * dx + x0
                yt[i_t:(i_t + N)] = t_arr * dy + y0
                v0_x = dx
                v0_y = dy
            i_t += N

            v0_norm = sqrt(v0_x * v0_x + v0_y * v0_y)
            if (v0_norm > 0.): # (infinity / infinity) won't happen here
                v0_x /= v0_norm; v0_y /= v0_norm
            else:
                v0_x = 0.; v0_y = 0.

        xt[i_t] = x1; yt[i_t] = y1

    return (xt, yt)


def calc_start_v(x: np_ndarray, y: np_ndarray,
    max_speed: float, normalize: bool = True
) -> (float, float):
    "Calculate a start velocity using the currently specified algorithm.\n" \
    "x, y, and max_speed are like the arguments with the same name in " \
    "uni_accel_interp().\n" \
    "x and y must have at least 3 elements; only the first three " \
    "elements are used.\n" \
    "Returns a start velocity in tuple form (v0_x, v0_y) if successful.\n" \
    "* normalize = True makes the returned velocity's norm be either " \
    "1 or 0.\nReturns (NaN, NaN) otherwise."
    x = np_ravel(x)
    y = np_ravel(y)
    if (x.size < 3 or y.size < 3):
        print("Error in calc_start_v(): Length of x or y is smaller than 3.",
              file = stderr)
        return (the_nan, the_nan)
    max_speed = abs(float(max_speed))
    if (not (max_speed > 1.)): # check NaN by the way
        print("Error in calc_start_v(): max_speed is too small: %G; "
              % max_speed + "it should be larger than 1.", file = stderr)
        return (the_nan, the_nan)
    dx1 = float(x[1] - x[0]); dy1 = float(y[1] - y[0])
    dx2 = float(x[2] - x[0]); dy2 = float(y[2] - y[0])
    (v0_x, v0_y) = _v0_from_3p(dx1, dy1, dx2, dy2, max_speed) # speed unit |v0|
    if (normalize):
        v0_norm = sqrt(v0_x * v0_x + v0_y * v0_y)
        if (v0_norm == the_inf):
            print("Warning in calc_start_v(): Norm of the calculated start "
                  "velocity (%G, %G) is infinity; the start velocity will be "
                  "returned as zero." % (v0_x, v0_y), file = stderr)
            v0_x = 0.; v0_y = 0.
        elif (v0_norm > 0.):
            v0_x /= v0_norm; v0_y /= v0_norm
        elif (v0_x != 0. or v0_y != 0.): # also dealing with NaN
            print("Warning in calc_start_v(): Norm of the calculated start "
                  "velocity (%G, %G) is too small or not a number; the start "
                  "velocity will be set to zero." % (v0_x, v0_y), file = stderr)
            v0_x = 0.; v0_y = 0.
        # 'else' is when (v0_x == 0. and v0_y == 0.)
    return (v0_x, v0_y)


def _v0_from_3p(dx1, dy1, dx2, dy2, Rv_0) -> (float, float):
    "Internal function to determine the default start velocity.\n" \
    "Note that no argument check is done.\n" \
    "Returns an unnormalised start velocity in tuple form (v0_x, v0_y).\n"
    # dx := x - x0, dy := y - y0
    # Rv_0 := (max speed at the second point) / |v0|; 'R' means 'ratio'
    global _1st_v0_choice
    dot11 = dx1 * dx1 + dy1 * dy1
    dot22 = dx2 * dx2 + dy2 * dy2
    if (isfinite(dot11) and isfinite(dot22)):
        dot12 = dx1 * dx2 + dy1 * dy2
        cross12 = dx1 * dy2 - dy1 * dx2
        if (_1st_v0_choice == 1): # Choice 1 'prefer_x'
            d_factor = .5 * sqrt(2.)
            d01 = min(abs(dx1), d_factor * (abs(dy1) + abs(dx1)))
            d12 = max(abs(dx2 - dx1),
                      d_factor * (abs(dy2 - dy1) + abs(dx2 - dx1)))
            if (d01 > 0.): # this ensures dot11 is non-zero as well
                if (dot22 > 0.):
                    cross12_n = abs(cross12) / (sqrt(dot11) * sqrt(dot22))
                        # double sqrt() to reduce overflow
                else:
                    cross12_n = 0.
                b_sq = 1. + d12 / d01 * (1. + Rv_0 * cross12_n)
                v0_x = b_sq * dx1 - dx2
                v0_y = b_sq * dy1 - dy2
            else:
                v0_x = 0.; v0_y = 0.
            # Parameter b^2 >= 1: v0 // (b^2 * Δx1 - Δx2)
            # Consider a uniformly accelerated motion through (x0, y0)
            # -> (x1, y1) -> (x2, y2), then b = [(t2 - t0) / (t1 - t0)].
        elif (_1st_v0_choice == 2): # Choice 2 'sketch'
            if (dot11 > 0.):
                if (dot22 > 0.):
                    cross12_n = abs(cross12) / (sqrt(dot11) * sqrt(dot22))
                        # double sqrt() to reduce overflow
                else:
                    cross12_n = 0.
                l12_sq = (dx2 - dx1)**2. + (dy2 - dy1)**2.
                b_sq = 1. + sqrt(l12_sq / dot11) * sqrt(1. + Rv_0 * cross12_n)
                b_sq *= b_sq
                v0_x = b_sq * dx1 - dx2
                v0_y = b_sq * dy1 - dy2
            else:
                v0_x = 0.; v0_y = 0.
        else: # Choice 0 'parabola'
            Rv_0_sq = Rv_0 * Rv_0
            if (Rv_0_sq != the_inf and dot11 > 0. and cross12 != 0.):
                (s, u) = _solve_s(dot11, dot12, cross12, Rv_0_sq)
                # Parameter s within (-1, 1): s := v0 × Δx1 / |v0||Δx1|
                # u := (1 - s^2)^0.5
                v0_x = u * dx1 + s * dy1
                v0_y = u * dy1 - s * dx1
            else:
                v0_x = dx1; v0_y = dy1
    else:
        v0_x = 0.; v0_y = 0.
    return (v0_x, v0_y)


def _solve_s(dot11, dot12, cross12, Rv_sq, v0 = 1.) -> (float, float):
    "Internal function to solve the start velocity's s parameter for " \
    "choice 'parabola'.\n" \
    "Note that all arguments are assumed to be physical values.\n" \
    "Returns a tuple of 2 float numbers (s, u):\n" \
    " -- s is within (-1, 1); sgn(s) = -sgn(cross12)\n" \
    " -- u = (1 - s**2)**0.5"
    ''' Let dx := x - x0, dt := t - t0, then:
        dot11 = dx1 * dx1 + dy1 * dy1
        dot12 = dx1 * dx2 + dy1 * dy2
        cross12 = dx1 * dy2 - dy1 * dx2
        Rv_sq := (max speed)^2 / v0^2 > 1 # 'Ratio of v squared'
        s := v0 × dx1_vector / (|v0| * (dot11)^0.5)
        u := v0 ∙ dx1_vector / (|v0| * (dot11)^0.5)
      * Fixing the three points (and v0), a root should definity exists within
        a certain |v1| (which is fixed to max_speed in the calculation to
        minimise (t1 - t0)) interval, but the interval is often exclusive of
        the chosen max_speed (sometimes exlusive of all max_speed > |v0|):
        Firstly it's necessary to note how the root is determined: given s, then
        calculate a dt1 value that matches the max_speed requirement, denoted
        as dt1_1 -> use dt1_1 and the cross product relation (including cross12)
        between dt1 & dt2 to calculate a value of dt2, denoted as dt2_2 -> use
        dt1_1 & dt2_2 and the dot product relation (including dot12) between dt1
        & dt2 to calculate a value of dt1, denoted as dt1_3 -> compare dt1_1 and
        dt1_3 (the 'f' in the code below is (dt1_1 - dt1_3) / dt2_2), and if
        they match (f = 0), s is the root.
        Ignoring cases where dot11 = 0 or cross12 = 0, there are some chosen
        constraints apart from |v1| = max_speed > |v0|:
            dt1 > 0, dt2 > dt1, u >= 0.
            (sgn(s) = sgn(cross12) is a natural result of the equations.)
        Applying, all these constraints, the f = 0 equation is often unsolvable
        within the domain (of definition) of s, and when this happens:
        Since dt1_1 > 0 and dt2_2 > dt1_1 (could be infinity, so the code below
        uses the reciprocal of it ('t2_inv')) is guaranteed by the equations
        and the domain of s, and u >= 0 is also ensured by replacing u by
        (1 - s^2)^0.5 in the equations, the only two constraints left possibly
        unfulfilled is dt1_3 > 0 and dt2_2 > dt1_3. It is therefore natural to
        choose some discriminants to judge dt1_3. Firstly, dt1_3 can be complex,
        and if not so, dt1_3 > 0 is guaranteed by the equations. Secondly,
        dt1_3 / dt2_2 can be larger than 1. These two cases lead to the 'disc' &
        'Rt_iter' variable in the code. The full constraints introduce an actual
        domain of definition of s, and I take the edge of the domain as the
        'best s value we can get' if no root exists: here the edge refers to the
        non-zero end of the actual domain of s, since s = 0 always satisfy the
        constraints (at s = 0 with non-0 cross12, dt1_1 = dt1_3 > 0, dt2_2 is
        infinity).
      * There are some extra statements and variables in the code if one follows
        what is wrote in the preceding comment. They are the legacy from some
        previous discriminants used. I keep them as they are in case of changing
        discriminants in the future, and they're also kind of 'integrated' with
        the whole code, which makes it unwise to change them especially as they
        shouldn't use much performance.
    '''
    # No need to consider dot11 = 0 or cross12 = 0: they should have been \
    # separated from possible input beforehand.
    # 'C' means 'constant':
    dot11_rt = sqrt(dot11)
    C_4xx_vv = 4. * dot12 / (v0 * v0)
    C_4x_v = 4. * dot11_rt / v0
    C_xx_vx = dot12 / (v0 * dot11_rt)
    C_tD_Rv = 2. * dot11_rt / (v0 * (Rv_sq - 1.))
    C_d = 4. * cross12 / (dot11_rt * v0)

    scan_step_power = 8 # scan step length = 2^(-scan_step_power)
    base_steps = 30
        # the step length is 2^(-base_steps) after the first loop
        # 2^(-30) is round 9.3E-10
    fine_steps = 10
        # append some steps (max = fine_steps) if the result s is small
        # 2^(-40) is round 9.1E-13

    sgn_s = 1. if (C_d > 0.) else -1.
    ds = sgn_s # the current step
    s_min = .5**(base_steps + 1) * sgn_s
    s_max = sgn_s - s_min
    avoid_bad = True
        # whether to retreat when encountering bad discriminant points

    # SPANNING DOMAIN (Calculation #3/3, see #1 for the original) -->
    ds_scan = .5**scan_step_power * sgn_s
    s = np_arange(0., ds + ds_scan, ds_scan)
    s[0] = s_min; s[-1] = s_max
    s_sq = s * s
    u = np_sqrt(1. - s_sq)
    t1 = C_tD_Rv * (np_sqrt(Rv_sq - s_sq) - u)
    st1 = s * t1
    t2_inv = 2. * s / (st1 + sgn_s * np_sqrt(st1 * (st1 + C_d)))
    disc = 1. - s_sq + t2_inv * (
        C_4xx_vv * t2_inv - C_4x_v * u
    )
    disc[disc < 0.] = the_nan
    Rt_iter = (np_sqrt(disc) - u) / (2. * (C_xx_vx * t2_inv - u))
    iter_bad_test = np_logical_not(Rt_iter < 1.)
    test0_i = iter_bad_test.argmax()
    good0_i = 0
    if (test0_i == 0 and iter_bad_test[0]):
        test0_i = iter_bad_test.argmin()
        if (test0_i == 0): good0_i = None
        else: good0_i = test0_i
        avoid_bad = False
    f = t1 * t2_inv - Rt_iter

    if (good0_i is not None):
        s_init = s[0] # s[0] == s_min
        f_init = the_nan if (iter_bad_test[0]) else f[0]
        over0_i = (f[1:] * f[0:-1] < 0.).argmax()
        if (f[over0_i] * f[over0_i + 1] < 0.):
            init_i = (good0_i - 1) if (good0_i) else 0
            s_init = s[init_i]
            f_init = f[init_i]
            ds = s[over0_i] - s_init
                # if the scan found a root, \
                # the root is guaranteed to be found later
            avoid_bad = False
        s_good0 = s[good0_i]
        f_good0 = f[good0_i]
            # the known first good discriminant point
    else: # when no good disc. section between s_min & s_max
        return (s[0], u[0]) # s[0] should be s_min
    if (avoid_bad and test0_i != 0):
        ds = s[test0_i]
            # if the scan found a bad disc. section (not near 0), \
            # the section is guaranteed to be found later (if wanted, \
            # i.e. if a root or 'a good disc. section near 0' is not found)
            # Bad disc. section at both ends is almost certainly to be found.

    s = s_init; f = f_init
    if (f == 0.): return (s, sqrt(1. - s * s))
    del iter_bad_test

    # THE FIRST LOOP (Calculation #1/3) -->
    for n in range(base_steps):
        f0 = f # f0 is 'the previous f(s)'
        ds *= .5
        s += ds

        ''' NOTE: Changes here should apply to other 3 places as well. '''
        s_sq = s * s
        u = sqrt(1. - s_sq) # u := √(1 - s^2)
        t1 = C_tD_Rv * (sqrt(Rv_sq - s_sq) - u)
        st1 = s * t1
        t2_inv = 2. * s / (st1 + sgn_s * sqrt(st1 * (st1 + C_d)))
        disc = 1. - s_sq + t2_inv * (
            C_4xx_vv * t2_inv - C_4x_v * u
        ) # discriminant to test if the iterated t1 (t1') is real
        if (disc >= 0.):
            Rt_iter = (sqrt(disc) - u) / (2. * (C_xx_vx * t2_inv - u))
                # Ratio of t1' to t2; it should be definitely larger than 0
                # the constraints also require it to be smaller than 1
            if (Rt_iter < 1.):
                f = t1 * t2_inv - Rt_iter
                if (f * f0 < 0.): ds = -ds
                elif (f == 0.): return (s, u)
            elif (avoid_bad):
                # copy of the 'elif' statement outside below
                s -= ds
        elif (avoid_bad):
            # unless the scan didn't found a bad disc. section near 0, \
            # and it didn't found a root either
            s -= ds

    if (not isfinite(f)):
        print("Warning for uni_accel_interp: In _solve_s(): Abnormal "
              "'f' value (%G) after search. With: s = %G, t1 = %G, "
              "t2_inv = %G, disc = %G, Rt_iter = %G. The solution will be "
              "returned as 0." % (f, s, t1, t2_inv, disc, Rt_iter),
              file = stderr)
        s = 0.; u = 1.
        return (s, u)

    # avoid_bad is False if the scan found a root or \
    # 'a bad disc. section near 0'.
    if (s == s_max):
        # if no root is found or the scan missed the root \
        # and 'NO bad disc. section near the end is found' or \
        # 'the whole domain is a bad disc. section (already returned)'
        if (s_good0 == s_min):
            # if there's a good disc. section at the start, \
            # the first loop can almost certainly find an existing root, \
            # so nothing can be done; choose to return s (= s_max) here
            return (s, u)
        else:
            # if the good disc. section is away from 0, the first loop \
            # might missed a root before s_good0, and if no root, \
            # it's time to refine the s_good0 result
            s = s_good0 + s_min # s_max + s_min == sgn_s
            f = f_good0
            more_steps = int(.094 / max(abs(s), .093 / fine_steps)) \
                + base_steps - scan_step_power
            ds = -(ds_scan + s_min)
            avoid_bad = True
    else:
        # conversely, it's time to refine the result
        more_steps = int(.094 / max(abs(s), .093 / fine_steps))
            # steps to be appended

    # THE APPENDED LOOP (Calculation #2/3, see #1 for the original) -->
    for n in range(more_steps):
        # copy of the previous (i.e. the first) loop
        f0 = f
        ds *= .5
        s += ds

        s_sq = s * s
        u = sqrt(1. - s_sq)
        t1 = C_tD_Rv * (sqrt(Rv_sq - s_sq) - u)
        st1 = s * t1
        t2_inv = 2. * s / (st1 + sgn_s * sqrt(st1 * (st1 + C_d)))
        disc = 1. - s_sq + t2_inv * (
            C_4xx_vv * t2_inv - C_4x_v * u
        )
        if (disc >= 0.):
            Rt_iter = (sqrt(disc) - u) / (2. * (C_xx_vx * t2_inv - u))
            if (Rt_iter < 1.):
                f = t1 * t2_inv - Rt_iter
                if (f * f0 < 0.): ds = -ds
                elif (f == 0.): return (s, u)
            elif (avoid_bad):
                s -= ds
        elif (avoid_bad):
            s -= ds

    if (disc <= 0. and avoid_bad):
        u = sqrt(1. - s * s)
    return (s, u)


def choose_start_v(choice_name: Union[str, int, None] = None) -> int:
    "Sets the algorithm used to calculate the default start velocity.\n" \
    "choice_name (case-insensitive) = 'parabola', 'prefer_x', 'sketch'\n" \
    "-- or as int type = 0, 1, 2. The default built-in value is 0.\n" \
    "Returns the current choice of type int if set successfully or " \
    "received None.\nReturns None otherwise."
    global _1st_v0_choice
    if (choice_name is not None):
        choice_is_valid = False
        global _choice_str_tuple
        if (type(choice_name) == str):
            choice_name = choice_name.lower()
            for i_choice in range(len(_choice_str_tuple)):
                if (choice_name == _choice_str_tuple[i_choice]):
                    choice_name = i_choice
                    choice_is_valid = True
                    break
        else:
            choice_name = int(choice_name)
            if (choice_name >= 0 and choice_name < len(_choice_str_tuple)):
                choice_is_valid = True
        if (choice_is_valid):
            _1st_v0_choice = choice_name
            return _1st_v0_choice
        else:
            print("Error in choose_start_v(): Invalid choice_name value '"
                  + str(choice_name) + "'.", file = stderr)
            return None
    else:
        return _1st_v0_choice
