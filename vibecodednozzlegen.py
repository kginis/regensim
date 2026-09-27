"""Generate a chamber and chart-based Rao-style bell nozzle as x_m,radius_m CSV.

The chart angles below are approximate readings of NASA SP-125, figure 4-14
(expansion area ratio 5..50, 60..100 percent length). For area ratios 2..5,
the first chart interval is linearly extrapolated; those angles are estimates
outside the published chart. The published construction
is a 1.5 Rt upstream throat arc, a 0.382 Rt downstream throat arc, and a
parabolic approximation from N to the exit. It is tangent-continuous (C1).

SMOOTH_TRANSITIONS adds short C2 blends at the throat and at N. This is a
locally modified Rao-style contour, not the unmodified published construction.
Set it to False when an exact arc-plus-parabola construction is required.

Reference: NASA SP-125, Design of Liquid Propellant Rocket Engines, section
4.3, figures 4-13 and 4-14, and equation 4-7.
https://ntrs.nasa.gov/citations/19710019929
"""

import csv
import math
from pathlib import Path


# User inputs
OUTPUT_FILENAME = "nozzle.csv"
THROAT_RADIUS_M = 0.020
CONTRACTION_RATIO = 5.0
EXPANSION_RATIO = 4.0  # Set as low as 2.0 for a labeled extrapolation.
L_STAR_M = 0.6
CHAMBER_RADIUS_M = None

CONVERGING_HALF_ANGLE_DEG = 25.0
CHAMBER_FILLET_RADIUS_M = 0.045
UPSTREAM_THROAT_FILLET_RT = 1.5
DOWNSTREAM_THROAT_FILLET_RT = 0.382
BELL_LENGTH_FRACTION = 1.0

# Approximate digitization of NASA SP-125 figure 4-14. Each row holds
# (expansion ratio, theta_N in degrees, theta_E in degrees). The published chart
# has finite line thickness, so these are design estimates, not exact values.
RAO_CHART = {
    0.60: ((5, 28.4, 20.0), (10, 32.0, 17.3), (14, 33.5, 15.5),
           (20, 35.3, 14.8), (30, 37.1, 13.9), (40, 38.2, 13.2),
           (50, 39.0, 12.8)),
    0.70: ((5, 25.3, 16.5), (10, 28.4, 13.5), (14, 30.0, 12.5),
           (20, 31.6, 11.5), (30, 33.5, 10.7), (40, 34.7, 10.4),
           (50, 35.5, 10.2)),
    0.80: ((5, 23.5, 13.0), (10, 25.8, 11.0), (14, 27.4, 9.8),
           (20, 28.9, 9.4), (30, 30.8, 8.7), (40, 32.0, 8.2),
           (50, 33.0, 8.0)),
    0.90: ((5, 21.6, 10.3), (10, 24.0, 8.8), (14, 25.1, 7.8),
           (20, 27.0, 7.2), (30, 29.0, 6.6), (40, 30.3, 6.3),
           (50, 31.4, 6.1)),
    1.00: ((5, 20.3, 8.0), (10, 22.7, 7.1), (14, 24.0, 6.2),
           (20, 25.6, 5.4), (30, 27.6, 5.1), (40, 28.9, 5.0),
           (50, 30.0, 5.0)),
}

# Optional local curvature blends. Values are half-widths in throat radii.
SMOOTH_TRANSITIONS = True
THROAT_BLEND_HALF_WIDTH_RT = 0.025
BELL_BLEND_HALF_WIDTH_RT = 0.080

CHAMBER_POINTS = 40
CHAMBER_FILLET_POINTS = 20
CONVERGING_POINTS = 40
UPSTREAM_THROAT_FILLET_POINTS = 30
DOWNSTREAM_THROAT_FILLET_POINTS = 30
DIVERGING_POINTS = 70
BLEND_POINTS = 25


def linspace(start, stop, count):
    if count < 2:
        raise ValueError("Each contour section needs at least two points.")
    return [start + (stop - start) * i / (count - 1) for i in range(count)]


def append_unique(destination, new_points):
    for x, radius in new_points:
        if destination and math.isclose(x, destination[-1][0], abs_tol=1e-12):
            if not math.isclose(radius, destination[-1][1], abs_tol=1e-10):
                raise ValueError("Section endpoints do not meet.")
            destination[-1] = (x, radius)
        else:
            destination.append((x, radius))


def volume_of_revolution(points):
    return sum(
        math.pi * 0.5 * (r0 * r0 + r1 * r1) * (x1 - x0)
        for (x0, r0), (x1, r1) in zip(points, points[1:])
    )


def interpolate_row(row, expansion_ratio):
    if not 2.0 <= expansion_ratio <= row[-1][0]:
        raise ValueError("Use 2 <= EXPANSION_RATIO <= 50; values below 5 are extrapolated.")
    if expansion_ratio < row[0][0]:
        left, right = row[0], row[1]
        f = (expansion_ratio - left[0]) / (right[0] - left[0])
        return (left[1] + f * (right[1] - left[1]),
                left[2] + f * (right[2] - left[2]))
    for left, right in zip(row, row[1:]):
        if left[0] <= expansion_ratio <= right[0]:
            f = (expansion_ratio - left[0]) / (right[0] - left[0])
            return (left[1] + f * (right[1] - left[1]),
                    left[2] + f * (right[2] - left[2]))
    return row[-1][1:]


def rao_angles(expansion_ratio, length_fraction):
    fractions = tuple(RAO_CHART)
    if not fractions[0] <= length_fraction <= fractions[-1]:
        raise ValueError("Chart interpolation requires 0.60 <= BELL_LENGTH_FRACTION <= 1.00.")
    for lower, upper in zip(fractions, fractions[1:]):
        if lower <= length_fraction <= upper:
            a = interpolate_row(RAO_CHART[lower], expansion_ratio)
            b = interpolate_row(RAO_CHART[upper], expansion_ratio)
            f = (length_fraction - lower) / (upper - lower)
            return a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1])
    return interpolate_row(RAO_CHART[fractions[-1]], expansion_ratio)


def circle_at_x(x, fillet_radius, throat_radius):
    """Return radius, first derivative, and second derivative of a throat arc."""
    root = math.sqrt(fillet_radius * fillet_radius - x * x)
    return (throat_radius + fillet_radius - root,
            x / root,
            fillet_radius * fillet_radius / root**3)


def quintic_blend(x, x0, x1, start, end):
    """Interpolate r, dr/dx and d2r/dx2 at both ends."""
    h = x1 - x0
    u = (x - x0) / h
    y0, m0, k0 = start
    y1, m1, k1 = end
    a0, a1, a2 = y0, h * m0, h * h * k0 / 2
    b0 = y1 - a0 - a1 - a2
    b1 = h * m1 - a1 - 2 * a2
    b2 = h * h * k1 - 2 * a2
    a3 = 10 * b0 - 4 * b1 + b2 / 2
    a4 = -15 * b0 + 7 * b1 - b2
    a5 = 6 * b0 - 3 * b1 + b2 / 2
    return a0 + u * (a1 + u * (a2 + u * (a3 + u * (a4 + u * a5))))


def converging_geometry(chamber_radius, throat_radius, dense=False):
    theta = math.radians(CONVERGING_HALF_ANGLE_DEG)
    chamber_fillet = CHAMBER_FILLET_RADIUS_M
    throat_fillet = UPSTREAM_THROAT_FILLET_RT * throat_radius
    if not 0 < theta < math.pi / 2:
        raise ValueError("CONVERGING_HALF_ANGLE_DEG must be between 0 and 90.")
    if chamber_fillet <= 0 or throat_fillet <= 0:
        raise ValueError("Fillet radii must be positive.")
    x_throat_tangent = -throat_fillet * math.sin(theta)
    r_throat_tangent = throat_radius + throat_fillet * (1 - math.cos(theta))
    r_chamber_tangent = chamber_radius - chamber_fillet * (1 - math.cos(theta))
    if r_chamber_tangent <= r_throat_tangent:
        raise ValueError("Chamber and throat fillets overlap.")
    x_chamber_tangent = x_throat_tangent - (
        (r_chamber_tangent - r_throat_tangent) / math.tan(theta))
    x_chamber_end = x_chamber_tangent - chamber_fillet * math.sin(theta)
    n_arc = 2001 if dense else CHAMBER_FILLET_POINTS
    n_line = 2001 if dense else CONVERGING_POINTS
    n_throat = 2001 if dense else UPSTREAM_THROAT_FILLET_POINTS
    n_blend = 1001 if dense else BLEND_POINTS

    chamber_arc = [
        (x_chamber_end + chamber_fillet * math.sin(a),
         chamber_radius - chamber_fillet * (1 - math.cos(a)))
        for a in linspace(0, theta, n_arc)
    ]
    straight = [
        (x, r_chamber_tangent - math.tan(theta) * (x - x_chamber_tangent))
        for x in linspace(x_chamber_tangent, x_throat_tangent, n_line)
    ]
    d = THROAT_BLEND_HALF_WIDTH_RT * throat_radius if SMOOTH_TRANSITIONS else 0.0
    if d < 0 or d >= -x_throat_tangent:
        raise ValueError("THROAT_BLEND_HALF_WIDTH_RT is too large.")
    throat_arc = [
        (x, circle_at_x(x, throat_fillet, throat_radius)[0])
        for x in linspace(x_throat_tangent, -d, n_throat)
    ]
    left_blend = []
    if SMOOTH_TRANSITIONS:
        common_curvature = 1 / math.sqrt(
            throat_fillet * DOWNSTREAM_THROAT_FILLET_RT * throat_radius)
        for x in linspace(-d, 0.0, n_blend):
            left_blend.append((x, quintic_blend(
                x, -d, 0.0,
                circle_at_x(-d, throat_fillet, throat_radius),
                (throat_radius, 0.0, common_curvature))))

    points = []
    for section in (chamber_arc, straight, throat_arc, left_blend):
        append_unique(points, section)
    return points, x_chamber_end


def bezier_state(t, n, q, e):
    omt = 1 - t
    x = omt * omt * n[0] + 2 * omt * t * q[0] + t * t * e[0]
    r = omt * omt * n[1] + 2 * omt * t * q[1] + t * t * e[1]
    dx = 2 * (omt * (q[0] - n[0]) + t * (e[0] - q[0]))
    dr = 2 * (omt * (q[1] - n[1]) + t * (e[1] - q[1]))
    ddx = 2 * (e[0] - 2 * q[0] + n[0])
    ddr = 2 * (e[1] - 2 * q[1] + n[1])
    return x, r, dr / dx, (ddr * dx - dr * ddx) / dx**3


def bezier_at_x(x, n, q, e):
    lo, hi = 0.0, 1.0
    for _ in range(55):
        t = (lo + hi) / 2
        if bezier_state(t, n, q, e)[0] < x:
            lo = t
        else:
            hi = t
    return bezier_state((lo + hi) / 2, n, q, e)


def divergent_geometry(throat_radius, exit_radius, theta_n_deg, theta_e_deg):
    downstream_fillet = DOWNSTREAM_THROAT_FILLET_RT * throat_radius
    theta_n, theta_e = map(math.radians, (theta_n_deg, theta_e_deg))
    if downstream_fillet <= 0 or not 0 < theta_e < theta_n < math.pi / 2:
        raise ValueError("Downstream fillet and chart angles are invalid.")
    x_n = downstream_fillet * math.sin(theta_n)
    r_n = throat_radius + downstream_fillet * (1 - math.cos(theta_n))
    if exit_radius <= r_n:
        raise ValueError("Expansion ratio is too small for the divergent fillet.")

    cone_angle = math.radians(15)
    # NASA SP-125, equation 4-7: include the 0.382 Rt arc correction.
    cone_length = (exit_radius - throat_radius
                   + downstream_fillet * (1 / math.cos(cone_angle) - 1)
                   ) / math.tan(cone_angle)
    x_exit = BELL_LENGTH_FRACTION * cone_length
    if x_exit <= x_n:
        raise ValueError("The bell ends inside the downstream throat fillet.")
    slope_n, slope_e = math.tan(theta_n), math.tan(theta_e)
    x_q = (exit_radius - r_n + slope_n * x_n - slope_e * x_exit
           ) / (slope_n - slope_e)
    r_q = r_n + slope_n * (x_q - x_n)
    if not x_n < x_q < x_exit:
        raise ValueError("Bell tangent intersection is outside the nozzle.")
    n, q, e = (x_n, r_n), (x_q, r_q), (x_exit, exit_radius)

    if not SMOOTH_TRANSITIONS:
        arc = [(x, circle_at_x(x, downstream_fillet, throat_radius)[0])
               for x in linspace(0.0, x_n, DOWNSTREAM_THROAT_FILLET_POINTS)]
        bell = [(bezier_state(t, n, q, e)[0], bezier_state(t, n, q, e)[1])
                for t in linspace(0.0, 1.0, DIVERGING_POINTS)]
        points = []
        append_unique(points, arc)
        append_unique(points, bell)
        return points, x_n

    d = THROAT_BLEND_HALF_WIDTH_RT * throat_radius
    w = BELL_BLEND_HALF_WIDTH_RT * throat_radius
    if d <= 0 or w <= 0 or d >= x_n - w or x_n + w >= x_exit:
        raise ValueError("Curvature blends overlap or extend beyond the bell.")
    upstream_fillet = UPSTREAM_THROAT_FILLET_RT * throat_radius
    common_curvature = 1 / math.sqrt(upstream_fillet * downstream_fillet)
    right_blend = [
        (x, quintic_blend(x, 0.0, d,
                          (throat_radius, 0.0, common_curvature),
                          circle_at_x(d, downstream_fillet, throat_radius)))
        for x in linspace(0.0, d, BLEND_POINTS)
    ]
    arc = [
        (x, circle_at_x(x, downstream_fillet, throat_radius)[0])
        for x in linspace(d, x_n - w, DOWNSTREAM_THROAT_FILLET_POINTS)
    ]
    start = circle_at_x(x_n - w, downstream_fillet, throat_radius)
    inflection = (r_n, slope_n, 0.0)
    _, end_r, end_slope, end_second = bezier_at_x(x_n + w, n, q, e)
    bell_blend_left = [
        (x, quintic_blend(x, x_n - w, x_n, start, inflection))
        for x in linspace(x_n - w, x_n, BLEND_POINTS)
    ]
    bell_blend_right = [
        (x, quintic_blend(x, x_n, x_n + w, inflection,
                          (end_r, end_slope, end_second)))
        for x in linspace(x_n, x_n + w, BLEND_POINTS)
    ]
    bell = [
        (x, bezier_at_x(x, n, q, e)[1])
        for x in linspace(x_n + w, x_exit, DIVERGING_POINTS)
    ]
    points = []
    for section in (right_blend, arc, bell_blend_left, bell_blend_right, bell):
        append_unique(points, section)
    return points, x_n


def build_contour():
    rt = THROAT_RADIUS_M
    if rt <= 0 or CONTRACTION_RATIO <= 1 or EXPANSION_RATIO <= 1 or L_STAR_M <= 0:
        raise ValueError("Throat radius and L* must be positive; area ratios must exceed 1.")
    rc = rt * math.sqrt(CONTRACTION_RATIO) if CHAMBER_RADIUS_M is None else CHAMBER_RADIUS_M
    if rc <= rt:
        raise ValueError("CHAMBER_RADIUS_M must exceed THROAT_RADIUS_M.")
    re = rt * math.sqrt(EXPANSION_RATIO)
    theta_n, theta_e = rao_angles(EXPANSION_RATIO, BELL_LENGTH_FRACTION)
    converging, x_chamber_end = converging_geometry(rc, rt)
    dense_converging, _ = converging_geometry(rc, rt, dense=True)
    target_volume = L_STAR_M * math.pi * rt * rt
    converging_volume = volume_of_revolution(dense_converging)
    cylindrical_volume = target_volume - converging_volume
    if cylindrical_volume <= 0:
        raise ValueError(f"L_STAR_M must exceed {converging_volume / (math.pi * rt * rt):.6g} m.")
    chamber_length = cylindrical_volume / (math.pi * rc * rc)
    x_injector = x_chamber_end - chamber_length
    chamber = [(x, rc) for x in linspace(x_injector, x_chamber_end, CHAMBER_POINTS)]
    divergent, x_n = divergent_geometry(rt, re, theta_n, theta_e)
    points = []
    for section in (chamber, converging, divergent):
        append_unique(points, section)
    if any(b[0] <= a[0] for a, b in zip(points, points[1:])):
        raise ValueError("Contour x coordinates are not strictly increasing.")
    if min(r for _, r in points) < rt - 1e-10:
        raise ValueError("Contour dips below the requested throat radius.")
    chamber_points = [point for point in points if point[0] <= 0]
    return {
        "points": points, "chamber_radius": rc, "exit_radius": re,
        "chamber_length": chamber_length, "x_injector": x_injector,
        "x_exit": points[-1][0], "x_n": x_n,
        "theta_n_deg": theta_n, "theta_e_deg": theta_e,
        "target_chamber_volume": target_volume,
        "calculated_chamber_volume": volume_of_revolution(chamber_points),
    }


def main():
    result = build_contour()
    output_path = Path(__file__).resolve().with_name(OUTPUT_FILENAME)
    with output_path.open("w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(("x_m", "radius_m"))
        writer.writerows(result["points"])
    actual_l_star = result["calculated_chamber_volume"] / (math.pi * THROAT_RADIUS_M**2)
    print(f"Nozzle CSV: {output_path}")
    print(f"Construction: {'locally curvature-smoothed Rao-style' if SMOOTH_TRANSITIONS else 'standard arc + parabola'}")
    print(f"Chart angles: theta_N={result['theta_n_deg']:.2f} deg, theta_E={result['theta_e_deg']:.2f} deg")
    print(f"Angle source: {'linear extrapolation below chart range' if EXPANSION_RATIO < 5 else 'approximate NASA SP-125 chart readings'}")
    print(f"Expansion ratio: {EXPANSION_RATIO:.3f}")
    print(f"Bell length fraction: {BELL_LENGTH_FRACTION:.3f}")
    print(f"Points: {len(result['points'])}")
    print(f"Throat: (0, {THROAT_RADIUS_M:.6f}) m")
    print(f"Nominal N x: {result['x_n']:.6f} m")
    print(f"Exit: ({result['x_exit']:.6f}, {result['exit_radius']:.6f}) m")
    print(f"Chamber radius: {result['chamber_radius']:.6f} m")
    print(f"Cylindrical chamber length: {result['chamber_length']:.6f} m")
    print(f"Requested L*: {L_STAR_M:.6f} m; calculated L*: {actual_l_star:.6f} m")


if __name__ == "__main__":
    main()
