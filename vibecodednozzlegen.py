"""Generate an axisymmetric rocket chamber/nozzle contour as x_m,radius_m CSV.

The throat is placed at x = 0. Chamber coordinates are negative and nozzle
coordinates are positive. Only the Python standard library is required.
"""

import csv
import math
from pathlib import Path


# =============================================================================
# USER-EDITABLE INPUTS
# =============================================================================
OUTPUT_FILENAME = "nozzle.csv"

# Primary sizing
THROAT_RADIUS_M = 0.018
CONTRACTION_RATIO = 6.0
EXPANSION_RATIO = 4.5
L_STAR_M = 0.65
CHAMBER_RADIUS_M = None

# Chamber and converging section
CONVERGING_HALF_ANGLE_DEG = 35.0
CHAMBER_FILLET_RADIUS_M = 0.025
UPSTREAM_THROAT_FILLET_RT = 1.5

# 100% length Rao-style bell approximation
DOWNSTREAM_THROAT_FILLET_RT = 0.382
BELL_LENGTH_FRACTION = 1.00
BELL_INITIAL_ANGLE_DEG = 30.0
BELL_EXIT_ANGLE_DEG = 10.0

# Output resolution
CHAMBER_POINTS = 40
CHAMBER_FILLET_POINTS = 20
CONVERGING_POINTS = 40
UPSTREAM_THROAT_FILLET_POINTS = 30
DOWNSTREAM_THROAT_FILLET_POINTS = 30
DIVERGING_POINTS = 70

# =============================================================================
# GEOMETRY HELPERS
# =============================================================================

def linspace(start, stop, count):
    """Return count evenly spaced values, including both endpoints."""
    if count < 2:
        raise ValueError("Each contour section needs at least two points.")
    step = (stop - start) / (count - 1)
    return [start + step * i for i in range(count)]


def append_unique(destination, new_points):
    """Append points while removing duplicated section endpoints."""
    for x, radius in new_points:
        if destination and math.isclose(
            x, destination[-1][0], rel_tol=0.0, abs_tol=1.0e-12
        ):
            destination[-1] = (x, radius)
        else:
            destination.append((x, radius))


def volume_of_revolution(points):
    """Calculate pi*integral(radius^2 dx) with trapezoidal integration."""
    volume = 0.0
    for (x_1, r_1), (x_2, r_2) in zip(points, points[1:]):
        volume += math.pi * 0.5 * (r_1**2 + r_2**2) * (x_2 - x_1)
    return volume


def converging_geometry(chamber_radius, throat_radius):
    """Build the tangent chamber fillet, cone, and upstream throat fillet."""
    theta = math.radians(CONVERGING_HALF_ANGLE_DEG)
    chamber_fillet = CHAMBER_FILLET_RADIUS_M
    throat_fillet = UPSTREAM_THROAT_FILLET_RT * throat_radius

    if not 0.0 < theta < math.pi / 2.0:
        raise ValueError("CONVERGING_HALF_ANGLE_DEG must be between 0 and 90.")
    if chamber_fillet <= 0.0 or throat_fillet <= 0.0:
        raise ValueError("Fillet radii must be positive.")

    # The throat is fixed at (0, throat_radius). The upstream throat fillet is
    # tangent to the throat and to the straight converging line.
    x_throat_tangent = -throat_fillet * math.sin(theta)
    r_throat_tangent = (
        throat_radius + throat_fillet * (1.0 - math.cos(theta))
    )

    # Radius at the downstream end of the chamber/convergent fillet.
    r_chamber_tangent = (
        chamber_radius - chamber_fillet * (1.0 - math.cos(theta))
    )

    if r_chamber_tangent <= r_throat_tangent:
        raise ValueError(
            "The chamber and throat fillets overlap. Reduce a fillet radius, "
            "increase contraction ratio, or change the converging angle."
        )

    # Locate the chamber fillet so that the intervening line has slope
    # -tan(theta) and is tangent to both circular arcs.
    x_chamber_tangent = x_throat_tangent - (
        (r_chamber_tangent - r_throat_tangent) / math.tan(theta)
    )
    x_chamber_end = x_chamber_tangent - chamber_fillet * math.sin(theta)

    # Chamber-to-convergent circular fillet.
    chamber_arc = []
    for angle in linspace(0.0, theta, CHAMBER_FILLET_POINTS):
        x = x_chamber_end + chamber_fillet * math.sin(angle)
        radius = chamber_radius - chamber_fillet * (1.0 - math.cos(angle))
        chamber_arc.append((x, radius))

    # Straight converging section.
    straight = []
    for x in linspace(x_chamber_tangent, x_throat_tangent, CONVERGING_POINTS):
        radius = r_chamber_tangent - math.tan(theta) * (
            x - x_chamber_tangent
        )
        straight.append((x, radius))

    # Upstream throat circular fillet.
    throat_arc = []
    for angle in linspace(theta, 0.0, UPSTREAM_THROAT_FILLET_POINTS):
        x = -throat_fillet * math.sin(angle)
        radius = throat_radius + throat_fillet * (1.0 - math.cos(angle))
        throat_arc.append((x, radius))

    points = []
    append_unique(points, chamber_arc)
    append_unique(points, straight)
    append_unique(points, throat_arc)

    return points, x_chamber_end


def dense_converging_geometry(chamber_radius, throat_radius):
    """High-resolution convergent used only for accurate volume integration."""
    global CHAMBER_FILLET_POINTS
    global CONVERGING_POINTS
    global UPSTREAM_THROAT_FILLET_POINTS

    original_counts = (
        CHAMBER_FILLET_POINTS,
        CONVERGING_POINTS,
        UPSTREAM_THROAT_FILLET_POINTS,
    )

    try:
        CHAMBER_FILLET_POINTS = 2001
        CONVERGING_POINTS = 2001
        UPSTREAM_THROAT_FILLET_POINTS = 2001
        return converging_geometry(chamber_radius, throat_radius)
    finally:
        (
            CHAMBER_FILLET_POINTS,
            CONVERGING_POINTS,
            UPSTREAM_THROAT_FILLET_POINTS,
        ) = original_counts


def divergent_geometry(throat_radius, exit_radius):
    """Build the downstream throat fillet and Rao-style bell divergent."""
    downstream_fillet = DOWNSTREAM_THROAT_FILLET_RT * throat_radius
    initial_angle = math.radians(BELL_INITIAL_ANGLE_DEG)

    if downstream_fillet <= 0.0:
        raise ValueError("The downstream throat fillet radius must be positive.")
    if not 0.0 < initial_angle < math.pi / 2.0:
        raise ValueError("The initial divergent angle must be between 0 and 90.")

    # Circular throat expansion fillet, tangent to a horizontal line at the
    # throat and to the divergent at initial_angle.
    throat_arc = []
    for angle in linspace(0.0, initial_angle, DOWNSTREAM_THROAT_FILLET_POINTS):
        x = downstream_fillet * math.sin(angle)
        radius = throat_radius + downstream_fillet * (1.0 - math.cos(angle))
        throat_arc.append((x, radius))

    x_n, r_n = throat_arc[-1]

    if exit_radius <= r_n:
        raise ValueError("Expansion ratio is too small for the divergent fillet.")

    exit_angle = math.radians(BELL_EXIT_ANGLE_DEG)
    if not 0.0 < exit_angle < initial_angle:
        raise ValueError(
            "BELL_EXIT_ANGLE_DEG must be positive and smaller than "
            "BELL_INITIAL_ANGLE_DEG."
        )
    if not 0.0 < BELL_LENGTH_FRACTION <= 1.0:
        raise ValueError("BELL_LENGTH_FRACTION must be in the range (0, 1].")

    reference_conical_length = (
        (exit_radius - throat_radius) / math.tan(math.radians(15.0))
    )
    x_exit = BELL_LENGTH_FRACTION * reference_conical_length

    if x_exit <= x_n:
        raise ValueError("The selected bell length ends inside the throat fillet.")

    # The control point is the intersection of the initial and exit tangent
    # lines. A quadratic Bezier through it gives a smooth Rao-style bell.
    initial_slope = math.tan(initial_angle)
    exit_slope = math.tan(exit_angle)

    x_control = (
        exit_radius
        - r_n
        + initial_slope * x_n
        - exit_slope * x_exit
    ) / (initial_slope - exit_slope)
    r_control = r_n + initial_slope * (x_control - x_n)

    if not x_n < x_control < x_exit:
        raise ValueError(
            "Bell tangent intersection lies outside the nozzle. Adjust "
            "the bell length or divergent angles."
        )

    divergent = []
    for t in linspace(0.0, 1.0, DIVERGING_POINTS):
        one_minus_t = 1.0 - t
        x = (
            one_minus_t**2 * x_n
            + 2.0 * one_minus_t * t * x_control
            + t**2 * x_exit
        )
        radius = (
            one_minus_t**2 * r_n
            + 2.0 * one_minus_t * t * r_control
            + t**2 * exit_radius
        )
        divergent.append((x, radius))

    points = []
    append_unique(points, throat_arc)
    append_unique(points, divergent)
    return points


def build_contour():
    throat_radius = THROAT_RADIUS_M
    if throat_radius <= 0.0:
        raise ValueError("THROAT_RADIUS_M must be positive.")
    if CONTRACTION_RATIO <= 1.0:
        raise ValueError("CONTRACTION_RATIO must be greater than 1.")
    if EXPANSION_RATIO <= 1.0:
        raise ValueError("EXPANSION_RATIO must be greater than 1.")
    if L_STAR_M <= 0.0:
        raise ValueError("L_STAR_M must be positive.")

    if CHAMBER_RADIUS_M is None:
        chamber_radius = throat_radius * math.sqrt(CONTRACTION_RATIO)
    else:
        chamber_radius = CHAMBER_RADIUS_M

    exit_radius = throat_radius * math.sqrt(EXPANSION_RATIO)
    throat_area = math.pi * throat_radius**2
    target_chamber_volume = L_STAR_M * throat_area

    converging, x_chamber_end = converging_geometry(
        chamber_radius, throat_radius
    )
    dense_converging, _ = dense_converging_geometry(
        chamber_radius, throat_radius
    )
    converging_volume = volume_of_revolution(dense_converging)

    cylindrical_volume = target_chamber_volume - converging_volume
    if cylindrical_volume <= 0.0:
        minimum_l_star = converging_volume / throat_area
        raise ValueError(
            f"L_STAR_M is too small for the selected convergent. It must be "
            f"greater than {minimum_l_star:.6g} m."
        )

    chamber_area = math.pi * chamber_radius**2
    chamber_length = cylindrical_volume / chamber_area
    x_injector = x_chamber_end - chamber_length

    chamber = [
        (x, chamber_radius)
        for x in linspace(x_injector, x_chamber_end, CHAMBER_POINTS)
    ]

    divergent = divergent_geometry(throat_radius, exit_radius)

    contour = []
    append_unique(contour, chamber)
    append_unique(contour, converging)
    append_unique(contour, divergent)

    return {
        "points": contour,
        "chamber_radius": chamber_radius,
        "exit_radius": exit_radius,
        "chamber_length": chamber_length,
        "x_injector": x_injector,
        "x_exit": contour[-1][0],
        "target_chamber_volume": target_chamber_volume,
        "calculated_chamber_volume": volume_of_revolution(
            [point for point in contour if point[0] <= 0.0]
        ),
    }


def main():
    result = build_contour()
    output_path = Path(__file__).resolve().with_name(OUTPUT_FILENAME)

    with output_path.open("w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["x_m", "radius_m"])
        writer.writerows(result["points"])

    actual_l_star = (
        result["calculated_chamber_volume"]
        / (math.pi * THROAT_RADIUS_M**2)
    )

    print(f"Nozzle CSV created at: {output_path}")
    print(f"Points: {len(result['points'])}")
    print(f"Chamber radius: {result['chamber_radius']:.6f} m")
    print(f"Throat radius: {THROAT_RADIUS_M:.6f} m")
    print(f"Exit radius: {result['exit_radius']:.6f} m")
    print(f"Cylindrical chamber length: {result['chamber_length']:.6f} m")
    print(f"Injector x: {result['x_injector']:.6f} m")
    print(f"Exit x: {result['x_exit']:.6f} m")
    print(f"Requested L*: {L_STAR_M:.6f} m")
    print(f"Calculated L*: {actual_l_star:.6f} m")


if __name__ == "__main__":
    main()