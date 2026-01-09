#!/usr/bin/env python3
"""
Analytical solution for heading rotation angle.

EQUATION:
=========
Given:
  - e_x: initial body x-axis (after elevation/azimuth rotation, before heading)
  - k: rotation axis (normalized tether direction)
  - h: target heading
  - wind_norm: wind direction (normalized)

The heading components vary with rotation angle θ as:
  heading_y(θ) = A1*sin(θ) + B1*cos(θ) + C1
  heading_z(θ) = A2*sin(θ) + B2*cos(θ) + C2

Where coefficients are extracted by sampling at θ = 0, π/2, π:
  C1 = (hy(0) + hy(π)) / 2
  B1 = hy(0) - C1
  A1 = hy(π/2) - C1
  (same for A2, B2, C2 using hz)

The heading equation atan2(hy, hz) = h is equivalent to:
  hy*cos(h) - hz*sin(h) = 0

Which gives: A*sin(θ) + B*cos(θ) + C = 0
  where: A = A1*cos(h) - A2*sin(h)
         B = B1*cos(h) - B2*sin(h)
         C = C1*cos(h) - C2*sin(h)

SOLUTION:
  θ = atan2(A, B) ± acos(-C / sqrt(A² + B²))

Pick the solution that gives the correct heading (two candidates, one is correct).
"""

import numpy as np
from numpy.linalg import norm

def rodrigues(v, k, theta):
    k = k / norm(k)
    return v * np.cos(theta) + np.cross(k, v) * np.sin(theta) + k * np.dot(k, v) * (1 - np.cos(theta))

def rotate_around_y(v, angle):
    c, s = np.cos(angle), np.sin(angle)
    R = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    return R @ v

def rotate_around_z(v, angle):
    c, s = np.cos(angle), np.sin(angle)
    R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    return R @ v

def calc_heading(e_x, wind_norm):
    minus_ex = -e_x
    proj_on_wind = np.dot(minus_ex, wind_norm) * wind_norm
    e_x_perp = minus_ex - proj_on_wind
    wind_cross_z = np.array([wind_norm[1], -wind_norm[0], 0])
    heading_y = np.dot(e_x_perp, wind_cross_z)
    heading_z = e_x_perp[2]
    return np.arctan2(heading_y, heading_z)

def wrap_to_pi(angle):
    return (angle + np.pi) % (2 * np.pi) - np.pi

def get_heading_components(e_x, k, theta, wind_norm):
    e_x_rot = rodrigues(e_x, k, theta)
    minus_ex = -e_x_rot
    proj_on_wind = np.dot(minus_ex, wind_norm) * wind_norm
    e_x_perp = minus_ex - proj_on_wind
    wind_cross_z = np.array([wind_norm[1], -wind_norm[0], 0])
    hy = np.dot(e_x_perp, wind_cross_z)
    hz = e_x_perp[2]
    return hy, hz

def solve_heading_analytical(e_x, k, target_h, wind_norm):
    """Analytical solution for heading rotation angle."""
    k = k / norm(k)

    # Extract coefficients by sampling
    hy_0, hz_0 = get_heading_components(e_x, k, 0, wind_norm)
    hy_90, hz_90 = get_heading_components(e_x, k, np.pi/2, wind_norm)
    hy_180, hz_180 = get_heading_components(e_x, k, np.pi, wind_norm)

    C1 = (hy_0 + hy_180) / 2
    B1 = hy_0 - C1
    A1 = hy_90 - C1

    C2 = (hz_0 + hz_180) / 2
    B2 = hz_0 - C2
    A2 = hz_90 - C2

    ch = np.cos(target_h)
    sh = np.sin(target_h)

    A = A1 * ch - A2 * sh
    B = B1 * ch - B2 * sh
    C = C1 * ch - C2 * sh

    r = np.sqrt(A**2 + B**2)

    if r < 1e-10:
        return 0.0

    base_angle = np.arctan2(A, B)
    arg = np.clip(-C / r, -1, 1)
    delta = np.arccos(arg)

    theta1 = base_angle - delta
    return theta1

# Test configuration
wind_norm = np.array([1.0, 0.0, 0.0])

for elev_deg, azim_deg in [(10, 0), (45, 0), (45, 20), (45, -30), (30, 45), (-50, 0)]:
    print("=" * 60)
    print(f"ELEVATION = {elev_deg}°, AZIMUTH = {azim_deg}°")
    print("=" * 60)

    elevation = np.radians(elev_deg)
    azimuth = np.radians(azim_deg)

    # e_x: rotate [0,0,-1] by -elevation around y, then by -azimuth around z
    e_x = rotate_around_y(np.array([0.0, 0.0, -1.0]), -elevation)
    e_x = rotate_around_z(e_x, -azimuth)

    # k: tether direction at given elevation and azimuth
    k = np.array([np.cos(elevation) * np.cos(azimuth),
                  np.cos(elevation) * np.sin(azimuth),
                  np.sin(elevation)])

    print(f"e_x = {e_x}")
    print(f"k = {k}")
    print()

    all_pass = True
    for target_deg in [0, 30, 45, 60, 90, 120, 150, 180, -30, -60, -90, -120, -150]:
        target_h = np.radians(target_deg)
        theta_solved = solve_heading_analytical(e_x, k, target_h, wind_norm)
        e_x_rot = rodrigues(e_x, k, theta_solved)
        actual_h = calc_heading(e_x_rot, wind_norm)
        error = wrap_to_pi(actual_h - target_h)

        status = "✓" if abs(error) < 0.01 else "✗"
        if abs(error) >= 0.01:
            all_pass = False
        print(f"{status} Target: {target_deg:7.1f}°, θ: {np.degrees(theta_solved):8.2f}°, "
              f"actual: {np.degrees(actual_h):8.2f}°, error: {np.degrees(error):8.4f}°")

    if all_pass:
        print(f"\n✓ All tests passed for elevation {elev_deg}°!")
    else:
        print(f"\n✗ Some tests failed for elevation {elev_deg}° (singularity)")
    print()

print("=" * 60)
print("ANALYTICAL SOLUTION SUMMARY")
print("=" * 60)
print("""
Given:
  e_x = body x-axis after elevation/azimuth rotation (before heading)
  k   = rotation axis (normalized tether/wing direction)
  h   = target heading

Extract coefficients by sampling heading components at θ = 0, π/2, π:
  hy(θ), hz(θ) = heading_y and heading_z components

  C1 = (hy(0) + hy(π)) / 2,  B1 = hy(0) - C1,  A1 = hy(π/2) - C1
  C2 = (hz(0) + hz(π)) / 2,  B2 = hz(0) - C2,  A2 = hz(π/2) - C2

Combine for target heading h:
  A = A1·cos(h) - A2·sin(h)
  B = B1·cos(h) - B2·sin(h)
  C = C1·cos(h) - C2·sin(h)
  r = √(A² + B²)

Solution:
  θ = atan2(A, B) - acos(-C/r)   OR   θ = atan2(A, B) + acos(-C/r)

  (verify which solution gives correct heading)
""")
