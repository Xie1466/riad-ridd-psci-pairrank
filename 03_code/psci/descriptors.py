"""Exact array-level assembly rules for the frozen PSCI representation.

These functions consume already extracted descriptor arrays. They do not read
or process trajectories. Atom/residue selections and GROMACS provenance are
specified in the bundled exact feature dictionary and documentation.
"""

from __future__ import annotations

import numpy as np

PRIMARY_WINDOW_NS = (0.0, 75.0)
STRIDE_PS = 100
CONTACT_CUTOFF_NM = 0.45
MODEL_FEATURE_ORDER = ("D", "C", "F1", "F4", "F5", "L")


def window_mean(time_ns, values, start=0.0, end=75.0):
    time = np.asarray(time_ns, float)
    value = np.asarray(values, float)
    use = (time >= start) & (time <= end) & np.isfinite(value)
    if not np.any(use):
        raise ValueError("no finite descriptor values in requested inclusive window")
    return float(np.mean(value[use]))


def d_normalized(com_distance_nm, reference_rg1_nm, reference_rg2_nm):
    denominator = (float(reference_rg1_nm) + float(reference_rg2_nm)) / 2.0
    if denominator <= 0:
        raise ValueError("reference coordinate Rg must be positive")
    return np.asarray(com_distance_nm, float) / denominator


def c_density(contact_count, n_heavy_1, n_heavy_2):
    denominator = np.sqrt(int(n_heavy_1) * int(n_heavy_2))
    if denominator <= 0:
        raise ValueError("heavy-atom group sizes must be positive")
    return np.asarray(contact_count, float) / denominator


def component_frame0_ratio(component_1, component_2):
    a, b = np.asarray(component_1, float), np.asarray(component_2, float)
    if a[0] <= 0 or b[0] <= 0:
        raise ValueError("frame-0 reference must be positive")
    return 0.5 * (a / a[0] + b / b[0])


def orientation_occupancy(cosine_1, cosine_2, threshold=0.5):
    return ((np.asarray(cosine_1, float) >= threshold) &
            (np.asarray(cosine_2, float) >= threshold)).astype(float)


def f2_terminal_active(distance_1_nm, distance_2_nm, reference_rg1_nm, reference_rg2_nm):
    return 0.5 * (np.asarray(distance_1_nm, float) / float(reference_rg1_nm) +
                  np.asarray(distance_2_nm, float) / float(reference_rg2_nm))


def f4_nonself_occupancy(min_distance_1_nm, min_distance_2_nm, cutoff_nm=0.45):
    return ((np.asarray(min_distance_1_nm, float) <= cutoff_nm) |
            (np.asarray(min_distance_2_nm, float) <= cutoff_nm)).astype(float)


def f5_terminal_constraint(per_residue_rmsf_nm):
    values = np.asarray(per_residue_rmsf_nm, float)
    if not len(values) or not np.isfinite(values).all():
        raise ValueError("finite terminal RMSF values are required")
    return float(1.0 / (1.0 + np.mean(values)))


def linker_extension(end_to_end_1, contour_1, end_to_end_2, contour_2):
    c1, c2 = np.asarray(contour_1, float), np.asarray(contour_2, float)
    if np.any(c1 <= 0) or np.any(c2 <= 0):
        raise ValueError("instantaneous linker contours must be positive")
    return 0.5 * (np.asarray(end_to_end_1, float) / c1 +
                  np.asarray(end_to_end_2, float) / c2)


def f3_native_interface(*_args, **_kwargs):
    raise NotImplementedError("F3 is SCIENTIFICALLY_NOT_IDENTIFIABLE for Internal16")


def q_topology(*_args, **_kwargs):
    raise NotImplementedError("Q_topology is SCIENTIFICALLY_NOT_IDENTIFIABLE for Internal16")
