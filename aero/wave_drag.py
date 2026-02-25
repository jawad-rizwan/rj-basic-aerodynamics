"""
Transonic and supersonic wave drag estimation.
Raymer Chapter 12, Sections 12.5.9 and 12.5.10
"""

import numpy as np


def mdd_wing(t_c, sweep_qc_deg, cl_design=0.0, supercritical=False):
    """
    Wing drag-divergence Mach number (Boeing definition).
    Raymer Eq. 12.46:
        M_DD = M_DD(L=0) * L_FDD - 0.05 * CL_design

    M_DD(L=0) is read from Fig. 12.29 (approximated here).
    L_FDD is the lift adjustment from Fig. 12.30 (approximated here).

    Parameters
    ----------
    t_c : float
        Wing thickness-to-chord ratio.
    sweep_qc_deg : float
        Quarter-chord sweep (degrees).
    cl_design : float
        Design lift coefficient.
    supercritical : bool
        If True, use kA = 0.95 instead of 0.87.

    Returns
    -------
    float : M_DD (Boeing definition, ~20 counts drag rise)
    """
    # Approximate M_DD at zero lift from Raymer Fig. 12.29
    # Based on curve fits for various t/c vs sweep
    sweep_rad = np.radians(sweep_qc_deg)
    cos_sweep = np.cos(sweep_rad)

    # Korn equation approximation (widely used, consistent with Fig 12.29)
    # M_DD ≈ (kA / cos Λ) - (t/c) / cos^2(Λ) - CL / (10 * cos^3(Λ))
    # kA ≈ 0.87 for conventional airfoils, 0.95 for supercritical
    # The supercritical benefit is captured entirely by kA; t/c is NOT reduced.
    ka = 0.95 if supercritical else 0.87
    mdd = ka / cos_sweep - t_c / cos_sweep**2 - cl_design / (10.0 * cos_sweep**3)

    return mdd


def mdd_body(Ln, d):
    """
    Body (fuselage) drag-divergence Mach number.
    From Raymer Fig. 12.31, approximated.

    Parameters
    ----------
    Ln : float
        Length from nose to where cross-section becomes constant (ft).
    d : float
        Body diameter at that location (ft).

    Returns
    -------
    float : M_DD for the body
    """
    fineness = 2.0 * Ln / d
    # Approximation from Fig. 12.31
    if fineness >= 14:
        return 0.98
    elif fineness >= 10:
        return 0.92 + (fineness - 10) * 0.015
    elif fineness >= 6:
        return 0.82 + (fineness - 6) * 0.025
    else:
        return 0.70 + (fineness - 2) * 0.03


def wave_drag_sears_haack(A_max, length):
    """
    Sears-Haack body minimum wave drag (D/q).
    Raymer Eq. 12.44: (D/q)_wave = 9*pi/2 * (A_max / l)^2 * A_max

    (Actually Eq. 12.44 gives D/q = (9*pi*A_max^2) / (2*l^2),
    but dimensionally this simplifies to D/q in area units.)

    Parameters
    ----------
    A_max : float
        Maximum cross-sectional area (ft^2).
    length : float
        Total aircraft length, minus constant-section portions (ft).

    Returns
    -------
    float : (D/q) Sears-Haack wave drag (ft^2)
    """
    # Raymer Eq. 12.44
    return 9.0 * np.pi / 2.0 * (A_max / length) ** 2 * A_max


def wave_drag_supersonic(A_max, length, E_wd, mach, sweep_le_deg):
    """
    Supersonic wave drag estimation using Sears-Haack correlation.
    Raymer Eq. 12.45:
        (D/q)_wave = E_WD * [1 - 0.386*(M-1.2)^0.57] *
                     (1 - pi*Λ_LE^0.77 / 100) * (D/q)_Sears-Haack

    Parameters
    ----------
    A_max : float
        Maximum cross-sectional area (ft^2). Subtract inlet capture area.
    length : float
        Effective aircraft length (ft). Subtract constant-section portions.
    E_wd : float
        Wave-drag efficiency factor.
        1.0 = perfect Sears-Haack
        1.2 = very clean blended design
        1.8-2.2 = typical supersonic fighter/SST
        2.5-3.0 = poor supersonic design
    mach : float
        Freestream Mach (> 1.0).
    sweep_le_deg : float
        Leading-edge sweep (degrees).

    Returns
    -------
    float : (D/q) wave drag (ft^2)
    """
    dq_sh = wave_drag_sears_haack(A_max, length)

    mach_term = 1.0 - 0.386 * (mach - 1.2) ** 0.57
    sweep_term = 1.0 - np.pi * sweep_le_deg**0.77 / 100.0

    # Raymer Eq. 12.45
    return E_wd * mach_term * sweep_term * dq_sh


def drag_rise_transonic(cd0_subsonic, mdd, cd_wave_m105, mach_array):
    """
    Construct transonic drag rise curve using Raymer Fig. 12.32 method.

    Points defined:
        E (M_cr) = M_DD - 0.08 : drag rise = 0
        D (M_DD) : drag rise = 0.002
        C (M=1.0) : drag rise ≈ 0.5 * cd_wave_m105
        B (M=1.05): drag rise = cd_wave_m105
        A (M≥1.2) : from supersonic wave drag

    Parameters
    ----------
    cd0_subsonic : float
        Subsonic CD0 (below M_cr).
    mdd : float
        Drag-divergence Mach number (Boeing definition).
    cd_wave_m105 : float
        Wave drag coefficient at Mach 1.05.
    mach_array : array-like
        Mach numbers to evaluate.

    Returns
    -------
    np.ndarray : Total CD0 at each Mach number
    """
    mach = np.asarray(mach_array, dtype=float)
    cd = np.full_like(mach, cd0_subsonic)

    m_cr = mdd - 0.08

    for i, m in enumerate(mach):
        if m <= m_cr:
            cd[i] = cd0_subsonic
        elif m <= mdd:
            # Smooth rise from 0 to 0.002 counts
            frac = (m - m_cr) / (mdd - m_cr)
            cd[i] = cd0_subsonic + 0.002 * frac**2
        elif m <= 1.0:
            # Rise from MDD to Mach 1.0
            frac = (m - mdd) / (1.0 - mdd)
            cd_at_1 = 0.5 * cd_wave_m105
            cd[i] = cd0_subsonic + 0.002 + (cd_at_1 - 0.002) * frac
        elif m <= 1.05:
            # Rise from Mach 1.0 to 1.05
            frac = (m - 1.0) / 0.05
            cd_at_1 = 0.5 * cd_wave_m105
            cd[i] = cd0_subsonic + cd_at_1 + (cd_wave_m105 - cd_at_1) * frac
        else:
            cd[i] = cd0_subsonic + cd_wave_m105

    return cd
