"""
Drag-due-to-lift (induced drag) calculations.
Raymer Chapter 12, Section 12.6
"""

import numpy as np


def oswald_e(A, sweep_le_deg):
    """
    Oswald span efficiency factor.
    Raymer Eq. 12.48 (straight wing) and Eq. 12.49 (swept wing, 30-60 deg LE sweep).
    Linearly interpolated for 0-30 deg sweep.

    Parameters
    ----------
    A : float
        Wing aspect ratio (use effective A if winglets present).
    sweep_le_deg : float
        Leading-edge sweep angle (degrees).

    Returns
    -------
    float : Oswald efficiency factor e
    """
    # Raymer Eq. 12.48 — straight-wing aircraft
    e_straight = 1.78 * (1.0 - 0.045 * A**0.68) - 0.64

    if sweep_le_deg <= 0:
        return max(e_straight, 0.1)

    # Raymer Eq. 12.49 — swept-wing aircraft (30 < Λ_LE < 60 deg)
    sweep_le_rad = np.radians(sweep_le_deg)
    e_swept = (4.61 * (1.0 - 0.045 * A**0.68)
               * np.cos(sweep_le_rad)**0.15 - 3.1)

    if sweep_le_deg >= 30:
        return max(e_swept, 0.1)

    # Linear interpolation between 0 and 30 deg
    frac = sweep_le_deg / 30.0
    e = e_straight * (1.0 - frac) + e_swept * frac
    return max(e, 0.1)


def k_factor(A, e):
    """
    Drag-due-to-lift factor K from Oswald efficiency.
    Raymer Eq. 12.47: K = 1 / (pi * A * e)

    Parameters
    ----------
    A : float
        Wing aspect ratio (use effective A if winglets present).
    e : float
        Oswald span efficiency factor.

    Returns
    -------
    float : K factor
    """
    return 1.0 / (np.pi * A * e)


def k_100pct(A):
    """
    K value at 100% leading-edge suction (ideal induced drag only).
    K_100 = 1 / (pi * A)

    Parameters
    ----------
    A : float
        Aspect ratio.

    Returns
    -------
    float
    """
    return 1.0 / (np.pi * A)


def k_0pct(cl_alpha):
    """
    K value at 0% leading-edge suction.
    Raymer Eq. 12.56: K_0 = 1 / CL_alpha (per radian)

    Parameters
    ----------
    cl_alpha : float
        Lift-curve slope (per radian).

    Returns
    -------
    float
    """
    return 1.0 / cl_alpha


def k_factor_leading_edge_suction(A, cl_alpha, S_pct):
    """
    K factor using the leading-edge suction method.
    Raymer Eq. 12.57: K = S * K_100 + (1 - S) * K_0

    This method is superior to the simple Oswald method as it accounts
    for viscous separation changes with lift coefficient and Mach effects.

    Parameters
    ----------
    A : float
        Wing aspect ratio.
    cl_alpha : float
        Lift-curve slope (per radian) at the flight condition.
    S_pct : float
        Leading-edge suction parameter (0.0 to 1.0).
        Typically 0.85-0.95 for subsonic cruise of a well-designed wing.
        See Raymer Fig. 12.38 for values vs CL and design CL.

    Returns
    -------
    float : K factor
    """
    k100 = k_100pct(A)
    k0 = k_0pct(cl_alpha)
    # Raymer Eq. 12.57
    return S_pct * k100 + (1.0 - S_pct) * k0


def k_supersonic(A, mach, sweep_le_rad):
    """
    Supersonic drag-due-to-lift factor.
    Raymer Eq. 12.51:
        K = A*(M^2 - 1) * cos(Λ_LE) / (4*A*sqrt(M^2-1) - 2)

    Parameters
    ----------
    A : float
        Aspect ratio.
    mach : float
        Freestream Mach (must be > 1.0).
    sweep_le_rad : float
        Leading-edge sweep (radians).

    Returns
    -------
    float : K factor (supersonic)
    """
    beta = np.sqrt(mach**2 - 1.0)
    num = A * (mach**2 - 1.0) * np.cos(sweep_le_rad)
    den = 4.0 * A * beta - 2.0
    return num / den


def ground_effect_factor(h, b):
    """
    Ground effect reduction factor on K.
    Raymer Eq. 12.60: K_eff/K = 33*(h/b)^1.5 / (1 + 33*(h/b)^1.5)

    Parameters
    ----------
    h : float
        Wing height above ground at quarter chord (ft).
    b : float
        Wing span (ft).

    Returns
    -------
    float : Multiplier on K (< 1.0 means drag reduction)
    """
    hb = h / b
    val = 33.0 * hb**1.5
    return val / (1.0 + val)
