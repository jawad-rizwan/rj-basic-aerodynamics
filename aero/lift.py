"""
Lift calculations: lift-curve slope, maximum lift.
Raymer Chapter 12, Section 12.4
"""

import numpy as np


def cl_alpha_subsonic(A, mach, sweep_max_t, eta=0.95,
                      s_exposed_ratio=1.0, F=1.0):
    """
    Subsonic wing lift-curve slope (per radian).
    Raymer Eq. 12.6:
        CLa = (2*pi*A) / (2 + sqrt(4 + (A*beta/eta)^2 * (1 + tan^2(Λ_mt)/beta^2))) * (S_exp/S_ref) * F

    where beta^2 = 1 - M^2  (Eq. 12.7)
    and   eta = Cl_alpha / (2*pi/beta)  (Eq. 12.8), typically ~0.95

    Parameters
    ----------
    A : float
        Wing aspect ratio.
    mach : float
        Freestream Mach number (must be < 1.0).
    sweep_max_t : float
        Sweep of maximum-thickness line (radians).
    eta : float
        Airfoil efficiency, typically 0.95. (Eq. 12.8)
    s_exposed_ratio : float
        S_exposed / S_ref. Ratio of exposed planform to reference area.
    F : float
        Fuselage lift factor. (Eq. 12.9)

    Returns
    -------
    float : CLa (per radian)
    """
    beta_sq = 1.0 - mach**2
    if beta_sq <= 0:
        raise ValueError(f"Mach {mach} >= 1.0; use cl_alpha_supersonic instead.")

    beta = np.sqrt(beta_sq)
    AB_eta = A * beta / eta

    inner = 4.0 + AB_eta**2 * (1.0 + np.tan(sweep_max_t)**2 / beta_sq)
    cl_a = (2.0 * np.pi * A) / (2.0 + np.sqrt(inner))
    cl_a *= s_exposed_ratio * F

    return cl_a


def fuselage_lift_factor(d, b):
    """
    Fuselage lift factor accounting for lift spill-over.
    Raymer Eq. 12.9: F = 1.07 * (1 + d/b)^2

    Parameters
    ----------
    d : float
        Fuselage diameter (ft).
    b : float
        Wing span (ft).

    Returns
    -------
    float : F
    """
    return 1.07 * (1.0 + d / b) ** 2


def effective_aspect_ratio_winglet(A, h, b):
    """
    Effective aspect ratio with winglets.
    Raymer Eq. 12.11: A_eff = A * (1 + h/b)^2

    Parameters
    ----------
    A : float
        Geometric aspect ratio.
    h : float
        Winglet height (ft).
    b : float
        Wing span (ft).

    Returns
    -------
    float : Effective aspect ratio
    """
    return A * (1.0 + h / b) ** 2


def cl_alpha_supersonic(mach):
    """
    Supersonic 2-D theoretical lift-curve slope.
    Raymer Eq. 12.12: CLa = 4 / beta

    where beta = sqrt(M^2 - 1)  (Eq. 12.13)

    Parameters
    ----------
    mach : float
        Freestream Mach number (must be > 1.0).

    Returns
    -------
    float : CLa (per radian)
    """
    if mach <= 1.0:
        raise ValueError(f"Mach {mach} <= 1.0; use cl_alpha_subsonic instead.")
    beta = np.sqrt(mach**2 - 1.0)
    return 4.0 / beta


def cl_max_clean(cl_max_airfoil, sweep_qc):
    """
    Maximum lift coefficient for a clean wing (no flaps).
    Raymer Eq. 12.15: CLmax = 0.9 * Cl_max_airfoil * cos(Λ_0.25c)

    Valid for moderate sweep, high-aspect-ratio wings.

    Parameters
    ----------
    cl_max_airfoil : float
        2-D airfoil maximum lift coefficient (from airfoil data at similar Re).
    sweep_qc : float
        Quarter-chord sweep angle (radians).

    Returns
    -------
    float : Wing CLmax (clean)
    """
    return 0.9 * cl_max_airfoil * np.cos(sweep_qc)


def cl_max_flaps(cl_max_clean, delta_cl_max_airfoil, s_flapped_ratio,
                 sweep_hl, delta_cl_max_le=0.0):
    """
    Maximum lift with high-lift devices deployed.
    Raymer Eq. 12.21:
        ΔCLmax = 0.9 * ΔCl_max * (S_flapped/S_ref) * cos(Λ_HL)

    Parameters
    ----------
    cl_max_clean : float
        Clean wing CLmax.
    delta_cl_max_airfoil : float
        2-D airfoil lift increment from flap/slat.
        Raymer Table 12.2:
            Plain/split: 0.9
            Slotted:     1.3
            Fowler:      1.3 * c'/c
            Double slot: 1.6 * c'/c
            Triple slot: 1.9 * c'/c
            Slat:        0.4 * c'/c
    s_flapped_ratio : float
        S_flapped / S_ref (ratio of flapped wing area to reference area).
    sweep_hl : float
        Hinge-line sweep angle (radians).
    delta_cl_max_le : float
        Additional CLmax from leading-edge devices.

    Returns
    -------
    float : Total wing CLmax with high-lift devices
    """
    # Raymer Eq. 12.21
    delta_cl = 0.9 * delta_cl_max_airfoil * s_flapped_ratio * np.cos(sweep_hl)
    return cl_max_clean + delta_cl + delta_cl_max_le
