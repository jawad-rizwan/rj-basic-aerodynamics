"""
Component form factors for parasite drag estimation.
Raymer Chapter 12, Section 12.5.4
"""

import numpy as np


def ff_wing(x_c_max, t_c, mach, sweep_max_t):
    """
    Form factor for wing, tail, strut, or pylon.
    Raymer Eq. 12.30:
        FF = [1 + 0.6/(x/c)_m * (t/c) + 100*(t/c)^4] * [1.34 * M^0.18 * (cos Λ_m)^0.28]

    Parameters
    ----------
    x_c_max : float
        Chordwise location of max thickness (fraction of chord).
        ~0.3 for low-speed airfoils, ~0.5 for high-speed (supercritical).
    t_c : float
        Thickness-to-chord ratio.
    mach : float
        Freestream Mach number.
    sweep_max_t : float
        Sweep of the maximum-thickness line (radians).

    Returns
    -------
    float : Form factor FF
    """
    term1 = 1.0 + 0.6 / x_c_max * t_c + 100.0 * t_c**4
    term2 = 1.34 * mach**0.18 * np.cos(sweep_max_t) ** 0.28
    return term1 * term2


def ff_fuselage(length, d):
    """
    Form factor for fuselage or smooth canopy.
    Raymer Eq. 12.31: FF = 0.9 + 5/f + f/400

    where f = l / d (fineness ratio), Eq. 12.33.

    A form factor of no less than 1.05 is recommended.

    Parameters
    ----------
    length : float
        Fuselage length (ft).
    d : float
        Fuselage equivalent diameter (ft).
        For non-circular: d = sqrt(4*A_max / pi).

    Returns
    -------
    float : Form factor FF
    """
    # Raymer Eq. 12.33
    f = length / d
    ff = 0.9 + 5.0 / f + f / 400.0
    return max(ff, 1.05)


def ff_nacelle(length, d):
    """
    Form factor for nacelle or smooth external store.
    Raymer Eq. 12.32: FF = 1 + 0.35 / f

    where f = l / d (fineness ratio), Eq. 12.33.

    Parameters
    ----------
    length : float
        Nacelle length (ft).
    d : float
        Nacelle equivalent diameter (ft).

    Returns
    -------
    float : Form factor FF
    """
    # Raymer Eq. 12.33
    f = length / d
    return 1.0 + 0.35 / f


def ff_tail_with_hinge(x_c_max, t_c, mach, sweep_max_t):
    """
    Form factor for tail surface with hinged rudder or elevator.
    Raymer Section 12.5.4: ~10% higher than Eq. 12.30 due to
    gap drag between tail and control surface.

    Parameters
    ----------
    (same as ff_wing)

    Returns
    -------
    float : Form factor FF (with 10% increment on pressure portion)
    """
    ff_base = ff_wing(x_c_max, t_c, mach, sweep_max_t)
    # Apply 10% increment to the pressure-caused portion (above 1.0)
    return 1.0 + 1.10 * (ff_base - 1.0)
