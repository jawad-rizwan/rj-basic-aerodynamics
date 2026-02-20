"""
Flat-plate skin-friction coefficient calculations.
Raymer Chapter 12, Section 12.5.3
"""

import numpy as np


def reynolds_number(rho, V, length, mu):
    """
    Reynolds number.
    Raymer Eq. 12.25: R = rho * V * l / mu

    Parameters
    ----------
    rho : float
        Air density (slug/ft^3).
    V : float
        Freestream velocity (ft/s).
    length : float
        Characteristic length (ft). Chord for wing/tail, total length for fuselage.
    mu : float
        Dynamic viscosity (slug/(ft·s)).

    Returns
    -------
    float : Reynolds number
    """
    return rho * V * length / mu


def cutoff_reynolds(length, k, mach=0.0):
    """
    Cutoff Reynolds number for surface roughness effects.
    Raymer Eq. 12.28 (subsonic) and Eq. 12.29 (transonic/supersonic).

    Parameters
    ----------
    length : float
        Characteristic length (ft).
    k : float
        Skin roughness height (ft). See Raymer Table 12.5:
            Smooth molded composite: 0.17e-5 ft
            Polished sheet metal:    0.50e-5 ft
            Production sheet metal:  1.33e-5 ft
            Smooth paint:            2.08e-5 ft
            Camouflage paint:        3.33e-5 ft

    mach : float
        Freestream Mach number.

    Returns
    -------
    float : Cutoff Reynolds number
    """
    if mach < 0.9:
        # Raymer Eq. 12.28
        return 38.21 * (length / k) ** 1.053
    else:
        # Raymer Eq. 12.29
        return 44.62 * (length / k) ** 1.053 * mach ** 1.16


def cf_laminar(Re):
    """
    Laminar flat-plate skin-friction coefficient.
    Raymer Eq. 12.26: Cf = 1.328 / sqrt(R)

    Parameters
    ----------
    Re : float
        Reynolds number.

    Returns
    -------
    float : Laminar Cf
    """
    return 1.328 / np.sqrt(Re)


def cf_turbulent(Re, mach=0.0):
    """
    Turbulent flat-plate skin-friction coefficient with compressibility correction.
    Raymer Eq. 12.27: Cf = 0.455 / (log10(R))^2.58 / (1 + 0.144*M^2)^0.65

    Parameters
    ----------
    Re : float
        Reynolds number.
    mach : float
        Freestream Mach number.

    Returns
    -------
    float : Turbulent Cf
    """
    return 0.455 / (np.log10(Re) ** 2.58 * (1.0 + 0.144 * mach**2) ** 0.65)


def cf_avg(Re, mach, pct_laminar, k=None, length=None):
    """
    Weighted-average skin-friction coefficient combining laminar and turbulent regions.
    Applies cutoff Reynolds number correction for roughness if k and length provided.

    Parameters
    ----------
    Re : float
        Actual Reynolds number.
    mach : float
        Freestream Mach number.
    pct_laminar : float
        Fraction of wetted area with laminar flow (0.0 to 1.0).
    k : float, optional
        Skin roughness height (ft).
    length : float, optional
        Characteristic length (ft).

    Returns
    -------
    float : Weighted-average Cf
    """
    Re_eff = Re
    if k is not None and length is not None:
        Re_cutoff = cutoff_reynolds(length, k, mach)
        Re_eff = min(Re, Re_cutoff)

    cf_lam = cf_laminar(Re_eff)
    cf_turb = cf_turbulent(Re_eff, mach)

    return pct_laminar * cf_lam + (1.0 - pct_laminar) * cf_turb
