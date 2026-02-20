"""
Standard atmosphere model for aerodynamic calculations.
Uses the 1976 US Standard Atmosphere (ISA).
"""

import numpy as np


def isa_properties(altitude_ft):
    """
    Compute ISA standard atmosphere properties.

    Parameters
    ----------
    altitude_ft : float
        Geometric altitude in feet.

    Returns
    -------
    dict with keys:
        rho     : density (slug/ft^3)
        T       : temperature (R)
        P       : pressure (lb/ft^2)
        a       : speed of sound (ft/s)
        mu      : dynamic viscosity (slug/(ft·s))
        rho_mu  : rho/mu ratio for Reynolds number (1/(ft·s))
    """
    h = altitude_ft

    # Sea-level constants
    rho_sl = 0.002377   # slug/ft^3
    T_sl = 518.67       # R (59°F)
    P_sl = 2116.22      # lb/ft^2
    a_sl = 1116.45      # ft/s
    mu_sl = 3.737e-7    # slug/(ft·s)

    if h <= 36089:
        # Troposphere: linear lapse rate
        T = T_sl - 0.00356616 * h
        P = P_sl * (T / T_sl) ** 5.2561
    elif h <= 65617:
        # Lower stratosphere: isothermal
        T = 389.97  # R
        P = P_sl * 0.22336 * np.exp(-4.80614e-5 * (h - 36089))
    elif h <= 104987:
        # Upper stratosphere: inversion
        T = 389.97 + 0.00054864 * (h - 65617)
        P = P_sl * 0.02541 * (T / 389.97) ** -34.1632
    else:
        raise ValueError(f"Altitude {h} ft exceeds model range (0-104987 ft)")

    rho = P / (1716.49 * T)
    a = np.sqrt(1.4 * 1716.49 * T)

    # Sutherland's law for dynamic viscosity
    T_R = T
    mu = mu_sl * (T_R / T_sl) ** 1.5 * (T_sl + 198.72) / (T_R + 198.72)

    return {
        "rho": rho,
        "T": T,
        "P": P,
        "a": a,
        "mu": mu,
        "rho_mu": rho / mu,
    }


def dynamic_pressure(mach, altitude_ft):
    """Compute dynamic pressure q = 0.5 * rho * V^2 = 0.7 * P * M^2."""
    atm = isa_properties(altitude_ft)
    V = mach * atm["a"]
    return 0.5 * atm["rho"] * V**2
