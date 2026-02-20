"""
Parasite (zero-lift) drag estimation methods.
Raymer Chapter 12, Sections 12.5.1 and 12.5.2
"""

import numpy as np
from .skin_friction import reynolds_number, cf_avg
from .form_factors import ff_wing, ff_fuselage, ff_nacelle, ff_tail_with_hinge

# Raymer Table 12.3: Equivalent skin-friction coefficients
CFE_TABLE = {
    "bomber": 0.0030,
    "civil_transport": 0.0026,
    "military_cargo": 0.0035,
    "air_force_fighter": 0.0035,
    "navy_fighter": 0.0040,
    "clean_supersonic": 0.0025,
    "light_single": 0.0055,
    "light_twin": 0.0045,
    "prop_seaplane": 0.0065,
    "jet_seaplane": 0.0040,
    "uav_prop": 0.0110,
}


def cd0_equivalent_cfe(s_wet, s_ref, aircraft_type="civil_transport"):
    """
    Quick CD0 estimate using equivalent skin-friction method.
    Raymer Eq. 12.23: CD0 = Cfe * (S_wet / S_ref)

    Parameters
    ----------
    s_wet : float
        Total aircraft wetted area (ft^2).
    s_ref : float
        Wing reference area (ft^2).
    aircraft_type : str
        Key from CFE_TABLE.

    Returns
    -------
    float : CD0 estimate
    """
    cfe = CFE_TABLE[aircraft_type]
    return cfe * s_wet / s_ref


def cd0_component_buildup(components, s_ref, mach, altitude_ft,
                          cd_misc=0.0, leak_pct=0.03):
    """
    Subsonic parasite drag via component buildup method.
    Raymer Eq. 12.24:
        CD0 = Σ(Cf_c * FF_c * Q_c * S_wet_c) / S_ref + CD_misc + CD_L&P

    Parameters
    ----------
    components : list of dict
        Each dict contains:
            name       : str   - component name
            s_wet      : float - wetted area (ft^2)
            length     : float - characteristic length (ft)
            ff         : float - form factor (pre-computed)
            Q          : float - interference factor (from Raymer Table 12.6)
            pct_laminar: float - fraction of laminar flow (0-1)
            k          : float - surface roughness (ft), optional
    s_ref : float
        Wing reference area (ft^2).
    mach : float
        Freestream Mach number.
    altitude_ft : float
        Altitude (ft).
    cd_misc : float
        Miscellaneous drag coefficient (flaps, upsweep, base, etc.).
    leak_pct : float
        Leakage and protuberance drag as fraction of parasite drag.
        Raymer Table 12.9: 0.02-0.05 for jet transports.

    Returns
    -------
    dict with:
        cd0_total  : float - total CD0
        cd0_skin   : float - skin friction + form + interference
        cd_misc    : float - miscellaneous
        cd_leak    : float - leakage and protuberance
        breakdown  : list  - per-component CD0 values
    """
    from .atmosphere import isa_properties

    atm = isa_properties(altitude_ft)
    V = mach * atm["a"]
    rho = atm["rho"]
    mu = atm["mu"]

    cd0_skin = 0.0
    breakdown = []

    for comp in components:
        Re = reynolds_number(rho, V, comp["length"], mu)

        k = comp.get("k", None)
        cf = cf_avg(Re, mach, comp["pct_laminar"], k=k, length=comp["length"])

        cd_comp = cf * comp["ff"] * comp["Q"] * comp["s_wet"] / s_ref
        cd0_skin += cd_comp

        breakdown.append({
            "name": comp["name"],
            "Re": Re,
            "Cf": cf,
            "FF": comp["ff"],
            "Q": comp["Q"],
            "S_wet": comp["s_wet"],
            "CD0": cd_comp,
        })

    cd_leak = leak_pct * cd0_skin
    cd0_total = cd0_skin + cd_misc + cd_leak

    return {
        "cd0_total": cd0_total,
        "cd0_skin": cd0_skin,
        "cd_misc": cd_misc,
        "cd_leak": cd_leak,
        "breakdown": breakdown,
    }
