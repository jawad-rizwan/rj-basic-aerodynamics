"""
ZRJ100 — 100-seat high-wing regional jet geometry and configuration.
All units Imperial (ft, lbs, deg).
"""

import numpy as np

AIRCRAFT = {
    "name": "ZRJ100",
    "description": "100-seat high-wing regional jet",

    # --- Wing (shared geometry) ---
    "S_ref": 792.47,              # ft^2, wing reference area (trapezoidal)
    "AR": 7.8,                    # aspect ratio
    "b": 78.62,                   # ft, wing span (from spreadsheet)
    "MAC": 10.93,                 # ft, mean aerodynamic chord (from spreadsheet)
    "taper": 0.33,                # taper ratio (from AVRO RJ)
    "t_c_wing": 0.12,             # thickness-to-chord ratio (NASA SC(3)-0712B, 12%)
    "x_c_max_wing": 0.37,         # chordwise location of max thickness
    "sweep_qc_deg": 22.9,         # quarter-chord sweep (deg)
    "sweep_le_deg": 26.0,         # leading-edge sweep (deg)
    "sweep_mt_deg": 21.4,         # max-thickness line sweep (deg)
    "S_exposed": 649.7,           # ft^2, exposed planform (minus fuselage cover)
    "winglet_h": 0.0,             # ft, no winglets (anhedral tips)

    # --- Fuselage ---
    "fuse_d": 9.83,               # ft, max fuselage diameter
    "fuse_length": 108.2,         # ft, fuselage length

    # --- Weights ---
    "MTOW": 81624.0,              # lbs

    # --- Horizontal tail (T-tail) ---
    "S_htail": 190.58,            # ft^2, htail reference area
    "S_htail_exposed": 174.8,     # ft^2 (exposed span 25.76 ft, minus structural box)
    "MAC_htail": 7.38,            # ft
    "t_c_htail": 0.12,            # *** UPDATE *** (airfoil TBD, using 12% placeholder)
    "x_c_max_htail": 0.37,        # *** UPDATE *** (airfoil TBD, using 37% placeholder)
    "sweep_mt_htail_deg": 23.8,   # *** UPDATE *** max-thickness sweep — recompute when airfoil chosen

    # --- Vertical tail (enhanced area) ---
    "S_vtail": 146.35,            # ft^2 (enhanced area)
    "S_vtail_exposed": 124.1,     # ft^2 (exposed span 13.73 ft)
    "MAC_vtail": 11.03,           # ft
    "t_c_vtail": 0.12,            # *** UPDATE *** (airfoil TBD, using 12% placeholder)
    "x_c_max_vtail": 0.37,        # *** UPDATE *** (airfoil TBD, using 37% placeholder)
    "sweep_mt_vtail_deg": 38.5,   # *** UPDATE *** max-thickness sweep — recompute when airfoil chosen

    # --- Nacelles (x2) ---
    "nacelle_length": 13.46,      # ft, per nacelle
    "nacelle_d": 6.17,            # ft, diameter
    "n_nacelles": 2,

    # --- Flight condition ---
    "cruise_mach": 0.78,
    "cruise_alt": 41000.0,        # ft

    # --- Airfoil data ---
    "cl_max_airfoil": 2.22,       # *** UPDATE *** 2D airfoil CLmax — verify from test data (typical SC ~1.6-1.8)

    # --- Surface finish ---
    "k_surface": 1.33e-5,         # ft, production sheet metal (Raymer Table 12.5)

    # --- Interference factors Q — Raymer Table 12.6 ---
    "Q_wing": 1.00,               # high wing, well-filleted (Raymer Table 12.6)
    "Q_fuse": 1.00,               # fuselage
    "Q_htail": 1.05,              # conventional tail group
    "Q_vtail": 1.05,              # conventional tail group
    "Q_nac": 1.30,                # nacelle ~1 diameter away from fuselage

    # --- Laminar flow percentages — Raymer Table 12.4 ---
    "lam_wing": 0.10,
    "lam_fuse": 0.05,
    "lam_tail": 0.05,
    "lam_nac": 0.0,

    # --- Leakage & protuberance — Raymer Table 12.9 ---
    "leak_pct": 0.03,

    # --- Leading-edge suction ---
    "S_suction": 0.90,            # typical for well-designed wing at design CL

    # --- High-lift devices ---
    "delta_cl_te_factor": 1.3 * 1.25,   # Fowler flap, c'/c = 1.25
    "te_flap_takeoff_fraction": 0.70,    # partial deflection for takeoff
    "s_flapped_ratio": 0.70,             # ~70% of wing span has flaps
    "sweep_hl_deg": 20.0,                # flap hinge line sweep
    "delta_cl_le_factor": 0.4 * 1.10,   # slat, c'/c = 1.10
    "s_slat_ratio": 0.75,
}
