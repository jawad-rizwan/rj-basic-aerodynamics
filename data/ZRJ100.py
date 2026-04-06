"""
ZRJ100 — 100-seat high-wing regional jet geometry and configuration.
All units Imperial (ft, lbs, deg).
"""

import numpy as np

AIRCRAFT = {
    "name": "ZRJ100",
    "description": "100-seat high-wing regional jet",

    # --- Wing (shared geometry) ---
    "S_ref": 1016.58,             # ft^2, wing reference area (trapezoidal)
    "AR": 7.8,                    # aspect ratio
    "b": 89.05,                   # ft, wing span
    "MAC": 12.38,                 # ft, mean aerodynamic chord
    "taper": 0.33,                # taper ratio
    "t_c_wing": 0.123,            # chord-weighted avg t/c (root 0.14, mid 0.12, tip 0.10)
    "x_c_max_wing": 0.37,         # chordwise location of max thickness
    "sweep_qc_deg": 22.9,         # quarter-chord sweep (deg)
    "sweep_le_deg": 26.0,         # leading-edge sweep (deg)
    "sweep_mt_deg": 21.4,         # max-thickness line sweep (deg)
    "S_exposed": 843.5,           # ft^2, exposed planform (minus fuselage cover)
    "winglet_h": 6.9741,           # ft, scimitar winglet (5.67 upper + 1.30 lower)

    # --- Fuselage ---
    "fuse_d": 10.5,               # ft, max fuselage diameter
    "fuse_length": 110.2,         # ft, fuselage length

    # --- Weights ---
    "MTOW": 93255.0,              # lbs

    # --- Horizontal tail (T-tail) ---
    "S_htail": 276.92,            # ft^2, htail reference area
    "S_htail_exposed": 260.2,     # ft^2, *** UPDATE *** estimate — verify structural box width
    "MAC_htail": 8.90,            # ft
    "t_c_htail": 0.10,            # NASA SC(2)-0010 airfoil
    "x_c_max_htail": 0.37,        # chordwise location of max thickness
    "taper_htail": 0.4,            # htail taper ratio
    "sweep_mt_htail_deg": 23.7,   # max-thickness line sweep (deg)
    "cl_max_htail_airfoil": 1.0,  # *** UPDATE *** 2D cl_max for SC(2)-0010 — verify from XFLR5 at tail Re

    # --- Vertical tail ---
    "S_vtail": 200.8,             # ft^2
    "S_vtail_exposed": 172.5,     # ft^2, *** UPDATE *** estimate — verify structural carrythrough
    "MAC_vtail": 14.46,           # ft
    "t_c_vtail": 0.12,            # NASA SC(2)-0012 airfoil
    "x_c_max_vtail": 0.37,        # chordwise location of max thickness
    "taper_vtail": 0.6,            # vtail taper ratio
    "AR_vtail": 1.0,               # vtail aspect ratio
    "sweep_mt_vtail_deg": 31.8,   # max-thickness line sweep (deg)
    "cl_max_vtail_airfoil": 1.1,  # *** UPDATE *** 2D cl_max for SC(2)-0012 — verify from XFLR5 at tail Re

    # --- Nacelles (x2) ---
    "nacelle_length": 13.45833333,  # ft, per nacelle
    "nacelle_d": 6.166666667,       # ft, diameter
    "n_nacelles": 2,

    # --- Belly fairing (landing gear fairing) ---
    "fairing_s_wet_delta": 37.412,  # ft^2, incremental wetted area vs no-fairing baseline
    "fairing_length": 20.0,         # ft, 240 in
    "fairing_d_eq": 3.008,          # ft, equivalent diameter from 108.55 in width x 12 in depth
    "Q_fairing": 1.05,              # well-blended external fairing
    "lam_fairing": 0.0,             # assume fully turbulent
    "k_fairing": 1.33e-5,           # ft, production sheet metal

    # --- Flight condition ---
    "cruise_mach": 0.78,
    "cruise_alt": 35000.0,        # ft

    # --- Airfoil data ---
    "cl_max_airfoil": 1.40,       # 2D cl_max (tip-limited, NASA SC(2)-0710 at Re ~9M)
    "alpha_0L_deg": -2.0,         # deg, *** UPDATE *** zero-lift AoA — estimate for supercritical, verify from XFLR5
    "cm_0_airfoil": -0.09,        # *** UPDATE *** 2D Cm_0 at aero center — estimate for supercritical, verify from XFLR5

    # --- Surface finish ---
    "k_composite": 0.50e-5,       # ft, polished composite (Raymer Table 12.5) — wing, tail, nacelles
    "k_metal": 1.33e-5,           # ft, production sheet metal (Raymer Table 12.5) — fuselage

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
    "delta_cl_te_factor": 1.6 * 1.25,   # double slotted flap, c'/c = 1.25 *** UPDATE *** verify c'/c from flap design
    "te_flap_takeoff_fraction": 0.70,    # *** UPDATE *** partial deflection for takeoff
    "s_flapped_ratio": 0.70,             # *** UPDATE *** ~70% of wing span has flaps — verify from planform
    "sweep_hl_deg": 20.0,                # *** UPDATE *** flap hinge line sweep — verify from planform
    "delta_cl_le_factor": 0.4 * 1.10,   # slotted slat, c'/c = 1.10 *** UPDATE *** verify c'/c from slat design
    "s_slat_ratio": 0.75,               # *** UPDATE *** ~75% of wing span has slats — verify from planform
}
