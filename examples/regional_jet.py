"""
Regional Jet Aerodynamic Analysis — high-wing regional jet
Using Raymer Chapter 12 methods.

This script computes CD0, Oswald e, K, CLmax, L/D, MDD, and required
cruise thrust for a high-wing regional jet using the component buildup method.

Aircraft data is loaded from dedicated files in the data/ directory:
  - data/ZRJ70.py   (76-seat variant)
  - data/ZRJ100.py  (100-seat variant)

Outputs are intended to feed into Chapter 6 refined sizing
(rj-mission-sizing).
"""

import sys
import os
import importlib
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aero.atmosphere import isa_properties, dynamic_pressure
from aero.skin_friction import reynolds_number
from aero.form_factors import ff_wing, ff_fuselage, ff_nacelle, ff_tail_with_hinge
from aero.parasite_drag import cd0_component_buildup, cd0_equivalent_cfe
from aero.lift import (
    cl_alpha_subsonic, fuselage_lift_factor,
    cl_max_clean, cl_max_flaps, effective_aspect_ratio_winglet,
)
from aero.induced_drag import oswald_e, k_factor, k_factor_leading_edge_suction
from aero.wave_drag import mdd_wing
from aero.drag_polar import DragPolar


def analyse(ac):
    """Run full Raymer Ch.12 aerodynamic analysis for one aircraft configuration."""

    name = ac["name"]

    # --- Unpack geometry ---
    S_ref        = ac["S_ref"]
    AR           = ac["AR"]
    b            = ac["b"]
    MAC          = ac["MAC"]
    t_c_wing     = ac["t_c_wing"]
    x_c_max_wing = ac["x_c_max_wing"]
    sweep_qc_rad = np.radians(ac["sweep_qc_deg"])
    sweep_le_deg = ac["sweep_le_deg"]
    sweep_mt_rad = np.radians(ac["sweep_mt_deg"])
    S_exposed    = ac["S_exposed"]
    winglet_h    = ac["winglet_h"]

    fuse_d       = ac["fuse_d"]
    fuse_length  = ac["fuse_length"]
    MTOW         = ac["MTOW"]
    W_cruise_approx = 0.92 * MTOW

    S_htail_exposed = ac["S_htail_exposed"]
    MAC_htail    = ac["MAC_htail"]
    t_c_htail    = ac["t_c_htail"]
    x_c_max_htail = ac["x_c_max_htail"]
    sweep_mt_htail_rad = np.radians(ac["sweep_mt_htail_deg"])

    S_vtail_exposed = ac["S_vtail_exposed"]
    MAC_vtail    = ac["MAC_vtail"]
    t_c_vtail    = ac["t_c_vtail"]
    x_c_max_vtail = ac["x_c_max_vtail"]
    sweep_mt_vtail_rad = np.radians(ac["sweep_mt_vtail_deg"])

    nacelle_length = ac["nacelle_length"]
    nacelle_d    = ac["nacelle_d"]
    n_nacelles   = ac["n_nacelles"]
    fairing_s_wet = ac.get("fairing_s_wet_delta", 0.0)
    fairing_length = ac.get("fairing_length", 0.0)
    fairing_d_eq = ac.get("fairing_d_eq", 0.0)

    cruise_mach  = ac["cruise_mach"]
    cruise_alt   = ac["cruise_alt"]
    atm          = isa_properties(cruise_alt)
    V_cruise     = cruise_mach * atm["a"]

    cl_max_airfoil = ac["cl_max_airfoil"]
    k_composite  = ac["k_composite"]
    k_metal      = ac["k_metal"]

    Q_wing       = ac["Q_wing"]
    Q_fuse       = ac["Q_fuse"]
    Q_htail      = ac["Q_htail"]
    Q_vtail      = ac["Q_vtail"]
    Q_nac        = ac["Q_nac"]
    Q_fairing    = ac.get("Q_fairing", 1.0)

    lam_wing     = ac["lam_wing"]
    lam_fuse     = ac["lam_fuse"]
    lam_tail     = ac["lam_tail"]
    lam_nac      = ac["lam_nac"]
    lam_fairing  = ac.get("lam_fairing", 0.0)
    leak_pct     = ac["leak_pct"]
    S_suction    = ac["S_suction"]
    k_fairing    = ac.get("k_fairing", k_metal)

    delta_cl_te  = ac["delta_cl_te_factor"]
    te_to_frac   = ac["te_flap_takeoff_fraction"]
    s_flapped_ratio = ac["s_flapped_ratio"]
    sweep_hl_rad = np.radians(ac["sweep_hl_deg"])
    delta_cl_le  = ac["delta_cl_le_factor"]
    s_slat_ratio = ac["s_slat_ratio"]

    # =========================================================================
    #  HEADER
    # =========================================================================
    print(f"\n{'#'*65}")
    print(f"  {name} — {ac['description']}")
    print(f"  (fuse = {fuse_length:.1f} ft, MTOW = {MTOW:.0f} lbs)")
    print(f"{'#'*65}")

    # =========================================================================
    #  1. LIFT-CURVE SLOPE
    # =========================================================================
    AR_eff = effective_aspect_ratio_winglet(AR, winglet_h, b)
    F = fuselage_lift_factor(fuse_d, b)
    s_exp_ratio = S_exposed / S_ref
    sf_product = min(s_exp_ratio * F, 0.98)

    CLa = cl_alpha_subsonic(
        A=AR_eff, mach=cruise_mach, sweep_max_t=sweep_mt_rad,
        eta=0.95, s_exposed_ratio=sf_product, F=1.0,
    )

    print(f"\n--- LIFT ---")
    print(f"  Wing span              = {b:.1f} ft")
    print(f"  MAC                    = {MAC:.2f} ft")
    print(f"  AR (geometric)         = {AR:.1f}")
    print(f"  AR (effective, winglet) = {AR_eff:.2f}")
    print(f"  Fuselage lift factor F = {F:.3f}")
    print(f"  (S_exp/S_ref)*F        = {s_exp_ratio * F:.3f}")
    print(f"  CLa (cruise, M={cruise_mach})   = {CLa:.4f} /rad = {CLa * np.pi/180:.5f} /deg")

    # =========================================================================
    #  2. WETTED AREAS & FORM FACTORS
    # =========================================================================
    S_wet_wing  = 2.0 * S_exposed * (1.0 + 0.25 * t_c_wing)
    S_wet_htail = 2.0 * S_htail_exposed * (1.0 + 0.25 * t_c_htail)
    S_wet_vtail = 2.0 * S_vtail_exposed * (1.0 + 0.25 * t_c_vtail)
    S_wet_nac   = np.pi * nacelle_d * nacelle_length * 0.8

    fuse_f = fuse_length / fuse_d
    S_wet_fuse = (np.pi * fuse_d * fuse_length
                  * (1.0 - 2.0/fuse_f)**(2.0/3.0)
                  * (1.0 + 1.0/fuse_f**2))

    S_wet_total = (S_wet_wing + S_wet_fuse + S_wet_htail
                   + S_wet_vtail + n_nacelles * S_wet_nac)

    FF_wing  = ff_wing(x_c_max_wing, t_c_wing, cruise_mach, sweep_mt_rad)
    FF_fuse  = ff_fuselage(fuse_length, fuse_d)
    FF_htail = ff_tail_with_hinge(x_c_max_htail, t_c_htail, cruise_mach, sweep_mt_htail_rad)
    FF_vtail = ff_tail_with_hinge(x_c_max_vtail, t_c_vtail, cruise_mach, sweep_mt_vtail_rad)
    FF_nac   = ff_nacelle(nacelle_length, nacelle_d)
    FF_fairing = ff_nacelle(fairing_length, fairing_d_eq) if fairing_s_wet > 0.0 else None

    print(f"\n--- WETTED AREAS ---")
    print(f"  Wing          = {S_wet_wing:.0f} ft^2")
    print(f"  Fuselage      = {S_wet_fuse:.0f} ft^2")
    print(f"  H-tail        = {S_wet_htail:.0f} ft^2")
    print(f"  V-tail        = {S_wet_vtail:.0f} ft^2")
    print(f"  Nacelles (x{n_nacelles}) = {n_nacelles * S_wet_nac:.0f} ft^2")
    if fairing_s_wet > 0.0:
        print(f"  Belly fairing = {fairing_s_wet:.0f} ft^2")
        S_wet_total += fairing_s_wet
    print(f"  TOTAL         = {S_wet_total:.0f} ft^2")
    print(f"  S_wet / S_ref = {S_wet_total / S_ref:.2f}")

    print(f"\n--- FORM FACTORS ---")
    print(f"  Wing (Eq. 12.30)     = {FF_wing:.4f}")
    print(f"  Fuselage (Eq. 12.31) = {FF_fuse:.4f}")
    print(f"  H-tail (Eq. 12.30+)  = {FF_htail:.4f}")
    print(f"  V-tail (Eq. 12.30+)  = {FF_vtail:.4f}")
    print(f"  Nacelle (Eq. 12.32)  = {FF_nac:.4f}")
    if FF_fairing is not None:
        print(f"  Fairing (Eq. 12.32)  = {FF_fairing:.4f}")

    # =========================================================================
    #  3. CD0 COMPONENT BUILDUP (Eq. 12.24)
    # =========================================================================
    components = [
        {"name": "Wing",       "s_wet": S_wet_wing,    "length": MAC,
         "ff": FF_wing,  "Q": Q_wing,  "pct_laminar": lam_wing, "k": k_composite},
        {"name": "Fuselage",   "s_wet": S_wet_fuse,    "length": fuse_length,
         "ff": FF_fuse,  "Q": Q_fuse,  "pct_laminar": lam_fuse, "k": k_metal},
        {"name": "H-tail",     "s_wet": S_wet_htail,   "length": MAC_htail,
         "ff": FF_htail, "Q": Q_htail, "pct_laminar": lam_tail, "k": k_composite},
        {"name": "V-tail",     "s_wet": S_wet_vtail,   "length": MAC_vtail,
         "ff": FF_vtail, "Q": Q_vtail, "pct_laminar": lam_tail, "k": k_composite},
    ]
    for i in range(n_nacelles):
        side = "L" if i == 0 else "R"
        components.append(
            {"name": f"Nacelle {side}", "s_wet": S_wet_nac, "length": nacelle_length,
             "ff": FF_nac, "Q": Q_nac, "pct_laminar": lam_nac, "k": k_composite})
    if fairing_s_wet > 0.0:
        components.append(
            {"name": "Fairing", "s_wet": fairing_s_wet, "length": fairing_length,
             "ff": FF_fairing, "Q": Q_fairing, "pct_laminar": lam_fairing, "k": k_fairing})

    result = cd0_component_buildup(
        components, S_ref, cruise_mach, cruise_alt,
        cd_misc=0.0, leak_pct=leak_pct,
    )

    print(f"\n--- CD0 COMPONENT BUILDUP (Eq. 12.24) ---")
    print(f"  {'Component':<12s} {'Re':>10s} {'Cf':>9s} {'FF':>7s} {'Q':>5s} {'Swet':>7s} {'CD0':>9s}")
    print(f"  {'-'*60}")
    for c in result["breakdown"]:
        print(f"  {c['name']:<12s} {c['Re']:10.0f} {c['Cf']:9.6f} {c['FF']:7.4f} {c['Q']:5.2f} {c['S_wet']:7.0f} {c['CD0']:9.6f}")
    print(f"  {'-'*60}")
    print(f"  {'Skin friction + FF + Q':<40s} = {result['cd0_skin']:.6f}")
    print(f"  {'Leakage & protuberance':<40s} = {result['cd_leak']:.6f}  ({leak_pct*100:.0f}%)")
    print(f"  {'Miscellaneous':<40s} = {result['cd_misc']:.6f}")
    print(f"  {'TOTAL CD0':<40s} = {result['cd0_total']:.6f}")

    cd0_check = cd0_equivalent_cfe(S_wet_total, S_ref, "civil_transport")
    print(f"\n  Cross-check: Cfe method (Eq. 12.23)    = {cd0_check:.6f}")

    # =========================================================================
    #  4. DRAG DUE TO LIFT
    # =========================================================================
    e = oswald_e(AR_eff, sweep_le_deg)
    K = k_factor(AR_eff, e)
    K_les = k_factor_leading_edge_suction(AR_eff, CLa, S_suction)
    e_from_les = 1.0 / (np.pi * AR_eff * K_les)

    print(f"\n--- DRAG DUE TO LIFT ---")
    print(f"  Oswald e (Eq. 12.48/49)    = {e:.4f}  (conservative)")
    print(f"  K (Oswald, Eq. 12.47)      = {K:.5f}")
    print(f"  K (LE suction, Eq. 12.57)  = {K_les:.5f}  (S = {S_suction})")
    print(f"  e equiv (LE suction)       = {e_from_les:.4f}  (more realistic)")

    # =========================================================================
    #  5. DRAG POLAR & L/D
    # =========================================================================
    cd0 = result["cd0_total"]
    polar_oswald = DragPolar(cd0, K)
    polar_les = DragPolar(cd0, K_les)
    polar = polar_les

    ld_max, cl_ldmax = polar.ld_max()
    ld_max_oswald, _ = polar_oswald.ld_max()

    q_cruise = dynamic_pressure(cruise_mach, cruise_alt)
    cl_cruise = polar.cl_for_cruise(W_cruise_approx, q_cruise, S_ref)
    ld_cruise = polar.ld(cl_cruise)
    cd_cruise = polar.cd(cl_cruise)

    print(f"\n--- DRAG POLAR (using LE suction K) ---")
    polar.summary()
    print(f"  (L/D)_max (Oswald K)  = {ld_max_oswald:.2f}  (conservative)")
    print(f"\n  Cruise condition (M={cruise_mach}, {cruise_alt:.0f} ft):")
    print(f"    q_cruise       = {q_cruise:.1f} lb/ft^2")
    print(f"    V_cruise       = {V_cruise:.0f} ft/s = {V_cruise/1.688:.0f} kt")
    print(f"    CL_cruise      = {cl_cruise:.4f}  (W={W_cruise_approx:.0f} lbs)")
    print(f"    CD_cruise      = {cd_cruise:.6f}")
    print(f"    L/D_cruise     = {ld_cruise:.2f}")

    # --- Required thrust for level unaccelerated cruise ---
    D_cruise = q_cruise * S_ref * cd_cruise          # lbs
    T_required = D_cruise                             # T = D in steady cruise
    T_per_engine = T_required / n_nacelles            # lbs per engine
    T_over_W = T_required / W_cruise_approx

    print(f"\n  Cruise thrust required (T = D):")
    print(f"    D_cruise       = {D_cruise:.0f} lbs")
    print(f"    T_required     = {T_required:.0f} lbs  (total)")
    print(f"    T per engine   = {T_per_engine:.0f} lbs  (x{n_nacelles})")
    print(f"    T/W_cruise     = {T_over_W:.4f}")

    # =========================================================================
    #  6. MAXIMUM LIFT
    # =========================================================================
    cl_max = cl_max_clean(cl_max_airfoil, sweep_qc_rad)

    cl_max_to = cl_max_flaps(cl_max, delta_cl_te * te_to_frac, s_flapped_ratio,
                              sweep_hl_rad, delta_cl_max_le=0.0)

    cl_max_land = cl_max_flaps(cl_max, delta_cl_te, s_flapped_ratio, sweep_hl_rad,
                                delta_cl_max_le=0.9 * delta_cl_le * s_slat_ratio * np.cos(sweep_hl_rad))

    print(f"\n--- MAXIMUM LIFT ---")
    print(f"  CL_max (clean, Eq. 12.15)  = {cl_max:.3f}")
    print(f"  CL_max (takeoff flaps)     = {cl_max_to:.3f}")
    print(f"  CL_max (landing, full)     = {cl_max_land:.3f}")

    # =========================================================================
    #  7. TRANSONIC
    # =========================================================================
    m_dd = mdd_wing(t_c_wing, ac["sweep_qc_deg"], cl_design=cl_cruise,
                     supercritical=True)

    print(f"\n--- TRANSONIC ---")
    print(f"  M_DD (Boeing, Eq. 12.46)   = {m_dd:.3f}")
    print(f"  M_cr (approx)              = {m_dd - 0.08:.3f}")

    # =========================================================================
    #  8. STABILITY & CONTROL VALUES
    # =========================================================================
    alpha_0L = ac["alpha_0L_deg"]
    cm_0_airfoil = ac["cm_0_airfoil"]

    # Wing angles of attack
    alpha_cl_max = cl_max / (CLa * np.pi/180) + alpha_0L      # deg
    alpha_cl_cruise = cl_cruise / (CLa * np.pi/180) + alpha_0L  # deg

    # Wing Cm_0 about AC (3D correction for swept wing)
    cos_qc = np.cos(sweep_qc_rad)
    cm_0_ac = cm_0_airfoil * (AR * cos_qc**2) / (AR + 2*cos_qc)

    # Horizontal tail: compute AR from taper, MAC, S
    taper_ht = ac["taper_htail"]
    S_htail = ac["S_htail"]
    c_root_ht = MAC_htail / ((2.0/3.0) * (1 + taper_ht + taper_ht**2) / (1 + taper_ht))
    b_ht = 2.0 * S_htail / (c_root_ht * (1 + taper_ht))
    AR_htail = b_ht**2 / S_htail
    s_exp_ratio_ht = S_htail_exposed / S_htail

    # Htail sweep of quarter-chord (derive from max-thickness sweep)
    sweep_qc_ht_rad = np.arctan(
        np.tan(sweep_mt_htail_rad)
        + (4.0/AR_htail) * (x_c_max_htail - 0.25) * (1 - taper_ht) / (1 + taper_ht)
    )

    CLa_htail = cl_alpha_subsonic(
        A=AR_htail, mach=cruise_mach, sweep_max_t=sweep_mt_htail_rad,
        eta=0.95, s_exposed_ratio=s_exp_ratio_ht, F=1.0,
    )

    cl_max_htail_airfoil = ac["cl_max_htail_airfoil"]
    CLmax_htail = 0.9 * cl_max_htail_airfoil * np.cos(sweep_qc_ht_rad)
    alpha_CLmax_htail = CLmax_htail / (CLa_htail * np.pi/180) + 0.0  # symmetric airfoil, alpha_0L ~ 0

    # Vertical tail: AR given directly
    AR_vtail = ac["AR_vtail"]
    taper_vt = ac["taper_vtail"]
    S_vtail = ac["S_vtail"]
    s_exp_ratio_vt = S_vtail_exposed / S_vtail

    sweep_qc_vt_rad = np.arctan(
        np.tan(sweep_mt_vtail_rad)
        + (4.0/AR_vtail) * (x_c_max_vtail - 0.25) * (1 - taper_vt) / (1 + taper_vt)
    )

    CLa_vtail = cl_alpha_subsonic(
        A=AR_vtail, mach=cruise_mach, sweep_max_t=sweep_mt_vtail_rad,
        eta=0.95, s_exposed_ratio=s_exp_ratio_vt, F=1.0,
    )

    cl_max_vtail_airfoil = ac["cl_max_vtail_airfoil"]
    CLmax_vtail = 0.9 * cl_max_vtail_airfoil * np.cos(sweep_qc_vt_rad)
    alpha_CLmax_vtail = CLmax_vtail / (CLa_vtail * np.pi/180) + 0.0  # symmetric airfoil

    print(f"\n--- STABILITY & CONTROL VALUES ---")
    print(f"  (Values marked * are estimates — see *** UPDATE *** in data files)")
    print(f"\n  WING:")
    print(f"    a (dCL/da)    = {CLa:.4f} /rad = {CLa*np.pi/180:.5f} /deg")
    print(f"    Cm_0_ac       = {cm_0_ac:.4f} *")
    print(f"    CLmax (clean) = {cl_max:.3f}")
    print(f"    a_CLmax       = {alpha_cl_max:.1f} deg *")
    print(f"    CL_cruise     = {cl_cruise:.4f}")
    print(f"    a_CL_cruise   = {alpha_cl_cruise:.2f} deg *")
    print(f"    alpha_0L      = {alpha_0L:.1f} deg *")
    print(f"\n  HORIZONTAL TAIL (AR = {AR_htail:.2f}):")
    print(f"    a_t (dCL/da)  = {CLa_htail:.4f} /rad = {CLa_htail*np.pi/180:.5f} /deg")
    print(f"    CLmax_tail    = {CLmax_htail:.3f} *")
    print(f"    a_CLmax_tail  = {alpha_CLmax_htail:.1f} deg *")
    print(f"\n  VERTICAL TAIL (AR = {AR_vtail:.2f}):")
    print(f"    a_t (dCL/da)  = {CLa_vtail:.4f} /rad = {CLa_vtail*np.pi/180:.5f} /deg")
    print(f"    CLmax_tail    = {CLmax_vtail:.3f} *")
    print(f"    a_CLmax_tail  = {alpha_CLmax_vtail:.1f} deg *")
    print(f"\n  NOT COMPUTABLE (need CG position, tail arm, downwash):")
    print(f"    dCm/da (wing, htail, vtail)")
    print(f"    CL_cruise_tail, a_cruise_tail (need trim analysis)")

    # =========================================================================
    #  SUMMARY
    # =========================================================================
    print(f"\n{'='*65}")
    print(f"  {name} — VALUES FOR MISSION SIZING")
    print(f"{'='*65}")
    print(f"  CD0          = {cd0:.5f}")
    print(f"  e (Oswald)   = {e:.4f}  (conservative)")
    print(f"  e (LE suct.) = {e_from_les:.4f}  (recommended)")
    print(f"  K (Oswald)   = {K:.5f}")
    print(f"  K (LE suct.) = {K_les:.5f}  (recommended)")
    print(f"  (L/D)_max    = {ld_max:.2f}  (LE suction K)")
    print(f"  CL_cruise    = {cl_cruise:.4f}")
    print(f"  CL_max clean = {cl_max:.3f}")
    print(f"  CL_max T/O   = {cl_max_to:.3f}")
    print(f"  CL_max Land  = {cl_max_land:.3f}")
    print(f"  M_DD         = {m_dd:.3f}")
    print(f"  CLa          = {CLa:.4f} /rad")
    print(f"{'='*65}")


# =============================================================================
#  MAIN — load aircraft data files and run analysis
# =============================================================================
if __name__ == "__main__":
    from data.ZRJ70 import AIRCRAFT as ZRJ70
    from data.ZRJ100 import AIRCRAFT as ZRJ100

    print("=" * 65)
    print("  REGIONAL JET AERODYNAMIC ANALYSIS (Raymer Ch.12 Methods)")
    print("=" * 65)

    for ac in [ZRJ70, ZRJ100]:
        analyse(ac)
