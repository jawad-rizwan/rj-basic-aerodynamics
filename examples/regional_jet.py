"""
Regional Jet Aerodynamic Analysis — CRJ700-class baseline
Using Raymer Chapter 12 methods.

This script computes CD0, Oswald e, K, CLmax, L/D, and MDD
for a CRJ700-class regional jet using the component buildup method.

Outputs are intended to feed into Chapter 6 refined sizing
(rj-mission-sizing).
"""

import sys
import os
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


# =============================================================================
#  AIRCRAFT GEOMETRY — CRJ700-class regional jet (PLACEHOLDER DATA)
#  Replace every value marked *** UPDATE *** with your actual design geometry.
# =============================================================================

# --- Wing ---
S_ref = 520.0         # ft^2, wing reference area (trapezoidal)       *** UPDATE ***
AR = 7.8              # aspect ratio                                  *** UPDATE ***
b = np.sqrt(AR * S_ref)  # span (ft)
MAC = S_ref / b       # mean aerodynamic chord approximation (ft)
taper = 0.30          # taper ratio                                   *** UPDATE ***
t_c_wing = 0.11       # thickness-to-chord ratio                      *** UPDATE ***
x_c_max_wing = 0.37   # chordwise location of max thickness           *** UPDATE ***
sweep_qc_deg = 26.0   # quarter-chord sweep (deg)                     *** UPDATE ***
sweep_qc_rad = np.radians(sweep_qc_deg)
sweep_le_deg = 29.0   # leading-edge sweep (deg)                      *** UPDATE ***
sweep_le_rad = np.radians(sweep_le_deg)
# Approximate max-thickness sweep from quarter-chord sweep
sweep_mt_deg = 23.0   # max-thickness line sweep (deg)                *** UPDATE ***
sweep_mt_rad = np.radians(sweep_mt_deg)
S_exposed = 460.0     # ft^2, exposed planform (minus fuselage cover) *** UPDATE ***
winglet_h = 3.5       # ft, winglet height                            *** UPDATE ***

# --- Fuselage ---
fuse_length = 106.0   # ft, total fuselage length                     *** UPDATE ***
fuse_d = 8.83         # ft, equivalent diameter                       *** UPDATE ***
fuse_Amax = np.pi / 4 * fuse_d**2  # ft^2, max cross-section

# --- Horizontal tail ---
S_htail = 120.0       # ft^2, htail reference area                    *** UPDATE ***
S_htail_exposed = 108.0  # ft^2                                       *** UPDATE ***
MAC_htail = 5.5       # ft                                            *** UPDATE ***
t_c_htail = 0.09      #                                               *** UPDATE ***
x_c_max_htail = 0.35  #                                               *** UPDATE ***
sweep_mt_htail_rad = np.radians(28.0)  # deg                          *** UPDATE ***

# --- Vertical tail ---
S_vtail = 100.0       # ft^2                                          *** UPDATE ***
S_vtail_exposed = 95.0  # ft^2                                        *** UPDATE ***
MAC_vtail = 7.0       # ft                                            *** UPDATE ***
t_c_vtail = 0.09      #                                               *** UPDATE ***
x_c_max_vtail = 0.35  #                                               *** UPDATE ***
sweep_mt_vtail_rad = np.radians(35.0)  # deg                          *** UPDATE ***

# --- Nacelles (x2, CF34-8C5 class) ---
nacelle_length = 10.5  # ft, per nacelle                              *** UPDATE ***
nacelle_d = 4.2        # ft, diameter                                 *** UPDATE ***
n_nacelles = 2         #                                              *** UPDATE ***

# --- Flight condition ---
cruise_mach = 0.78     #                                              *** UPDATE ***
cruise_alt = 41000.0   # ft                                           *** UPDATE ***
atm = isa_properties(cruise_alt)
V_cruise = cruise_mach * atm["a"]  # ft/s

# --- Airfoil data ---
cl_max_airfoil = 1.6   # 2D airfoil CLmax (supercritical, ~Re 10M)   *** UPDATE ***

# --- Surface finish ---
k_metal = 1.33e-5      # ft, production sheet metal (Raymer Table 12.5)  *** UPDATE ***
k_composite = 0.50e-5  # ft, polished/smooth composite                  *** UPDATE ***

# =============================================================================
#  1. LIFT-CURVE SLOPE
# =============================================================================
print("=" * 65)
print("  REGIONAL JET AERODYNAMIC ANALYSIS (Raymer Ch.12 Methods)")
print("=" * 65)

# Effective aspect ratio with winglets — Raymer Eq. 12.11
AR_eff = effective_aspect_ratio_winglet(AR, winglet_h, b)

# Fuselage lift factor — Raymer Eq. 12.9
F = fuselage_lift_factor(fuse_d, b)
s_exp_ratio = S_exposed / S_ref
# Clamp (S_exp/S_ref)*F to max 0.98 per Raymer recommendation
sf_product = min(s_exp_ratio * F, 0.98)

# Subsonic lift-curve slope — Raymer Eq. 12.6
# Use clamped (S_exp/S_ref)*F product per Raymer recommendation
CLa = cl_alpha_subsonic(
    A=AR_eff,
    mach=cruise_mach,
    sweep_max_t=sweep_mt_rad,
    eta=0.95,  # Raymer Eq. 12.8
    s_exposed_ratio=sf_product,  # clamped (S_exp/S_ref)*F
    F=1.0,  # already folded into sf_product
)

print(f"\n--- LIFT ---")
print(f"  Wing span              = {b:.1f} ft")
print(f"  MAC                    = {MAC:.2f} ft")
print(f"  AR (geometric)         = {AR:.1f}")
print(f"  AR (effective, winglet) = {AR_eff:.2f}")
print(f"  Fuselage lift factor F = {F:.3f}")
print(f"  (S_exp/S_ref)*F        = {s_exp_ratio * F:.3f}")
print(f"  CLa (cruise, M={cruise_mach})   = {CLa:.4f} /rad = {CLa * np.pi/180:.5f} /deg")

# =============================================================================
#  2. CD0 — COMPONENT BUILDUP METHOD (Raymer Eq. 12.24)
# =============================================================================

# Wetted areas (approximate for CRJ700-class)
# Wing: S_wet ≈ 2 * S_exposed * (1 + 0.25 * t/c) for finite thickness
S_wet_wing = 2.0 * S_exposed * (1.0 + 0.25 * t_c_wing)

# Fuselage: S_wet ≈ pi * d * L * (1 - 2/f)^(2/3) * (1 + 1/f^2), simplified
fuse_f = fuse_length / fuse_d
S_wet_fuse = np.pi * fuse_d * fuse_length * (1.0 - 2.0/fuse_f)**(2.0/3.0) * (1.0 + 1.0/fuse_f**2)

# Horizontal tail
S_wet_htail = 2.0 * S_htail_exposed * (1.0 + 0.25 * t_c_htail)

# Vertical tail
S_wet_vtail = 2.0 * S_vtail_exposed * (1.0 + 0.25 * t_c_vtail)

# Nacelles (each)
S_wet_nac = np.pi * nacelle_d * nacelle_length * 0.8  # ~80% due to pylon/cutout

S_wet_total = S_wet_wing + S_wet_fuse + S_wet_htail + S_wet_vtail + n_nacelles * S_wet_nac

print(f"\n--- WETTED AREAS ---")
print(f"  Wing          = {S_wet_wing:.0f} ft^2")
print(f"  Fuselage      = {S_wet_fuse:.0f} ft^2")
print(f"  H-tail        = {S_wet_htail:.0f} ft^2")
print(f"  V-tail        = {S_wet_vtail:.0f} ft^2")
print(f"  Nacelles (x2) = {n_nacelles * S_wet_nac:.0f} ft^2")
print(f"  TOTAL         = {S_wet_total:.0f} ft^2")
print(f"  S_wet / S_ref = {S_wet_total / S_ref:.2f}")

# Form factors
FF_wing = ff_wing(x_c_max_wing, t_c_wing, cruise_mach, sweep_mt_rad)
FF_fuse = ff_fuselage(fuse_length, fuse_d)
FF_htail = ff_tail_with_hinge(x_c_max_htail, t_c_htail, cruise_mach, sweep_mt_htail_rad)
FF_vtail = ff_tail_with_hinge(x_c_max_vtail, t_c_vtail, cruise_mach, sweep_mt_vtail_rad)
FF_nac = ff_nacelle(nacelle_length, nacelle_d)

print(f"\n--- FORM FACTORS ---")
print(f"  Wing (Eq. 12.30)     = {FF_wing:.4f}")
print(f"  Fuselage (Eq. 12.31) = {FF_fuse:.4f}")
print(f"  H-tail (Eq. 12.30+)  = {FF_htail:.4f}")
print(f"  V-tail (Eq. 12.30+)  = {FF_vtail:.4f}")
print(f"  Nacelle (Eq. 12.32)  = {FF_nac:.4f}")

# Interference factors Q — Raymer Table 12.6
Q_wing = 1.05    # high/mid wing, well-filleted                      *** UPDATE ***
Q_fuse = 1.00    # fuselage                                          *** UPDATE ***
Q_htail = 1.05   # conventional tail group                           *** UPDATE ***
Q_vtail = 1.05   # conventional tail group                           *** UPDATE ***
Q_nac = 1.30     # nacelle ~1 diameter away from fuselage            *** UPDATE ***

# Laminar flow percentages — Raymer Table 12.4 (civil jet, classic metal)
lam_wing = 0.10  #                                                   *** UPDATE ***
lam_fuse = 0.05  #                                                   *** UPDATE ***
lam_tail = 0.05  #                                                   *** UPDATE ***
lam_nac = 0.0    #                                                   *** UPDATE ***

# Build component list for Eq. 12.24
components = [
    {"name": "Wing",       "s_wet": S_wet_wing,         "length": MAC,
     "ff": FF_wing,  "Q": Q_wing,  "pct_laminar": lam_wing, "k": k_metal},
    {"name": "Fuselage",   "s_wet": S_wet_fuse,         "length": fuse_length,
     "ff": FF_fuse,  "Q": Q_fuse,  "pct_laminar": lam_fuse, "k": k_metal},
    {"name": "H-tail",     "s_wet": S_wet_htail,        "length": MAC_htail,
     "ff": FF_htail, "Q": Q_htail, "pct_laminar": lam_tail, "k": k_metal},
    {"name": "V-tail",     "s_wet": S_wet_vtail,        "length": MAC_vtail,
     "ff": FF_vtail, "Q": Q_vtail, "pct_laminar": lam_tail, "k": k_metal},
    {"name": "Nacelle L",  "s_wet": S_wet_nac,          "length": nacelle_length,
     "ff": FF_nac,   "Q": Q_nac,   "pct_laminar": lam_nac,  "k": k_metal},
    {"name": "Nacelle R",  "s_wet": S_wet_nac,          "length": nacelle_length,
     "ff": FF_nac,   "Q": Q_nac,   "pct_laminar": lam_nac,  "k": k_metal},
]

# Leakage & protuberance — Raymer Table 12.9: 2-5% for jet transports
leak_pct = 0.03       #                                               *** UPDATE ***

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

# Quick check vs equivalent Cfe method — Raymer Eq. 12.23
cd0_check = cd0_equivalent_cfe(S_wet_total, S_ref, "civil_transport")
print(f"\n  Cross-check: Cfe method (Eq. 12.23)    = {cd0_check:.6f}")

# =============================================================================
#  3. OSWALD EFFICIENCY & K FACTOR
# =============================================================================
e = oswald_e(AR_eff, sweep_le_deg)
K = k_factor(AR_eff, e)

# Leading-edge suction method — Raymer Eq. 12.57
S_suction = 0.90  # typical for well-designed wing at design CL      *** UPDATE ***
K_les = k_factor_leading_edge_suction(AR_eff, CLa, S_suction)

# NOTE: Raymer Eqs. 12.48/12.49 are known to be conservative.
# Raymer says "a good wing should see an e of no less than 0.9" at design CL.
# The LE suction method (Eq. 12.57) typically gives more realistic results.
e_from_les = 1.0 / (np.pi * AR_eff * K_les)

print(f"\n--- DRAG DUE TO LIFT ---")
print(f"  Oswald e (Eq. 12.48/49)    = {e:.4f}  (conservative)")
print(f"  K (Oswald, Eq. 12.47)      = {K:.5f}")
print(f"  K (LE suction, Eq. 12.57)  = {K_les:.5f}  (S = {S_suction})")
print(f"  e equiv (LE suction)       = {e_from_les:.4f}  (more realistic)")

# =============================================================================
#  4. DRAG POLAR & L/D
# =============================================================================
cd0 = result["cd0_total"]

# Build polar with both K methods
polar_oswald = DragPolar(cd0, K)
polar_les = DragPolar(cd0, K_les)  # more realistic for well-designed wing
polar = polar_les  # use LE suction as primary

ld_max, cl_ldmax = polar.ld_max()
ld_max_oswald, _ = polar_oswald.ld_max()

# Cruise CL estimate (need weight estimate)
W_cruise_approx = 65000.0  # lbs, mid-mission weight guess            *** UPDATE ***
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

# =============================================================================
#  5. MAXIMUM LIFT
# =============================================================================
cl_max = cl_max_clean(cl_max_airfoil, sweep_qc_rad)

# With trailing-edge flaps (slotted Fowler, c'/c ≈ 1.25)
# Raymer Table 12.2: Fowler ΔCl_max ≈ 1.3 * 1.25 = 1.625
delta_cl_te = 1.3 * 1.25   # Fowler flap                             *** UPDATE ***
s_flapped_ratio = 0.70     # ~70% of wing span has flaps              *** UPDATE ***
sweep_hl_rad = np.radians(20.0)  # flap hinge line sweep              *** UPDATE ***

# With leading-edge slat
# Raymer Table 12.2: slat ΔCl_max ≈ 0.4 * c'/c
delta_cl_le = 0.4 * 1.10   # slat with c'/c ≈ 1.10                   *** UPDATE ***
s_slat_ratio = 0.75        #                                          *** UPDATE ***

cl_max_to = cl_max_flaps(cl_max, delta_cl_te * 0.70, s_flapped_ratio, sweep_hl_rad,
                          delta_cl_max_le=0.0)  # no slats for takeoff typically

cl_max_land = cl_max_flaps(cl_max, delta_cl_te, s_flapped_ratio, sweep_hl_rad,
                            delta_cl_max_le=0.9 * delta_cl_le * s_slat_ratio * np.cos(sweep_hl_rad))

print(f"\n--- MAXIMUM LIFT ---")
print(f"  CL_max (clean, Eq. 12.15)  = {cl_max:.3f}")
print(f"  CL_max (takeoff flaps)     = {cl_max_to:.3f}")
print(f"  CL_max (landing, full)     = {cl_max_land:.3f}")

# =============================================================================
#  6. DRAG DIVERGENCE MACH
# =============================================================================
m_dd = mdd_wing(t_c_wing, sweep_qc_deg, cl_design=cl_cruise,
                 supercritical=True)

print(f"\n--- TRANSONIC ---")
print(f"  M_DD (Boeing, Eq. 12.46)   = {m_dd:.3f}")
print(f"  M_cr (approx)              = {m_dd - 0.08:.3f}")

# =============================================================================
#  SUMMARY — values to feed into rj-mission-sizing
# =============================================================================
print(f"\n{'='*65}")
print(f"  VALUES FOR MISSION SIZING (replace *** UPDATE *** placeholders)")
print(f"{'='*65}")
print(f"  CD0          = {cd0:.5f}")
print(f"  e (Oswald)   = {e:.4f}  (conservative)")
print(f"  e (LE suct.) = {e_from_les:.4f}  (recommended)")
print(f"  K (Oswald)   = {K:.5f}")
print(f"  K (LE suct.) = {K_les:.5f}  (recommended)")
print(f"  (L/D)_max    = {ld_max:.2f}  (LE suction K)")
print(f"  CL_max clean = {cl_max:.3f}")
print(f"  CL_max T/O   = {cl_max_to:.3f}")
print(f"  CL_max Land  = {cl_max_land:.3f}")
print(f"  M_DD         = {m_dd:.3f}")
print(f"  CLa          = {CLa:.4f} /rad")
print(f"{'='*65}")
