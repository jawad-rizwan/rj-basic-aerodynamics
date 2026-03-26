"""
ZRJ70 — Cruise L/D and (L/D)_max vs Altitude.

Sweeps altitude from 30,000 to 50,000 ft in 500-ft increments,
recomputing CD0 (via component buildup) and the drag polar at each
altitude.  Plots cruise L/D (solid blue) and (L/D)_max (dashed orange).

Usage (from project root):
    python3 examples/plot_ld_vs_alt.py
"""

import sys
import os

# Ensure project root is on the path so aero/ and data/ are importable.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from data.ZRJ70 import AIRCRAFT

from aero.atmosphere import isa_properties, dynamic_pressure
from aero.form_factors import ff_wing, ff_fuselage, ff_nacelle, ff_tail_with_hinge
from aero.parasite_drag import cd0_component_buildup
from aero.lift import (
    cl_alpha_subsonic, fuselage_lift_factor, effective_aspect_ratio_winglet,
)
from aero.induced_drag import oswald_e, k_factor_leading_edge_suction
from aero.drag_polar import DragPolar

# ── Unpack aircraft data ─────────────────────────────────────────────────
ac = AIRCRAFT

S_ref        = ac["S_ref"]
AR           = ac["AR"]
b            = ac["b"]
MAC          = ac["MAC"]
t_c_wing     = ac["t_c_wing"]
x_c_max_wing = ac["x_c_max_wing"]
sweep_qc_deg = ac["sweep_qc_deg"]
sweep_le_deg = ac["sweep_le_deg"]
sweep_mt_rad = np.radians(ac["sweep_mt_deg"])
S_exposed    = ac["S_exposed"]
winglet_h    = ac["winglet_h"]

fuse_d       = ac["fuse_d"]
fuse_length  = ac["fuse_length"]
MTOW         = ac["MTOW"]
W_cruise     = 0.92 * MTOW

S_htail_exposed = ac["S_htail_exposed"]
MAC_htail       = ac["MAC_htail"]
t_c_htail       = ac["t_c_htail"]
x_c_max_htail   = ac["x_c_max_htail"]
sweep_mt_htail_rad = np.radians(ac["sweep_mt_htail_deg"])

S_vtail_exposed = ac["S_vtail_exposed"]
MAC_vtail       = ac["MAC_vtail"]
t_c_vtail       = ac["t_c_vtail"]
x_c_max_vtail   = ac["x_c_max_vtail"]
sweep_mt_vtail_rad = np.radians(ac["sweep_mt_vtail_deg"])

nacelle_length = ac["nacelle_length"]
nacelle_d      = ac["nacelle_d"]
n_nacelles     = ac["n_nacelles"]

cruise_mach  = ac["cruise_mach"]

k_composite  = ac["k_composite"]
k_metal      = ac["k_metal"]

Q_wing   = ac["Q_wing"]
Q_fuse   = ac["Q_fuse"]
Q_htail  = ac["Q_htail"]
Q_vtail  = ac["Q_vtail"]
Q_nac    = ac["Q_nac"]

lam_wing = ac["lam_wing"]
lam_fuse = ac["lam_fuse"]
lam_tail = ac["lam_tail"]
lam_nac  = ac["lam_nac"]
leak_pct = ac["leak_pct"]
S_suction = ac["S_suction"]

# ── Geometry that does not change with altitude ──────────────────────────
AR_eff = effective_aspect_ratio_winglet(AR, winglet_h, b)
F_lift = fuselage_lift_factor(fuse_d, b)
sf_product = min((S_exposed / S_ref) * F_lift, 0.98)

# Wetted areas (altitude-independent)
S_wet_wing  = 2.0 * S_exposed * (1.0 + 0.25 * t_c_wing)
S_wet_htail = 2.0 * S_htail_exposed * (1.0 + 0.25 * t_c_htail)
S_wet_vtail = 2.0 * S_vtail_exposed * (1.0 + 0.25 * t_c_vtail)
S_wet_nac   = np.pi * nacelle_d * nacelle_length * 0.8

fuse_f     = fuse_length / fuse_d
S_wet_fuse = (np.pi * fuse_d * fuse_length
              * (1.0 - 2.0 / fuse_f) ** (2.0 / 3.0)
              * (1.0 + 1.0 / fuse_f ** 2))

# Form factors (altitude-independent at fixed Mach)
FF_wing  = ff_wing(x_c_max_wing, t_c_wing, cruise_mach, sweep_mt_rad)
FF_fuse  = ff_fuselage(fuse_length, fuse_d)
FF_htail = ff_tail_with_hinge(x_c_max_htail, t_c_htail, cruise_mach, sweep_mt_htail_rad)
FF_vtail = ff_tail_with_hinge(x_c_max_vtail, t_c_vtail, cruise_mach, sweep_mt_vtail_rad)
FF_nac   = ff_nacelle(nacelle_length, nacelle_d)


# ── Altitude sweep ───────────────────────────────────────────────────────
altitudes = np.arange(30000, 50500, 500)  # 30,000 to 50,000 in 500-ft steps
ld_cruise_arr = np.empty_like(altitudes, dtype=float)
ld_max_arr    = np.empty_like(altitudes, dtype=float)

for idx, alt in enumerate(altitudes):
    # Lift-curve slope varies with Mach (fixed here) but Re changes via atm;
    # CLa itself only depends on Mach/geometry, so compute once outside loop
    # would be fine. But we keep it inside for clarity.
    CLa = cl_alpha_subsonic(
        AR_eff, cruise_mach, sweep_mt_rad,
        eta=0.95, s_exposed_ratio=sf_product, F=1.0,
    )

    # Induced-drag K via LE suction method
    K_les = k_factor_leading_edge_suction(AR_eff, CLa, S_suction)

    # Component buildup at this altitude
    components = [
        {"name": "Wing",      "s_wet": S_wet_wing,  "length": MAC,
         "ff": FF_wing,  "Q": Q_wing,  "pct_laminar": lam_wing, "k": k_composite},
        {"name": "Fuselage",  "s_wet": S_wet_fuse,  "length": fuse_length,
         "ff": FF_fuse,  "Q": Q_fuse,  "pct_laminar": lam_fuse, "k": k_metal},
        {"name": "H-tail",    "s_wet": S_wet_htail, "length": MAC_htail,
         "ff": FF_htail, "Q": Q_htail, "pct_laminar": lam_tail, "k": k_composite},
        {"name": "V-tail",    "s_wet": S_wet_vtail, "length": MAC_vtail,
         "ff": FF_vtail, "Q": Q_vtail, "pct_laminar": lam_tail, "k": k_composite},
    ]
    for i in range(n_nacelles):
        side = "L" if i == 0 else "R"
        components.append(
            {"name": f"Nacelle {side}", "s_wet": S_wet_nac, "length": nacelle_length,
             "ff": FF_nac, "Q": Q_nac, "pct_laminar": lam_nac, "k": k_composite})

    result = cd0_component_buildup(
        components, S_ref, cruise_mach, float(alt),
        cd_misc=0.0, leak_pct=leak_pct,
    )
    cd0 = result["cd0_total"]

    # Drag polar at this altitude
    polar = DragPolar(cd0, K_les)

    # (L/D)_max
    ld_max_val, _ = polar.ld_max()
    ld_max_arr[idx] = ld_max_val

    # Cruise L/D (CL set by weight / q / S_ref)
    q = dynamic_pressure(cruise_mach, float(alt))
    cl_cruise = polar.cl_for_cruise(W_cruise, q, S_ref)
    ld_cruise_arr[idx] = polar.ld(cl_cruise)


# ── Find optimum altitude ────────────────────────────────────────────────
opt_idx     = np.argmax(ld_cruise_arr)
opt_alt     = altitudes[opt_idx]
opt_ld      = ld_cruise_arr[opt_idx]

baseline_alt = 35000
baseline_idx = np.argmin(np.abs(altitudes - baseline_alt))
baseline_ld  = ld_cruise_arr[baseline_idx]


# ── Plot ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 6))

ax.plot(altitudes / 1000, ld_cruise_arr,
        color="#2563eb", linewidth=2.2, label="Cruise $L/D$")
ax.plot(altitudes / 1000, ld_max_arr,
        color="#ea580c", linewidth=2.0, linestyle="--", label="$(L/D)_{max}$")

# Mark optimal altitude
ax.plot(opt_alt / 1000, opt_ld, "o",
        color="#2563eb", markersize=10, zorder=5)
ax.annotate(
    f"Optimum: {opt_alt/1000:.1f} kft\n$L/D$ = {opt_ld:.2f}",
    xy=(opt_alt / 1000, opt_ld),
    xytext=(30, -30), textcoords="offset points", fontsize=10,
    color="#2563eb",
    arrowprops=dict(arrowstyle="->", color="#2563eb", lw=1.3),
)

# Mark 35 kft baseline
ax.plot(baseline_alt / 1000, baseline_ld, "s",
        color="#16a34a", markersize=9, zorder=5)
ax.annotate(
    f"Baseline: {baseline_alt/1000:.0f} kft\n$L/D$ = {baseline_ld:.2f}",
    xy=(baseline_alt / 1000, baseline_ld),
    xytext=(-30, 30), textcoords="offset points", fontsize=10,
    color="#16a34a",
    arrowprops=dict(arrowstyle="->", color="#16a34a", lw=1.3),
)

ax.set_xlabel("Altitude (kft)", fontsize=12)
ax.set_ylabel("$L/D$", fontsize=12)
ax.set_title(r"ZRJ70 — Cruise $L/D$ vs Altitude ($M$ = 0.78)",
             fontsize=13, fontweight="bold")
ax.legend(fontsize=11, loc="best")
ax.grid(True, alpha=0.4)

fig.tight_layout()

out_path = os.path.join(os.path.dirname(__file__), "charts", "ld_vs_altitude_ZRJ70.png")
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Saved: {out_path}")
