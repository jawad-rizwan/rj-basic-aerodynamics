"""
Aerodynamic charts for the regional jet analysis.
Generates 6 figures from the computed data in regional_jet.py.

Usage:
    python3 examples/plot_aero.py          # show interactive plots
    python3 examples/plot_aero.py --save   # save PNGs to examples/charts/
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from aero.atmosphere import isa_properties, dynamic_pressure
from aero.form_factors import ff_wing, ff_fuselage, ff_nacelle, ff_tail_with_hinge
from aero.parasite_drag import cd0_component_buildup
from aero.lift import cl_alpha_subsonic, fuselage_lift_factor, cl_max_clean, cl_max_flaps, effective_aspect_ratio_winglet
from aero.induced_drag import oswald_e, k_factor, k_factor_leading_edge_suction
from aero.wave_drag import mdd_wing, drag_rise_transonic
from aero.drag_polar import DragPolar

from data.ZRJ70 import AIRCRAFT as ZRJ70
from data.ZRJ100 import AIRCRAFT as ZRJ100

if "--save" in sys.argv:
    matplotlib.use("Agg")


# ── Style ────────────────────────────────────────────────────────────────
plt.style.use("seaborn-v0_8-whitegrid")
COLORS = {"ZRJ70": "#2563eb", "ZRJ100": "#dc2626"}  # blue / red

VARIANTS = {"ZRJ70": ZRJ70, "ZRJ100": ZRJ100}


# ── Per-variant computations ─────────────────────────────────────────────
variant_data = {}
for vname, ac in VARIANTS.items():
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
    W_cruise     = 0.92 * MTOW
    cruise_mach  = ac["cruise_mach"]
    cruise_alt   = ac["cruise_alt"]
    k_composite  = ac["k_composite"]
    k_metal      = ac["k_metal"]
    leak_pct     = ac["leak_pct"]
    cl_max_airfoil = ac["cl_max_airfoil"]
    S_suction    = ac["S_suction"]

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

    Q_wing   = ac["Q_wing"]
    Q_fuse   = ac["Q_fuse"]
    Q_htail  = ac["Q_htail"]
    Q_vtail  = ac["Q_vtail"]
    Q_nac    = ac["Q_nac"]
    lam_wing = ac["lam_wing"]
    lam_fuse = ac["lam_fuse"]
    lam_tail = ac["lam_tail"]
    lam_nac  = ac["lam_nac"]

    delta_cl_te = ac["delta_cl_te_factor"]
    te_to_frac  = ac["te_flap_takeoff_fraction"]
    s_flapped   = ac["s_flapped_ratio"]
    sweep_hl    = np.radians(ac["sweep_hl_deg"])
    delta_cl_le = ac["delta_cl_le_factor"]
    s_slat      = ac["s_slat_ratio"]

    AR_eff = effective_aspect_ratio_winglet(AR, winglet_h, b)
    F_lift = fuselage_lift_factor(fuse_d, b)
    sf_product = min((S_exposed / S_ref) * F_lift, 0.98)

    CLa = cl_alpha_subsonic(AR_eff, cruise_mach, sweep_mt_rad, eta=0.95,
                            s_exposed_ratio=sf_product, F=1.0)

    S_wet_wing  = 2.0 * S_exposed * (1.0 + 0.25 * t_c_wing)
    S_wet_htail = 2.0 * S_htail_exposed * (1.0 + 0.25 * t_c_htail)
    S_wet_vtail = 2.0 * S_vtail_exposed * (1.0 + 0.25 * t_c_vtail)
    S_wet_nac   = np.pi * nacelle_d * nacelle_length * 0.8

    fuse_f     = fuse_length / fuse_d
    S_wet_fuse = (np.pi * fuse_d * fuse_length
                  * (1.0 - 2.0/fuse_f)**(2.0/3.0)
                  * (1.0 + 1.0/fuse_f**2))

    FF_wing  = ff_wing(x_c_max_wing, t_c_wing, cruise_mach, sweep_mt_rad)
    FF_fuse  = ff_fuselage(fuse_length, fuse_d)
    FF_htail = ff_tail_with_hinge(x_c_max_htail, t_c_htail, cruise_mach, sweep_mt_htail_rad)
    FF_vtail = ff_tail_with_hinge(x_c_max_vtail, t_c_vtail, cruise_mach, sweep_mt_vtail_rad)
    FF_nac   = ff_nacelle(nacelle_length, nacelle_d)

    e_oswald = oswald_e(AR_eff, sweep_le_deg)
    K_oswald = k_factor(AR_eff, e_oswald)
    K_les    = k_factor_leading_edge_suction(AR_eff, CLa, S_suction)
    e_les    = 1.0 / (np.pi * AR_eff * K_les)

    cl_max_val = cl_max_clean(cl_max_airfoil, sweep_qc_rad)
    cl_max_to  = cl_max_flaps(cl_max_val, delta_cl_te * te_to_frac, s_flapped, sweep_hl)
    cl_max_land = cl_max_flaps(cl_max_val, delta_cl_te, s_flapped, sweep_hl,
                               delta_cl_max_le=0.9 * delta_cl_le * s_slat * np.cos(sweep_hl))

    q_cruise = dynamic_pressure(cruise_mach, cruise_alt)

    components = [
        {"name": "Wing",      "s_wet": S_wet_wing,  "length": MAC,
         "ff": FF_wing, "Q": Q_wing, "pct_laminar": lam_wing, "k": k_composite},
        {"name": "Fuselage",  "s_wet": S_wet_fuse,   "length": fuse_length,
         "ff": FF_fuse, "Q": Q_fuse, "pct_laminar": lam_fuse, "k": k_metal},
        {"name": "H-tail",    "s_wet": S_wet_htail,  "length": MAC_htail,
         "ff": FF_htail,"Q": Q_htail,"pct_laminar": lam_tail, "k": k_composite},
        {"name": "V-tail",    "s_wet": S_wet_vtail,  "length": MAC_vtail,
         "ff": FF_vtail,"Q": Q_vtail,"pct_laminar": lam_tail, "k": k_composite},
        {"name": "Nacelle L", "s_wet": S_wet_nac,    "length": nacelle_length,
         "ff": FF_nac,  "Q": Q_nac,  "pct_laminar": lam_nac,  "k": k_composite},
        {"name": "Nacelle R", "s_wet": S_wet_nac,    "length": nacelle_length,
         "ff": FF_nac,  "Q": Q_nac,  "pct_laminar": lam_nac,  "k": k_composite},
    ]

    result = cd0_component_buildup(components, S_ref, cruise_mach, cruise_alt,
                                    cd_misc=0.0, leak_pct=leak_pct)

    cd0 = result["cd0_total"]
    polar = DragPolar(cd0, K_les)
    ld_max, cl_ldmax = polar.ld_max()
    cl_cruise = polar.cl_for_cruise(W_cruise, q_cruise, S_ref)
    m_dd = mdd_wing(t_c_wing, ac["sweep_qc_deg"], cl_design=cl_cruise, supercritical=True)

    # Component CD0 breakdown (merge nacelles)
    breakdown = result["breakdown"]
    comp_cd0 = {}
    for c in breakdown:
        key = c["name"]
        if key.startswith("Nacelle"):
            comp_cd0["Nacelles"] = comp_cd0.get("Nacelles", 0) + c["CD0"]
        else:
            comp_cd0[key] = c["CD0"]
    comp_cd0["Leak + Misc"] = result["cd_leak"] + result["cd_misc"]

    S_wet_total = S_wet_wing + S_wet_fuse + S_wet_htail + S_wet_vtail + n_nacelles * S_wet_nac

    variant_data[vname] = {
        "cd0": cd0, "polar": polar, "K": K_les, "e": e_les,
        "ld_max": ld_max, "cl_ldmax": cl_ldmax,
        "cl_cruise": cl_cruise, "W_cruise": W_cruise,
        "m_dd": m_dd, "comp_cd0": comp_cd0,
        "S_wet_fuse": S_wet_fuse, "S_wet_total": S_wet_total,
        "S_wet_wing": S_wet_wing, "S_wet_htail": S_wet_htail,
        "S_wet_vtail": S_wet_vtail, "S_wet_nac": S_wet_nac,
        "n_nacelles": n_nacelles,
        "cruise_mach": cruise_mach, "cruise_alt": cruise_alt,
        "cl_max": cl_max_val, "cl_max_to": cl_max_to, "cl_max_land": cl_max_land,
    }


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 1 — Drag Polar  (CD vs CL)
# ═══════════════════════════════════════════════════════════════════════
fig1, ax1 = plt.subplots(figsize=(8, 6))
CL = np.linspace(0, 1.4, 300)

for vname, col in [("ZRJ70", COLORS["ZRJ70"]), ("ZRJ100", COLORS["ZRJ100"])]:
    vd = variant_data[vname]
    CD = vd["polar"].cd(CL)
    ax1.plot(CD, CL, color=col, linewidth=2, label=vname)
    # Cruise point
    ax1.plot(vd["polar"].cd(vd["cl_cruise"]), vd["cl_cruise"],
             'o', color=col, markersize=8, zorder=5)
    ax1.annotate(f'Cruise\nCL={vd["cl_cruise"]:.3f}',
                 xy=(vd["polar"].cd(vd["cl_cruise"]), vd["cl_cruise"]),
                 xytext=(25, -10), textcoords='offset points', fontsize=8,
                 color=col, arrowprops=dict(arrowstyle='->', color=col, lw=1))
    # L/D max point
    ax1.plot(vd["polar"].cd(vd["cl_ldmax"]), vd["cl_ldmax"],
             's', color=col, markersize=7, zorder=5)

ax1.set_xlabel("$C_D$", fontsize=12)
ax1.set_ylabel("$C_L$", fontsize=12)
ax1.set_title("Drag Polar — ZRJ Regional Jet (Raymer Ch. 12)", fontsize=13, fontweight='bold')
ax1.legend(fontsize=10)
ax1.set_xlim(0, 0.14)
ax1.set_ylim(0, 1.4)
ax1.xaxis.set_minor_locator(MultipleLocator(0.005))
fig1.tight_layout()


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 2 — L/D vs CL
# ═══════════════════════════════════════════════════════════════════════
fig2, ax2 = plt.subplots(figsize=(8, 6))
CL_ld = np.linspace(0.05, 1.4, 300)

ldmax_offsets  = {"ZRJ70": (20, 15),  "ZRJ100": (20, -25)}
cruise_offsets = {"ZRJ70": (-80, -30), "ZRJ100": (-80, 15)}

for vname, col in [("ZRJ70", COLORS["ZRJ70"]), ("ZRJ100", COLORS["ZRJ100"])]:
    vd = variant_data[vname]
    LD = vd["polar"].ld(CL_ld)
    ax2.plot(CL_ld, LD, color=col, linewidth=2, label=vname)
    # Mark L/D max
    ax2.plot(vd["cl_ldmax"], vd["ld_max"], 's', color=col, markersize=8, zorder=5)
    ax2.annotate(f'(L/D)max = {vd["ld_max"]:.1f}',
                 xy=(vd["cl_ldmax"], vd["ld_max"]),
                 xytext=ldmax_offsets[vname], textcoords='offset points', fontsize=9,
                 color=col, arrowprops=dict(arrowstyle='->', color=col, lw=1))
    # Mark cruise
    ld_cr = vd["polar"].ld(vd["cl_cruise"])
    ax2.plot(vd["cl_cruise"], ld_cr, 'o', color=col, markersize=8, zorder=5)
    ax2.annotate(f'Cruise L/D = {ld_cr:.1f}',
                 xy=(vd["cl_cruise"], ld_cr),
                 xytext=cruise_offsets[vname], textcoords='offset points', fontsize=9,
                 color=col, arrowprops=dict(arrowstyle='->', color=col, lw=1))

ax2.set_xlabel("$C_L$", fontsize=12)
ax2.set_ylabel("$L/D$", fontsize=12)
ax2.set_title("Lift-to-Drag Ratio vs $C_L$", fontsize=13, fontweight='bold')
ax2.legend(fontsize=10)
ax2.set_xlim(0, 1.4)
ax2.set_ylim(0, 18)
fig2.tight_layout()


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 3 — CD0 Component Breakdown  (stacked bar)
# ═══════════════════════════════════════════════════════════════════════
fig3, ax3 = plt.subplots(figsize=(8, 5))
comp_names = ["Wing", "Fuselage", "H-tail", "V-tail", "Nacelles", "Leak + Misc"]
comp_colors = ["#3b82f6", "#ef4444", "#22c55e", "#f59e0b", "#8b5cf6", "#94a3b8"]
x_pos = np.arange(len(VARIANTS))
bar_w = 0.5

for i, vname in enumerate(VARIANTS):
    bottom = 0
    vd = variant_data[vname]
    for j, cn in enumerate(comp_names):
        val = vd["comp_cd0"].get(cn, 0) * 10000  # counts
        ax3.bar(x_pos[i], val, bar_w, bottom=bottom, color=comp_colors[j],
                label=cn if i == 0 else None, edgecolor='white', linewidth=0.5)
        if val > 4:  # only label bars wide enough
            ax3.text(x_pos[i], bottom + val/2, f'{val:.1f}',
                     ha='center', va='center', fontsize=8, color='white', fontweight='bold')
        bottom += val

ax3.set_xticks(x_pos)
ax3.set_xticklabels([f"{n}\nCD0 = {variant_data[n]['cd0']:.5f}" for n in VARIANTS], fontsize=10)
ax3.set_ylabel("$C_{D0}$ (drag counts)", fontsize=12)
ax3.set_title("Parasite Drag Breakdown — Component Buildup", fontsize=13, fontweight='bold')
ax3.legend(loc='upper right', fontsize=9)
fig3.tight_layout()


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 4 — CLmax Comparison  (grouped bar)
# ═══════════════════════════════════════════════════════════════════════
fig4, ax4 = plt.subplots(figsize=(7, 5))
# Use ZRJ70 values (shared wing, so identical for both)
vd0 = variant_data["ZRJ70"]
configs = ["Clean", "Takeoff\n(dbl-slot flaps)", "Landing\n(dbl-slot + slats)"]
clmax_vals = [vd0["cl_max"], vd0["cl_max_to"], vd0["cl_max_land"]]
bar_colors = ["#64748b", "#f97316", "#dc2626"]

bars = ax4.bar(configs, clmax_vals, width=0.5, color=bar_colors, edgecolor='white', linewidth=1.5)
for bar, val in zip(bars, clmax_vals):
    ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.04,
             f'{val:.3f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

ax4.set_ylabel("$C_{L,max}$", fontsize=12)
ax4.set_title("Maximum Lift Coefficients", fontsize=13, fontweight='bold')
ax4.set_ylim(0, 3.2)
fig4.tight_layout()


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 5 — Transonic Drag Rise  (CD0 vs Mach)
# ═══════════════════════════════════════════════════════════════════════
fig5, ax5 = plt.subplots(figsize=(8, 5))
M_arr = np.linspace(0.5, 1.05, 400)

for vname, col in [("ZRJ70", COLORS["ZRJ70"]), ("ZRJ100", COLORS["ZRJ100"])]:
    vd = variant_data[vname]
    cd_wave_105 = 0.012  # representative wave drag at M=1.05
    cd0_vs_M = drag_rise_transonic(vd["cd0"], vd["m_dd"], cd_wave_105, M_arr)
    ax5.plot(M_arr, cd0_vs_M * 10000, color=col, linewidth=2, label=vname)
    # MDD line
    ax5.axvline(vd["m_dd"], color=col, linestyle='--', linewidth=1, alpha=0.7)
    ax5.text(vd["m_dd"] + 0.005, vd["cd0"] * 10000 + 25,
             f'$M_{{DD}}$ = {vd["m_dd"]:.3f}', fontsize=9, color=col, rotation=90, va='bottom')

# Mark cruise Mach
cruise_mach = ZRJ70["cruise_mach"]
ax5.axvline(cruise_mach, color='#16a34a', linestyle=':', linewidth=2, alpha=0.8)
ax5.text(cruise_mach - 0.005, ax5.get_ylim()[0] + 5, f'Cruise\nM={cruise_mach}',
         fontsize=9, color='#16a34a', ha='right', va='bottom')

ax5.set_xlabel("Mach Number", fontsize=12)
ax5.set_ylabel("$C_{D0}$ (drag counts)", fontsize=12)
ax5.set_title("Transonic Drag Rise (Raymer Fig. 12.32 Method)", fontsize=13, fontweight='bold')
ax5.set_xlim(0.5, 1.05)
ax5.legend(fontsize=10)
fig5.tight_layout()


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 6 — Wetted Area Breakdown  (pie chart)
# ═══════════════════════════════════════════════════════════════════════
fig6, axes6 = plt.subplots(1, 2, figsize=(11, 5))

for idx, vname in enumerate(VARIANTS):
    vd = variant_data[vname]
    areas = [vd["S_wet_wing"], vd["S_wet_fuse"], vd["S_wet_htail"],
             vd["S_wet_vtail"], vd["n_nacelles"] * vd["S_wet_nac"]]
    labels_pie = ["Wing", "Fuselage", "H-tail", "V-tail", "Nacelles"]
    pie_colors = ["#3b82f6", "#ef4444", "#22c55e", "#f59e0b", "#8b5cf6"]

    wedges, texts, autotexts = axes6[idx].pie(
        areas, labels=labels_pie, colors=pie_colors, autopct='%1.1f%%',
        startangle=90, pctdistance=0.75, textprops={'fontsize': 9})
    for at in autotexts:
        at.set_fontweight('bold')
        at.set_fontsize(9)
    axes6[idx].set_title(f"{vname}\n$S_{{wet,total}}$ = {vd['S_wet_total']:.0f} ft$^2$",
                          fontsize=11, fontweight='bold')

fig6.suptitle("Wetted Area Breakdown", fontsize=13, fontweight='bold', y=1.02)
fig6.tight_layout()


# ═══════════════════════════════════════════════════════════════════════
#  Save or show
# ═══════════════════════════════════════════════════════════════════════
if "--save" in sys.argv:
    chart_dir = os.path.join(os.path.dirname(__file__), "charts")
    os.makedirs(chart_dir, exist_ok=True)
    for i, fig in enumerate([fig1, fig2, fig3, fig4, fig5, fig6], 1):
        path = os.path.join(chart_dir, f"fig{i}.png")
        fig.savefig(path, dpi=200, bbox_inches='tight')
        print(f"  Saved {path}")
    print("Done.")
else:
    plt.show()
