# rj-basic-aerodynamics

Conceptual-level aerodynamic analysis for a high-wing regional jet design,
implementing the methods from **Raymer, _Aircraft Design: A Conceptual
Approach_ (7th Ed.), Chapter 12 — Aerodynamics**.

The script computes CD0, Oswald span efficiency, induced-drag K factor,
lift-curve slope, CLmax, drag-divergence Mach number, and required cruise
thrust using the component buildup method. Both ZRJ70 (76-seat) and ZRJ100 (100-seat) variants are
analysed (shared wing/tail/nacelle geometry, different fuselage lengths).

Outputs are intended to feed into Chapter 6 refined mission sizing
(see [`rj-mission-sizing`](https://github.com/jawad-rizwan/rj-mission-sizing)).

## Repository Structure

```
rj-basic-aerodynamics/
├── aero/                    # Core library modules
│   ├── atmosphere.py        # ISA standard atmosphere model
│   ├── skin_friction.py     # Flat-plate Cf (Eq. 12.25–12.29)
│   ├── form_factors.py      # Component form factors (Eq. 12.30–12.35)
│   ├── parasite_drag.py     # CD0: component buildup & Cfe method (Eq. 12.23–12.24)
│   ├── lift.py              # CLa, CLmax clean & with flaps (Eq. 12.6–12.22)
│   ├── induced_drag.py      # Oswald e, K factor, LE suction (Eq. 12.47–12.57)
│   ├── wave_drag.py         # MDD, wave drag (Eq. 12.42–12.46)
│   └── drag_polar.py        # DragPolar class (Eq. 12.4–12.5)
├── data/                    # Aircraft configuration data
│   ├── ZRJ70.py             # 76-seat variant geometry & parameters
│   └── ZRJ100.py            # 100-seat variant geometry & parameters
├── examples/
│   └── regional_jet.py      # Full analysis for both fuselage variants
└── requirements.txt
```

## Quick Start

```bash
pip install -r requirements.txt
python3 examples/regional_jet.py
```

## Aircraft Configuration

- **Wing**: S_ref = 1016.58 ft², AR = 7.8 (effective 9.07 with winglets), taper = 0.33, supercritical airfoils (root NASA SC(2)-0714, mid SC(3)-0712, tip SC(2)-0710), 26° LE sweep, scimitar winglets (5.67 ft upper + 1.30 ft lower = 6.97 ft total), composite construction
- **Fuselage**: 10.5 ft diameter, 96.7 ft (ZRJ70) / 110.2 ft (ZRJ100), aluminium
- **Horizontal tail**: T-tail, S = 276.92 ft², NASA SC(2)-0010 airfoil, composite
- **Vertical tail**: S = 200.8 ft², NASA SC(2)-0012 airfoil, composite
- **Nacelles**: 13.46 ft length, 6.17 ft diameter, x2, composite
- **Cruise**: Mach 0.78 at 35,000 ft

## Module Reference

Each module corresponds to a section of Raymer Chapter 12. Every function
includes the Raymer equation number in a comment.

| Module | Key Functions | Raymer Equations |
|--------|--------------|-----------------|
| `atmosphere.py` | `isa_properties`, `dynamic_pressure` | ISA model |
| `skin_friction.py` | `cf_laminar`, `cf_turbulent`, `cf_avg` | Eq. 12.25–12.29 |
| `form_factors.py` | `ff_wing`, `ff_fuselage`, `ff_nacelle` | Eq. 12.30–12.35 |
| `parasite_drag.py` | `cd0_component_buildup`, `cd0_equivalent_cfe` | Eq. 12.23–12.24 |
| `lift.py` | `cl_alpha_subsonic`, `cl_max_clean`, `cl_max_flaps` | Eq. 12.6–12.22 |
| `induced_drag.py` | `oswald_e`, `k_factor`, `k_factor_leading_edge_suction` | Eq. 12.47–12.57 |
| `wave_drag.py` | `mdd_wing`, `wave_drag_sears_haack` | Eq. 12.42–12.46 |
| `drag_polar.py` | `DragPolar` class | Eq. 12.4–12.5 |

## Sample Output

| Parameter | ZRJ70 | ZRJ100 | Notes |
|-----------|-------|--------|-------|
| CD0 | 0.01843 | 0.01910 | Component buildup (Eq. 12.24) |
| e (Oswald) | 0.557 | 0.557 | Eq. 12.48/49 — conservative |
| e (LE suction) | 0.727 | 0.727 | Eq. 12.57 — recommended |
| K (LE suction) | 0.0483 | 0.0483 | Recommended for sizing |
| (L/D)_max | 16.76 | 16.46 | With LE suction K |
| L/D cruise | 13.60 | 13.94 | M 0.78 at 35,000 ft |
| CL cruise | 0.3166 | 0.3477 | W_cruise / (q * S_ref) |
| CLa | 5.99 /rad | 5.99 /rad | At cruise Mach |
| CLmax clean | 1.161 | 1.161 | Eq. 12.15 |
| CLmax takeoff | 1.990 | 1.990 | Double slotted flaps, no slats |
| CLmax landing | 2.624 | 2.624 | Double slotted flaps + slats |
| MDD | 0.846 | 0.842 | Korn equation (Eq. 12.46) |
| T_required | 5,017 lbs | 5,376 lbs | Total cruise thrust (T = D) |
| T per engine | 2,509 lbs | 2,688 lbs | x2 engines |
| T/W_cruise | 0.0735 | 0.0717 | Thrust-to-weight at cruise |

## Stability & Control Values

Values marked with \* are estimates — update `alpha_0L_deg`, `cm_0_airfoil`,
`cl_max_htail_airfoil`, and `cl_max_vtail_airfoil` in the data files with
XFLR5 results for improved accuracy.

### Wing

| Parameter | ZRJ70 | ZRJ100 | Notes |
|-----------|-------|--------|-------|
| a (dCL/dα) | 5.9874 /rad | 5.9874 /rad | Eq. 12.6 at cruise Mach |
| Cm_0_ac | −0.0618 \* | −0.0618 \* | 3D correction of 2D cm_0 |
| CLmax (clean) | 1.161 | 1.161 | Eq. 12.15 |
| α_CLmax | 9.1° \* | 9.1° \* | Depends on α_0L |
| CL_cruise | 0.3166 | 0.3477 | W_cruise / (q × S_ref) |
| α_CL_cruise | 1.03° \* | 1.33° \* | Depends on α_0L |

### Horizontal Tail (AR = 3.94)

| Parameter | Value | Notes |
|-----------|-------|-------|
| a_t (dCL/dα) | 4.0450 /rad | Eq. 12.6, symmetric SC(2)-0010 |
| CLmax_tail | 0.808 \* | Needs cl_max from XFLR5 |
| α_CLmax_tail | 11.4° \* | Symmetric airfoil (α_0L ≈ 0) |

### Vertical Tail (AR = 1.00)

| Parameter | Value | Notes |
|-----------|-------|-------|
| a_t (dCL/dα) | 1.2838 /rad | Eq. 12.6, symmetric SC(2)-0012 |
| CLmax_tail | 0.796 \* | Needs cl_max from XFLR5 |
| α_CLmax_tail | 35.5° \* | High due to very low AR |

### Not Yet Computable

The following require CG position, tail moment arms, and downwash gradient
(dε/dα), which are beyond the scope of Raymer Ch. 12 aerodynamics:

- dCm/dα (wing, htail, vtail contributions)
- CL_cruise and α_cruise for tail surfaces (need trim analysis)

## Notes

- The Oswald e from Eq. 12.48/12.49 is known to be conservative. Raymer
  himself notes that a well-designed wing should achieve e ~ 0.85–0.90 at
  design CL. The leading-edge suction method (Eq. 12.57) gives more
  realistic results and is the recommended approach.
- The `DragPolar` class supports both uncambered (Eq. 12.4) and cambered
  (Eq. 12.5) polars via the optional `cl_min_drag` parameter.
- All computations use Imperial units (ft, lbs, slugs, Rankine) to match
  Raymer's conventions.
- Surface finish uses polished composite (k = 0.50e-5 ft) for wing, tail,
  and nacelles, and production sheet metal (k = 1.33e-5 ft) for the
  aluminium fuselage.

## Requirements

- Python 3.8+
- NumPy >= 1.20
