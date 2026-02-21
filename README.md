# rj-basic-aerodynamics

Conceptual-level aerodynamic analysis for regional jet design, implementing
the methods from **Raymer, _Aircraft Design: A Conceptual Approach_ (7th Ed.),
Chapter 12 — Aerodynamics**.

These scripts compute CD0, Oswald span efficiency, induced-drag K factor,
lift-curve slope, CLmax, and drag-divergence Mach number using the component
buildup method. The outputs are intended to feed into Chapter 6 refined mission
sizing (see [`rj-mission-sizing`](https://github.com/jawad-rizwan/rj-mission-sizing)).

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
├── examples/
│   └── regional_jet.py      # Full example with CRJ700-class placeholder data
└── requirements.txt
```

## Quick Start

```bash
pip install -r requirements.txt
python3 examples/regional_jet.py
```

## Using Your Own Aircraft Data

The example script (`examples/regional_jet.py`) is pre-loaded with CRJ700-class
geometry as placeholder data. Every value you need to replace is marked with
`*** UPDATE ***` in the comments. Search for that tag to find all of them:

```bash
grep "UPDATE" examples/regional_jet.py
```

The placeholders cover:
- **Wing**: S_ref, AR, taper, t/c, sweep angles, S_exposed, winglet height
- **Fuselage**: length, diameter
- **Horizontal tail**: S, S_exposed, MAC, t/c, x/c_max, sweep
- **Vertical tail**: S, S_exposed, MAC, t/c, x/c_max, sweep
- **Nacelles**: length, diameter, count
- **Flight condition**: cruise Mach, cruise altitude
- **Airfoil data**: 2D CLmax
- **Surface finish**: roughness height k
- **Tuning parameters**: interference factors Q, laminar %, leakage %, LE suction
- **High-lift devices**: flap/slat deltas, flapped span ratios, hinge-line sweep
- **Weight**: mid-mission cruise weight estimate

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

## Sample Output (CRJ700-class placeholders)

| Parameter | Value | Notes |
|-----------|-------|-------|
| CD0 | 0.02541 | Component buildup (Eq. 12.24) |
| e (Oswald) | 0.543 | Eq. 12.48/49 — conservative |
| e (LE suction) | 0.730 | Eq. 12.57 — recommended |
| K (LE suction) | 0.0502 | Recommended for sizing |
| (L/D)_max | 14.0 | With LE suction K |
| CLa | 5.88 /rad | At cruise Mach |
| CLmax clean | 1.294 | Eq. 12.15 |
| CLmax landing | 2.535 | With flaps + slats |
| MDD | 0.867 | Korn equation (Eq. 12.46) |

## Notes

- The Oswald e from Eq. 12.48/12.49 is known to be conservative. Raymer
  himself notes that a well-designed wing should achieve e ~ 0.85–0.90 at
  design CL. The leading-edge suction method (Eq. 12.57) gives more
  realistic results and is the recommended approach.
- The `DragPolar` class supports both uncambered (Eq. 12.4) and cambered
  (Eq. 12.5) polars via the optional `cl_min_drag` parameter.
- All computations use Imperial units (ft, lbs, slugs, Rankine) to match
  Raymer's conventions.

## Requirements

- Python 3.8+
- NumPy >= 1.20
