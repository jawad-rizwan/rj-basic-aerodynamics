# rj-basic-aerodynamics

Conceptual-level aerodynamic analysis for a high-wing regional jet design
(AN-148/AVRO RJ class), implementing the methods from **Raymer, _Aircraft
Design: A Conceptual Approach_ (7th Ed.), Chapter 12 — Aerodynamics**.

The script computes CD0, Oswald span efficiency, induced-drag K factor,
lift-curve slope, CLmax, and drag-divergence Mach number using the component
buildup method. Both 76-seat and 100-seat fuselage variants are analysed
(shared wing/tail/nacelle geometry, different fuselage lengths).

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

- **Wing**: S_ref = 792.47 ft², AR = 7.8, taper = 0.33, NASA SC(3)-0712B airfoil (12% t/c), 26° LE sweep, no winglets (anhedral tips)
- **Fuselage**: 9.83 ft diameter, 96.44 ft (76-seat) / 108.2 ft (100-seat)
- **Horizontal tail**: T-tail, S = 190.58 ft², NASA SC(2)-0012 airfoil
- **Vertical tail**: S = 146.35 ft² (enhanced), NASA SC(2)-0012 airfoil
- **Nacelles**: 13.46 ft length, 6.17 ft diameter, x2
- **Cruise**: Mach 0.78 at 41,000 ft

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

| Parameter | 76-seat | 100-seat | Notes |
|-----------|---------|----------|-------|
| CD0 | 0.02113 | 0.02185 | Component buildup (Eq. 12.24) |
| e (Oswald) | 0.639 | 0.639 | Eq. 12.48/49 — conservative |
| e (LE suction) | 0.753 | 0.753 | Eq. 12.57 — recommended |
| K (LE suction) | 0.0542 | 0.0542 | Recommended for sizing |
| (L/D)_max | 14.78 | 14.53 | With LE suction K |
| L/D cruise | 14.40 | 14.39 | At cruise Mach / altitude |
| CL cruise | 0.4974 | 0.5507 | W_cruise / (q * S_ref) |
| CLa | 5.73 /rad | 5.73 /rad | At cruise Mach |
| CLmax clean | 1.841 | 1.841 | Eq. 12.15 |
| CLmax takeoff | 2.514 | 2.514 | Fowler flaps, no slats |
| CLmax landing | 3.082 | 3.082 | Fowler flaps + slats |
| MDD | 0.826 | 0.819 | Korn equation (Eq. 12.46) |

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
