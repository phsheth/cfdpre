<h1 align="center">
<img src="https://raw.githubusercontent.com/phsheth/cfdpre/refs/heads/main/cfdprelogo.png" width="300">
</h1><br>


What is CFDPre?
----------------------

CFDPre is an open-source collection of object-oriented software tools for
calculating boundary layer mesh dimensions for Computational Fluid Dynamics simulations.
Among other things, it can be used to:

* Calculate First Layer Thickness, Growth Ratio and Final Layer Thickess


Installation
----------------------

[![Install](https://img.shields.io/pypi/v/cfdpre?label=CFDPre)](
https://pypi.org/project/cfdpre/) [![PyPI Downloads](https://img.shields.io/pypi/dm/cfdpre?label=PyPI%20Downloads)](
https://pypistats.org/packages/cfdpre)

In your command line, within Python environment:
```pythom
    pip install cfdpre
```
- The Python module can also be installed using pip on Windows, macOS, and Linux.


Usage
----------------------

`yhgrcalc` supports two flow types:

- **`flow_type='internal'`** (default) — fully-developed pipe/duct flow. `massflow_kgpersec`
  and `hydraulicdia_mm` are used to derive velocity and Reynolds number. Skin friction uses
  the laminar (`Re < 2300`, `16/Re`) or the **Haaland (1983)** explicit turbulent correlation,
  valid across the full turbulent range for smooth and rough pipes (`roughness_mm`) — unlike
  Blasius, which is only accurate below `Re ≈ 1e5`. The prism-stack thickness defaults to the
  pipe radius, but you can set it explicitly with `bl_thickness_mm` or as a fraction of the
  radius with `bl_thickness_fraction` (filling the entire radius with prisms is rarely what
  you actually mesh).
- **`flow_type='external'`** — flow over a flat plate / aerodynamic surface.
  `hydraulicdia_mm` is interpreted as a characteristic length (e.g. chord length), and
  `flow_velocity_mpersec` is **required** (the mass-flow/area velocity is meaningless for
  external flow). Skin friction and boundary layer thickness use the laminar (Blasius) or
  turbulent (Schlichting) flat-plate correlations; `bl_thickness_mm` can override the
  correlation-based thickness.

If the computed growth ratio exceeds ~1.3 (or the boundary layer is too thin to resolve
with the requested `target_yplus`/`num_layers`), a `UserWarning` is raised so you can
adjust your inputs.

```python
from cfdpre import yhgrcalc
yhgrcalc('Air', 50, 10, 2.5, 125, 1, 8)
```

Output:
```python
{'fluid': 'Air',
 'temperature [C]': 50,
 'pressure [bar]': 10,
 'massflow [kg/sec]': 2.5,
 'hydraulicdia [mm]': 125,
 'target yplus': 1,
 'number of layers': 8,
 'flow type': 'internal',
 'roughness [mm]': 0.0,
 'dynvisc [N-sec/m^2]': 1.9762497305390764e-05,
 'thermal conductivity [W/m-k]': 0.028357331300649127,
 'specific heat [cp] [J/kg-k]': 1019.3146170790077,
 'density [kg/m^3]': 10.792698589669245,
 'kinematic viscosity [m^2/s]': 1.8310987878701066e-06,
 'flow velocity [m/sec]': 18.875569021507275,
 'reynolds number': 1288541.1444310248,
 'prandtl number': 0.7103701741111368,
 'skin friction coefficient [cf]': 0.002776948990231969,
 'wall shear stress [tau_wall]': 5.339100066909961,
 'boundary layer thickness [delta99] [m]': 0.0625,
 'height of cell centroid from wall [yp] [m]': 2.6034112015544278e-06,
 'first layer height [yh] [m]': 5.2068224031088555e-06,
 'Growth Ratio': 3.6553638704281375,
 'Final Layer Thickness [m]': 0.04540326342524223}
```

> Note: with `target_yplus=1` and only 8 layers at this high Reynolds number (~1.29e6),
> the growth ratio comes out to ~3.66, well above the recommended ~1.3 maximum — this
> example raises a `UserWarning` accordingly. In practice you'd either use far more layers
> (20-40+) or restrict the prism stack to part of the radius. For example, spanning only
> 30% of the radius:
>
> ```python
> yhgrcalc('Air', 50, 10, 2.5, 125, 1, 8, bl_thickness_fraction=0.3)  # GR ~3.04
> ```

External flow example (`flow_velocity_mpersec` is required):
```python
yhgrcalc('Air', 25, 1.01325, 2.5, 1000, 1, 15, flow_type='external', flow_velocity_mpersec=30)
```

Documentation
----------------------

In progress - not yet made!

- **Project Home Page:** https://cfdpre.github.io/ [under construction]
- **Users Group:** https://groups.google.com/g/cfdpre
- **Source code:** https://github.com/phsheth/cfdpre
- **PyPI Page:** https://pypi.org/project/cfdpre/



Call for Contributors
----------------------

The CFDPre project welcomes your expertise and enthusiasm! Better to discuss on the users group before starting to contribute!


Project Log
----------------------
January 2025:
1. Created Library

June 2026:
1. Added `flow_type` ('internal'/'external') support. Internal flow uses laminar `16/Re`
   and the Haaland (1983) turbulent correlation (valid to high Re, smooth/rough via
   `roughness_mm`); external flow uses Blasius/Schlichting flat-plate correlations and now
   requires `flow_velocity_mpersec`. Prism-stack thickness is selectable via `bl_thickness_mm`
   or `bl_thickness_fraction`. Added growth-ratio (>1.3) and boundary-layer sanity warnings.


Project RoadMap:
----------------------

1. Documentation for existing functionality.
2. Include example data within library.
3. Object-oriented API refactor + input validation.
4. Generate `snappyHexMeshDict`-ready `addLayersControls` snippets.
5. Bundled examples, unit tests, CI, and docs site.
6. Parametric/batch sweeps.






