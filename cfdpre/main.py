# main.py


import warnings

import numpy as np
import CoolProp.CoolProp as CP



def calculate_f(r, N, delta99, yH):

    """Calculates the function f(r) for root finding."""

    return r**N - r * (delta99 / yH) + (delta99 / yH - 1)



def calculate_df_dr(r, N, delta99, yH):

    """Calculates the derivative of f(r) with respect to r."""

    return N * r**(N-1) - (delta99 / yH)



def newton_raphson(N, delta99, yH, initial_guess=5, tolerance=1e-12, max_iterations=1000):

    """Finds the growth ratio using the Newton-Raphson method."""

    r = initial_guess

    for _ in range(max_iterations):

        f_r = calculate_f(r, N, delta99, yH)

        df_dr = calculate_df_dr(r, N, delta99, yH)

        r_new = r - f_r / df_dr
       

        if np.abs(r_new - r) < tolerance:

            return r_new

        r = r_new

    return r




def yhgrcalc(fluid, temperature_c, pressure_bar, massflow_kgpersec, hydraulicdia_mm,
             target_yplus, num_layers, flow_type='internal', flow_velocity_mpersec=None,
             bl_thickness_mm=None, bl_thickness_fraction=None, roughness_mm=0.0):

    """Calculate boundary layer mesh sizing (first layer height, growth ratio, final layer thickness).

    Two flow types are supported:

    - ``'internal'`` (default): fully-developed pipe/duct flow. ``massflow_kgpersec`` and
      ``hydraulicdia_mm`` are used to derive the flow velocity and Reynolds number. The skin
      friction coefficient uses the laminar (``Re < 2300``, ``16/Re``) or the Haaland (1983)
      explicit turbulent correlation, which is valid across the full turbulent range for both
      smooth and rough pipes (unlike the Blasius ``0.079*Re^-0.25`` form, which is only valid
      for ``Re < ~1e5``). The prism-stack thickness to be spanned defaults to the pipe radius,
      but can be set explicitly with ``bl_thickness_mm`` or as a fraction of the radius with
      ``bl_thickness_fraction`` (filling the whole radius is rarely what you mesh in practice).
    - ``'external'``: flow over a flat plate / aerodynamic surface. ``hydraulicdia_mm`` is
      interpreted as the characteristic length (e.g. chord or plate length) along the flow
      direction. ``flow_velocity_mpersec`` is **required** (the mass-flow/area velocity is
      meaningless for external flow). The skin friction coefficient and boundary layer
      thickness use the laminar (``Re < 5e5``, Blasius flat plate) or turbulent (Schlichting)
      flat-plate correlations; ``bl_thickness_mm`` can override the correlation-based thickness.

    Args:
        fluid (str): CoolProp fluid name (e.g. ``'Air'``, ``'Water'``).
        temperature_c (float): Static temperature in degrees Celsius.
        pressure_bar (float): Static pressure in bar (absolute).
        massflow_kgpersec (float): Mass flow rate in kg/s (used for internal flow).
        hydraulicdia_mm (float): Hydraulic diameter (internal flow) or characteristic
            length (external flow) in mm.
        target_yplus (float): Target y+ value at the first cell centroid.
        num_layers (int): Number of inflation/prism layers.
        flow_type (str): ``'internal'`` (default) or ``'external'``.
        flow_velocity_mpersec (float, optional): Free-stream velocity in m/s. Required for
            ``flow_type='external'``; ignored for internal flow (derived from mass flow).
        bl_thickness_mm (float, optional): Explicit total prism-stack thickness in mm to be
            spanned by the layers. Overrides the default (pipe radius for internal,
            correlation BL thickness for external).
        bl_thickness_fraction (float, optional): Internal flow only. Prism-stack thickness as
            a fraction of the pipe radius (0 < f <= 1), e.g. 0.3 for 30% of the radius.
            Ignored if ``bl_thickness_mm`` is given.
        roughness_mm (float): Absolute wall roughness in mm for the Haaland turbulent
            internal correlation. Default 0.0 (hydraulically smooth).

    Returns:
        dict: Fluid properties, flow conditions, and boundary layer mesh sizing results,
        including ``'Growth Ratio'`` and ``'Final Layer Thickness [m]'``.

    Warns:
        UserWarning: If the boundary layer thickness is not larger than the first layer
            height (mesh sizing inputs are inconsistent), or if the resulting growth ratio
            exceeds the commonly recommended maximum of ~1.3.

    Raises:
        ValueError: If ``flow_type`` is invalid, if ``flow_type='external'`` is requested
            without ``flow_velocity_mpersec``, or if ``bl_thickness_fraction`` is outside (0, 1].

    Example:
        >>> from cfdpre import yhgrcalc
        >>> yhgrcalc('Air', 50, 10, 2.5, 125, 1, 8)
    """

    if flow_type not in ('internal', 'external'):
        raise ValueError("flow_type must be 'internal' or 'external', got %r" % (flow_type,))

    if flow_type == 'external' and flow_velocity_mpersec is None:
        raise ValueError(
            "flow_type='external' requires flow_velocity_mpersec (free-stream velocity in "
            "m/s); deriving velocity from mass flow and a pipe-area assumption is meaningless "
            "for external flow."
        )

    if bl_thickness_fraction is not None and not (0 < bl_thickness_fraction <= 1):
        raise ValueError(
            "bl_thickness_fraction must be in the interval (0, 1], got %r" % (bl_thickness_fraction,)
        )

    pressure_pa = pressure_bar * 1e5
    hydraulicdia_m = hydraulicdia_mm/1000
    temperature_k = temperature_c + 273.15

    dynvisc_pas = CP.PropsSI('V', 'T', temperature_k, 'P', pressure_pa, fluid)
    dynvisc_nsm2 = dynvisc_pas

    thermal_conductivity_wpermk = CP.PropsSI('L', 'T', temperature_k, 'P', pressure_pa, fluid)
    specific_heat_cp_jperkgk = CP.PropsSI('C', 'T', temperature_k, 'P', pressure_pa, fluid)
    density_kgperm3 =  CP.PropsSI('D', 'T', temperature_k, 'P', pressure_pa, fluid)
    kinevisc_m2s = dynvisc_nsm2 / density_kgperm3
    volflowrate_m3persec = massflow_kgpersec / density_kgperm3

    if flow_type == 'external':
        flowvelocity_mpersec = flow_velocity_mpersec
    else:
        flowvelocity_mpersec = volflowrate_m3persec / ((np.pi/4) * np.power(hydraulicdia_m, 2))

    reynolds = flowvelocity_mpersec * hydraulicdia_m / kinevisc_m2s
    prandtl = specific_heat_cp_jperkgk * dynvisc_nsm2 / thermal_conductivity_wpermk

    if flow_type == 'internal':
        radius_m = hydraulicdia_m / 2
        if reynolds < 2300:
            cf = 16 / reynolds  # laminar pipe flow, Fanning friction factor (= Hagen-Poiseuille)
        else:
            # Haaland (1983) explicit approximation to Colebrook. Valid across the full
            # turbulent range for smooth and rough pipes, avoiding Blasius's ~Re<1e5 ceiling.
            rough_ratio = (roughness_mm / 1000) / hydraulicdia_m
            f_darcy = (1.0 / (-1.8 * np.log10((6.9 / reynolds) + (rough_ratio / 3.7) ** 1.11))) ** 2
            cf = f_darcy / 4.0  # convert Darcy -> Fanning (skin friction coefficient)
        # prism-stack thickness to span: explicit > fraction-of-radius > full radius
        if bl_thickness_mm is not None:
            delta99 = bl_thickness_mm / 1000
        elif bl_thickness_fraction is not None:
            delta99 = bl_thickness_fraction * radius_m
        else:
            delta99 = radius_m
    else:
        if reynolds < 5e5:
            cf = 0.664 / np.sqrt(reynolds)  # Blasius laminar flat plate
            delta99 = 4.91 * hydraulicdia_m / np.sqrt(reynolds)
        else:
            cf = np.power(((2 * np.log10(reynolds) - 0.65)), (-2.3))  # Schlichting turbulent flat plate
            delta99 = 0.38 * hydraulicdia_m * reynolds**(-1/5)
        if bl_thickness_mm is not None:
            delta99 = bl_thickness_mm / 1000

    tau_wall = 0.5 * density_kgperm3 * np.square(flowvelocity_mpersec) * cf # wall shear stress
    u_tau = np.sqrt(tau_wall/density_kgperm3) #friction velocity
    yp_m = (target_yplus * dynvisc_nsm2)  /(u_tau * density_kgperm3)
    yh_m = yp_m * 2

    if delta99 <= yh_m:
        warnings.warn(
            "Boundary layer thickness (%.3e m) is not larger than the first layer height "
            "(%.3e m). Reduce target_yplus/num_layers or check the input flow conditions "
            "and geometry." % (delta99, yh_m),
            UserWarning,
        )
        growth_ratio = float('nan')
        final_layer_thickness_m = float('nan')
    else:
        growth_ratio = newton_raphson(num_layers, delta99, yh_m)
        final_layer_thickness_m = yh_m * growth_ratio**(num_layers - 1)
        if growth_ratio > 1.3:
            warnings.warn(
                "Computed growth ratio (%.3f) exceeds the commonly recommended maximum of "
                "~1.3 for boundary layer mesh quality. Consider increasing num_layers or "
                "revisiting target_yplus." % growth_ratio,
                UserWarning,
            )

    result = {
        'fluid' : fluid,
        'temperature [C]' : temperature_c,
        'pressure [bar]' : pressure_bar,
        'massflow [kg/sec]' : massflow_kgpersec,
        'hydraulicdia [mm]' : hydraulicdia_mm,
        'target yplus' : target_yplus,
        'number of layers' : num_layers,
        'flow type' : flow_type,
        'roughness [mm]' : roughness_mm,
        'dynvisc [N-sec/m^2]' : dynvisc_nsm2,
        'thermal conductivity [W/m-k]' : thermal_conductivity_wpermk,
        'specific heat [cp] [J/kg-k]' : specific_heat_cp_jperkgk,
        'density [kg/m^3]' : density_kgperm3,
        'kinematic viscosity [m^2/s]' : kinevisc_m2s,
        'flow velocity [m/sec]' : flowvelocity_mpersec,
        'reynolds number' : reynolds,
        'prandtl number' : prandtl,
        'skin friction coefficient [cf]' : cf,
        'wall shear stress [tau_wall]' : tau_wall,
        'boundary layer thickness [delta99] [m]' : delta99,
        'height of cell centroid from wall [yp] [m]' : yp_m,
        'first layer height [yh] [m]' : yh_m,
        'Growth Ratio' : growth_ratio,
        'Final Layer Thickness [m]' : final_layer_thickness_m
    }

    return result



    
