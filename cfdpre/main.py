# main.py


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




def yhgrcalc(fluid, temperature_c, pressure_bar, massflow_kgpersec, hydraulicdia_mm, target_yplus, num_layers):
    
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
    flowvelocity_mpersec = volflowrate_m3persec / ((np.pi/4) * np.power(hydraulicdia_m, 2))
    reynolds = flowvelocity_mpersec * hydraulicdia_m / kinevisc_m2s
    prandtl = specific_heat_cp_jperkgk * dynvisc_nsm2 / thermal_conductivity_wpermk
    cf = np.power(((2 * np.log10(reynolds)- 0.65)), (-2.3)) # skin friction coefficient
    tau_wall = 0.5 * density_kgperm3 * np.square(flowvelocity_mpersec) * cf # wall shear stress
    u_tau = np.sqrt(tau_wall/density_kgperm3) #friction velocity
    yp_m = (target_yplus * dynvisc_nsm2)  /(u_tau * density_kgperm3)
    yh_m = yp_m * 2   

    if reynolds < 5e5:
        delta99 = 4.91 * hydraulicdia_m / np.sqrt(reynolds)
    else:
        delta99 = 0.38 * hydraulicdia_m * reynolds**(-1/5)

    growth_ratio = newton_raphson(num_layers, delta99, yh_m)
    final_layer_thickness_m = yh_m * growth_ratio**(num_layers - 1)
    

    result = {
        'fluid' : fluid,
        'temperature [C]' : temperature_c,
        'pressure [bar]' : pressure_bar,
        'massflow [kg/sec]' : massflow_kgpersec,
        'hydraulicdia [mm]' : hydraulicdia_mm,
        'target yplus' : target_yplus,
        'number of layers' : num_layers,
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
        'height of cell centroid from wall [yp] [m]' : yp_m,
        'first layer height [yh] [m]' : yh_m,
        'Growth Ratio' : growth_ratio,
        'Final Layer Thickness [m]' : final_layer_thickness_m
    }
    
    return result



    
