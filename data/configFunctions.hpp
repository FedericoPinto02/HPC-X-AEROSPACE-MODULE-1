#pragma once

#include <cmath>
#include <algorithm>

namespace ConfigFuncs {
    // --------------------------------------------- Plane Couette Poiseuille ------------------------------------
    // ------------------------------------------------------------------------------------------------
    // --- Physical Parameters ---
    constexpr double nu = 1.0;          // Kinematic viscosity
    
    // --- Geometry Definitions  ---
    constexpr double L_x = 1.0;
    constexpr double U_max = 5.0;
    constexpr double G = 100.0;

    // --- Brinkman Penalization Parameters ---
    constexpr double K_fluid = 1e15;    // High permeability (Fluid region)

    // Permeability field calculation
    inline double calc_k(double x, double y, double z) {
        return K_fluid;
    }

    // ------------------------------------
    // Boundary Conditions (BCs)
    // ------------------------------------

    // Inlet velocity (U) - Analytical Poiseuille profile
    inline double bcu_func(double x, double y, double z, double t) {
        double couette = std::sin(t) * y * U_max / L_x;
        double poiseuille = - std::sin(t) * y / nu * G * (L_x - y) * 0.5;
        return couette + poiseuille ; 
    }

    inline double bcv_func(double, double, double, double) { return 0.0; }
    inline double bcw_func(double, double, double, double) { return 0.0; }

    // ------------------------------------
    // Source terms (Driving Force)
    // ------------------------------------

    // Momentum source term in X (Pressure gradient substitute)
    inline double fx_func(double x, double y, double z, double t) {
        return std::sin(t) * G;
    }

    inline double fy_func(double, double, double, double) { return 0.0; }
    inline double fz_func(double, double, double, double) { return 0.0; }

    // ------------------------------------
    // Initial Conditions
    // ------------------------------------
     inline double u_init_func(double x, double y, double z, double t = 0) {
        return bcu_func(x, y, z, t);
    }

    inline double v_init_func(double x, double y, double z, double t = 0) {
        return bcv_func(x, y, z, t);
    }

    inline double w_init_func(double x, double y, double z, double t = 0) {
        return bcw_func(x, y, z, t);
    }
     inline double p_init_func(double, double, double, double) { return 10.0; }

    // Export permeability field
    inline double k_func(double x, double y, double z, double /*t*/ = 0) {
        return calc_k(x, y, z);
    }
}