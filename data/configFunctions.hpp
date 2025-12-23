#pragma once

#include <cmath>
#include <algorithm>

namespace ConfigFuncs {

    // --- Physical Parameters ---
    constexpr double nu = 1.0;          // Kinematic viscosity
    
    // --- Geometry Definitions (Straight tube in a 0.32^3 domain) ---
    constexpr double L_x = 1.0;
    constexpr double yc = L_x / 2;         // Cylinder center Y
    constexpr double zc = L_x / 2;         // Cylinder center Z
    constexpr double R_vessel = L_x / 4;   // Vessel radius
    
    // --- Brinkman Penalization Parameters ---
    constexpr double K_fluid = 1e20;    // High permeability (Fluid region)
    constexpr double K_solid = 1e-20;    // Low permeability (Solid region)
    constexpr double epsilon = 1e-10;   // Interface smoothing width
    constexpr double G_force = 10.0;    // Imposed pressure gradient (Body force)

    // Helper: Squared distance from the cylinder axis
    inline double get_dist_to_vessel_sq(double y, double z) {
        return std::pow(y - yc, 2) + std::pow(z - zc, 2);
    }

    // Permeability field calculation
    inline double calc_k(double x, double y, double z) {
        double dist = std::sqrt(get_dist_to_vessel_sq(y, z));
        // Sigmoid transition from solid to fluid
        double indicator = 0.5 * (1.0 - std::tanh((dist - R_vessel) / epsilon));
        return K_solid + indicator * (K_fluid - K_solid);
    }

    // ------------------------------------
    // Boundary Conditions (BCs)
    // ------------------------------------

    // Inlet velocity (U) - Analytical Poiseuille profile
    inline double bcu_func(double x, double y, double z, double t) {
        double dist_sq = get_dist_to_vessel_sq(y, z);
        double R2 = R_vessel * R_vessel;

        if (dist_sq < R2) {
            // Time signal synced with the source term (cosine)
            double time_signal = std::cos(t); 
            
            // Poiseuille solution: u(r) = G/(4*nu) * (R^2 - r^2)
            double poiseuille = (G_force / (4.0 * nu)) * (R2 - dist_sq);
            
            return poiseuille * time_signal;
        }

        return 0.0; // No-slip condition
    }

    inline double bcv_func(double, double, double, double) { return 0.0; }
    inline double bcw_func(double, double, double, double) { return 0.0; }

    // ------------------------------------
    // Source terms (Driving Force)
    // ------------------------------------

    // Momentum source term in X (Pressure gradient substitute)
    inline double fx_func(double x, double y, double z, double t) {
        double dist_sq = get_dist_to_vessel_sq(y, z);
        double R_eff_sq = std::pow(R_vessel + epsilon, 2);

        // Apply force G only within the fluidic channel
        if (dist_sq < R_eff_sq) {
            return G_force * std::cos(t); 
        }
        return 0.0;
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