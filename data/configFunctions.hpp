#pragma once

#include <cmath>
#include <algorithm>

namespace ConfigFuncs {

    // Physical parameters
    constexpr double nu = 1.0;
    
    // Geometry definitions (Straight tube in a 0.32^3 domain)
    constexpr double L_x = 0.32;
    constexpr double yc = 0.16; // Center Y
    constexpr double zc = 0.16; // Center Z
    constexpr double R_vessel = 0.05; // Slightly larger radius for the simple test
    
    // Penalization parameters
    // K_solid represents "zero porosity" (solid)
    // K_fluid represents "high porosity" (fluid)
    constexpr double K_fluid = 1e10;
    constexpr double K_solid = 1e-8;
    constexpr double epsilon = 0.005; // Interface smoothness

    // Helper: Distance to a straight cylinder along X-axis
    inline double get_dist_to_vessel(double /*x*/, double y, double z) {
        // Distance from the centerline (yc, zc)
        return std::sqrt(std::pow(y - yc, 2) + std::pow(z - zc, 2));
    }

    // Brinkman Permeability field
    inline double calc_k(double x, double y, double z) {
        double dist = get_dist_to_vessel(x, y, z);
        // Sigmoid transition: 1 inside, 0 outside
        double indicator = 0.5 * (1.0 - std::tanh((dist - R_vessel) / epsilon));
        return K_solid + indicator * (K_fluid - K_solid);
    }

    // ------------------------------------
    // Boundary Conditions (BCs)
    // ------------------------------------

    inline double bcu_func(double x, double y, double z, double t) {
        double dist_sq = std::pow(y - yc, 2) + std::pow(z - zc, 2);
        double R2 = R_vessel * R_vessel;
        
        // Apply flow only within the cylinder radius at Inlet (x=0) and Outlet (x=Lx)
        if (x < 1e-5 || x > L_x - 1e-5) {
            if (dist_sq < R2) {
                // Time-dependent parabolic profile
                return 10.0 * (1.0 - dist_sq / R2) * std::sin(t);
            }
        }
        
        return 0.0; // No-slip on all other boundaries
    }

    inline double bcv_func(double /*x*/, double /*y*/, double /*z*/, double /*t*/) {
        return 0.0;
    }

    inline double bcw_func(double /*x*/, double /*y*/, double /*z*/, double /*t*/) {
        return 0.0;
    }

    // ------------------------------------
    // Source terms & Initial Conditions
    // ------------------------------------
    inline double fx_func(double x, double y, double z, double t) { return 0.0; }
    inline double fy_func(double x, double y, double z, double t) { return 0.0; }
    inline double fz_func(double x, double y, double z, double t) { return 0.0; }

    inline double u_init_func(double x, double y, double z, double t = 0) { return 0.0; }
    inline double v_init_func(double x, double y, double z, double t = 0) { return 0.0; }
    inline double w_init_func(double x, double y, double z, double t = 0) { return 0.0; }
    inline double p_init_func(double x, double y, double z, double t = 0) { return 0.0; }

    // Physics Fields Export
    inline double k_func(double x, double y, double z, double /*t*/ = 0) {
        return calc_k(x, y, z);
    }
}