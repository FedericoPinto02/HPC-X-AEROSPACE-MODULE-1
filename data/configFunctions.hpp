#pragma once

#include <cmath> // For std::sin, std::cos

/**
 * @namespace ConfigFuncs
 * @brief Contains analytical functions for initial conditions, boundary conditions, and forcing terms.
 * * These functions define a "Manufactured Solution" used to validate the numerical accuracy 
 * of the solver. The default setup represents a divergence-free velocity field.
 */
namespace ConfigFuncs {

    // Physical parameters
    constexpr double nu = 6.0;
    constexpr double Re = 1.0;

    // Helper: Permeability field k(x, y, z)
    inline double calc_k(double x, double y, double z) {
        return 10.0 * (2.0 + std::cos(x) * std::cos(y) * std::cos(z));
    }

    // ------------------------------------
    // Boundary Conditions (BCs) & Exact Solution
    // ------------------------------------

    // u = sin(x)cos(t+y)sin(z)
    inline double bcu_func(double x, double y, double z, double t) {
        return std::sin(x) * std::cos(t + y) * std::sin(z);
    }

    // v = cos(x)sin(t+y)sin(z)
    inline double bcv_func(double x, double y, double z, double t) {
        return std::cos(x) * std::sin(t + y) * std::sin(z);
    }

    // w = 2cos(x)cos(t+y)cos(z)
    inline double bcw_func(double x, double y, double z, double t) {
        return 2.0 * std::cos(x) * std::cos(t + y) * std::cos(z);
    }

    // ------------------------------------
    // Forces
    // f = dt_u - nu*lapl_u + (nu/k)*u + grad_p
    // Broken down into: Transient + Viscous + Drag + Pressure
    // ------------------------------------

    inline double fx_func(double x, double y, double z, double t) {
        double k = calc_k(x, y, z);
        
        // Terms
        double transient = -std::sin(x) * std::sin(t + y) * std::sin(z);
        double viscous   = 3.0 * nu * std::sin(x) * std::cos(t + y) * std::sin(z);
        double drag      = (nu / k) * std::sin(x) * std::cos(t + y) * std::sin(z);
        double pressure  = -(3.0 / Re) * std::sin(x) * std::cos(t + y) * std::cos(z);

        return transient + viscous + drag + pressure;
    }

    inline double fy_func(double x, double y, double z, double t) {
        double k = calc_k(x, y, z);

        // Terms
        double transient = std::cos(x) * std::cos(t + y) * std::sin(z);
        double viscous   = 3.0 * nu * std::cos(x) * std::sin(t + y) * std::sin(z);
        double drag      = (nu / k) * std::cos(x) * std::sin(t + y) * std::sin(z);
        double pressure  = -(3.0 / Re) * std::cos(x) * std::sin(t + y) * std::cos(z);

        return transient + viscous + drag + pressure;
    }

    inline double fz_func(double x, double y, double z, double t) {
        double k = calc_k(x, y, z);

        // Terms (Note: w has a factor of 2, so viscous is 3*nu*2 = 6*nu)
        double transient = -2.0 * std::cos(x) * std::sin(t + y) * std::cos(z);
        double viscous   = 6.0 * nu * std::cos(x) * std::cos(t + y) * std::cos(z); 
        double drag      = (2.0 * nu / k) * std::cos(x) * std::cos(t + y) * std::cos(z);
        double pressure  = -(3.0 / Re) * std::cos(x) * std::cos(t + y) * std::sin(z);

        return transient + viscous + drag + pressure;
    }

    // ------------------------------------
    // Initial Conditions (ICs)
    // t = 0
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

    // p = (3/Re) cos(x)cos(t+y)cos(z)
    inline double p_init_func(double x, double y, double z, double t = 0) {
        return (3.0 / Re) * std::cos(x) * std::cos(t + y) * std::cos(z);
    }

    // ------------------------------------
    // Physics Fields
    // ------------------------------------

    // Permeability
    inline double k_func(double x, double y, double z, double /*t*/ = 0) {
        return calc_k(x, y, z);
    }

}