#pragma once

#include <cmath> // For std::sin, std::cos, std::exp, M_PI

namespace ConfigFuncs {

    // ------------------------------------
    // Boundary Conditions (BCs)
    // f(x, y, z, t)
    // ------------------------------------

    // Boundary Condition U (u_expr)
    inline double bcu_func(double x, double y, double z, double t) {
        return std::sin(t) * std::sin(x) * std::sin(y) * std::sin(z);
    }

    // Boundary Condition V (v_expr)
    inline double bcv_func(double x, double y, double z, double t) {
        return std::sin(t) * std::cos(x) * std::cos(y) * std::cos(z);
    }

    // Boundary Condition W (w_expr)
    inline double bcw_func(double x, double y, double z, double t) {
        return std::sin(t) * std::cos(x) * std::sin(y) * (std::sin(z) + std::cos(z));
    }

    // ------------------------------------
    // Forces
    // f(t, x, y, z)
    // ------------------------------------

    // Force Fx (fx_expr)
    inline double fx_func(double x, double y, double z, double t) {
        double term1 = std::cos(t) * std::sin(x) * std::sin(y) * std::sin(z);
        double term2 = std::sin(t) * ((21.0 + 6e-20) * std::sin(x) * std::sin(y) * std::sin(z) -
                                      3.0 * std::sin(x) * std::sin(y) * std::cos(z));
        return term1 + term2;
    }

    // Force Fy (fy_expr)
    inline double fy_func(double x, double y, double z, double t) {
        double term1 = std::cos(t) * std::cos(x) * std::cos(y) * std::cos(z);
        double term2 = std::sin(t) * ((21.0 + 6e-20) * std::cos(x) * std::cos(y) * std::cos(z) -
                                      3.0 * std::cos(x) * std::cos(y) * std::sin(z));
        return term1 + term2;
    }

    // Force Fz (fz_expr)
    inline double fz_func(double x, double y, double z, double t) {
        double term1 = std::cos(t) * std::cos(x) * std::sin(y) * (std::sin(z) + std::cos(z));
        double term2 = std::sin(t) * ((15.0 + 6e-20) * std::cos(x) * std::sin(y) * (std::sin(z) + std::cos(z)));
        return term1 + term2;
    }

    // ------------------------------------
    // Initial Conditions (ICs)
    // f(x, y, z, t) - Note: t will always be 0.0 in setup
    // ------------------------------------

    // Initial Velocity U (u_expr)
    inline double u_init_func(double x, double y, double z, double t = 0) {
        return std::sin(t) * std::sin(x) * std::sin(y) * std::sin(z);
    }

    // Initial Velocity V (v_expr)
    inline double v_init_func(double x, double y, double z, double t = 0) {
        return std::sin(t) * std::cos(x) * std::cos(y) * std::cos(z);
    }

    // Initial Velocity W (w_expr)
    inline double w_init_func(double x, double y, double z, double t = 0) {
        return std::sin(t) * std::cos(x) * std::sin(y) * (std::sin(z) + std::cos(z));
    }

    // Initial Pressure (p_expr)
    inline double p_init_func(double x, double y, double z, double t = 0) {
        return -std::sin(t) * 3.0 * std::cos(x) * std::sin(y) * (std::sin(z) - std::cos(z));
    }

    // ------------------------------------
    // Physics
    // f(x, y, z, t) - Note: t will always be 0.0 if time-independent
    // ------------------------------------

    // Permeability K (k_expr)
    inline double k_func(double /*x*/, double /*y*/, double /*z*/, double /*t*/ = 0) {
        return 1e20;
    }

}