#!/usr/bin/env python3

"""
This script performs a spatial convergence analysis for a Method of
Manufactured Solutions (MMS) test.

It reads multiple VTK files, each from a simulation with a different
grid resolution. For each file, it computes the L2 RMS error
of:
1. Velocity components (u, v, w)
2. Pressure (p) - Both L2 and L1 norms
3. Pressure Gradient (gradP_x, gradP_y, gradP_z)
4. Velocity Divergence (divU) - Analytical solution is 0.

This version is adapted for staggered grid solvers. It assumes the VTK output
contains vector fields for Velocity and Pressure Gradient.
"""

# --- 0. IMPORTS ---
import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import os
import sys
import re 
from matplotlib import ticker

# --- 1. USER CONFIGURATION ---

# Physical Constants (Must match C++ config)
RE = 1.0

# List of simulation results to analyze.
# NOTE: This list is now IGNORED if run_convergence_study_unified.py is used.
SIMULATIONS = [
    {"nx": 30, "file": "../output/n30simulation_output_0030.vtk"},
    {"nx": 20, "file": "../output/n20simulation_output_0030.vtk"},
]

DOMAIN_LENGTH_X = 6.0 # Assumed, check your domain size
DT = 0.001 

# Field names in the VTK file
VELOCITY_FIELD_NAME = "velocity"
PRESSURE_FIELD_NAME = "pressure"
GRAD_PRESSURE_FIELD_NAME = "gradP" 
DIV_VELOCITY_FIELD_NAME = "divU"  # Field name for Divergence


# Helper function to get step number from filename
def get_step_from_filename(filepath: str) -> int:
    match = re.search(r'_(\d+)\.vtk$', os.path.basename(filepath))
    if match:
        return int(match.group(1))
    else:
        raise ValueError(f"Could not parse step number from filename: {filepath}.")


# --- 2. ANALYTICAL SOLUTION FUNCTIONS ---

# --- Velocity (u, v, w) ---
# u = sin(x)cos(t+y)sin(z)
def compute_analytical_u_at_x_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    x = points[:, 0] + h / 2.0  
    y = points[:, 1]
    z = points[:, 2]
    u_ana = np.sin(x) * np.cos(t + y) * np.sin(z)
    return u_ana

# v = cos(x)sin(t+y)sin(z)
def compute_analytical_v_at_y_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1] + h / 2.0  
    z = points[:, 2]
    v_ana = np.cos(x) * np.sin(t + y) * np.sin(z)
    return v_ana

# w = 2cos(x)cos(t+y)cos(z)
def compute_analytical_w_at_z_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2] + h / 2.0  
    w_ana = 2.0 * np.cos(x) * np.cos(t + y) * np.cos(z)
    return w_ana

# --- Pressure (p) ---
# p = (3/Re) * cos(x)cos(t+y)cos(z)
def compute_analytical_p_at_centers(points: np.ndarray, t: float) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    p_ana = (3.0 / RE) * np.cos(x) * np.cos(t + y) * np.cos(z)
    return p_ana

# --- Pressure Gradient (gradP) ---
# Evaluated at velocity points (staggered faces)
# nabla p = -3/Re * [sin(x)cos(t+y)cos(z), cos(x)sin(t+y)cos(z), cos(x)cos(t+y)sin(z)]^T

def compute_analytical_gradP_x_at_x_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    # dp/dx = - (3/Re) * sin(x) * cos(t+y) * cos(z)
    x = points[:, 0] + h / 2.0
    y = points[:, 1]
    z = points[:, 2]
    dp_dx = -(3.0 / RE) * np.sin(x) * np.cos(t + y) * np.cos(z)
    return dp_dx

def compute_analytical_gradP_y_at_y_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    # dp/dy = - (3/Re) * cos(x) * sin(t+y) * cos(z)
    x = points[:, 0]
    y = points[:, 1] + h / 2.0
    z = points[:, 2]
    dp_dy = -(3.0 / RE) * np.cos(x) * np.sin(t + y) * np.cos(z)
    return dp_dy

def compute_analytical_gradP_z_at_z_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    # dp/dz = - (3/Re) * cos(x) * cos(t+y) * sin(z)
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2] + h / 2.0  
    dp_dz = -(3.0 / RE) * np.cos(x) * np.cos(t + y) * np.sin(z)
    return dp_dz


# --- 3. ERROR CALCULATION HELPERS ---
def calculate_l2_rms(data: np.ndarray) -> float:
    data_sq = data**2
    if data.ndim == 2:
        data_sq = np.sum(data_sq, axis=1)
    mean_squared_value = np.mean(data_sq)
    return np.sqrt(mean_squared_value + 1e-16)

def calculate_l2_rms_error(field_numerical: np.ndarray, field_analytical: np.ndarray) -> float:
    error_field = field_numerical - field_analytical
    rms_error = calculate_l2_rms(error_field)
    return rms_error

def calculate_l1_mean_error(field_numerical: np.ndarray, field_analytical: np.ndarray) -> float:
    error_field = field_numerical - field_analytical
    # L1 norm (Mean Absolute Error)
    l1_error = np.mean(np.abs(error_field))
    return l1_error

# --- 4. PLOTTING FUNCTION ---

def plot_convergence(x_values: list, error_data: dict, title: str, save_filename: str, x_label: str, expected_order: int):
    x_arr = np.array(x_values)

    if len(x_arr) == 0:
        print("No data to plot.", file=sys.stderr)
        return

    sort_idx = np.argsort(x_arr)
    x_arr = x_arr[sort_idx]

    plottable_error_data = {}
    for label, errors in error_data.items():
        if len(errors) == len(x_arr):
            plottable_error_data[label] = np.array(errors)[sort_idx]
        else:
            print(f"WARNING: Skipping '{label}' (dim mismatch)", file=sys.stderr)

    if not plottable_error_data:
        return

    fig, ax = plt.subplots(figsize=(10, 7))
    markers = ['o', 's', '^', 'D', 'v', '*']

    # Plotting Data
    for i, (label, errors) in enumerate(plottable_error_data.items()):
        ax.plot(x_arr, errors, marker=markers[i % len(markers)], markersize=8, linestyle='-', label=label)

    # --- Reference Lines ---
    if len(x_arr) > 1:
        x_ref = x_arr[-1]
        err_ref = list(plottable_error_data.values())[0][-1]
        x_line = np.array([x_arr.min(), x_arr.max()])

        C_order = err_ref / (x_ref**expected_order)
        ax.plot(x_line, C_order * (x_line**expected_order), 'k--', label=f'O($x^{expected_order}$) Reference') 

        C_1 = err_ref / (x_ref**1)
        ax.plot(x_line, C_1 * (x_line**1), 'k:', label='O($x$) Reference') 

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(x_label, fontsize=14) 
    ax.set_ylabel("Error", fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.grid(True, which="both", linestyle=':', linewidth=0.5)
    ax.legend(fontsize=12)

    plt.tight_layout()
    print(f"Saving plot to {save_filename}...")
    plt.savefig(save_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)

# --- 5. MAIN ANALYSIS LOOP ---

def run_analysis(custom_simulations=None):
    if custom_simulations is not None:
        sims_to_analyze = custom_simulations
    else:
        sims_to_analyze = SIMULATIONS
        for s in sims_to_analyze:
            s['dt'] = DT
            s['mode'] = 'SPATIAL_ONLY'

    if not sims_to_analyze:
        print("No data to analyze.")
        return

    study_mode = sims_to_analyze[0].get('mode', 'SPATIAL_ONLY')
    print(f"--- Starting Convergence Analysis (Mode: {study_mode}) ---")
    
    # Lists to store results
    nx_values = []
    dt_values = []
    h_values = []
    
    # Error containers
    errors_u = []; errors_v = []; errors_w = []
    errors_p = []
    errors_p_l1 = [] # New container for Pressure L1 error
    errors_gp_x = []; errors_gp_y = []; errors_gp_z = []
    errors_div_u = []

    # Sort logic
    if study_mode == "TEMPORAL_ONLY":
        sorted_simulations = sorted(sims_to_analyze, key=lambda s: s['dt'])
    else:
        sorted_simulations = sorted(sims_to_analyze, key=lambda s: s['nx'])
    
    # Update header
    print(f"{'Nx':>4} | {'h':>10} | {'L2_Err_U':>10} | {'L2_Err_P':>10} | {'L1_Err_P':>10} | {'L2_Err_GradP':>12} | {'L2_Err_DivU':>12}")
    print("-" * 90)

    for sim in sorted_simulations:
        nx = sim['nx']
        dt = sim['dt']
        filepath = sim['file']

        try:
            step_number = get_step_from_filename(filepath)
            current_time = step_number * dt
        except ValueError as e:
            print(f"ERROR: {e}", file=sys.stderr); continue

        h = DOMAIN_LENGTH_X / (nx - 0.5)
        nx_values.append(nx); dt_values.append(dt); h_values.append(h)

        try:
            grid = pv.read(filepath)
            points = grid.points
            
            # --- 1. Velocity Error ---
            vel_numerical = grid.point_data[VELOCITY_FIELD_NAME]
            u_ana = compute_analytical_u_at_x_faces(points, current_time, h)
            v_ana = compute_analytical_v_at_y_faces(points, current_time, h)
            w_ana = compute_analytical_w_at_z_faces(points, current_time, h)

            errors_u.append(calculate_l2_rms_error(vel_numerical[:, 0], u_ana))
            errors_v.append(calculate_l2_rms_error(vel_numerical[:, 1], v_ana))
            errors_w.append(calculate_l2_rms_error(vel_numerical[:, 2], w_ana))

            # --- 2. Pressure Error ---
            p_numerical = grid.point_data[PRESSURE_FIELD_NAME]
            p_ana = compute_analytical_p_at_centers(points, current_time)
            
            # L2 Error
            errors_p.append(calculate_l2_rms_error(p_numerical, p_ana))
            # L1 Error (New)
            errors_p_l1.append(calculate_l1_mean_error(p_numerical, p_ana))

            # --- 3. Pressure Gradient Error ---
            if GRAD_PRESSURE_FIELD_NAME in grid.point_data:
                gradp_numerical = grid.point_data[GRAD_PRESSURE_FIELD_NAME]
                
                gp_x_ana = compute_analytical_gradP_x_at_x_faces(points, current_time, h)
                gp_y_ana = compute_analytical_gradP_y_at_y_faces(points, current_time, h)
                gp_z_ana = compute_analytical_gradP_z_at_z_faces(points, current_time, h)

                e_gpx = calculate_l2_rms_error(gradp_numerical[:, 0], gp_x_ana)
                e_gpy = calculate_l2_rms_error(gradp_numerical[:, 1], gp_y_ana)
                e_gpz = calculate_l2_rms_error(gradp_numerical[:, 2], gp_z_ana)

                errors_gp_x.append(e_gpx)
                errors_gp_y.append(e_gpy)
                errors_gp_z.append(e_gpz)
                
                total_grad_err = np.sqrt(e_gpx**2 + e_gpy**2 + e_gpz**2)
            else:
                total_grad_err = 0.0
                if len(h_values) == 1: 
                    print(f"Warning: Field '{GRAD_PRESSURE_FIELD_NAME}' not found.", file=sys.stderr)

            # --- 4. Divergence Error ---
            if DIV_VELOCITY_FIELD_NAME in grid.point_data:
                div_u_num = grid.point_data[DIV_VELOCITY_FIELD_NAME]
                div_u_ana = np.zeros_like(div_u_num)
                err_div = calculate_l2_rms_error(div_u_num, div_u_ana)
                errors_div_u.append(err_div)
            else:
                errors_div_u.append(0.0)
                if len(h_values) == 1:
                    print(f"Warning: Field '{DIV_VELOCITY_FIELD_NAME}' not found.", file=sys.stderr)

        except Exception as e:
            print(f"ERROR processing file {filepath}: {e}", file=sys.stderr); continue

        # Print summary line
        print(f"{nx:>4} | {h:>10.4e} | {errors_u[-1]:>10.4e} | {errors_p[-1]:>10.4e} | {errors_p_l1[-1]:>10.4e} | {total_grad_err:>12.4e} | {errors_div_u[-1]:>12.4e}")

    print("-" * 90)
    print("Analysis complete. Saving plots...\n")

    # --- 6. PLOTTING ---
    if not nx_values:
        return

    if study_mode == "TEMPORAL_ONLY":
        x_values = dt_values; x_label = 'Time Step $\Delta t$ (s)'; expected_order = 2 
        filename_suffix = "_temporal_only"
    else:
        x_values = h_values; x_label = 'Grid Size $h$ (m)'; expected_order = 2 
        filename_suffix = f"_{study_mode.lower()}"
    
    # Plot 1: Velocity
    plot_convergence(x_values, 
                     {'u': errors_u, 'v': errors_v, 'w': errors_w},
                     f'Convergence: Velocity ({study_mode})',
                     f"../results/convergence_velocity{filename_suffix}.png",
                     x_label, expected_order)

    # Plot 2: Pressure (L2)
    plot_convergence(x_values, 
                     {'p (L2)': errors_p},
                     f'Convergence: Pressure L2 ({study_mode})',
                     f"../results/convergence_pressure{filename_suffix}.png",
                     x_label, expected_order)

    # Plot 3: Pressure (L1 - New)
    plot_convergence(x_values, 
                     {'p (L1)': errors_p_l1},
                     f'Convergence: Pressure L1 ({study_mode})',
                     f"../results/convergence_pressure_l1{filename_suffix}.png",
                     x_label, expected_order)

    # Plot 4: Pressure Gradient
    if errors_gp_x:
        plot_convergence(x_values, 
                        {'GradP_x': errors_gp_x, 'GradP_y': errors_gp_y, 'GradP_z': errors_gp_z},
                        f'Convergence: Pressure Gradient ({study_mode})',
                        f"../results/convergence_gradP{filename_suffix}.png",
                        x_label, expected_order)
        
    # Plot 5: Divergence
    if errors_div_u and any(e > 0 for e in errors_div_u):
        plot_convergence(x_values,
                         {'Div U': errors_div_u},
                         f'Convergence: Divergence ({study_mode})',
                         f"../results/convergence_divU{filename_suffix}.png",
                         x_label, expected_order)

    print("\nAll plots saved successfully.")

if __name__ == "__main__":
    run_analysis()