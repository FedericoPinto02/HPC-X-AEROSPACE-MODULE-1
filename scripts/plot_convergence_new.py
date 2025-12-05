#!/usr/bin/env python3

"""
This script performs a spatial convergence analysis for a Method of
Manufactured Solutions (MMS) test.

It reads multiple VTK files, each from a simulation with a different
grid resolution. For each file, it computes the L2 RMS error
of velocity components (u, v, w) and pressure (p) against a
known analytical solution.

This version is adapted for staggered grid solvers where the VTK output
is co-located (e.g., at cell centers). It compares the numerical
co-located data against the analytical solution evaluated at the
appropriate staggered locations (e.g., cell faces).

Script modified from a previous version to:
1. Use only L2 RMS (absolute) errors.
2. Add O(h) and O(h^2) reference lines.
3. Translate all comments and outputs to English.
4. Remove all divergence calculations.
5. Save plots to high-definition image files instead of showing them.
"""

# --- 0. IMPORTS ---
import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import os
import sys
import re 
# Importiamo la libreria ticker per la formattazione
from matplotlib import ticker

# --- 1. USER CONFIGURATION (Mantenute le costanti per compatibilità) ---

# List of simulation results to analyze.
# NOTE: This list is now IGNORED if run_convergence_study_unified.py is used.
SIMULATIONS = [
{"nx": 30, "file": "../output/n30simulation_output_0030.vtk"},
{"nx": 20, "file": "../output/n20simulation_output_0030.vtk"},
# ... (altre simulazioni rimosse)
]

DOMAIN_LENGTH_X = 6
# DT is now superseded by the DT in the custom_simulations list
DT = 0.001 

VELOCITY_FIELD_NAME = "velocity"
PRESSURE_FIELD_NAME = "pressure"


# Helper function to get step number from filename
def get_step_from_filename(filepath: str) -> int:
    """
    [... Funzione get_step_from_filename non modificata ...]
    """
    # This regex looks for one or more digits (\d+) right before ".vtk"
    # It assumes the format is like "_123.vtk" or "abc_123.vtk"
    match = re.search(r'_(\d+)\.vtk$', os.path.basename(filepath))

    if match:
        return int(match.group(1))
    else:
        # If no match, raise an error
        raise ValueError(f"Could not parse step number from filename: {filepath}. "
                         "Expected format like '..._123.vtk'")


# --- 2. ANALYTICAL SOLUTION FUNCTIONS (STAGGERED) ---
# [Le funzioni analitiche u, v, w, p non sono modificate]
def compute_analytical_u_at_x_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    x = points[:, 0] + h / 2.0  
    y = points[:, 1]
    z = points[:, 2]
    u_ana = np.sin(t) * np.sin(x) * np.sin(y) * np.sin(z)
    return u_ana

def compute_analytical_v_at_y_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1] + h / 2.0  
    z = points[:, 2]
    v_ana = np.sin(t) * np.cos(x) * np.cos(y) * np.cos(z)
    return v_ana

def compute_analytical_w_at_z_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2] + h / 2.0  
    w_ana = np.sin(t) * np.cos(x) * np.sin(y) * (np.sin(z) + np.cos(z))
    return w_ana

def compute_analytical_p_at_centers(points: np.ndarray, t: float) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    p_ana = -np.sin(t) * 3.0 * np.cos(x) * np.sin(y) * (np.sin(z) - np.cos(z))
    return p_ana

# --- 3. ERROR CALCULATION HELPERS (Non modificati) ---
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

# --- 4. PLOTTING FUNCTION MODIFIED ---

def plot_convergence(x_values: list, error_data: dict, title: str, save_filename: str, x_label: str, expected_order: int):
    """
    Plots the convergence data using generic x-axis values and a dynamic label/order.
    """
    x_arr = np.array(x_values)

    if len(x_arr) == 0:
        print("No data to plot.", file=sys.stderr)
        return

    # Sort data based on x-axis values (Ascending)
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
        # Use the rightmost point (largest h or largest dt) as the anchor
        x_ref = x_arr[-1]
        err_ref = list(plottable_error_data.values())[0][-1]

        # Define the reference line span (min x-value to max x-value)
        x_line = np.array([x_arr.min(), x_arr.max()])

        # O(x^p) Reference (using 'expected_order' for the expected power)
        C_order = err_ref / (x_ref**expected_order)
        ax.plot(x_line, C_order * (x_line**expected_order), 'k--', label=f'O($x^{expected_order}$) Reference') 

        # O(x) Reference (Always useful for context)
        C_1 = err_ref / (x_ref**1)
        ax.plot(x_line, C_1 * (x_line**1), 'k:', label='O($x$) Reference') 

    # --- Formatting ---
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(x_label, fontsize=14) 
    ax.set_ylabel("L2 RMS Error", fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.grid(True, which="both", linestyle=':', linewidth=0.5)
    ax.legend(fontsize=12)

    plt.tight_layout()

    print(f"Saving plot to {save_filename}...")
    plt.savefig(save_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)

# --- 5. MAIN ANALYSIS LOOP MODIFIED ---

def run_analysis(custom_simulations=None):
    """
    Main function to loop through simulations, compute errors, and plot.
    Now supports different convergence modes.
    """

    if custom_simulations is not None:
        sims_to_analyze = custom_simulations
    else:
        # If run directly, use the global list and assume SPATIAL_ONLY
        sims_to_analyze = SIMULATIONS
        # Add default mode/dt for direct run compatibility (optional)
        for s in sims_to_analyze:
            s['dt'] = DT
            s['mode'] = 'SPATIAL_ONLY'

    if not sims_to_analyze:
        print("No data to analyze.")
        return

    # Determine mode from the first element
    study_mode = sims_to_analyze[0].get('mode', 'SPATIAL_ONLY')
    print(f"--- Starting Convergence Analysis (Mode: {study_mode}) ---")
    
    # Lists to store results
    nx_values = []
    dt_values = []
    h_values = []
    errors_u = []
    errors_v = []
    errors_w = []
    errors_p = []

    # Ensure simulations are sorted for a smooth plot
    if study_mode == "TEMPORAL_ONLY":
        # Sort by dt (smallest first)
        sorted_simulations = sorted(sims_to_analyze, key=lambda s: s['dt'])
    else:
        # Sort by nx (smallest first) which means largest h first
        sorted_simulations = sorted(sims_to_analyze, key=lambda s: s['nx'])
    
    # Print table header
    print(f"{'Nx':>4} | {'h':>10} | {'dt':>10} | {'Time':>10} | {'L2_Err_u':>10} | {'L2_Err_v':>10} | {'L2_Err_w':>10} | {'L2_Err_p':>10}")
    print("-" * 85)

    for sim in sorted_simulations:
        nx = sim['nx']
        dt = sim['dt']
        filepath = sim['file']

        # Calculate CURRENT_TIME and h
        try:
            step_number = get_step_from_filename(filepath)
            current_time = step_number * dt
        except ValueError as e:
            print(f"ERROR: {e}", file=sys.stderr); print("Skipping.", file=sys.stderr); continue

        h = DOMAIN_LENGTH_X / (nx - 0.5)
        
        # Store for plotting later
        nx_values.append(nx)
        dt_values.append(dt)
        h_values.append(h)

        # --- Load Data and Compute Errors (Logic remains the same) ---
        try:
            grid = pv.read(filepath)
            points = grid.points
            
            # 1. Velocity Error
            vel_numerical = grid.point_data[VELOCITY_FIELD_NAME]
            u_analytical = compute_analytical_u_at_x_faces(points, current_time, h)
            v_analytical = compute_analytical_v_at_y_faces(points, current_time, h)
            w_analytical = compute_analytical_w_at_z_faces(points, current_time, h)

            err_u = calculate_l2_rms_error(vel_numerical[:, 0], u_analytical)
            err_v = calculate_l2_rms_error(vel_numerical[:, 1], v_analytical)
            err_w = calculate_l2_rms_error(vel_numerical[:, 2], w_analytical)

            errors_u.append(err_u); errors_v.append(err_v); errors_w.append(err_w)

            # 2. Pressure Error
            p_numerical = grid.point_data[PRESSURE_FIELD_NAME]
            p_analytical = compute_analytical_p_at_centers(points, current_time)
            err_p = calculate_l2_rms_error(p_numerical, p_analytical)
            errors_p.append(err_p)

        except Exception as e:
            print(f"ERROR processing file {filepath}: {e}", file=sys.stderr); continue

        # --- Print Results to Terminal ---
        print(f"{nx:>4} | {h:>10.4e} | {dt:>10.4e} | {current_time:>10.4f} | {err_u:>10.4e} | {err_v:>10.4e} | {err_w:>10.4e} | {err_p:>10.4e}")

    print("-" * 85)
    print("Analysis complete. Saving plots...\n")

    # --- 6. DETERMINE PLOT PARAMETERS AND PLOT RESULTS ---

    if not nx_values:
        print("No data was successfully processed. Exiting.", file=sys.stderr)
        return

    # Logic to select X-axis and Expected Order
    if study_mode == "TEMPORAL_ONLY":
        x_values = dt_values
        x_label = 'Time Step $\Delta t$ (s)'
        # CN is O(dt^2)
        expected_order = 2 
        filename_suffix = "_temporal_only"
    else: # SPATIAL_ONLY or SPACE_TIME_COUPLED
        x_values = h_values
        x_label = 'Grid Size $h$ (m)'
        # DF is O(h^2)
        expected_order = 2 
        filename_suffix = f"_{study_mode.lower()}"
    
    # Plot 1: Velocity
    vel_error_dict = {'u-velocity': errors_u, 'v-velocity': errors_v, 'w-velocity': errors_w}
    plot_convergence(x_values,
                     vel_error_dict,
                     f'Convergence of Velocity Components ({study_mode})',
                     save_filename=f"../results/convergence_velocity{filename_suffix}.png",
                     x_label=x_label,
                     expected_order=expected_order)

    # Plot 2: Pressure
    p_error_dict = {'Pressure': errors_p}
    plot_convergence(x_values,
                     p_error_dict,
                     f'Convergence of Pressure ({study_mode})',
                     save_filename=f"../results/convergence_pressure{filename_suffix}.png",
                     x_label=x_label,
                     expected_order=expected_order)

    print("\nAll plots saved successfully.")


if __name__ == "__main__":
    # If run directly without arguments, it defaults to the old behavior (Spatial Only)
    run_analysis()