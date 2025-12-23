#!/usr/bin/env python3

"""
Script for Spatial Convergence Analysis (Method of Manufactured Solutions).

This script reads VTK files from simulations with varying grid resolutions
and computes the L2 and L1 error norms against an analytical solution.

Changes in this version:
- Fixed SyntaxWarning regarding invalid escape sequences in LaTeX strings.
- Reverted to standard sans-serif fonts (cleaner look).
- High-contrast colors and explicit log-log grids.
- Dynamic LaTeX titles indicating scaling (e.g., dt ~ dx).
- Calculates and prints Convergence Order table.
- UPDATED: Thinner lines and new color palette.
"""

# --- 0. IMPORTS ---
import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import re 
from matplotlib.ticker import LogLocator, NullFormatter

# --- 1. PLOTTING STYLE CONFIGURATION (CLEAN & SHARP) ---
plt.rcParams.update({
    'font.family': 'sans-serif',     # Clean standard font
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 15,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'lines.linewidth': 1.5,          # UPDATED: Thinner lines (was 2.5)
    'lines.markersize': 8,           # Slightly smaller markers
    'figure.figsize': (10, 8),       # Larger figure
    'savefig.dpi': 400,              # Very high resolution
    'axes.grid': True,
    'grid.alpha': 0.5,
})

# --- 2. USER CONFIGURATION ---

# Physical Constants
RE = 1.0

# List of simulation results.
SIMULATIONS = [
    {"nx": 30, "file": "../output/n30simulation_output_0030.vtk"},
    {"nx": 20, "file": "../output/n20simulation_output_0030.vtk"},
]

DOMAIN_LENGTH_X = 1.0 
DT = 0.001            

# Field names in VTK
VELOCITY_FIELD_NAME = "velocity"
PRESSURE_FIELD_NAME = "pressure"

# --- 3. HELPER FUNCTIONS ---

def get_step_from_filename(filepath: str) -> int:
    filename = os.path.basename(filepath)

    # 1. Try matching new format: ..._RANK_STEP.vtk
    # We capture two groups of digits. Group 1 is RANK, Group 2 is STEP.
    match_rank = re.search(r'_(\d+)_(\d+)\.vtk$', filename)
    if match_rank:
        return int(match_rank.group(2)) # Return the STEP

    # 2. Fallback to legacy format: ..._STEP.vtk
    match_legacy = re.search(r'_(\d+)\.vtk$', filename)
    if match_legacy:
        return int(match_legacy.group(1))

    raise ValueError(f"Could not parse step number from filename: {filepath}.")

def compute_convergence_orders(x_values: list, errors: list) -> tuple:
    x = np.array(x_values)
    e = np.array(errors)
    if len(x) < 2:
        return 0.0, []
    
    # Global Order (Linear Regression)
    coeffs = np.polyfit(np.log(x), np.log(e), 1)
    global_order = coeffs[0]

    # Local Orders
    local_orders = []
    for i in range(1, len(x)):
        p = np.log(e[i] / e[i-1]) / np.log(x[i] / x[i-1])
        local_orders.append(p)
        
    return global_order, local_orders

# --- 4. ANALYTICAL SOLUTION FUNCTIONS ---

def compute_analytical_u_at_x_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    x = points[:, 0] + h / 2.0 
    y = points[:, 1]
    z = points[:, 2]
    U_max = 5.0
    L = 1.0
    nu = 1.0
    G = 100.0

    return np.sin(t) * U_max / L * y - y / nu * G *(L - y) *0.5 * np.sin(t)

def compute_analytical_v_at_y_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1] + h / 2.0  
    z = points[:, 2]
    return 0

def compute_analytical_w_at_z_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2] + h / 2.0  
    return 0

def compute_analytical_p_at_centers(points: np.ndarray, t: float) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    return 10

# --- 5. ERROR CALCULATION ---

def calculate_l2_rms(data: np.ndarray) -> float:
    data_sq = data**2
    if data.ndim == 2:
        data_sq = np.sum(data_sq, axis=1)
    return np.sqrt(np.mean(data_sq) + 1e-16)

def calculate_l2_rms_error(field_numerical: np.ndarray, field_analytical: np.ndarray) -> float:
    return calculate_l2_rms(field_numerical - field_analytical)

def calculate_l1_mean_error(field_numerical: np.ndarray, field_analytical: np.ndarray) -> float:
    return np.mean(np.abs(field_numerical - field_analytical))

# --- 6. PLOTTING FUNCTION (IMPROVED) ---

def plot_clean_convergence(x_values: list, error_data: dict, title_prefix: str, save_filename: str, x_label: str, expected_order: int, mode_label: str):
    """
    Generates a clean, high-contrast log-log plot.
    """
    x_arr = np.array(x_values)
    if len(x_arr) == 0:
        return

    # Sort
    sort_idx = np.argsort(x_arr)
    x_arr = x_arr[sort_idx]

    # Prepare Data
    plottable_error_data = {}
    for label, errors in error_data.items():
        if len(errors) == len(x_arr):
            plottable_error_data[label] = np.array(errors)[sort_idx]

    if not plottable_error_data:
        return

    fig, ax = plt.subplots()

    # UPDATED: Standard Tab10 Palette (Blue, Orange, Green, Red, Purple)
    # Distinct and classic
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'] 
    markers = ['o', 's', '^', 'D', 'v'] 
    
    # Plot Data
    for i, (label, errors) in enumerate(plottable_error_data.items()):
        ax.loglog(x_arr, errors, 
                  marker=markers[i % len(markers)], 
                  linestyle='-', 
                  color=colors[i % len(colors)],
                  label=label, alpha=0.9)

    # Reference Lines
    if len(x_arr) > 1:
        x_ref = x_arr[-1]
        err_ref = list(plottable_error_data.values())[0][-1]
        x_line = np.geomspace(x_arr.min(), x_arr.max(), 100)

        # Order N
        C_order = err_ref / (x_ref**expected_order)
        y_order = C_order * (x_line**expected_order)
        # Fix: Use raw f-string (rf) to handle LaTeX backslashes correctly
        ax.loglog(x_line, y_order, 'k--', linewidth=1.5, label=rf'Ideal Slope {expected_order}: $\mathcal{{O}}(x^{{{expected_order}}})$') 

        # Order 1 (if expected is not 1)
        if expected_order != 1:
            C_1 = err_ref / (x_ref**1)
            y_1 = C_1 * (x_line**1)
            # Fix: Use raw string (r) to handle LaTeX backslashes correctly
            ax.loglog(x_line, y_1, 'k:', linewidth=1.5, label=r'Reference Slope 1: $\mathcal{O}(x)$')

    # Grid - Crucial for Log-Log
    ax.grid(True, which='major', linestyle='-', linewidth=0.7, color='gray', alpha=0.5)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.5, color='gray', alpha=0.3)
    
    # Title Logic
    full_title = f"{title_prefix}\n{mode_label}"
    ax.set_title(full_title)
    
    ax.set_xlabel(x_label)
    ax.set_ylabel("RMS Error ($L_2$ Norm)")
    ax.legend(loc='best', framealpha=0.9, edgecolor='gray')

    print(f"Saving high-quality plot to {save_filename}...")
    plt.tight_layout()
    plt.savefig(save_filename, bbox_inches='tight')
    plt.close(fig)

# --- 7. MAIN LOOP ---

def run_analysis(custom_simulations=None):
    if custom_simulations is not None:
        sims_to_analyze = custom_simulations
    else:
        sims_to_analyze = SIMULATIONS
        for s in sims_to_analyze:
            s['dt'] = DT
            s['mode'] = 'SPATIAL_ONLY' # Change this if running coupled tests

    if not sims_to_analyze:
        print("No data.")
        return

    # Check variation to determine label
    dts = [s['dt'] for s in sims_to_analyze]
    nxs = [s['nx'] for s in sims_to_analyze]
    dt_varies = len(set(dts)) > 1
    nx_varies = len(set(nxs)) > 1
    
    # Determine Label Mode
    if dt_varies and nx_varies:
        study_mode = "COUPLED"
        mode_label_latex = r"Scaling: $\Delta t \sim \Delta x$"
        filename_suffix = "_coupled"
    elif dt_varies:
        study_mode = "TEMPORAL_ONLY"
        mode_label_latex = r"Temporal Refinement ($\Delta x$ fixed)"
        filename_suffix = "_temporal"
    else:
        study_mode = "SPATIAL_ONLY"
        mode_label_latex = r"Spatial Refinement ($\Delta t$ fixed)"
        filename_suffix = "_spatial"

    print(f"--- Starting Analysis (Detected Mode: {study_mode}) ---")
    
    nx_values, dt_values, h_values = [], [], []
    errors_u, errors_v, errors_w = [], [], []
    errors_p, errors_p_l1 = [], []

    # Sort based on the primary variable
    if study_mode == "TEMPORAL_ONLY":
        sorted_simulations = sorted(sims_to_analyze, key=lambda s: s['dt'])
    else:
        sorted_simulations = sorted(sims_to_analyze, key=lambda s: s['nx'])
    
    print(f"{'Nx':>4} | {'h':>10} | {'L2_Err_U':>12} | {'L2_Err_V':>12} | {'L2_Err_W':>12} | {'L2_Err_P':>12} | {'L1_Err_P':>12}")
    print("-" * 95)

    for sim in sorted_simulations:
        nx = sim['nx']; dt = sim['dt']; filepath = sim['file']
        try:
            step = get_step_from_filename(filepath)
            time = step * dt
        except ValueError as e:
            print(f"Skip {filepath}: {e}", file=sys.stderr); continue

        h = DOMAIN_LENGTH_X / (nx - 0.5)
        nx_values.append(nx); dt_values.append(dt); h_values.append(h)

        try:
            grid = pv.read(filepath)
            pts = grid.points
            
            # Velocity
            if VELOCITY_FIELD_NAME in grid.point_data:
                vel = grid.point_data[VELOCITY_FIELD_NAME]
                u_a = compute_analytical_u_at_x_faces(pts, time, h)
                v_a = compute_analytical_v_at_y_faces(pts, time, h)
                w_a = compute_analytical_w_at_z_faces(pts, time, h)
                errors_u.append(calculate_l2_rms_error(vel[:, 0], u_a))
                errors_v.append(calculate_l2_rms_error(vel[:, 1], v_a))
                errors_w.append(calculate_l2_rms_error(vel[:, 2], w_a))
            
            # Pressure
            if PRESSURE_FIELD_NAME in grid.point_data:
                p_num = grid.point_data[PRESSURE_FIELD_NAME]
                p_a = compute_analytical_p_at_centers(pts, time)
                errors_p.append(calculate_l2_rms_error(p_num, p_a))
                errors_p_l1.append(calculate_l1_mean_error(p_num, p_a))

        except Exception as e:
            print(f"Error {filepath}: {e}", file=sys.stderr); continue

        print(f"{nx:>4} | {h:>10.4e} | {errors_u[-1]:>12.4e} | {errors_v[-1]:>12.4e} | {errors_w[-1]:>12.4e} | {errors_p[-1]:>12.4e} | {errors_p_l1[-1]:>12.4e}")

    print("-" * 95)

    # --- ORDER ANALYSIS ---
    x_analysis = dt_values if study_mode == "TEMPORAL_ONLY" else h_values
    
    print("\n--- Convergence Order Analysis ---")
    print(f"{'Variable':<20} | {'Global Order':<12} | {'Local Orders':<40}")
    print("-" * 90)

    for name, errs in [("Vel U", errors_u), ("Vel V", errors_v), ("Vel W", errors_w), ("Press L2", errors_p), ("Press L1", errors_p_l1)]:
        if errs:
            g_ord, l_ords = compute_convergence_orders(x_analysis, errs)
            l_str = ", ".join([f"{o:.2f}" for o in l_ords]) if l_ords else "N/A"
            print(f"{name:<20} | {g_ord:>12.2f} | {l_str}")

    print("-" * 90)

    # --- PLOTTING ---
    if not nx_values: return

    if study_mode == "TEMPORAL_ONLY":
        x_plot = dt_values; x_lbl = r'Time Step $\Delta t$'
    else:
        x_plot = h_values; x_lbl = r'Grid Size $h$'
    
    # 1. Velocity
    plot_clean_convergence(x_plot, 
                     {'U-Component': errors_u, 'V-Component': errors_v, 'W-Component': errors_w},
                     'Velocity Convergence',
                     f"../results/convergence_velocity{filename_suffix}.png",
                     x_lbl, 2, mode_label_latex)

    # 2. Pressure
    plot_clean_convergence(x_plot, 
                     {'Pressure ($L_2$)': errors_p, 'Pressure ($L_1$)': errors_p_l1},
                     'Pressure Convergence',
                     f"../results/convergence_pressure{filename_suffix}.png",
                     x_lbl, 2, mode_label_latex)

    print("\nHigh-quality plots saved.")

if __name__ == "__main__":
    run_analysis()