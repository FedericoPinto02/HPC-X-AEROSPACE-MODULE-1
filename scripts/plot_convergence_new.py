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

DOMAIN_LENGTH_X = 6.0 
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
    return np.sin(x) * np.cos(t + y) * np.sin(z)

def compute_analytical_v_at_y_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1] + h / 2.0  
    z = points[:, 2]
    return np.cos(x) * np.sin(t + y) * np.sin(z)

def compute_analytical_w_at_z_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2] + h / 2.0  
    return 2.0 * np.cos(x) * np.cos(t + y) * np.cos(z)

def compute_analytical_p_at_centers(points: np.ndarray, t: float) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    return (3.0 / RE) * np.cos(x) * np.cos(t + y) * np.cos(z)

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
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_clean_convergence(x_values: list, error_data: dict, title_prefix: str, save_filename: str, x_label: str, expected_order: int, mode_label: str):
    # Configurazione stile MathText professionale
    plt.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "font.size": 11,
        "axes.labelsize": 12,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "lines.linewidth": 0.8,
        "lines.markersize": 5,
        "axes.linewidth": 0.6,
        "figure.autolayout": True
    })

    x_arr = np.array(x_values)
    if len(x_arr) == 0:
        return

    # Ordina i valori x (dal più piccolo al più grande)
    sort_idx = np.argsort(x_arr)
    x_arr = x_arr[sort_idx]

    plottable_error_data = {}
    for label, errors in error_data.items():
        if len(errors) == len(x_arr):
            plottable_error_data[label] = np.array(errors)[sort_idx]

    if not plottable_error_data:
        return

    fig, ax = plt.subplots(figsize=(6, 4.5))

    colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628']
    markers = ['o', 's', '^', 'D', 'v', 'p', '*']

    # Plot dei dati
    for i, (label, errors) in enumerate(plottable_error_data.items()):
        ax.loglog(x_arr, errors,
                  marker=markers[i % len(markers)],
                  linestyle='-',
                  color=colors[i % len(colors)],
                  label=label,
                  alpha=0.9)

    # Linee di Riferimento
    if len(x_arr) > 1:
        # Ancoraggio al punto più coarse (ultimo elemento dell'array ordinato)
        x_anchor = x_arr[-1]
        
        # Prende l'errore della prima serie disponibile a quel punto come ancora
        y_anchor = list(plottable_error_data.values())[0][-1]
        
        # Estensione della linea su tutto il range
        x_line = np.geomspace(x_arr.min(), x_arr.max(), 100)

        # Calcolo pendenza ideale (Order N) partendo dall'ancora
        # Formula: y = y_anchor * (x / x_anchor)^order
        y_ideal = y_anchor * (x_line / x_anchor)**expected_order
        ax.loglog(x_line, y_ideal, 'k--', linewidth=0.8,
                  label=rf"Ideal $\mathcal{{O}}(x^{{{expected_order}}})$")

        # Calcolo pendenza riferimento (Order 1) partendo dalla STESSA ancora
        if expected_order != 1:
            y_one = y_anchor * (x_line / x_anchor)**1
            ax.loglog(x_line, y_one, 'k:', linewidth=0.8,
                      label=r"Ref $\mathcal{O}(x)$")

    # Griglia e Assi
    ax.grid(True, which='major', linestyle='-', linewidth=0.4, color='#bbbbbb')
    ax.grid(True, which='minor', linestyle=':', linewidth=0.2, color='#dddddd')

    # Etichette (x_label usata direttamente per evitare conflitti LaTeX)
    ax.set_xlabel(x_label)
    ax.set_ylabel(r"$\|\epsilon\|_2$")
    
    # Nessun titolo impostato come richiesto
    
    ax.legend(loc='best', frameon=True, fancybox=False, edgecolor='black', framealpha=1)

    # Salvataggio
    base, ext = os.path.splitext(save_filename)
    if ext.lower() not in ['.pdf', '.svg', '.eps']:
        save_filename = base + ".pdf"
    
    plt.savefig(save_filename, format='pdf', bbox_inches='tight', dpi=300)
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
    
    print(f"{'Nx':>4} | {'h':>10}| {'dt':>10} | {'L2_Err_U':>12} | {'L2_Err_V':>12} | {'L2_Err_W':>12} | {'L2_Err_P':>12} | {'L1_Err_P':>12}")
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

        print(f"{nx:>4} | {h:>10.4e} | {dt:>10.4e} | {errors_u[-1]:>12.4e} | {errors_v[-1]:>12.4e} | {errors_w[-1]:>12.4e} | {errors_p[-1]:>12.4e} | {errors_p_l1[-1]:>12.4e}")

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
        x_plot = dt_values; x_lbl = r' $\Delta t$'
    else:
        x_plot = h_values; x_lbl = r'$h$'
    
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